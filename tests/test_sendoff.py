# mypy: allow-untyped-decorators
"""Test sendoff parsing of metadata on generated SDFs."""
import copy
import errno
import io
import os
import pickle
import threading
import typing
from itertools import chain
from pathlib import Path
from typing import Generator, Iterable, Iterator, cast

import pytest

from sendoff.ctable import (
    IndicesDuplicateError,
    IndicesMismatchError,
    IndicesOutOfOrderError,
)
from sendoff.sdblock import Pathy, SDBlock, parse_sdf, write_sdf


def test_single_mol_read(single_mol_sdf: Pathy) -> None:
    """An sdf of a single molecule contains one molecule.

    Args:
        single_mol_sdf: pytest fixture of a Path to the sdf
    """
    mols = list(parse_sdf(single_mol_sdf))
    assert len(mols) == 1


def test_double_mol_read(double_mol_sdf: Pathy) -> None:
    """An sdf of two molecule contains two molecules.

    Args:
        double_mol_sdf: pytest fixture of a Path to the sdf
    """
    mols = list(parse_sdf(double_mol_sdf))
    assert len(mols) == 2


@pytest.mark.parametrize("trailing_newline", ["", "\n"])
def test_unterminated_final_record_is_returned(
    tmp_path: Path, trailing_newline: str
) -> None:
    """Return a final block without a delimiter, with or without a newline.

    Args:
        tmp_path: Temporary directory for the SDF fixture.
        trailing_newline: Newline appended after the unterminated value.
    """
    sdf = tmp_path / "unterminated.sdf"
    sdf.write_text(
        "Closed\n\nM  END\n$$$$\n"
        "Title\n\nM  END\n> <Record>\nValue" + trailing_newline
    )

    blocks = list(parse_sdf(sdf))

    assert [block.title for block in blocks] == ["Closed", "Title"]
    assert blocks[-1].records() == [("Record", "Value")]


def test_empty_input_has_no_blocks() -> None:
    """Do not synthesize a block for empty input."""
    assert list(SDBlock.from_lines([])) == []


def test_bytes_lines_preserve_invalid_bytes() -> None:
    """Byte-backed blocks retain non-UTF-8 title and metadata octets."""
    raw = [
        b"\xff title \n",
        b"\n",
        b"M  END\n",
        b"> <field\xfe>\n",
        b"value\xfd\n",
        b"\n",
        b"$$$$\n",
        b"second\n",
        b"\n",
        b"M  END\n",
    ]

    blocks = SDBlock.from_lines(raw)
    block = next(blocks)
    tail = next(blocks)

    assert block.title == b"\xff title"
    assert block.records() == [(b"field\xfe", b"value\xfd")]
    assert tail.title == b"second"
    output = io.BytesIO()
    block.write(output)
    assert b"\xff" in output.getvalue()
    assert b"field\xfe" in output.getvalue()
    assert b"value\xfd" in output.getvalue()
    reparsed = next(SDBlock.from_lines(output.getvalue().splitlines(keepends=True)))
    assert reparsed.title == block.title
    assert reparsed.records() == block.records()

    with pytest.raises(NotImplementedError, match="byte-backed"):
        block.ctable()  # type: ignore[misc]


def test_bytes_path_uses_binary_input(tmp_path: Path) -> None:
    """A bytes filesystem path yields bytes and preserves invalid UTF-8.

    Args:
        tmp_path: Temporary directory for the SDF fixture.
    """
    path = tmp_path / "bytes.sdf"
    path.write_bytes(b"\xff title\n\nM  END\n> <field>\nvalue\xfe\n\n$$$$\n")

    block = next(parse_sdf(os.fsencode(path)))

    assert block.title == b"\xff title"
    assert block.records() == [(b"field", b"value\xfe")]


def test_text_path_preserves_text_values(tmp_path: Path) -> None:
    """A str filesystem path continues to use text decoding.

    Args:
        tmp_path: Temporary directory for the SDF fixture.
    """
    path = tmp_path / "text.sdf"
    path.write_text("title é\n\nM  END\n> <field>\nvalue é\n\n$$$$\n")

    block = next(parse_sdf(path))

    assert block.title == "title é"
    assert block.records() == [("field", "value é")]


def test_from_lines_rejects_mixed_text_and_bytes() -> None:
    """A single line source cannot mix text and bytes."""
    mixed_lines = cast(Iterable[str], [b"title\n", "M  END\n"])
    with pytest.raises(TypeError, match="mixed str and bytes"):
        list(SDBlock.from_lines(mixed_lines))


def test_from_lines_yields_incrementally_and_retains_blocks() -> None:
    """Only consume through the next delimiter and retain already parsed data."""
    consumed: list[str] = []
    first_record = [
        "First\n",
        "\n",
        "M  END\n",
        "> <Record>\n",
        "one\n",
        "\n",
        "$$$$\n",
    ]
    second_record = ["Second\n", "\n", "M  END\n", "$$$$\n"]

    def lines() -> Iterator[str]:
        for line in first_record + second_record:
            consumed.append(line)
            yield line

    blocks = SDBlock.from_lines(lines())

    assert consumed == []
    assert iter(blocks) is blocks
    first = next(blocks)
    assert consumed == first_record

    second = next(blocks)

    assert first.title == "First"
    assert first.records() == [("Record", "one")]
    assert second.title == "Second"
    assert consumed == first_record + second_record


def test_from_lines_surfaces_source_errors_on_advance() -> None:
    """Yield completed blocks before surfacing a later source error."""

    def lines() -> Iterator[str]:
        yield "First\n"
        yield "\n"
        yield "M  END\n"
        yield "$$$$\n"
        raise OSError("read failed")

    blocks = SDBlock.from_lines(lines())

    assert next(blocks).title == "First"
    with pytest.raises(OSError, match="read failed"):
        next(blocks)
    with pytest.raises(StopIteration):
        next(blocks)


def test_from_lines_does_not_close_caller_owned_iterable() -> None:
    """Abandoning the parser leaves the caller's source open."""
    closed: list[bool] = []

    def lines() -> Generator[str, None, None]:
        try:
            yield "First\n"
            yield "\n"
            yield "M  END\n"
            yield "$$$$\n"
        finally:
            closed.append(True)

    source = lines()
    blocks = SDBlock.from_lines(source)
    assert next(blocks).title == "First"
    del blocks

    assert closed == []
    source.close()
    assert closed == [True]


def open_paths() -> set[str]:
    """Return the paths this process has open or mapped, from /proc/self.

    Returns:
        Resolved targets of the open file descriptors, and mapped files.
    """
    paths = set()
    for fd in os.listdir("/proc/self/fd"):
        try:
            paths.add(os.readlink(os.path.join("/proc/self/fd", fd)))
        except OSError:
            continue
    with open("/proc/self/maps") as maps:
        for line in maps:
            fields = line.split(maxsplit=5)
            if len(fields) == 6:
                paths.add(fields[5].rstrip("\n"))
    return paths


@pytest.mark.skipif(not os.path.isdir("/proc/self/fd"), reason="needs /proc")
def test_parse_sdf_opens_lazily_and_close_releases_the_file(tmp_path: Path) -> None:
    """The file opens on the first next() and close() releases it.

    Args:
        tmp_path: Temporary directory for the SDF fixture.
    """
    path = tmp_path / "two.sdf"
    path.write_text("First\n\nM  END\n$$$$\nSecond\n\nM  END\n$$$$\n")

    blocks = parse_sdf(path)
    assert str(path) not in open_paths()
    assert next(blocks).title == "First"
    assert str(path) in open_paths()
    blocks.close()
    assert str(path) not in open_paths()
    with pytest.raises(StopIteration):
        next(blocks)

    assert len(list(parse_sdf(path))) == 2
    assert str(path) not in open_paths()


@pytest.mark.skipif(not hasattr(os, "mkfifo"), reason="needs named pipes")
def test_parse_sdf_reads_files_that_do_not_map(tmp_path: Path) -> None:
    """Pipes and empty files, which are not mapped, are read instead.

    Args:
        tmp_path: Temporary directory for the pipe and the empty file.
    """
    pipe = tmp_path / "pipe.sdf"
    os.mkfifo(pipe)
    data = b"a\r\nM  END\r\n$$$$\r\nb\r$$$$\r"
    writer = threading.Thread(target=pipe.write_bytes, args=(data,), daemon=True)
    writer.start()
    assert [block.title for block in parse_sdf(pipe)] == ["a", "b"]
    writer.join()
    empty = tmp_path / "empty.sdf"
    empty.write_bytes(b"")
    assert list(parse_sdf(empty)) == []


def test_parse_sdf_raises_open_errors_on_first_next(tmp_path: Path) -> None:
    """Open errors surface lazily as the OSError subclass open() raises.

    Args:
        tmp_path: Temporary directory without the SDF.
    """
    missing = tmp_path / "missing.sdf"
    blocks = parse_sdf(missing)

    with pytest.raises(FileNotFoundError) as info:
        next(blocks)
    assert info.value.errno == errno.ENOENT
    assert info.value.filename == str(missing)
    with pytest.raises(IsADirectoryError):
        next(parse_sdf(os.fsencode(tmp_path)))
    with pytest.raises(TypeError):
        parse_sdf(1)  # type: ignore[call-overload]


def test_parse_sdf_text_paths_decode_utf8_on_access(tmp_path: Path) -> None:
    """A non-UTF-8 feed still splits on a str path; only bad values raise.

    Args:
        tmp_path: Temporary directory for the SDF fixture.
    """
    path = tmp_path / "latin1.sdf"
    path.write_bytes(
        b"first\n\nM  END\n> <F>\n\xa6 bad\n\n$$$$\nsecond\n\nM  END\n$$$$\n"
    )

    blocks = list(parse_sdf(path))

    assert [block.title for block in blocks] == ["first", "second"]
    with pytest.raises(UnicodeDecodeError):
        blocks[0].records()
    assert next(parse_sdf(os.fsencode(path))).records() == [(b"F", b"\xa6 bad")]


def test_from_lines_reads_file_objects_in_chunks() -> None:
    """File objects are read with read1 or read rather than line iteration."""

    class Reader:
        def __init__(self, data: bytes) -> None:
            self.data = data
            self.sizes: list[int] = []

        def read1(self, size: int) -> bytes:
            self.sizes.append(size)
            chunk, self.data = self.data[:size], self.data[size:]
            return chunk

        def __iter__(self) -> Iterator[bytes]:
            raise AssertionError("file objects are not iterated by line")

    reader = Reader(b"First\n\nM  END\n$$$$\nSecond\n\nM  END\n$$$$\n")
    assert [block.title for block in SDBlock.from_lines(reader)] == [
        b"First",
        b"Second",
    ]
    assert reader.sizes and min(reader.sizes) >= 1 << 16


def test_from_lines_rejects_non_string_reads() -> None:
    """A read method returning neither str nor bytes is a TypeError."""

    class Reader:
        def read(self, size: int) -> list[str]:
            return ["First\n"]

    with pytest.raises(TypeError, match="must read str or bytes"):
        next(SDBlock.from_lines(Reader()))  # type: ignore[arg-type]


def test_bytes_block_requires_binary_append_values() -> None:
    """Appending to a bytes block cannot introduce implicitly encoded text."""
    block = next(SDBlock.from_lines([b"title\n", b"\n", b"M  END\n"]))

    with pytest.raises(TypeError, match="record names must be bytes"):
        block.append_record("field", b"value")  # type: ignore[arg-type]
    with pytest.raises(TypeError, match="values must be bytes"):
        block.append_record(b"field", "value")  # type: ignore[arg-type]

    block.append_record(b"field", b"value")
    assert block.records() == [(b"field", b"value")]


def test_block_without_mend_round_trips() -> None:
    """A block without M  END writes one delimiter, so rereading is stable."""
    text = "broken\nx\n$$$$\nt2\nM  END\n> <A>\n1\n\n$$$$\n"
    for _ in range(2):
        blocks = list(SDBlock.from_lines(io.StringIO(text)))
        assert [(b.title, b.records()) for b in blocks] == [
            ("broken", []),
            ("t2", [("A", "1")]),
        ]
        assert blocks[0].molblock == "broken\nx\n"
        out = io.StringIO()
        for block in blocks:
            block.write(out)
        assert out.getvalue() == text
        text = out.getvalue()


def test_titled_mol_read(single_titled_mol_sdf: Pathy) -> None:
    """An sdf of a single titled molecule contains it with the expected title.

    Args:
        single_titled_mol_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock[str] = next(parse_sdf(single_titled_mol_sdf))
    assert mol.title == "Title"


@pytest.mark.xfail  # sendoff does not support SDF titles that start with $$$$
def test_first_delimiter_titled_mol_read(
    double_first_delimiter_titled_mol_sdf: Pathy,
) -> None:
    """SDF with 1st molecule delimiter titled, containing with expected titles.

    Args:
        double_first_delimiter_titled_mol_sdf: pytest fixture of a Path to sdf
    """
    supp = parse_sdf(double_first_delimiter_titled_mol_sdf)
    mol: SDBlock[str] = next(supp)
    other_mol: SDBlock[str] = next(supp)
    assert mol.title == "$$$$"
    assert other_mol.title == "Title"


def test_second_delimiter_titled_mol_read(
    double_second_delimiter_titled_mol_sdf: Pathy,
) -> None:
    """SDF with 2nd molecule delimiter titled, containing with expected titles.

    Args:
        double_second_delimiter_titled_mol_sdf: pytest fixture of a Path to sdf
    """
    supp = parse_sdf(double_second_delimiter_titled_mol_sdf)
    other_mol: SDBlock[str] = next(supp)
    mol: SDBlock[str] = next(supp)
    assert other_mol.title == "Title"
    assert mol.title == "$$$$"


def test_single_record_mol_read(single_record_mol_sdf: Pathy) -> None:
    """An sdf of a delimiter record valued molecule contains the expected data.

    Args:
        single_record_mol_sdf: pytest fixture of a Path to the sdf
    """
    mol = next(parse_sdf(single_record_mol_sdf))
    records = list(mol.records())
    assert ("Record", "Value") in records


def test_iter_records_preserves_order_duplicates_and_multiline_values() -> None:
    """Iterate records on demand without changing list-returning records()."""
    block = next(
        SDBlock.from_lines(
            [
                "Title\n",
                "\n",
                "M  END\n",
                "> <Record>\n",
                "first\n",
                "second\n",
                "\n",
                "> <Record>\n",
                "third\n",
                "\n",
                "$$$$\n",
            ]
        )
    )
    iterator = block.iter_records()
    assert iter(iterator) is iterator
    first = next(iterator)
    second = next(iterator)

    assert first == ("Record", "first\nsecond")
    assert second == ("Record", "third")
    assert first == ("Record", "first\nsecond")
    records = block.records()
    assert isinstance(records, list)
    assert records == [first, second]


def test_iter_records_defers_malformed_field_error() -> None:
    """A later malformed record does not prevent reading earlier fields."""
    block = next(
        SDBlock.from_lines(
            [
                "Title\n",
                "\n",
                "M  END\n",
                "> <Good>\n",
                "value\n",
                "\n",
                "> missing name\n",
                "other\n",
                "\n",
                "$$$$\n",
            ]
        )
    )
    iterator = block.iter_records()

    assert next(iterator) == ("Good", "value")
    with pytest.raises(ValueError, match="missing a name"):
        next(iterator)
    with pytest.raises(StopIteration):
        next(iterator)
    with pytest.raises(ValueError, match="missing a name"):
        block.records()


def test_iter_records_preserves_bytes_and_owns_its_block() -> None:
    """Byte record iteration remains lazy and retains source data safely."""
    block = next(
        SDBlock.from_lines(
            [
                b"Title\n",
                b"\n",
                b"M  END\n",
                b"> <Record>\n",
                b"first\xfe\n",
                b"second\n",
                b"\n",
                b"> <Record>\n",
                b"third\n",
                b"\n",
                b"$$$$\n",
            ]
        )
    )
    iterator = block.iter_records()
    del block

    record = next(iterator)

    assert record == (b"Record", b"first\xfe\nsecond")
    assert list(iterator) == [(b"Record", b"third")]


def test_get_and_to_dict_match_dict_of_records() -> None:
    """get() and to_dict() let later duplicate records win, like dict()."""
    text = "t\nM  END\n> <A>\n1\n\n> <B>\nx\ny\n\n> <A>\n2\n\n$$$$\n"
    block = next(SDBlock.from_lines(io.StringIO(text)))
    byte_block = next(SDBlock.from_lines(io.BytesIO(text.encode())))

    assert block.to_dict() == dict(block.records()) == {"A": "2", "B": "x\ny"}
    assert byte_block.to_dict() == {b"A": b"2", b"B": b"x\ny"}
    assert block.get("A") == "2"
    assert block.get("missing") is None
    assert block.get("missing", "default") == "default"
    assert byte_block.get(b"B") == b"x\ny"
    with pytest.raises(TypeError, match="must be bytes"):
        byte_block.get("A")  # type: ignore[call-overload]
    with pytest.raises(TypeError, match="must be str"):
        block.get(b"A")  # type: ignore[call-overload]


def test_record_mutators_keep_order_and_other_lines() -> None:
    """set_record edits in place; remove_record and replace_records rebuild."""
    text = "t\nM  END\n> <A>  (1)\n1\n\nnote\n\n> <B>\nb\n\n> <A>\n2\n\n$$$$\n"
    block = next(SDBlock.from_lines(io.StringIO(text)))

    block.set_record("A", "new")
    block.set_record("C", "c")
    assert block.records() == [("A", "new"), ("B", "b"), ("C", "c")]
    assert block.remove_record("B") == 1
    assert block.remove_record("B") == 0
    out = io.StringIO()
    block.write(out)
    assert out.getvalue() == ("t\nM  END\n> <A>\nnew\n\nnote\n\n> <C>\nc\n\n$$$$\n")

    block.replace_records({"X": "1", "Y": "2"})
    assert block.records() == [("X", "1"), ("Y", "2")]
    block.replace_records((name, value) for name, value in [("X", "a"), ("X", "b")])
    assert block.records() == [("X", "a"), ("X", "b")]
    block.replace_records([])
    assert block.records() == []
    assert block.molblock == "t\nM  END\n"


def test_record_mutators_check_types_and_headers() -> None:
    """Names and values must match the block's type; bad headers raise."""
    byte_block = next(SDBlock.from_lines([b"t\n", b"M  END\n", b"> <A>\n", b"1\n"]))
    byte_block.set_record(b"A", b"2")
    assert byte_block.to_dict() == {b"A": b"2"}
    with pytest.raises(TypeError, match="record names must be bytes"):
        byte_block.set_record("A", b"2")  # type: ignore[arg-type]
    with pytest.raises(TypeError, match="values must be bytes"):
        byte_block.replace_records({b"A": "2"})  # type: ignore[arg-type]
    with pytest.raises(TypeError, match="record names must be bytes"):
        byte_block.remove_record("A")  # type: ignore[arg-type]
    with pytest.raises(TypeError):
        byte_block.replace_records([b"AB"])  # type: ignore[list-item]

    bad = next(SDBlock.from_lines(["t\n", "M  END\n", "> no name\n", "1\n"]))
    with pytest.raises(ValueError, match="missing a name"):
        bad.set_record("A", "1")
    with pytest.raises(ValueError, match="missing a name"):
        bad.remove_record("A")
    bad.replace_records({"A": "1"})
    assert bad.records() == [("A", "1")]


def test_record_names_are_shared_between_blocks() -> None:
    """Repeated names reuse one object; evicted and long names stay exact."""

    def check(sdf: typing.AnyStr, names: list[typing.AnyStr]) -> None:
        first, second = SDBlock(sdf), SDBlock(sdf)
        assert [n for n, _ in first.records()] == names
        assert list(first.to_dict()) == [n for n, _ in second.records()]
        assert all(n + v[-1:] == v for n, v in second.iter_records())
        assert first.records()[-1][0] is second.records()[-1][0]

    names = [f"N{index}" for index in range(5000)] + ["L" * 65, "A"]
    sdf = "t\nM  END\n" + str.join("", [f"> <{name}>\n{name}v\n\n" for name in names])
    check(sdf, names)
    check(sdf.encode(), [name.encode() for name in names])


def test_record_names_follow_changing_orders() -> None:
    """Reordered, missing, repeated, and new names between blocks read exactly."""
    orders = [
        ["A", "B", "C", "D"],
        ["A", "C", "D"],
        ["D", "C", "B", "A"],
        ["A", "A", "E", "B"],
        [],
        ["B", "C", "D", "A"],
        ["A", "B", "C", "D"],
    ]
    sdf = str.join(
        "",
        [
            "t\nM  END\n"
            + str.join("", [f"> <{n}>\n{n}{i}\n\n" for n in order])
            + "$$$$\n"
            for i, order in enumerate(orders)
        ],
    )
    expected = [[(n, f"{n}{i}") for n in order] for i, order in enumerate(orders)]
    for _ in range(2):
        blocks = list(SDBlock.from_sdf(sdf))
        assert [block.records() for block in blocks] == expected
        assert [block.to_dict() for block in blocks] == [dict(e) for e in expected]
        assert [list(block.iter_records()) for block in blocks] == expected
    lazy = blocks[2].iter_records()
    assert next(lazy) == expected[2][0]
    assert blocks[3].to_dict() == dict(expected[3])
    assert list(lazy) == expected[2][1:]
    byte_blocks = list(SDBlock.from_sdf(sdf.encode()))
    assert [block.records() for block in byte_blocks] == [
        [(n.encode(), v.encode()) for n, v in e] for e in expected
    ]


def test_sdblock_is_generic_at_runtime() -> None:
    """SDBlock[str] and SDBlock[bytes] work in evaluated annotations."""
    alias = SDBlock[str]
    assert typing.get_origin(alias) is SDBlock
    assert typing.get_args(alias) == (str,)
    assert typing.get_args(SDBlock[bytes]) == (bytes,)


def test_mdl_is_the_ctab_text() -> None:
    """The mdl property is the CTAB through M  END as one value, without the title."""
    text = "title\r\n  prog\r\n\r\nM  END\r\n> <A>\r\n1\r\n"
    block = next(SDBlock.from_lines(io.StringIO(text, newline=None)))
    byte_block = next(SDBlock.from_lines(io.BytesIO(text.encode())))

    assert block.mdl == "  prog\n\nM  END\n"
    assert block.title + "\n" + block.mdl == block.molblock
    assert block.mdl.splitlines() == ["  prog", "", "M  END"]
    assert byte_block.mdl == b"  prog\r\n\r\nM  END\r\n"
    assert byte_block.mdl.splitlines() == [b"  prog", b"", b"M  END"]
    assert next(SDBlock.from_lines(["t\n", "M  END\n"])).mdl == "M  END\n"
    assert next(SDBlock.from_lines(["t\n"])).mdl == ""
    no_end = next(SDBlock.from_lines(["t\n", "  prog\n", "> <A>\n", "$$$$\n"]))
    assert no_end.mdl == "  prog\n> <A>\n"


def test_names_len_and_contains_skip_values() -> None:
    """names(), len() and in read record names without building values."""
    sdf = "t\nM  END\n> <A>\n1\n2\n\n> <B>\n\n> <A>\n3\n\n$$$$\n"
    (block,) = SDBlock.from_sdf(sdf)
    (byte_block,) = SDBlock.from_sdf(sdf.encode())

    assert block.names() == [name for name, _ in block.records()] == ["A", "B", "A"]
    assert byte_block.names() == [b"A", b"B", b"A"]
    assert len(block) == len(byte_block) == 3
    assert "A" in block and "B" in block and "C" not in block
    assert b"A" in byte_block and "" not in block
    with pytest.raises(TypeError, match="record names must be bytes"):
        assert "A" in byte_block  # type: ignore[operator]

    empty = SDBlock("t\n\nM  END\n")
    assert (empty.names(), len(empty), "A" in empty) == ([], 0, False)
    assert empty, "a block without records is still true"

    block.append_record("C", "4")
    block.remove_record("A")
    assert (block.names(), len(block), "A" in block) == (["B", "C"], 2, False)
    malformed = SDBlock("t\n\nM  END\n> <A>\n1\n\n> no name\n1\n\n$$$$\n")
    for call in (malformed.names, lambda: len(malformed), lambda: "A" in malformed):
        with pytest.raises(ValueError):
            call()


def test_constructor_parses_one_block() -> None:
    """SDBlock(sdf) parses a molblock or one record and appends records."""
    molblock = "title\n  prog\n\nM  END\n"
    block = SDBlock(molblock)
    assert (block.title, block.molblock, block.records()) == ("title", molblock, [])
    assert block.to_sdf() == molblock + "$$$$\n"

    with_data = SDBlock(molblock + "> <A>\n1\n\n$$$$\n\n", {"B": "2"})
    assert with_data.molblock == molblock
    assert with_data.records() == [("A", "1"), ("B", "2")]
    pairs = SDBlock(molblock.encode(), [(b"A", b"1"), (b"A", b"2")])
    assert pairs.records() == [(b"A", b"1"), (b"A", b"2")]

    assert SDBlock(molblock.replace("\n", "\r\n")).molblock == molblock
    assert SDBlock(molblock.encode().replace(b"\n", b"\r\n")).mdl.endswith(
        b"M  END\r\n"
    )
    empty = SDBlock("")
    assert (empty.title, empty.mdl, empty.records()) == ("", "", [])

    with pytest.raises(ValueError, match="one block"):
        SDBlock(molblock + "$$$$\n" + molblock)
    with pytest.raises(TypeError, match="str or bytes"):
        SDBlock(1)  # type: ignore[type-var]
    with pytest.raises(TypeError, match="record names must be str"):
        SDBlock(molblock, [(b"A", b"1")])  # type: ignore[type-var]


def test_blocks_compare_by_content_and_are_unhashable() -> None:
    """== compares kind and content, and blocks are mutable so unhashable."""
    sdf = "t\n  prog\n\nM  END\n> <A>\n1\n\n$$$$\n"
    block = SDBlock(sdf)
    assert block == SDBlock(sdf) and not block != SDBlock(sdf)
    assert block != SDBlock(sdf.encode())
    assert SDBlock(sdf.encode()) == SDBlock(sdf.encode())
    assert block != sdf and block != object()

    edited = SDBlock(sdf)
    edited.set_record("A", "2")
    assert edited != block
    edited.set_record("A", "1")
    assert edited == block

    # The same text parses differently when SD data follows an MDL
    # without M  END, so equal text alone does not make equal blocks.
    appended = SDBlock("t\n  prog\n")
    appended.append_record("A", "1")
    reparsed = SDBlock(appended.to_sdf())
    assert reparsed.to_sdf() == appended.to_sdf()
    assert (reparsed.records(), appended.records()) == ([], [("A", "1")])
    assert reparsed != appended

    with pytest.raises(TypeError, match="unhashable"):
        hash(block)


def _state_edge_blocks() -> list[SDBlock[str]]:
    """Build blocks that SDBlock(block.to_sdf()) does not rebuild exactly.

    Returns:
        Edited and unusual blocks
    """
    delimiter_title = SDBlock("t\n  prog\n\nM  END\n> <A>\n1\n\n")
    delimiter_title.title = "$$$$ not a delimiter"
    no_end = SDBlock("t\n  prog\n")
    no_end.append_record("A", "1")
    padded = SDBlock(" t \n  prog\n\nM  END\n")
    return [delimiter_title, no_end, padded, SDBlock("")]


@pytest.mark.parametrize("protocol", range(pickle.HIGHEST_PROTOCOL + 1))
def test_blocks_pickle_their_exact_state(protocol: int) -> None:
    """Pickling keeps the kind and the state, even where text would not.

    Args:
        protocol: pickle protocol
    """
    sdf = "t\r\n  prog\r\n\r\nM  END\r\n> <A>\r\n1\r\n\r\n$$$$\r\n"
    blocks: list[SDBlock[str] | SDBlock[bytes]] = [
        SDBlock(sdf),
        SDBlock(sdf.encode()),
        SDBlock("caf\u00e9\nM  END\n"),
        SDBlock(b"\xff\nM  END\n> <\xfe>\n\xfd\n\n"),
        *_state_edge_blocks(),
    ]
    for block in blocks:
        restored = pickle.loads(pickle.dumps(block, protocol))
        assert restored == block and restored is not block
        assert type(restored.title) is type(block.title)
        assert restored.to_sdf() == block.to_sdf()
        assert restored.records() == block.records()
    assert pickle.loads(pickle.dumps(blocks, protocol)) == blocks


def test_block_state_is_validated() -> None:
    """Unpickling rejects state that no block has instead of panicking."""
    from_state = SDBlock._from_state
    assert from_state(b"t\n$$$$\n", 2, 2, True) == SDBlock("t\n")
    for state in [
        (b"", 0, 0, True),
        (b"t\n", 2, 2, True),
        (b"t\n$$$$", 2, 2, False),
        (b"t\n$$$$\n", 1, 1, False),
        (b"t\nx\n$$$$\n", 2, 3, False),
        (b"t\nx\n$$$$\n", 4, 4, False),
        (b"t\nM  END\n$$$$\n", 9, 2, False),
        (b"t\nM  END\n$$$$\n", 2, 99, False),
        (b"\xff\nM  END\n$$$$\n", 2, 9, True),
    ]:
        with pytest.raises(ValueError, match="invalid SDBlock state"):
            from_state(*state)
    assert from_state(b"\xff\nM  END\n$$$$\n", 2, 9, False).title == b"\xff"


@pytest.mark.parametrize("copier", [copy.copy, copy.deepcopy])
def test_block_copies_are_independent(
    copier: typing.Callable[[SDBlock[str]], SDBlock[str]],
) -> None:
    """Copies, shallow or deep, are equal blocks that are edited separately.

    Args:
        copier: copy.copy or copy.deepcopy
    """
    for block in [SDBlock("t\nM  END\n> <A>\n1\n\n"), *_state_edge_blocks()]:
        duplicate = copier(block)
        assert duplicate == block and duplicate is not block
        duplicate.append_record("B", "2")
        duplicate.title = "other"
        assert duplicate != block and "B" not in block
    byte_block = SDBlock(b"t\nM  END\n")
    assert copier(byte_block) == byte_block  # type: ignore[arg-type]
    shared = SDBlock("t\n")
    first, second = copy.deepcopy([shared, shared])
    assert first is second and first is not shared


def test_block_repr_shows_the_title() -> None:
    """The repr shows the stripped title in the block's type."""
    assert repr(SDBlock(" aspirin \nM  END\n")) == "<SDBlock title='aspirin'>"
    assert repr(SDBlock(b"aspirin\nM  END\n")) == "<SDBlock title=b'aspirin'>"


def test_from_sdf_parses_like_a_file(tmp_path: Path) -> None:
    """from_sdf parses in-memory str with universal newlines, bytes as is.

    Args:
        tmp_path: pytest fixture of a temporary directory
    """
    # The title's CRLF straddles the 128 KiB chunk size, and the second
    # block's value spans several chunks.
    title = "x" * ((128 << 10) - 1)
    value = "v" * (1 << 20)
    text = f"{title}\r\nM  END\r> <A>\r\n1\n\n$$$$\r\nt\nM  END\n> <B>\n{value}\n$$$$\n"
    path = tmp_path / "input.sdf"
    path.write_text(text, newline="")
    data = text.encode()
    expected = [b.to_sdf() for b in SDBlock.from_sdf(text)]
    assert [b.to_sdf() for b in parse_sdf(path)] == expected
    with path.open() as handle:
        assert [b.to_sdf() for b in SDBlock.from_lines(handle)] == expected
    expected_bytes = [b.to_sdf() for b in SDBlock.from_sdf(data)]
    assert [b.to_sdf() for b in parse_sdf(os.fsencode(path))] == expected_bytes
    assert [b.to_sdf() for b in SDBlock.from_lines(io.BytesIO(data))] == (
        expected_bytes
    )
    assert [b.records() for b in SDBlock.from_sdf(text)] == [
        [("A", "1")],
        [("B", value)],
    ]
    assert list(SDBlock.from_sdf("")) == []
    with pytest.raises(TypeError, match="use SDBlock.from_sdf"):
        SDBlock.from_lines("text")
    with pytest.raises(TypeError, match="use SDBlock.from_sdf"):
        SDBlock.from_lines(b"bytes")  # type: ignore[arg-type]
    with pytest.raises(TypeError, match="str or bytes"):
        SDBlock.from_sdf(["lines"])  # type: ignore[call-overload]


def test_molblock_setter_replaces_the_structure_only() -> None:
    """Setting molblock swaps title and CTAB and keeps every record line."""
    text = "old\nsrc\nM  END\n> <A>  (1)\n1\n\nnote\n\n$$$$\n"
    block = next(SDBlock.from_lines(io.StringIO(text)))

    block.molblock = "new\r\n  prog\r\n\r\nM  END\r\n\r\n"
    assert block.title == "new"
    assert block.molblock == "new\n  prog\n\nM  END\n"
    out = io.StringIO()
    block.write(out)
    assert out.getvalue() == ("new\n  prog\n\nM  END\n> <A>  (1)\n1\n\nnote\n\n$$$$\n")
    for bad, message in [
        ("t\nx\n", "no M  END"),
        ("t\nM  END\n> <A>\n1\n", "after its M  END"),
        ("t\n$$$$\nM  END\n", r"\$\$\$\$ line"),
    ]:
        with pytest.raises(ValueError, match=message):
            block.molblock = bad
    with pytest.raises(TypeError, match="molblocks must be str"):
        block.molblock = b"t\nM  END\n"  # type: ignore[assignment]
    assert block.molblock == "new\n  prog\n\nM  END\n"

    byte_block = next(SDBlock.from_lines([b"t\n", b"M  END\n", b"> <A>\n", b"1\n"]))
    byte_block.molblock = b"b\r\nM  END\r\n"
    assert byte_block.molblock == b"b\nM  END\r\n"  # the title is stripped
    assert byte_block.title == b"b"
    assert byte_block.to_dict() == {b"A": b"1"}


def test_title_setter_replaces_only_the_title() -> None:
    """Setting the title keeps the CTAB and records and rejects line breaks."""
    block = next(SDBlock.from_lines(["old\n", "src\n", "M  END\n", "> <A>\n", "1\n"]))

    block.title = "a much longer title"
    assert block.title == "a much longer title"
    assert block.molblock == "a much longer title\nsrc\nM  END\n"
    assert block.records() == [("A", "1")]
    block.title = " padded "
    assert block.title == "padded"
    for bad in ["two\nlines", "cr\r"]:
        with pytest.raises(ValueError, match="single line"):
            block.title = bad
    with pytest.raises(TypeError, match="titles must be str"):
        block.title = b"bytes"  # type: ignore[assignment]
    assert block.title == "padded"


def test_text_strips_like_str_and_bytes_like_bytes() -> None:
    """Text blocks strip str.strip() whitespace; bytes blocks bytes.strip().

    A value line holding only Unicode whitespace is blank in text blocks, as
    with str.strip() in pure-Python sendoff 0.1.x, so it ends the record.
    """
    text_block = next(
        SDBlock.from_lines(
            [
                " \u2003Title\u3000 \n",
                "\n",
                "M  END\n",
                "> <Record>\u2003\n",
                "\u2003Text\x1f\n",
                "\n",
                "> <Blank>\n",
                "\u3000\n",
                "skipped\n",
                "\n",
                "$$$$\n",
            ]
        )
    )
    byte_block = next(
        SDBlock.from_lines(
            [
                b"\xc2\xa0Title\xc2\xa0\n",
                b"\n",
                b"M  END\n",
                b"> <Record>\n",
                b"\xc2\xa0Bytes\x1f\n",
                b"\n",
                b"$$$$\n",
            ]
        )
    )

    assert text_block.title == "Title"
    assert text_block.records() == [("Record", "Text"), ("Blank", "")]
    assert byte_block.title == b"\xc2\xa0Title\xc2\xa0"
    assert byte_block.records() == [(b"Record", b"\xc2\xa0Bytes\x1f")]


@pytest.mark.xfail  # sendoff does not support metadata vals starting with $$$$
def test_delimiter_record_mol_read(
    single_delimiter_record_mol_sdf: Pathy,
) -> None:
    """An sdf of a delimiter record valued molecule contains the expected data.

    Args:
        single_delimiter_record_mol_sdf: pytest fixture of a Path to the sdf
    """
    mol = next(parse_sdf(single_delimiter_record_mol_sdf))
    records = list(mol.records())
    assert ("Record", "$$$$") in records


def test_multiline_record_mol_read(
    single_multiline_record_mol_sdf: Pathy,
) -> None:
    """An sdf of a multiline record valued molecule contains the expected data.

    Args:
        single_multiline_record_mol_sdf: pytest fixture of a Path to the sdf
    """
    mol = next(parse_sdf(single_multiline_record_mol_sdf))
    records = list(mol.records())
    assert ("Record", "0\n1\n2\n3") in records


def test_empty_string_record_mol_read(
    single_empty_string_record_mol_sdf: Pathy,
) -> None:
    """An sdf of an empty string record valued molecule contains the expected data.

    Args:
        single_empty_string_record_mol_sdf: pytest fixture of a Path to the sdf
    """
    mol = next(parse_sdf(single_empty_string_record_mol_sdf))
    records = list(mol.records())
    assert ("Record", "") in records


def test_multiline_record_name_mol_read(
    single_multiline_record_name_mol_sdf: Pathy,
) -> None:
    """An sdf of a multiline record name molecule cannot contain the expected data.

    Args:
        single_multiline_record_name_mol_sdf: fixture of a Path to the sdf
    """
    mol = next(parse_sdf(single_multiline_record_name_mol_sdf))
    records = list(mol.records())
    assert "Rec\nord" not in map(str, chain(zip(*records)))


def test_single_mol_newline_write(single_mol_sdf: Pathy) -> None:
    """An sdf block gets written with as many lines as the input.

    Args:
        single_mol_sdf: pytest fixture of a Path to the sdf
    """
    mols = list(parse_sdf(single_mol_sdf))
    buffer = io.StringIO()
    for mol in mols:
        mol.write(buffer)
    assert len(buffer.getvalue().splitlines()) == len(
        open(single_mol_sdf).read().splitlines()
    )


def test_single_mol_newline_write_splitlines(single_mol_sdf: Pathy) -> None:
    """An sdf block gets written with as many lines as a splitlines input.

    Args:
        single_mol_sdf: pytest fixture of a Path to the sdf
    """
    mols = list(SDBlock.from_lines(open(single_mol_sdf).read().splitlines()))
    buffer = io.StringIO()
    for mol in mols:
        mol.write(buffer)
    assert len(buffer.getvalue().splitlines()) == len(
        open(single_mol_sdf).read().splitlines()
    )


def test_right_angle_bracket_record_name_mol_read(
    single_right_angle_bracket_record_name_mol_sdf: Pathy,
) -> None:
    """An sdf with a record name with a ">" is present in the read data as expected.

    Args:
        single_right_angle_bracket_record_name_mol_sdf: fixture of a Path to the sdf
    """
    mol: SDBlock[str] = next(parse_sdf(single_right_angle_bracket_record_name_mol_sdf))
    records = list(mol.records())
    assert records == [("Rec>ord", "Value")]


def test_0_atoms_mol_num_atoms(single_0_atoms_mol_sdf: Pathy) -> None:
    """An sdf block with 0 atoms written parses as 0 atoms.

    Args:
        single_0_atoms_mol_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock[str] = next(parse_sdf(single_0_atoms_mol_sdf))
    assert mol.num_atoms() == 0


def test_star_atom_mol_num_atoms(single_star_atom_mol_sdf: Pathy) -> None:
    """An sdf block with a single star atom written parses as 1 atom.

    Args:
        single_star_atom_mol_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock[str] = next(parse_sdf(single_star_atom_mol_sdf))
    assert mol.num_atoms() == 1


def test_999_atoms_mol_num_atoms(single_999_atoms_mol_sdf: Pathy) -> None:
    """An sdf block with 999 atoms written parses as 999 atoms.

    Args:
        single_999_atoms_mol_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock[str] = next(parse_sdf(single_999_atoms_mol_sdf))
    assert mol.num_atoms() == 999


def test_1001_atoms_mol_num_atoms(single_1001_atoms_mol_sdf: Pathy) -> None:
    """An sdf block with 1001 atoms written parses as 1001 atoms.

    Args:
        single_1001_atoms_mol_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock[str] = next(parse_sdf(single_1001_atoms_mol_sdf))
    assert mol.num_atoms() == 1001


def test_0_atoms_mol_num_bonds(single_0_atoms_mol_sdf: Pathy) -> None:
    """An sdf block with 0 atoms written parses as 0 bonds.

    Args:
        single_0_atoms_mol_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock[str] = next(parse_sdf(single_0_atoms_mol_sdf))
    assert mol.num_bonds() == 0


def test_star_atom_mol_num_bonds(single_star_atom_mol_sdf: Pathy) -> None:
    """An sdf block with a single star atom written parses as 0 bonds.

    Args:
        single_star_atom_mol_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock[str] = next(parse_sdf(single_star_atom_mol_sdf))
    assert mol.num_bonds() == 0


def test_999_atoms_mol_num_bonds(single_999_atoms_mol_sdf: Pathy) -> None:
    """An sdf block with 999 atoms written parses as 998 bonds.

    Args:
        single_999_atoms_mol_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock[str] = next(parse_sdf(single_999_atoms_mol_sdf))
    assert mol.num_bonds() == 998


def test_1001_atoms_mol_num_bonds(single_1001_atoms_mol_sdf: Pathy) -> None:
    """An sdf block with 1001 atoms written parses as 1000 bonds.

    Args:
        single_1001_atoms_mol_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock[str] = next(parse_sdf(single_1001_atoms_mol_sdf))
    assert mol.num_bonds() == 1000


def test_single_mol_v2000_sdf_valid_atom_indices_not_implemented(
    single_mol_v2000_sdf: Pathy,
) -> None:
    """A v2000 sdf block cannot yet be validated for atom indices.

    Args:
        single_mol_v2000_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock[str] = next(parse_sdf(single_mol_v2000_sdf))
    ctable = mol.ctable()
    with pytest.raises(NotImplementedError):
        assert ctable.valid_atom_indices()


def test_single_mol_v2000_sdf_valid_bond_indices_not_implemented(
    single_mol_v2000_sdf: Pathy,
) -> None:
    """A v2000 sdf block cannot yet be validated for bond indices.

    Args:
        single_mol_v2000_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock[str] = next(parse_sdf(single_mol_v2000_sdf))
    ctable = mol.ctable()
    with pytest.raises(NotImplementedError):
        assert ctable.valid_bond_indices()


def test_single_large_atom_index_valid_atom_indices(
    single_large_atom_index_v3000_sdf: Pathy,
) -> None:
    """An sdf block with large atom indices has valid indices.

    That means that the number of atoms matches the counts line,
    and each index is unique and integer.
    Does not check that the indices are 1-indexed and in order.

    Args:
        single_large_atom_index_v3000_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock[str] = next(parse_sdf(single_large_atom_index_v3000_sdf))
    ctable = mol.ctable()
    assert ctable.valid_atom_indices()


def test_single_large_atom_index_strict_invalid_atom_indices(
    single_large_atom_index_v3000_sdf: Pathy,
) -> None:
    """An sdf block with large atom indices has invalid indices, when strict.

    That means that the number of atoms matches the counts line,
    and each index is unique and integer.
    However, the indices are not 1-indexed and in order.

    Args:
        single_large_atom_index_v3000_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock[str] = next(parse_sdf(single_large_atom_index_v3000_sdf))
    ctable = mol.ctable()
    with pytest.raises(IndicesOutOfOrderError, match="atoms"):
        ctable.valid_atom_indices(strict=True)


def test_single_large_bond_index_valid_bond_indices(
    single_large_bond_index_v3000_sdf: Pathy,
) -> None:
    """An sdf block with large bond indices has valid indices.

    That means that the number of bonds matches the counts line,
    and each index is unique and integer.
    Does not check that the indices are 1-indexed and in order.

    Args:
        single_large_bond_index_v3000_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock[str] = next(parse_sdf(single_large_bond_index_v3000_sdf))
    ctable = mol.ctable()
    assert ctable.valid_bond_indices()


def test_single_large_bond_index_strict_invalid_bond_indices(
    single_large_bond_index_v3000_sdf: Pathy,
) -> None:
    """An sdf block with large bond indices has invalid indices, when strict.

    That means that the number of bonds matches the counts line,
    and each index is unique and integer.
    However, the indices are not 1-indexed and in order.

    Args:
        single_large_bond_index_v3000_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock[str] = next(parse_sdf(single_large_bond_index_v3000_sdf))
    ctable = mol.ctable()
    with pytest.raises(IndicesOutOfOrderError, match="bonds"):
        ctable.valid_bond_indices(strict=True)


def test_single_missing_atom_line_invalid_atom_indices(
    single_missing_atom_line_v3000_sdf: Pathy,
) -> None:
    """An sdf block with a missing atom line has invalid indices.

    That means that the number of atoms is fewer than the counts line.

    Args:
        single_missing_atom_line_v3000_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock[str] = next(parse_sdf(single_missing_atom_line_v3000_sdf))
    ctable = mol.ctable()
    with pytest.raises(IndicesMismatchError, match="fewer atom lines than count line"):
        ctable.valid_atom_indices()


def test_single_extra_atom_line_invalid_atom_indices(
    single_extra_atom_line_v3000_sdf: Pathy,
) -> None:
    """An sdf block with a extra atom line has invalid indices.

    That means that the number of atoms is more than the counts line.

    Args:
        single_extra_atom_line_v3000_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock[str] = next(parse_sdf(single_extra_atom_line_v3000_sdf))
    ctable = mol.ctable()
    with pytest.raises(IndicesMismatchError, match="more atom lines than count line"):
        ctable.valid_atom_indices()


def test_single_missing_bond_line_invalid_atom_indices(
    single_missing_bond_line_v3000_sdf: Pathy,
) -> None:
    """An sdf block with a missing bond line has invalid indices.

    That means that the number of bonds is fewer than the counts line.

    Args:
        single_missing_bond_line_v3000_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock[str] = next(parse_sdf(single_missing_bond_line_v3000_sdf))
    ctable = mol.ctable()
    with pytest.raises(IndicesMismatchError, match="fewer bond lines than count line"):
        ctable.valid_bond_indices()


def test_single_extra_bond_line_invalid_atom_indices(
    single_extra_bond_line_v3000_sdf: Pathy,
) -> None:
    """An sdf block with a extra bond line has invalid indices.

    That means that the number of bonds is more than the counts line.

    Args:
        single_extra_bond_line_v3000_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock[str] = next(parse_sdf(single_extra_bond_line_v3000_sdf))
    ctable = mol.ctable()
    with pytest.raises(IndicesMismatchError, match="more bond lines than count line"):
        ctable.valid_bond_indices()


def test_single_duplicate_atom_index_invalid_atom_indices(
    single_duplicate_atom_index_v3000_sdf: Pathy,
) -> None:
    """An sdf block with a duplicate atom index has invalid indices.

    Args:
        single_duplicate_atom_index_v3000_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock[str] = next(parse_sdf(single_duplicate_atom_index_v3000_sdf))
    ctable = mol.ctable()
    with pytest.raises(IndicesDuplicateError, match="atoms"):
        ctable.valid_atom_indices()


def test_single_duplicate_bond_index_invalid_bond_indices(
    single_duplicate_bond_index_v3000_sdf: Pathy,
) -> None:
    """An sdf block with a duplicate bond index has invalid indices.

    Args:
        single_duplicate_bond_index_v3000_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock[str] = next(parse_sdf(single_duplicate_bond_index_v3000_sdf))
    ctable = mol.ctable()
    with pytest.raises(IndicesDuplicateError, match="bonds"):
        ctable.valid_bond_indices()


def test_single_shuffled_atom_lines_valid_atom_indices(
    single_shuffled_atom_lines_v3000_sdf: Pathy,
) -> None:
    """An sdf block with shuffled atom lines has valid indices.

    That means that the number of atoms matches the counts line,
    and each index is unique and integer.
    Does not check that the indices are 1-indexed and in order.

    Args:
        single_shuffled_atom_lines_v3000_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock[str] = next(parse_sdf(single_shuffled_atom_lines_v3000_sdf))
    ctable = mol.ctable()
    assert ctable.valid_atom_indices()


def test_single_shuffled_atom_lines_strict_invalid_atom_indices(
    single_shuffled_atom_lines_v3000_sdf: Pathy,
) -> None:
    """An sdf block with shuffled atom lines has invalid indices, when strict.

    That means that the number of atoms matches the counts line,
    and each index is unique and integer.
    However, the indices are not 1-indexed and in order.

    Args:
        single_shuffled_atom_lines_v3000_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock[str] = next(parse_sdf(single_shuffled_atom_lines_v3000_sdf))
    ctable = mol.ctable()
    with pytest.raises(IndicesOutOfOrderError, match="atoms"):
        ctable.valid_atom_indices(strict=True)


def test_single_shuffled_bond_lines_valid_bond_indices(
    single_shuffled_bond_lines_v3000_sdf: Pathy,
) -> None:
    """An sdf block with shuffled bond lines has valid indices.

    That means that the number of bonds matches the counts line,
    and each index is unique and integer.
    Does not check that the indices are 1-indexed and in order.

    Args:
        single_shuffled_bond_lines_v3000_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock[str] = next(parse_sdf(single_shuffled_bond_lines_v3000_sdf))
    ctable = mol.ctable()
    assert ctable.valid_bond_indices()


def test_single_shuffled_bond_lines_strict_invalid_bond_indices(
    single_shuffled_bond_lines_v3000_sdf: Pathy,
) -> None:
    """An sdf block with shuffled bond lines has invalid indices, when strict.

    That means that the number of bonds matches the counts line,
    and each index is unique and integer.
    However, the indices are not 1-indexed and in order.

    Args:
        single_shuffled_bond_lines_v3000_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock[str] = next(parse_sdf(single_shuffled_bond_lines_v3000_sdf))
    ctable = mol.ctable()
    with pytest.raises(IndicesOutOfOrderError, match="bonds"):
        ctable.valid_bond_indices(strict=True)


def test_single_mol_v2000_sdf_renumber_indices_not_implemented(
    single_mol_v2000_sdf: Pathy,
) -> None:
    """A v2000 sdf block cannot yet have indices renumbered.

    Args:
        single_mol_v2000_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock[str] = next(parse_sdf(single_mol_v2000_sdf))
    with pytest.raises(NotImplementedError):
        mol.renumber_indices()


def test_single_large_atom_index_strict_valid_atom_indices_renumbered(
    single_large_atom_index_v3000_sdf: Pathy,
) -> None:
    """An sdf block with large atom indices has strict valid indices, after renumbering.

    That means that the number of atoms matches the counts line,
    and each index is unique and integer.
    The indices are 1-indexed and in order after renumbering.

    Args:
        single_large_atom_index_v3000_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock[str] = next(parse_sdf(single_large_atom_index_v3000_sdf))
    mol.renumber_indices()
    ctable = mol.ctable()
    assert ctable.valid_atom_indices(strict=True)


def test_single_large_bond_index_strict_valid_bond_indices_renumbered(
    single_large_bond_index_v3000_sdf: Pathy,
) -> None:
    """An sdf block with large bond indices has strict valid indices, after renumbering.

    That means that the number of atoms matches the counts line,
    and each index is unique and integer.
    The indices are 1-indexed and in order after renumbering.

    Args:
        single_large_bond_index_v3000_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock[str] = next(parse_sdf(single_large_bond_index_v3000_sdf))
    mol.renumber_indices()
    ctable = mol.ctable()
    assert ctable.valid_bond_indices(strict=True)


def test_single_missing_atom_line_valid_atom_indices_after_renumbering(
    single_missing_atom_line_v3000_sdf: Pathy,
) -> None:
    """An sdf block with a missing atom line has stric valid indices, after renumbering.

    Args:
        single_missing_atom_line_v3000_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock[str] = next(parse_sdf(single_missing_atom_line_v3000_sdf))
    mol.renumber_indices()
    ctable = mol.ctable()
    assert ctable.valid_atom_indices(strict=True)


def test_single_extra_atom_line_valid_atom_indices_after_renumbering(
    single_extra_atom_line_v3000_sdf: Pathy,
) -> None:
    """An sdf block with a extra atom line has strict valid indices, after renumbering.

    Args:
        single_extra_atom_line_v3000_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock[str] = next(parse_sdf(single_extra_atom_line_v3000_sdf))
    mol.renumber_indices()
    ctable = mol.ctable()
    assert ctable.valid_atom_indices(strict=True)


def test_single_missing_bond_line_valid_atom_indices_after_renumbering(
    single_missing_bond_line_v3000_sdf: Pathy,
) -> None:
    """An sdf block with a missing bond line has strict valid indices after renumbering.

    Args:
        single_missing_bond_line_v3000_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock[str] = next(parse_sdf(single_missing_bond_line_v3000_sdf))
    mol.renumber_indices()
    ctable = mol.ctable()
    assert ctable.valid_bond_indices(strict=True)


def test_single_extra_bond_line_valid_atom_indices_after_renumbering(
    single_extra_bond_line_v3000_sdf: Pathy,
) -> None:
    """An sdf block with a extra bond line has strict valid indices, after renumbering.

    Args:
        single_extra_bond_line_v3000_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock[str] = next(parse_sdf(single_extra_bond_line_v3000_sdf))
    mol.renumber_indices()
    ctable = mol.ctable()
    assert ctable.valid_bond_indices(strict=True)


def test_single_duplicate_atom_index_invalid_atom_indices_even_renumbering(
    single_duplicate_atom_index_v3000_sdf: Pathy,
) -> None:
    """An sdf block with a duplicate atom index in bonds is invalid in renumbering.

    Args:
        single_duplicate_atom_index_v3000_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock[str] = next(parse_sdf(single_duplicate_atom_index_v3000_sdf))
    with pytest.raises(IndicesDuplicateError, match="atom index mapping in bond"):
        mol.renumber_indices()


def test_single_duplicate_bond_index_valid_bond_indices_after_renumbering(
    single_duplicate_bond_index_v3000_sdf: Pathy,
) -> None:
    """An sdf block w/ duplicate bond index has strict valid indices after renumbering.

    Args:
        single_duplicate_bond_index_v3000_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock[str] = next(parse_sdf(single_duplicate_bond_index_v3000_sdf))
    mol.renumber_indices()
    ctable = mol.ctable()
    assert ctable.valid_bond_indices(strict=True)


def test_single_shuffled_atom_lines_strict_valid_atom_indices_after_renumbering(
    single_shuffled_atom_lines_v3000_sdf: Pathy,
) -> None:
    """An sdf block w/ shuffled atom lines has strict valid indices, after renumbering.

    That means that the number of atoms matches the counts line,
    and each index is unique and integer.
    After renumbering, the indices are 1-indexed and in order.

    Args:
        single_shuffled_atom_lines_v3000_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock[str] = next(parse_sdf(single_shuffled_atom_lines_v3000_sdf))
    mol.renumber_indices()
    ctable = mol.ctable()
    assert ctable.valid_atom_indices(strict=True)


def test_single_shuffled_bond_lines_strict_valid_bond_indices_after_renumbering(
    single_shuffled_bond_lines_v3000_sdf: Pathy,
) -> None:
    """An sdf block w/ shuffled bond lines has strict valid indices, after renumbering.

    That means that the number of bonds matches the counts line,
    and each index is unique and integer.
    After renumbering, the indices are 1-indexed and in order.

    Args:
        single_shuffled_bond_lines_v3000_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock[str] = next(parse_sdf(single_shuffled_bond_lines_v3000_sdf))
    mol.renumber_indices()
    ctable = mol.ctable()
    assert ctable.valid_bond_indices(strict=True)


def test_single_mol_v3000_sdf_idempotent_renumber(
    single_mol_v3000_sdf: Pathy,
) -> None:
    """A V3000 sdf block generated by RDKit, renumbered, is the same as the input.

    Args:
        single_mol_v3000_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock[str] = next(parse_sdf(single_mol_v3000_sdf))
    mol.renumber_indices()
    buffer = io.StringIO()
    mol.write(buffer)
    assert buffer.getvalue() == open(single_mol_v3000_sdf).read()


def test_single_mol_v3000_sdf_idempotent_renumber_with_splitlines(
    single_mol_v3000_sdf: Pathy,
) -> None:
    """A V3000 sdf block generated by RDKit, renumbered, is the same as the input.

    Args:
        single_mol_v3000_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock[str] = next(
        SDBlock.from_lines(open(single_mol_v3000_sdf).read().splitlines())
    )
    mol.renumber_indices()
    buffer = io.StringIO()
    mol.write(buffer)
    assert buffer.getvalue() == open(single_mol_v3000_sdf).read()


V3000_HEADER = "t\n  sendoff\n\n  0  0  0     0  0            999 V3000\n"


def v3000_block(
    atom: str = "M  V30 1 C 0 0 0 0",
    bond: str = "M  V30 1 1 1 2",
    counts: str = "M  V30 COUNTS 2 1 0 0 0",
) -> SDBlock[str]:
    """Build a two-atom V3000 block with one replaceable line of each kind.

    Args:
        atom: the first atom line
        bond: the only bond line
        counts: the V3000 counts line

    Returns:
        The parsed block.
    """
    text = (
        V3000_HEADER
        + f"M  V30 BEGIN CTAB\n{counts}\nM  V30 BEGIN ATOM\n{atom}\n"
        + "M  V30 2 O 0 0 0 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\n"
        + f"{bond}\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n$$$$\n"
    )
    return next(SDBlock.from_lines(text.splitlines()))


@pytest.mark.parametrize(
    "text, error",
    [
        ("t\nM  END\n$$$$\n", ValueError),
        (V3000_HEADER + "M  END\n$$$$\n", ValueError),
        ("t\n\n\n   \nM  END\n$$$$\n", IndexError),
        ("t\n\n\n  0  0  0     0  0            999 V9\nM  END\n$$$$\n", KeyError),
        ("t\n\n\n  x  0  0     0  0            999 V2000\nM  END\n$$$$\n", ValueError),
    ],
)
def test_malformed_ctable_header_raises_exception(
    text: str, error: type[Exception]
) -> None:
    """Malformed CTable headers raise ordinary exceptions, not Rust panics.

    Args:
        text: SDF text with a malformed CTable
        error: the expected exception type
    """
    block = next(SDBlock.from_lines(text.splitlines()))
    with pytest.raises(error):
        block.num_atoms()


def test_malformed_v3000_counts_raises_index_error() -> None:
    """A V3000 counts line without atom and bond fields raises IndexError."""
    with pytest.raises(IndexError):
        v3000_block(counts="M  V30 COUNTS 1").num_atoms()


@pytest.mark.parametrize(
    "atom, bond, error",
    [
        ("M  V30", "M  V30 1 1 1 2", IndexError),
        ("M  V30 x C 0 0 0 0", "M  V30 1 1 1 2", ValueError),
        ("M  V30 1 C 0 0 0 0", "M  V30", IndexError),
        ("M  V30 1 C 0 0 0 0", "M  V30 x 1 1 2", ValueError),
    ],
)
def test_malformed_index_lines_raise_exception(
    atom: str, bond: str, error: type[Exception]
) -> None:
    """Atom and bond lines without integer indices raise ordinary exceptions.

    Args:
        atom: the first atom line
        bond: the only bond line
        error: the expected exception type
    """
    ctable = v3000_block(atom=atom, bond=bond).ctable()
    with pytest.raises(error):
        ctable.valid_atom_indices()
        ctable.valid_bond_indices()


@pytest.mark.parametrize(
    "atom, bond, error",
    [
        ("M  V30", "M  V30 1 1 1 2", IndexError),
        ("M  V30 x C 0 0 0 0", "M  V30 1 1 1 2", ValueError),
        ("M  V30 1 C 0 0 0 0", "M  V30 1 1 1", IndexError),
        ("M  V30 1 C 0 0 0 0", "M  V30 1 1 1 x", ValueError),
        ("M  V30 1 C 0 0 0 0", "M  V30 1 1 1 9", IndexError),
        # A duplicated atom index outranks a missing one in the same bond.
        ("M  V30 2 C 0 0 0 0", "M  V30 1 1 9 2", IndicesDuplicateError),
        ("M  V30 2 C 0 0 0 0", "M  V30 1 1 2 9", IndicesDuplicateError),
    ],
)
def test_malformed_renumber_input_raises_exception(
    atom: str, bond: str, error: type[Exception]
) -> None:
    """Renumbering malformed or dangling lines raises and leaves the block intact.

    Args:
        atom: the first atom line
        bond: the only bond line
        error: the expected exception type
    """
    block = v3000_block(atom=atom, bond=bond)
    before = block.ctable().lines
    with pytest.raises(error):
        block.renumber_indices()
    assert block.ctable().lines == before


def test_renumber_slices_atom_lines_by_character() -> None:
    """Renumbering skips characters like Python slicing, even mid-UTF-8."""
    block = v3000_block(atom="a b 1 \u4e2d")
    block.renumber_indices()
    assert "M  V30 1" in block.ctable().lines


def test_ctable_lines_are_a_fresh_list() -> None:
    """Each CTable.lines access is a new list that does not alias the table."""
    table = v3000_block(atom="a b 1 C").ctable()
    lines = table.lines
    assert lines == table.lines and lines is not table.lines
    lines.clear()
    assert table.lines and table.lines[0] == "t"


def test_write_sdf_writes_what_block_write_would(tmp_path: Path) -> None:
    """write_sdf to a path or handle matches SDBlock.write, across chunks.

    Args:
        tmp_path: pytest fixture of a temporary directory
    """
    # Padded titles are stripped; 3 MiB of blocks span several chunks, and
    # non-ASCII text blocks go to a text handle on their own.
    text = str.join(
        "",
        [
            f"  t{i}  \r\nM  END\n> <A>\n{'é' * (i % 7 < 2)}{'v' * 1000}\n\n$$$$\n"
            for i in range(3000)
        ],
    )
    source = tmp_path / "input.sdf"
    source.write_text(text, newline="")
    out = tmp_path / "out.sdf"
    paths: list[Pathy] = [source, os.fsencode(source)]
    for path in paths:
        blocks = list(parse_sdf(path))
        handle = io.StringIO() if isinstance(path, Path) else io.BytesIO()
        expected = bytes.join(
            b"",
            [
                b.to_sdf().encode() if isinstance(b.title, str) else b.to_sdf()
                for b in blocks
            ],
        )
        assert write_sdf(out, blocks) == 3000
        assert out.read_bytes() == expected
        assert write_sdf(handle, iter(blocks)) == 3000
        written = handle.getvalue()
        assert (written.encode() if isinstance(written, str) else written) == expected
    assert write_sdf(str(out), []) == 0
    assert out.read_bytes() == b""


def test_write_sdf_rejects_bad_input(tmp_path: Path) -> None:
    """write_sdf reports bad paths, mixed handles, and undecodable text.

    Args:
        tmp_path: pytest fixture of a temporary directory
    """
    text, data = SDBlock("t\nM  END\n"), SDBlock(b"t\nM  END\n")
    with pytest.raises(TypeError, match="mixed str and bytes"):
        write_sdf(io.StringIO(), [text, data])  # type: ignore[list-item]
    assert write_sdf(tmp_path / "mixed.sdf", [text, data]) == 2
    with pytest.raises(TypeError):
        write_sdf(tmp_path / "x.sdf", ["t\nM  END\n"])  # type: ignore[list-item]
    with pytest.raises(FileNotFoundError) as error:
        write_sdf(tmp_path / "missing" / "x.sdf", [text])
    assert error.value.filename == str(tmp_path / "missing" / "x.sdf")

    source = tmp_path / "latin1.sdf"
    source.write_bytes(b"t\nM  END\n> <A>\n\xe9\n\n$$$$\n")
    dests: list[typing.Any] = [tmp_path / "bad.sdf", io.StringIO()]
    for dest in dests:
        with pytest.raises(UnicodeDecodeError):
            write_sdf(dest, parse_sdf(source))
