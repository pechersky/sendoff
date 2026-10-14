# mypy: allow-untyped-decorators
"""Test reading SDFs into blocks: parse_sdf, from_lines, from_sdf, decoding."""
import copy
import errno
import gc
import io
import os
import pickle
import signal
import subprocess
import sys
import threading
import typing
import weakref
from importlib import import_module
from itertools import chain
from pathlib import Path
from typing import Generator, Iterator

import pytest

from sendoff.sdblock import SDBlock, parse_sdf, write_sdf
from tests.fixtures import Pathy, sdf_of


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
    sdf.write_bytes(
        b"Closed\n\nM  END\n$$$$\n"
        b"Title\n\nM  END\n> <Record>\nValue" + trailing_newline.encode()
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
    assert b"\xff" in bytes(block)
    assert b"field\xfe" in bytes(block)
    assert b"value\xfd" in bytes(block)
    reparsed = next(SDBlock.from_sdf(bytes(block)))
    assert reparsed.title == block.title
    assert reparsed.records() == block.records()


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
    path.write_bytes("title é\n\nM  END\n> <field>\nvalue é\n\n$$$$\n".encode())

    block = next(parse_sdf(path))

    assert block.title == "title é"
    assert block.records() == [("field", "value é")]


class _Reads:
    """A file object whose `read` returns each of `chunks`, then ""."""

    def __init__(self, *chunks: object) -> None:
        self.chunks = iter(chunks)

    def read(self, size: int) -> object:
        return next(self.chunks, "")


@pytest.mark.parametrize(
    ("lines", "message"),
    [
        ([1], "SDBlock.from_lines lines must be str or bytes"),
        (["t\n", None], "SDBlock.from_lines lines must be str or bytes"),
        (
            [b"t\n", "M  END\n"],
            "SDBlock.from_lines does not accept mixed str and bytes lines",
        ),
        (_Reads(["t\n"]), "file objects must read str or bytes"),
        (_Reads("t\n", b"M  END\n"), "file objects must not mix str and bytes reads"),
    ],
    ids=["item", "later-item", "mixed-items", "read", "mixed-reads"],
)
def test_from_lines_rejects_items_as_blocks_are_produced(
    lines: typing.Any, message: str
) -> None:
    """Items or reads that are not homogeneous str or bytes raise on advance.

    Args:
        lines: the bad input
        message: the TypeError's message
    """
    blocks = SDBlock.from_lines(lines)
    with pytest.raises(TypeError, match=f"^{message}$"):
        next(blocks)


@pytest.mark.parametrize(
    ("lines", "message"),
    [
        (
            "t\n",
            "SDBlock.from_lines takes lines or a file object; use SDBlock.from_sdf",
        ),
        (
            b"t\n",
            "SDBlock.from_lines takes lines or a file object; use SDBlock.from_sdf",
        ),
        (5, "'int' object is not iterable"),
    ],
    ids=["str", "bytes", "non-iterable"],
)
def test_from_lines_rejects_arguments(lines: typing.Any, message: str) -> None:
    """A str, bytes, or non-iterable argument raises at the call.

    Args:
        lines: the bad argument
        message: the start of the TypeError's message
    """
    with pytest.raises(TypeError, match=f"^{message}"):
        SDBlock.from_lines(lines)


def test_from_lines_items_end_lines() -> None:
    """Each item ends a line, completed like its block's title line."""
    cases: dict[tuple[str, ...], list[str]] = {
        ("t", "M  END", "$$$$"): ["t\nM  END\n$$$$\n"],
        ("t\r\n", "M  END", "> <A>", "1", "", "$$$$"): [
            "t\r\nM  END\r\n> <A>\r\n1\r\n\r\n$$$$\r\n"
        ],
        ("t\r", "M  END"): ["t\r\nM  END\r\n$$$$\r\n"],
        ("t\n", "", "M  END\n"): ["t\n\nM  END\n$$$$\n"],
        ("t\nM  END\n$$$$\nu\nM  END\n$$$$\n",): [
            "t\nM  END\n$$$$\n",
            "u\nM  END\n$$$$\n",
        ],
        ("t\n", "M  END\n$$$$\nu\n", "M  END"): [
            "t\nM  END\n$$$$\n",
            "u\nM  END\n$$$$\n",
        ],
    }
    for lines, expected in cases.items():
        assert [sdf_of(block) for block in SDBlock.from_lines(lines)] == expected
        encoded = [line.encode() for line in lines]
        assert [sdf_of(block) for block in SDBlock.from_lines(encoded)] == [
            text.encode() for text in expected
        ]


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


class SelfParsingReader(io.StringIO):
    """A reader that keeps its own parser, making a reference cycle."""

    parser: typing.Any = None


class SelfParsingLines:
    """An iterable of lines that keeps its own parser, making a reference cycle."""

    parser: typing.Any = None

    def __iter__(self) -> Iterator[str]:
        """Yield one block's lines.

        Yields:
            The lines.
        """
        yield from ["t\n", "M  END\n", "$$$$\n"]


@pytest.mark.parametrize("source", [SelfParsingReader, SelfParsingLines])
def test_parsers_in_reference_cycles_are_collected(
    source: type[SelfParsingReader] | type[SelfParsingLines],
) -> None:
    """The garbage collector frees a source that refers to its own parser.

    Args:
        source: the class of source to parse
    """
    held = source()
    held.parser = SDBlock.from_lines(held)
    collected = weakref.ref(held)
    del held
    gc.collect()
    assert collected() is None


def open_paths() -> set[str]:
    """Return the paths this process has open or mapped, from /proc/self.

    Returns:
        Resolved targets of the open file descriptors, and mapped files.
    """
    paths = mapped_paths()
    for fd in os.listdir("/proc/self/fd"):
        try:
            paths.add(os.readlink(os.path.join("/proc/self/fd", fd)))
        except OSError:
            continue
    return paths


def mapped_paths() -> set[str]:
    """Return the files this process has mapped, from /proc/self/maps.

    Returns:
        Resolved paths of the mapped files.
    """
    with open("/proc/self/maps", encoding="utf-8") as maps:
        fields = (line.split(maxsplit=5) for line in maps)
        return {line[5].rstrip("\n") for line in fields if len(line) == 6}


@pytest.mark.skipif(not os.path.isdir("/proc/self/fd"), reason="needs /proc")
def test_parse_sdf_opens_lazily_and_close_releases_the_file(tmp_path: Path) -> None:
    """The file opens on the first next() and close() releases it.

    Args:
        tmp_path: Temporary directory for the SDF fixture.
    """
    path = tmp_path / "two.sdf"
    path.write_bytes(b"First\n\nM  END\n$$$$\nSecond\n\nM  END\n$$$$\n")

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


TWO = b"First\n\nM  END\n> <A>\n1\n\n$$$$\nSecond\n\nM  END\n> <A>\n2\n\n$$$$\n"


def borrowed_blocks(path: Path, binary: bool = False) -> list[SDBlock[typing.Any]]:
    """Parse ``path`` with ``borrow=True``, skipping unless the file maps.

    Args:
        path: an SDF on the temporary directory's filesystem
        binary: parse bytes blocks from a bytes path

    Returns:
        The blocks, which borrow the mapped file.
    """
    blocks = list(parse_sdf(os.fsencode(path) if binary else path, borrow=True))
    if not os.path.exists("/proc/self/maps") or str(path) not in mapped_paths():
        pytest.skip(f"{path.parent} is not mapped")
    return blocks


@pytest.mark.skipif(not os.path.exists("/proc/self/maps"), reason="needs /proc")
def test_borrowed_blocks_keep_the_file_mapped(tmp_path: Path) -> None:
    """Borrowed blocks hold the mapping until the last of them is freed.

    Args:
        tmp_path: Temporary directory for the SDF.
    """
    path = tmp_path / "two.sdf"
    path.write_bytes(TWO)
    copied = list(parse_sdf(path))
    assert str(path) not in mapped_paths()
    borrowed = borrowed_blocks(path)
    assert borrowed == copied
    del borrowed[0]
    gc.collect()
    assert str(path) in mapped_paths()
    del borrowed
    gc.collect()
    assert str(path) not in mapped_paths()


def test_editing_borrowed_blocks_copies_them(tmp_path: Path) -> None:
    """Edits and copies of borrowed blocks leave the file and others alone.

    Args:
        tmp_path: Temporary directory for the SDF.
    """
    path = tmp_path / "two.sdf"
    path.write_bytes(TWO)
    first, second = borrowed_blocks(path, binary=True)
    split = TWO.index(b"Second")
    view = memoryview(first)
    assert view.tobytes() == TWO[:split]
    with pytest.raises(BufferError, match="exports"):
        first.title = b"Edited"
    view.release()
    duplicate = copy.copy(second)
    first.title = b"Edited"
    first[b"B"] = b"x"
    duplicate[b"A"] = b"3"
    assert bytes(first).startswith(b"Edited\n\nM  END\n> <A>\n1\n\n> <B>\nx\n")
    assert second[b"A"] == b"2" and duplicate[b"A"] == b"3"
    assert pickle.loads(pickle.dumps(second)) == second
    assert bytes(second) == TWO[split:]
    assert path.read_bytes() == TWO


def test_borrowed_blocks_outlive_replacing_the_file(tmp_path: Path) -> None:
    """Blocks keep reading the file they mapped after a new one replaces it.

    Args:
        tmp_path: Temporary directory for the SDF.
    """
    path = tmp_path / "two.sdf"
    path.write_bytes(TWO)
    blocks = borrowed_blocks(path)
    write_sdf(tmp_path / "new.sdf", blocks[1:])
    os.replace(tmp_path / "new.sdf", path)
    assert [block.title for block in blocks] == ["First", "Second"]
    assert [block.title for block in parse_sdf(path, borrow=True)] == ["Second"]


def test_borrow_is_ignored_where_nothing_maps(tmp_path: Path) -> None:
    """File objects and lone-CR text are copied, and parse as without it.

    Args:
        tmp_path: Temporary directory for the SDF.
    """
    path = tmp_path / "mixed.sdf"
    path.write_bytes(TWO.replace(b"\n", b"\r", 3))
    with open(path, "rb") as file:
        assert list(parse_sdf(file, borrow=True)) == list(parse_sdf(os.fsencode(path)))
    assert list(parse_sdf(path, borrow=True)) == list(parse_sdf(path))


@pytest.mark.skipif(sys.platform != "linux", reason="only Linux maps files")
@pytest.mark.parametrize("borrow", [False, True])
def test_truncating_a_borrowed_file_raises_sigbus(tmp_path: Path, borrow: bool) -> None:
    """Truncating a file kills a process reading blocks that borrow it.

    A subprocess reads the blocks, truncates the file, and serializes them.
    It loads the extension these tests use.

    Args:
        tmp_path: Temporary directory for the SDF.
        borrow: whether the blocks borrow the mapped file
    """
    native = import_module("sendoff.native").__file__
    assert native is not None
    path = tmp_path / "big.sdf"
    path.write_bytes(TWO.replace(b"1", b"1" * 100_000))
    script = """if True:
        import importlib.util, sys

        spec = importlib.util.spec_from_file_location("sendoff.native", sys.argv[2])
        native = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(native)

        blocks = list(native.parse_sdf(sys.argv[1], borrow=sys.argv[3] == "1"))
        with open("/proc/self/maps", encoding="utf-8") as maps:
            print(sys.argv[1] in maps.read(), flush=True)
        open(sys.argv[1], "wb").close()
        print(sum(len(str(block)) for block in blocks))
    """
    result = subprocess.run(
        [sys.executable, "-c", script, str(path), native, str(int(borrow))],
        capture_output=True,
        encoding="utf-8",
        timeout=60,
    )
    if not borrow:
        assert (result.returncode, result.stdout) == (
            0,
            f"False\n{len(TWO) + 99_999}\n",
        )
    elif result.stdout != "True\n":
        pytest.skip(f"{tmp_path} is not mapped")
    else:
        assert result.returncode == -signal.SIGBUS


@pytest.mark.skipif(not hasattr(os, "mkfifo"), reason="needs named pipes")
def test_parse_sdf_waits_for_a_pipe_writer_without_the_gil(tmp_path: Path) -> None:
    """Opening a pipe releases the GIL, so a Python thread can open its end.

    A subprocess runs the reader so that a deadlock fails the timeout
    rather than hanging the tests. It loads the extension these tests use.

    Args:
        tmp_path: Temporary directory for the pipe.
    """
    native = import_module("sendoff.native").__file__
    assert native is not None
    pipe = tmp_path / "pipe.sdf"
    os.mkfifo(pipe)
    script = """if True:
        import importlib.util, sys, threading, time

        spec = importlib.util.spec_from_file_location("sendoff.native", sys.argv[2])
        native = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(native)

        def write():
            time.sleep(0.2)
            with open(sys.argv[1], "w") as pipe:
                pipe.write("t\\n\\nM  END\\n$$$$\\n")

        threading.Thread(target=write, daemon=True).start()
        print([block.title for block in native.parse_sdf(sys.argv[1])])
    """
    result = subprocess.run(
        [sys.executable, "-c", script, str(pipe), native],
        capture_output=True,
        check=True,
        encoding="utf-8",
        timeout=60,
    )
    assert result.stdout == "['t']\n"


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
    with pytest.raises(OSError) as expected:
        open(tmp_path, "rb")
    with pytest.raises(type(expected.value)) as raised:
        next(parse_sdf(os.fsencode(tmp_path)))
    assert raised.value.errno == expected.value.errno
    with pytest.raises(TypeError, match="path or a file object, not int"):
        parse_sdf(1)  # type: ignore[call-overload]


def test_parse_sdf_reads_file_objects_from_their_position(tmp_path: Path) -> None:
    """File objects give their own type's blocks, as from_lines reads them.

    Args:
        tmp_path: Temporary directory for the SDF fixture.
    """
    sdf = "skip\n\nM  END\n$$$$\nt\r\n\r\nM  END\r\n> <A>\r\n1\r\n\r\n$$$$\r\n"
    path = tmp_path / "in.sdf"
    path.write_bytes(sdf.encode())

    with path.open("rb") as handle:
        handle.readline()
        handle.read(len("\nM  END\n$$$$\n"))
        blocks = list(parse_sdf(handle))
        assert not handle.closed
    assert [sdf_of(block) for block in blocks] == [sdf.split("$$$$\n", 1)[1].encode()]
    assert blocks == list(SDBlock.from_lines(io.BytesIO(sdf.encode())))[1:]

    with path.open(encoding="utf-8") as handle:
        assert [block.to_dict() for block in parse_sdf(handle)] == [{}, {"A": "1"}]
    assert [block.title for block in parse_sdf(io.StringIO(sdf))] == ["skip", "t"]

    read, write = os.pipe()
    with os.fdopen(write, "wb") as writer:
        writer.write(sdf.encode())
    with os.fdopen(read, "rb") as reader:
        assert [block[b"A"] for block in parse_sdf(reader) if b"A" in block] == [b"1"]


# A str-path block with valid title but undecodable CTAB, value, and name,
# and one with an undecodable title.
LATIN1_BODY = (
    b"first\n\xa6 src\n\n  0  0  0  0  0  0  0  0  0  0999 V2000\nM  END\n"
    b"> <F>\n\xa6 bad\n\n> <\xa6>\nv\n\n$$$$\n"
)


LATIN1_TITLE = b"\xa6\n\n\n  0  0  0  0  0  0  0  0  0  0999 V2000\nM  END\n$$$$\n"


DECODE_HINT = (
    "parse_sdf(os.fsencode(path)) reads such files as bytes, and "
    "parse_sdf(open(path, encoding=...)) decodes another encoding."
)


def test_parse_sdf_text_paths_decode_utf8_on_access(tmp_path: Path) -> None:
    """A non-UTF-8 feed still splits on a str path; only bad text raises.

    Args:
        tmp_path: Temporary directory for the SDF fixture.
    """
    path = tmp_path / "latin1.sdf"
    path.write_bytes(LATIN1_BODY + b"second\n\nM  END\n$$$$\n")

    blocks = list(parse_sdf(path))

    assert [block.title for block in blocks] == ["first", "second"]
    assert blocks[1].records() == []
    assert next(parse_sdf(os.fsencode(path))).records()[0] == (b"F", b"\xa6 bad")


# What each SDBlock[str] accessor that decodes block bytes was reading, as
# its UnicodeDecodeError's note says, and the block it reads.
DECODING: dict[str, tuple[bytes, typing.Callable[[SDBlock[str]], object], str]] = {
    "title": (LATIN1_TITLE, lambda block: block.title, "Reading the title of"),
    "molblock": (LATIN1_BODY, lambda block: block.molblock, "Reading the molblock of"),
    "mdl": (LATIN1_BODY, lambda block: block.mdl, "Reading the MDL of"),
    "ctable": (LATIN1_BODY, lambda block: block.ctable(), "Reading the CTAB of"),
    "fix_counts": (
        LATIN1_BODY,
        lambda block: block.fix_counts(),
        "Reading the CTAB of",
    ),
    "renumber_indices": (
        LATIN1_BODY,
        lambda block: block.renumber_indices(),
        "Reading the CTAB of",
    ),
    "records": (LATIN1_BODY, lambda block: block.records(), "Reading record 'F' of"),
    "iter_records": (
        LATIN1_BODY,
        lambda block: list(block.iter_records()),
        "Reading record 'F' of",
    ),
    "get": (LATIN1_BODY, lambda block: block.get("F"), "Reading record 'F' of"),
    "getitem": (LATIN1_BODY, lambda block: block["F"], "Reading record 'F' of"),
    "iter": (LATIN1_BODY, lambda block: list(block), "Reading record '\ufffd' of"),
    "to_dict": (LATIN1_BODY, lambda block: block.to_dict(), "Reading record 'F' of"),
    "names": (LATIN1_BODY, lambda block: block.names(), "Reading record '\ufffd' of"),
    "str": (LATIN1_BODY, str, "Writing"),
}


# Public SDBlock names that read no block text: constructors, record edits
# that compare and write bytes, and counts parsed from ASCII digits.
NOT_DECODING = {
    "from_lines",
    "from_sdf",
    "append_record",
    "set_record",
    "remove_record",
    "replace_records",
    "num_atoms",
    "num_bonds",
}


def test_decoding_accessors_are_listed() -> None:
    """Every public SDBlock name is listed as decoding block text or not."""
    public = {name for name in dir(SDBlock) if not name.startswith("_")}
    assert public == DECODING.keys() - {"getitem", "iter", "str"} | NOT_DECODING


@pytest.mark.parametrize(
    ("data", "read", "action"),
    [
        *DECODING.values(),
        *(
            (
                LATIN1_BODY,
                lambda block, dest=dest: write_sdf(dest, [block]),
                "Writing block 0 (from 0),",
            )
            for dest in (io.StringIO(), os.devnull)
        ),
    ],
    ids=[*DECODING, "write_sdf-handle", "write_sdf-path"],
)
def test_decode_errors_name_the_block_and_the_fix(
    tmp_path: Path,
    data: bytes,
    read: typing.Callable[[SDBlock[str]], object],
    action: str,
) -> None:
    """Each decoding accessor's UnicodeDecodeError notes the block and fix.

    Args:
        tmp_path: Temporary directory for the SDF fixture.
        data: the non-UTF-8 SDF
        read: the accessor call
        action: what the note says was being done
    """
    path = tmp_path / "latin1.sdf"
    path.write_bytes(data)
    (block,) = parse_sdf(path)
    title = "\ufffd" if data is LATIN1_TITLE else "first"
    with pytest.raises(UnicodeDecodeError) as error:
        read(block)
    assert error.value.__notes__ == [
        f"{action} the SDBlock titled '{title}'. {DECODE_HINT}"
    ]


def test_byte_level_methods_work_on_undecodable_text(tmp_path: Path) -> None:
    """Methods that read no block text work on a str block that is not UTF-8.

    Args:
        tmp_path: Temporary directory for the SDF fixture.
    """
    path = tmp_path / "latin1.sdf"
    path.write_bytes(LATIN1_TITLE + LATIN1_BODY)
    bad_title, block = parse_sdf(path)
    assert repr(bad_title) == "<SDBlock title='\ufffd'>"
    assert (block.num_atoms, block.num_bonds) == (0, 0)
    assert len(block) == 2 and "F" in block and "G" not in block
    block.set_record("F", "ok")
    block.append_record("G", "1")
    assert block.remove_record("G") == 1
    assert block.__reduce__()[1][0].endswith(b"> <F>\nok\n\n> <\xa6>\nv\n\n$$$$\n")
    block.replace_records({"F": "new"})
    assert block.records() == [("F", "new")]
    with pytest.raises(UnicodeDecodeError):
        block.mdl


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
        out = str.join("", map(str, blocks))
        assert out == text
        text = out


def test_bytes_blocks_keep_crlf(tmp_path: Path) -> None:
    """Untouched CRLF bytes blocks write back byte for byte, and edits stay CRLF.

    Args:
        tmp_path: pytest fixture of a temporary directory
    """
    data = (
        b" t1 \r\n  prog\r\n\r\nM  END\r\n> <A>\r\n1\r\n\r\n$$$$\r\n"
        b"t2\r\nM  END\r\n> <B>\r\n\xe9\r\n\r\n$$$$\r\n"
    )
    source = tmp_path / "crlf.sdf"
    source.write_bytes(data)
    blocks = list(parse_sdf(os.fsencode(source)))
    out = io.BytesIO()
    assert write_sdf(out, blocks) == 2
    assert out.getvalue() == data
    assert blocks[0].molblock == b" t1 \r\n  prog\r\n\r\nM  END\r\n"

    blocks[1].title = b"new"
    blocks[1].append_record(b"C", b"3")
    assert sdf_of(blocks[1]) == (
        b"new\r\nM  END\r\n> <B>\r\n\xe9\r\n\r\n> <C>\r\n3\r\n\r\n$$$$\r\n"
    )
    assert sdf_of(SDBlock(b"t\r\nM  END\r\n")) == b"t\r\nM  END\r\n$$$$\r\n"

    # A file without a final newline still writes CRLF throughout.
    source.write_bytes(data.removesuffix(b"\r\n"))
    out = io.BytesIO()
    write_sdf(out, parse_sdf(os.fsencode(source)))
    assert out.getvalue() == data
    with open(source, "rb") as handle:
        lines = handle.readlines()
    assert bytes.join(b"", (sdf_of(b) for b in SDBlock.from_lines(lines))) == data


def test_titled_mol_read(single_titled_mol_sdf: Pathy) -> None:
    """An sdf of a single titled molecule contains it with the expected title.

    Args:
        single_titled_mol_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock[str] = next(parse_sdf(single_titled_mol_sdf))
    assert mol.title == "Title"


def test_first_delimiter_title_ends_a_record(
    double_first_delimiter_titled_mol_sdf: Path,
) -> None:
    """A $$$$ title line ends a record, as the CTfile spec requires.

    The spec forbids $$$$ in a title and ends a compound at any line that
    begins with $$$$, so the title line is an empty record of its own, and
    the rest of that molecule is another block.

    Args:
        double_first_delimiter_titled_mol_sdf: pytest fixture of a Path to sdf
    """
    blocks = list(parse_sdf(double_first_delimiter_titled_mol_sdf))
    assert [block.title for block in blocks][::2] == ["", "Title"]
    assert len(blocks) == 3
    assert sdf_of(blocks[0]) == "$$$$\n"
    assert_rereads_alike(blocks)


def test_second_delimiter_title_ends_a_record(
    double_second_delimiter_titled_mol_sdf: Path,
) -> None:
    """A $$$$ title line after a record is a block of its own.

    Args:
        double_second_delimiter_titled_mol_sdf: pytest fixture of a Path to sdf
    """
    blocks = list(parse_sdf(double_second_delimiter_titled_mol_sdf))
    assert [block.title for block in blocks][:2] == ["Title", ""]
    assert len(blocks) == 3
    assert sdf_of(blocks[1]) == "$$$$\n"
    assert_rereads_alike(blocks)


def assert_rereads_alike(blocks: list[SDBlock[str]]) -> None:
    """Check that the written blocks read back as the same blocks.

    Args:
        blocks: parsed blocks
    """
    written = [sdf_of(block) for block in blocks]
    reread = SDBlock.from_sdf(str.join("", written))
    assert [sdf_of(block) for block in reread] == written


def test_single_record_mol_read(single_record_mol_sdf: Pathy) -> None:
    """An sdf of a delimiter record valued molecule contains the expected data.

    Args:
        single_record_mol_sdf: pytest fixture of a Path to the sdf
    """
    mol = next(parse_sdf(single_record_mol_sdf))
    records = list(mol.records())
    assert ("Record", "Value") in records


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
    path.write_bytes(text.encode())
    data = text.encode()
    expected = [sdf_of(b) for b in SDBlock.from_sdf(text)]
    assert [sdf_of(b) for b in parse_sdf(path)] == expected
    with path.open(encoding="utf-8") as handle:
        assert [sdf_of(b) for b in SDBlock.from_lines(handle)] == expected
    expected_bytes = [sdf_of(b) for b in SDBlock.from_sdf(data)]
    assert [sdf_of(b) for b in parse_sdf(os.fsencode(path))] == expected_bytes
    assert [sdf_of(b) for b in SDBlock.from_lines(io.BytesIO(data))] == (expected_bytes)
    assert [b.records() for b in SDBlock.from_sdf(text)] == [
        [("A", "1")],
        [("B", value)],
    ]
    assert list(SDBlock.from_sdf("")) == []
    with pytest.raises(TypeError, match="str or bytes"):
        SDBlock.from_sdf(["lines"])  # type: ignore[call-overload]


def test_from_sdf_translates_many_crlf_lines() -> None:
    """A str with many CRLF lines, shorter and longer than 16 bytes, reads as LF."""
    lf = str.join("", [f"t{i}\nM  END\n> <A>\n{'v' * i}\n\n$$$$\n" for i in range(40)])
    crlf = lf.replace("\n", "\r\n")
    expected = [sdf_of(block) for block in SDBlock.from_sdf(lf)]
    assert [sdf_of(block) for block in SDBlock.from_sdf(crlf)] == expected
    assert str.join("", expected) == lf


def test_delimiter_value_ends_a_record(
    single_delimiter_record_mol_sdf: Path,
) -> None:
    """A $$$$ value line ends the record, as the CTfile spec requires.

    Args:
        single_delimiter_record_mol_sdf: pytest fixture of a Path to the sdf
    """
    blocks = list(parse_sdf(single_delimiter_record_mol_sdf))
    assert blocks[0].records() == [("Record", "")]
    assert str.join("", (sdf_of(block) for block in blocks)) == (
        single_delimiter_record_mol_sdf.read_text(encoding="utf-8")
    )


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


def test_from_lines_reads_read1_only_file_objects() -> None:
    """A file object with only read1, and no line iteration, is read."""

    class Reader:
        def __init__(self, data: bytes) -> None:
            self.data = data

        def read1(self, size: int, /) -> bytes:
            chunk, self.data = self.data[:size], self.data[size:]
            return chunk

    blocks = SDBlock.from_lines(Reader(b"First\n\nM  END\n$$$$\n"))
    assert [block.title for block in blocks] == [b"First"]
