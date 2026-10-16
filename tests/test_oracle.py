"""Differential tests: the native parser against the pure-Python reference.

Hypothesis generates SDF-shaped inputs (titles, MDL bodies with and without
``M  END``, record headers and malformed headers, multiline and blank values,
``$$$$`` variants at and after the start of a line, LF/CRLF/CR endings,
unterminated input) and every entry point must agree with
``tests/reference_sdblock.py`` block by block. The readers must also agree
with each other: mapped and streamed paths, file objects and their lines,
including input past the read chunk and the map's resident window, and
blocks borrowed from a mapped file.

The alphabet includes whitespace that only ``str.strip`` removes, such as
``\x1f`` and the ideographic space, which neither mode strips, but not
``str.splitlines`` breaks.
"""

from __future__ import annotations

import copy
import io
import os
import pickle
import tempfile
import threading
from typing import Any, Callable, Iterable, Union

import pytest
from hypothesis import HealthCheck, assume, given, settings
from hypothesis import strategies as st

from sendoff.sdblock import SDBlock, parse_sdf
from tests.fixtures import sdf_of
from tests.reference_sdblock import RefBlock, ref_from_lines

ALPHABET = list("abXZ09 \t<>$-_.()é中\x1f\xa0\u2003\u3000")
SPECIAL = [
    "",
    " ",
    "\t ",
    "\u3000",
    " \xa0\x1f",
    "M  END",
    "M  END  ",
    "M  CHG  1   2  -1",
    "$$$$",
    "$$$$ ",
    "$$$$x",
    "> <A>",
    "> <B>  (1)",
    ">  <A>",
    "  > <C>",
    "> <A><B>",
    "> <>",
    "> <é名>",
]
HEADERS = ["> <A>", "> <B>", "> <A>  (7)", ">  < spaced name >", "> <multi<angle>>"]
BAD_HEADERS = ["> no angle", ">  ", "> "]
ENDINGS = ["\n"] * 6 + ["\r\n"] * 3 + ["\r"]

free_text = st.text(alphabet=ALPHABET, max_size=10)
line_content = st.one_of(
    st.sampled_from(SPECIAL),
    free_text,
    st.builds(lambda a, b: a + b, st.sampled_from(SPECIAL), free_text),
    st.builds(lambda a, b: a + b, free_text, st.sampled_from(SPECIAL)),
)
value_line = st.one_of(free_text.filter(lambda s: s.strip() != ""), line_content)


@st.composite
def sdf_lines(draw: Callable[..., Any]) -> list[str]:
    """Draw newline-terminated SDF-shaped lines; the last may be unterminated.

    Args:
        draw: Hypothesis draw function

    Returns:
        Lines with their endings.
    """
    contents: list[str] = []
    for _ in range(draw(st.integers(0, 4))):
        contents.append(draw(line_content))
        contents.extend(draw(st.lists(line_content, max_size=4)))
        if draw(st.integers(0, 4)):
            contents.append("M  END")
        for _ in range(draw(st.integers(0, 3))):
            header = st.sampled_from(HEADERS + BAD_HEADERS * draw(st.integers(0, 1)))
            contents.append(draw(header))
            contents.extend(draw(st.lists(value_line, max_size=3)))
            if draw(st.integers(0, 3)):
                contents.append(draw(st.sampled_from(["", " ", "\t", "\u3000"])))
        if draw(st.integers(0, 5)):
            contents.append(draw(st.sampled_from(["$$$$", "$$$$", "$$$$ ", "$$$$x"])))
    contents.extend(draw(st.lists(line_content, max_size=2)))
    endings = [draw(st.sampled_from(ENDINGS)) for _ in contents]
    if endings and draw(st.booleans()):
        endings[-1] = ""
    return [content + ending for content, ending in zip(contents, endings)]


def outcome(call: Callable[[], Any]) -> Any:
    """Return a call's result, or a marker when it raises any Exception.

    Args:
        call: zero-argument callable

    Returns:
        The result, or ``"<error>"``; exception types may differ by design.
    """
    try:
        return call()
    except Exception:
        return "<error>"


def serialized(block: Any, text: bool) -> Any:
    """Return a block's SDF.

    Args:
        block: a native SDBlock
        text: whether the block holds str values

    Returns:
        str(block) for a str block, else bytes(block).
    """
    return str(block) if text else bytes(block)


def assert_blocks_match(
    native_blocks: Iterable[SDBlock[Any]], reference: list[RefBlock[Any]], text: bool
) -> None:
    """Compare block count, titles, records, writes, and appended writes.

    Args:
        native_blocks: blocks from the native parser
        reference: blocks from the reference parser
        text: whether the blocks hold str values
    """
    native = list(native_blocks)
    assert len(native) == len(reference)
    for index, (block, expected) in enumerate(zip(native, reference)):
        context = f"block {index}"
        assert block.title == expected.title, context
        assert outcome(lambda: list(block.records())) == outcome(expected.records), (
            context
        )
        assert outcome(lambda: list(block.iter_records())) == outcome(
            expected.records
        ), context
        written = serialized(block, text)
        assert written == expected.serialize(), context
        assert sdf_of(block) == written, context
        if cast_like(written, "\r") not in written:  # SDBlock() translates str \r
            assert serialized(SDBlock(written), text) == written, context
        assert block.molblock == expected.molblock(), context
        newline: Any = "\n" if text else b"\n"
        assert block.mdl == expected.molblock().partition(newline)[2], context
        names = outcome(lambda: [name for name, _ in expected.records()])
        assert outcome(block.names) == names, context
        count = len(names) if isinstance(names, list) else names
        assert outcome(lambda: len(block)) == count, context
        expected_dict = outcome(lambda: dict(expected.records()))
        assert outcome(block.to_dict) == expected_dict, context
        for field in ["A", "B", " spaced name ", "missing"]:
            key = field if text else field.encode()
            assert outcome(lambda: block.get(key, 0)) == outcome(
                lambda: expected_dict.get(key, 0)
            ), context
            assert outcome(lambda: key in block) == outcome(
                lambda: key in expected_dict.keys()
            ), context
        name, value = ("Z", "v\u00e9") if text else (b"Z", "v\u00e9".encode())
        for header in ("padded", "copy"):
            styled, model = copy.copy(block), copy.deepcopy(expected)
            styled.append_record(name, value, header=header)
            model.append_record(name, value, header)
            assert serialized(styled, text) == model.serialize(), context
        block.append_record(name, value)
        expected.append_record(name, value)
        assert serialized(block, text) == expected.serialize(), context
        assert outcome(lambda: list(block.records())) == outcome(expected.records), (
            context
        )
        assert_mutations_match(block, text)


def cast_like(like: Any, value: str) -> Any:
    """Return ``value`` as str or bytes to match ``like``.

    Args:
        like: a str or bytes value
        value: the text to convert

    Returns:
        ``value`` in the type of ``like``.
    """
    return value if isinstance(like, str) else value.encode()


def reparsed_records(block: SDBlock[Any], text: bool) -> Any:
    """Return the records the reference parses from what ``block`` writes.

    Args:
        block: a native SDBlock
        text: whether the block holds str values

    Returns:
        The reference parser's records for the written block.
    """
    written = serialized(block, text)
    mdl = block.molblock.split(cast_like(written, "\n"))[1:]
    if not any(line.startswith(cast_like(written, "M  END")) for line in mdl):
        return block.records()  # records after MDL without M  END read as MDL
    handle: Any = io.StringIO(written) if text else io.BytesIO(written)
    (reference,) = ref_from_lines(handle)
    return reference.records()


def assert_mutations_match(block: SDBlock[Any], text: bool) -> None:
    """Check record and title mutators against a list model of records.

    Args:
        block: a native SDBlock, which is mutated
        text: whether the block holds str values
    """

    def cast(value: str) -> Any:
        return value if text else value.encode()

    before = outcome(block.records)
    name, other, value = cast("A"), cast(" spaced name "), cast("new\u00e9 v")
    if before == "<error>":
        assert outcome(lambda: block.set_record(name, value)) == "<error>"
        assert outcome(lambda: block.remove_record(name)) == "<error>"
    else:
        expected = [(n, v) for n, v in before if n != name]
        matches = [i for i, (n, _) in enumerate(before) if n == name]
        position = matches[0] if matches else len(before)
        expected.insert(position, (name, value))
        block.set_record(name, value)
        assert block.records() == expected
        assert reparsed_records(block, text) == expected
        removed = sum(n == other for n, _ in expected)
        assert block.remove_record(other) == removed
        expected = [(n, v) for n, v in expected if n != other]
        assert block.records() == expected
        assert reparsed_records(block, text) == expected

    title: Any = block.title
    molblock: Any = block.molblock
    replacement = [(cast("B"), cast("1")), (cast("B"), cast("x\ny")), (name, value)]
    block.replace_records(replacement)
    assert block.records() == replacement
    block.replace_records(dict(replacement))
    assert block.records() == [(cast("B"), cast("x\ny")), (name, value)]
    assert reparsed_records(block, text) == block.records()
    assert (block.title, block.molblock) == (title, molblock)

    line = cast(" \u3000new title\t" if text else " new title\t")
    block.title = line
    assert block.title == line.rstrip()
    old_line = molblock.partition(cast("\n"))[0].removesuffix(cast("\r"))
    assert block.molblock == line + molblock.removeprefix(old_line)
    with pytest.raises(ValueError, match="cannot start with"):
        block.title = cast("$$$$ new")
    assert block.title == line.rstrip()
    assert block.records() == [(cast("B"), cast("x\ny")), (name, value)]


class ChunkedReader:
    """A file-like source whose ``read`` returns a few items at a time."""

    def __init__(self, data: Union[str, bytes], sizes: list[int]) -> None:
        """Serve ``data`` in chunks cycling through ``sizes``.

        Args:
            data: the whole str or bytes content
            sizes: chunk sizes to cycle through
        """
        self.data = data
        self.sizes = sizes
        self.pos = 0
        self.calls = 0

    def read(self, size: int) -> Union[str, bytes]:
        """Return at most ``size`` items, and an empty chunk at the end.

        Args:
            size: the most items the caller wants

        Returns:
            The next chunk.
        """
        step = min(size, self.sizes[self.calls % len(self.sizes)])
        self.calls += 1
        start, end = self.pos, self.pos + step
        chunk = self.data[start:end]
        self.pos = end
        return chunk


examples = settings(
    max_examples=int(os.environ.get("SENDOFF_ORACLE_EXAMPLES", "300")),
    deadline=None,
    suppress_health_check=[HealthCheck.too_slow],
)


@examples
@given(sdf_lines())
def test_from_text_lines_matches_reference(lines: list[str]) -> None:
    """Text lines with their endings parse like the reference.

    Args:
        lines: generated lines with endings
    """
    assert_blocks_match(
        SDBlock.from_lines(lines), list(ref_from_lines(lines)), text=True
    )


@examples
@given(sdf_lines())
def test_from_bytes_lines_matches_reference(lines: list[str]) -> None:
    """Byte lines with their endings parse like the reference.

    Args:
        lines: generated lines with endings
    """
    encoded = [line.encode() for line in lines]
    assert_blocks_match(
        SDBlock.from_lines(encoded), list(ref_from_lines(encoded)), text=False
    )


@examples
@given(sdf_lines())
def test_from_splitlines_matches_reference(lines: list[str]) -> None:
    """Lines without endings, as from ``str.splitlines()``, parse alike.

    Args:
        lines: generated lines with endings
    """
    stripped = str.join("", lines).splitlines()
    assert_blocks_match(
        SDBlock.from_lines(stripped), list(ref_from_lines(stripped)), text=True
    )


@examples
@given(sdf_lines(), st.lists(st.integers(1, 9), min_size=1, max_size=5))
def test_file_objects_match_reference(lines: list[str], sizes: list[int]) -> None:
    """File objects are read in chunks of any size and parse like their lines.

    Args:
        lines: generated lines with endings
        sizes: chunk sizes for the reader to cycle through
    """
    text = str.join("", lines)
    data = text.encode()
    text_sources: list[Any] = [ChunkedReader(text, sizes), io.StringIO(text)]
    for text_source in text_sources:
        text_reference = list(ref_from_lines(io.StringIO(text)))
        assert_blocks_match(SDBlock.from_lines(text_source), text_reference, True)
    bytes_sources: list[Any] = [ChunkedReader(data, sizes), io.BytesIO(data)]
    for bytes_source in bytes_sources:
        bytes_reference = list(ref_from_lines(io.BytesIO(data)))
        assert_blocks_match(SDBlock.from_lines(bytes_source), bytes_reference, False)


@examples
@given(sdf_lines())
def test_from_sdf_matches_reference(lines: list[str]) -> None:
    """In-memory str uses universal newlines; bytes split on LF only.

    Args:
        lines: generated lines with endings
    """
    text = str.join("", lines)
    data = text.encode()
    text_reference = list(ref_from_lines(io.StringIO(text, newline=None)))
    assert_blocks_match(SDBlock.from_sdf(text), text_reference, text=True)
    bytes_reference = list(ref_from_lines(io.BytesIO(data)))
    assert_blocks_match(SDBlock.from_sdf(data), bytes_reference, text=False)


@pytest.mark.parametrize("borrow", [False, True])
@examples
@given(sdf_lines())
def test_parse_sdf_matches_reference(borrow: bool, lines: list[str]) -> None:
    """Text paths use universal newlines; bytes paths split on LF only.

    Blocks borrowed from the mapped file read and edit as copied ones do.

    Args:
        borrow: whether blocks borrow the mapped file
        lines: generated lines with endings
    """
    data = str.join("", lines).encode()
    with tempfile.TemporaryDirectory() as directory:
        path = os.path.join(directory, "input.sdf")
        with open(path, "wb") as handle:
            handle.write(data)
        with open(path, encoding="utf-8") as handle:
            text_reference = list(ref_from_lines(handle))
        assert_blocks_match(parse_sdf(path, borrow=borrow), text_reference, text=True)
        bytes_reference = list(ref_from_lines(io.BytesIO(data)))
        raw = parse_sdf(os.fsencode(path), borrow=borrow)
        assert_blocks_match(raw, bytes_reference, text=False)


# Sizes that change how readers buffer input: the chunk that file objects
# and unmapped files are read in (CHUNK in rust/src/python/sdblock.rs), and
# the input a map keeps resident behind its cursor (RESIDENT in
# rust/src/core/mapped.rs), past twice which it drops pages.
READ_CHUNK = 128 << 10
MAP_RESIDENT = 1 << 20


def snapshot(blocks: Iterable[SDBlock[Any]], text: bool) -> list[Any]:
    """Return what tells blocks apart: what they write, and their records.

    Args:
        blocks: native blocks
        text: whether the blocks hold str values

    Returns:
        One (written, records) pair per block.
    """
    return [(serialized(block, text), outcome(block.records)) for block in blocks]


def piped(path: Union[str, bytes], data: bytes) -> list[SDBlock[Any]]:
    """Parse ``data`` written through a named pipe, which is read, not mapped.

    Args:
        path: where to make the pipe; a bytes path parses bytes
        data: the SDF content

    Returns:
        The blocks ``parse_sdf`` produces.
    """
    os.mkfifo(path)

    def write() -> None:
        with open(path, "wb") as pipe:
            pipe.write(data)

    writer = threading.Thread(target=write, daemon=True)
    writer.start()
    try:
        return list(parse_sdf(path))
    finally:
        writer.join()
        os.remove(path)


def assert_readers_agree(text: str) -> None:
    """Every reader splits ``text``, saved as UTF-8, into the same blocks.

    Readers that decode text themselves translate universal newlines;
    ``from_lines`` keeps the endings it is given, so it agrees with them
    only through a file in default text mode. A file opened with
    ``newline=""`` splits its own lines at a lone "\\r" too, so
    ``from_lines(f)`` and ``from_lines(list(f))`` agree only without one.

    Args:
        text: SDF content
    """
    data = text.encode()
    pipes = hasattr(os, "mkfifo")
    with tempfile.TemporaryDirectory() as directory:
        path = os.path.join(directory, "input.sdf")
        with open(path, "wb") as out:
            out.write(data)
        pipe = os.path.join(directory, "pipe.sdf")

        def opened(**kwargs: Any) -> Any:
            return open(path, encoding="utf-8", **kwargs)

        with opened() as handle, opened() as lines:
            translated = [
                snapshot(parse_sdf(path), True),
                snapshot(parse_sdf(path, borrow=True), True),
                snapshot(SDBlock.from_sdf(text), True),
                snapshot(SDBlock.from_lines(io.StringIO(text, newline=None)), True),
                snapshot(SDBlock.from_lines(handle), True),
                snapshot(SDBlock.from_lines(list(lines)), True),
            ]
        if pipes:
            translated.append(snapshot(piped(pipe, data), True))
        assert all(blocks == translated[0] for blocks in translated)

        with opened(newline="") as handle, opened(newline="") as lines:
            kept = [
                snapshot(SDBlock.from_lines(handle), True),
                snapshot(SDBlock.from_lines(io.StringIO(text)), True),
            ]
            split_at_cr = snapshot(SDBlock.from_lines(list(lines)), True)
        assert kept[0] == kept[1]
        lone_cr = "\r" in text.replace("\r\n", "")
        if not lone_cr:
            assert split_at_cr == kept[0]
        if "\r" not in text:
            assert kept[0] == translated[0]

        with open(path, "rb") as handle, open(path, "rb") as lines:
            raw = [
                snapshot(parse_sdf(os.fsencode(path)), False),
                snapshot(parse_sdf(os.fsencode(path), borrow=True), False),
                snapshot(SDBlock.from_sdf(data), False),
                snapshot(SDBlock.from_lines(io.BytesIO(data)), False),
                snapshot(SDBlock.from_lines(handle), False),
                snapshot(SDBlock.from_lines(list(lines)), False),
            ]
        if pipes:
            raw.append(snapshot(piped(os.fsencode(pipe), data), False))
        assert all(blocks == raw[0] for blocks in raw)


@examples
@given(sdf_lines())
def test_readers_agree(lines: list[str]) -> None:
    """Paths, pipes, in-memory SDF, file objects, and lines split alike.

    Args:
        lines: generated lines with endings
    """
    assert_readers_agree(str.join("", lines))


@settings(
    max_examples=max(1, int(os.environ.get("SENDOFF_ORACLE_EXAMPLES", "300")) // 60),
    deadline=None,
    suppress_health_check=[HealthCheck.too_slow],
)
@given(sdf_lines(), st.integers(0, 512))
def test_readers_agree_past_buffer_sizes(lines: list[str], offset: int) -> None:
    """Readers agree when input crosses chunk and map window boundaries.

    A padding block puts the generated SDF ``offset`` characters before a
    multiple of the read chunk, past twice the map's resident window.

    Args:
        lines: generated lines with endings
        offset: where the generated SDF starts before the boundary
    """
    text = str.join("", lines)
    assume(text)
    head, tail = "pad\n", "\n$$$$\n"
    fill = 2 * MAP_RESIDENT + READ_CHUNK - offset - len(head) - len(tail)
    padding = head + ("x" * 79 + "\n") * (fill // 80) + "x" * (fill % 80) + tail
    assert (len(padding) + offset) % READ_CHUNK == 0
    assert_readers_agree(padding + text * 2)


def assert_round_trips(block: SDBlock[Any]) -> None:
    """Assert copies, pickles, and molblock reassignment keep ``block`` as documented.

    Args:
        block: a block a reader produced
    """
    original = copy.copy(block)
    text = bool(original.__reduce__()[1][3])
    name: Any = "Z" if text else b"Z"
    for duplicate in (
        copy.copy(block),
        copy.deepcopy(block),
        pickle.loads(pickle.dumps(block)),
    ):
        assert duplicate == block and duplicate is not block
        duplicate.append_record(name, name)
        assert block == original
    try:
        molblock = block.molblock
    except UnicodeDecodeError:
        return
    block.molblock = molblock
    data, _, mdl_end, _ = block.__reduce__()[1]
    old_data, _, old_mdl_end, _ = original.__reduce__()[1]
    assert data[mdl_end:] == old_data[old_mdl_end:]
    if text:
        if "\r" in molblock.replace("\r\n", ""):
            return  # translating a lone "\r" splits its line
        expected = molblock.replace("\r\n", "\n")
    else:
        expected = molblock
    assert block.molblock == expected
    assert block.title == original.title
    raw = old_data[:old_mdl_end]
    if (expected.encode() if text else expected) == raw:
        assert block == original
    reassigned = copy.copy(block)
    reassigned.molblock = block.molblock
    assert reassigned == block


@examples
@given(sdf_lines())
def test_blocks_round_trip(lines: list[str]) -> None:
    """Blocks from every reader survive copy, pickle, and molblock = molblock.

    Args:
        lines: generated lines with endings
    """
    text = str.join("", lines)
    blocks: list[SDBlock[Any]] = [
        *SDBlock.from_lines(lines),
        *SDBlock.from_sdf(text),
        *SDBlock.from_lines(io.StringIO(text)),
        *SDBlock.from_lines([line.encode() for line in lines]),
    ]
    with tempfile.TemporaryDirectory() as directory:
        path = os.path.join(directory, "latin1.sdf")
        with open(path, "wb") as out:
            out.write(text.encode().replace("\u00e9".encode(), b"\xe9"))
        blocks.extend(parse_sdf(path))
    for block in blocks:
        assert_round_trips(block)
