"""Differential tests: the native parser against the pure-Python reference.

Hypothesis generates SDF-shaped inputs (titles, MDL bodies with and without
``M  END``, record headers and malformed headers, multiline and blank values,
``$$$$`` variants, LF/CRLF/CR endings, unterminated input) and every entry
point must agree with ``tests/reference_sdblock.py`` block by block.

The alphabet includes whitespace that only ``str.strip`` removes, such as
``\x1f`` and the ideographic space, but not ``str.splitlines`` breaks.
"""

from __future__ import annotations

import io
import os
import tempfile
from typing import Any, Callable, Iterable, Union

from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

from sendoff.sdblock import SDBlock, parse_sdf
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
    """Return what ``block.write`` writes.

    Args:
        block: a native SDBlock
        text: whether the block holds str values

    Returns:
        The written str or bytes.
    """
    handle: io.StringIO | io.BytesIO = io.StringIO() if text else io.BytesIO()
    block.write(handle)
    return handle.getvalue()


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
        assert outcome(lambda: list(block.records())) == outcome(
            expected.records
        ), context
        assert outcome(lambda: list(block.iter_records())) == outcome(
            expected.records
        ), context
        written = serialized(block, text)
        assert written == expected.serialize(), context
        assert block.to_sdf() == written, context
        if not (
            block.title.startswith(cast_like(written, "$$$$"))  # reads as delimiter
            or cast_like(written, "\r") in written  # SDBlock() translates str \r
        ):
            assert serialized(SDBlock(written), text) == written, context
        assert block.molblock == expected.molblock(), context
        newline: Any = "\n" if text else b"\n"
        assert block.mdl == expected.molblock().partition(newline)[2], context
        for needle in ["M  END", "END\n", "a", "\u3000", "", "\n"]:
            key = needle if text else needle.encode()
            assert block.mdl_contains(key) == (key in block.mdl), context
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
        block.append_record(name, value)
        expected.append_record(name, value)
        assert serialized(block, text) == expected.serialize(), context
        assert outcome(lambda: list(block.records())) == outcome(
            expected.records
        ), context
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
    if block.title.startswith(cast_like(written, "$$$$")):
        return block.records()  # a written $$$$ title reads as a delimiter
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

    block.title = cast(" \u3000new title\t" if text else " new title\t")
    assert block.title == cast("new title")
    assert block.molblock == cast("new title") + molblock.removeprefix(title)
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
    stripped = "".join(lines).splitlines()
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
    text = "".join(lines)
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
    text = "".join(lines)
    data = text.encode()
    text_reference = list(ref_from_lines(io.StringIO(text, newline=None)))
    assert_blocks_match(SDBlock.from_sdf(text), text_reference, text=True)
    bytes_reference = list(ref_from_lines(io.BytesIO(data)))
    assert_blocks_match(SDBlock.from_sdf(data), bytes_reference, text=False)


@examples
@given(sdf_lines())
def test_parse_sdf_matches_reference(lines: list[str]) -> None:
    """Text paths use universal newlines; bytes paths split on LF only.

    Args:
        lines: generated lines with endings
    """
    data = "".join(lines).encode()
    with tempfile.TemporaryDirectory() as directory:
        path = os.path.join(directory, "input.sdf")
        with open(path, "wb") as handle:
            handle.write(data)
        with open(path, encoding="utf-8") as handle:
            text_reference = list(ref_from_lines(handle))
        assert_blocks_match(parse_sdf(path), text_reference, text=True)
        bytes_reference = list(ref_from_lines(io.BytesIO(data)))
        assert_blocks_match(parse_sdf(os.fsencode(path)), bytes_reference, text=False)
