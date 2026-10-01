"""Freeze raw SD block parsing, mutable state, and iterable/write boundaries."""

from __future__ import annotations

import gc
import io
import os
import weakref
from collections import deque
from pathlib import Path
from types import GeneratorType
from typing import Any, Callable, Iterable, Iterator, SupportsIndex, cast

import pytest

from sendoff.ctable import CTableFormat
from sendoff.sdblock import Pathy, SDBlock, parse_sdf
from tests.compat_literals import V2000, V3000


def test_literal_block_preserves_opaque_lines_and_declared_counts() -> None:
    """Only the title is stripped; chemistry bodies are not validated."""
    for text, expected_format in (
        (V2000, CTableFormat.V2000),
        (V3000, CTableFormat.V3000),
    ):
        raw = text.splitlines(keepends=True)
        block = SDBlock.from_block_lines(raw + ["> <key>\n", " value \n", "\n"])
        assert block.title == "literal title"
        assert block.mdl == deque(raw[1:])
        assert block.metadata == deque(["> <key>\n", " value \n", "\n"])
        assert block.ctable().format is expected_format
        assert (block.num_atoms(), block.num_bonds()) == (2, 1)
        assert list(block.records()) == [("key", "value")]
    malformed = SDBlock.from_block_lines(
        ["t", "source", "comment", " 99 -2 V2000", "not chemistry", "M  END"]
    )
    assert (malformed.num_atoms(), malformed.num_bonds()) == (99, -2)
    assert malformed.mdl[-2] == "not chemistry"


def test_section_generators_are_lazy_and_consume_only_through_prefix_markers() -> None:
    """Keep prefix recognition and the external iterator's remaining tail."""
    consumed: list[str] = []

    def source() -> Iterator[str]:
        for line in ["atom", "M  END extra", "value", "$$$$ extra", "tail"]:
            consumed.append(line)
            yield line

    lines = source()
    mdl = iter(SDBlock.parse_mdl(lines))
    assert isinstance(mdl, GeneratorType)
    assert consumed == []
    assert next(mdl) == "atom"
    assert consumed == ["atom"]
    assert next(mdl) == "M  END extra"
    assert list(mdl) == []
    metadata = iter(SDBlock.parse_metadata(lines))
    assert isinstance(metadata, GeneratorType)
    assert consumed == ["atom", "M  END extra"]
    assert list(metadata) == ["value"]
    assert next(lines) == "tail"
    assert list(lines) == []
    assert list(SDBlock.parse_mdl([" M  END", "tail"])) == [" M  END", "tail"]
    assert list(SDBlock.parse_metadata([" $$$$", "tail"])) == [" $$$$", "tail"]


def test_block_factory_is_eager_and_does_not_consume_past_metadata_delimiter() -> None:
    """One iterator is shared by the two parsers, even for repeatable inputs."""
    lines = iter(["  t  ", "raw", "M  END", "> <key>", "v", "$$$$", "tail"])
    block = SDBlock.from_block_lines(lines)
    assert block == SDBlock("t", deque(["raw", "M  END"]), deque(["> <key>", "v"]))
    assert next(lines) == "tail"
    repeatable = ["t", "M  END", "metadata", "$$$$", "tail"]
    assert SDBlock.from_block_lines(repeatable).metadata == deque(["metadata"])
    assert repeatable[-1] == "tail"
    with pytest.raises(StopIteration) as caught:
        SDBlock.from_block_lines([])
    assert caught.value.args == ()
    assert SDBlock.from_block_lines(["title"]) == SDBlock("title", deque(), deque())
    missing_end = SDBlock.from_block_lines(["t", "raw", "> <key>", "v", "$$$$"])
    assert missing_end.mdl == deque(["raw", "> <key>", "v", "$$$$"])
    assert missing_end.metadata == deque()


def test_stream_factory_is_lazy_and_drops_unterminated_tail() -> None:
    """Freeze sendoff-712.1 without fixing framing or reading ahead a block."""
    consumed: list[str] = []

    def source() -> Iterator[str]:
        for line in ["one", "M  END", "$$$$", "two", "M  END", "$$$$", "lost"]:
            consumed.append(line)
            yield line

    external = source()
    blocks = SDBlock.from_lines(external)
    assert isinstance(blocks, GeneratorType)
    assert iter(blocks) is blocks
    assert consumed == []
    assert next(blocks).title == "one"
    assert consumed == ["one", "M  END", "$$$$"]
    assert next(blocks).title == "two"
    assert consumed[-3:] == ["two", "M  END", "$$$$"]
    assert list(blocks) == []
    assert consumed[-1] == "lost"
    assert list(external) == []
    assert list(SDBlock.from_lines([])) == []
    assert list(SDBlock.from_lines(V2000.splitlines())) == []
    assert SDBlock.from_block_lines(V2000.splitlines()).num_atoms() == 2


def test_delimiter_prefixes_and_delimiter_valued_data_current_fragmentation() -> None:
    """Freeze literal manifestations of sendoff-712.5; retain existing xfails."""
    blocks = list(
        SDBlock.from_lines(["$$$$ title", "M  END", "> <key>", "$$$$ value", "$$$$"])
    )
    assert blocks == [
        SDBlock("$$$$ title", deque(), deque()),
        SDBlock("M  END", deque(["> <key>", "$$$$ value"]), deque()),
        SDBlock("$$$$", deque(), deque()),
    ]
    assert next(SDBlock.from_lines(["t", "M  END", "$$$$ suffix"])).title == "t"


def test_records_normalization_duplicates_empty_values_and_chunk_boundaries() -> None:
    """Records are ordered tuples, with stripped values but untrimmed inner names."""
    block = SDBlock(
        "t",
        deque(),
        deque(
            [
                "ignored\n",
                "> <hidden>\n",
                "lost\n",
                "\n",
                "  > prefix < key > suffix  \n",
                "  first  \n",
                "\tsecond\t\n",
                "\n",
                "> <key>\n",
                "\n",
                "> <key>\n",
                "last\n",
                "> <not-a-new-record>\n",
                "still value\n",
                "\n",
                "orphan continuation\n",
            ]
        ),
    )
    records = iter(block.records())
    assert isinstance(records, GeneratorType)
    assert list(records) == [
        (" key ", "first\nsecond"),
        ("key", ""),
        ("key", "last\n> <not-a-new-record>\nstill value"),
    ]
    assert list(records) == []
    assert list(block.records())[0] == (" key ", "first\nsecond")


def test_record_header_errors_are_lazy_and_specific() -> None:
    """Malformed field headers are neither eagerly rejected nor silently repaired."""
    block = SDBlock("t", deque(), deque(["> missing brackets", "value"]))
    records = iter(block.records())
    with pytest.raises(IndexError) as caught:
        next(records)
    assert caught.value.args == ("list index out of range",)
    assert list(records) == []
    block.metadata = deque(["> <name", "value"])
    assert list(block.records()) == [("name", "value")]
    block.metadata = deque(["> <>", "value"])
    assert list(block.records()) == [("", "value")]


def test_append_and_direct_metadata_edits_are_visible_to_new_parses() -> None:
    """Append stores exact raw strings; record iteration reflects mutable deques."""
    metadata: deque[str] = deque()
    block = SDBlock("t", deque(), metadata)
    records = iter(block.records())
    assert (
        cast(Callable[[str, str], object], block.append_record)(
            "key", "  first  \nsecond"
        )
        is None
    )
    assert block.metadata is metadata
    assert metadata == deque(["> <key>\n", "  first  \nsecond\n", "\n"])
    assert list(records) == [("key", "first  \nsecond")]
    metadata[1] = " replacement \n"
    metadata.extend(["> <key>", "another", ""])
    assert list(block.records()) == [("key", "replacement"), ("key", "another")]
    metadata.clear()
    assert list(block.records()) == []
    assert cast(Callable[[str, str], object], block.append_record)("empty", "") is None
    assert metadata == deque(["> <empty>\n", "\n", "\n"])
    assert list(block.records()) == [("empty", "")]


def test_active_record_iterator_rejects_container_mutation() -> None:
    """Mutating a deque after iteration starts keeps the legacy RuntimeError."""
    block = SDBlock("t", deque(), deque(["> <one>", "v", "", "> <two>", "v", ""]))
    records = iter(block.records())
    assert next(records) == ("one", "v")
    block.metadata.append("changed")
    with pytest.raises(RuntimeError, match="^deque mutated during iteration$"):
        next(records)


def test_mdl_edits_reparse_counts_without_sharing_ctable_deques() -> None:
    """Blocks reparse current mdl; a previously created CTable stays a snapshot."""
    block = SDBlock.from_block_lines(V2000.splitlines(keepends=True))
    old_table = block.ctable()
    new_table = block.ctable()
    assert old_table is not new_table
    assert old_table.lines is not new_table.lines
    assert old_table.lines is not block.mdl
    mdl = block.mdl
    mdl[2] = "  7  3 V2000\n"
    block.title = " changed title "
    assert (block.num_atoms(), block.num_bonds()) == (7, 3)
    assert block.ctable().title == "changed title"
    assert (old_table.num_atoms, old_table.num_bonds) == (2, 1)
    assert old_table.title == "literal title"
    block.mdl = deque(V3000.splitlines(keepends=True)[1:])
    assert block.ctable().format is CTableFormat.V3000
    assert (block.num_atoms(), block.num_bonds()) == (2, 1)
    assert mdl[2] == "  7  3 V2000\n"


@pytest.mark.parametrize(
    ("with_newlines", "expected"),
    [
        (True, "t\nraw\nalready\ncr\r\n\n> <key>\nv\n$$$$\n"),
        (False, "t\nrawalready\ncr\r> <key>v$$$$\n"),
    ],
)
def test_write_exact_raw_and_newline_modes(with_newlines: bool, expected: str) -> None:
    """Title and delimiter always get LF; all other raw line text is preserved.

    Args:
        with_newlines: Whether to append LF to lines lacking it.
        expected: Exact output for the selected write mode.
    """
    block = SDBlock(
        "t", deque(["raw", "already\n", "cr\r", ""]), deque(["> <key>", "v"])
    )
    output = io.StringIO()
    assert (
        cast(Callable[..., object], block.write)(output, with_newlines=with_newlines)
        is None
    )
    assert output.getvalue() == expected
    assert block.mdl == deque(["raw", "already\n", "cr\r", ""])
    block.title = "raw title\n"
    output = io.StringIO()
    SDBlock(block.title, deque(), deque()).write(cast(Any, output), with_newlines=False)
    assert output.getvalue() == "raw title\n\n$$$$\n"


def test_literal_write_and_reparse_preserve_data_except_title_trimming() -> None:
    """Assert whole output text, not merely line counts or toolkit equivalence."""
    for text in (V2000, V3000):
        block = SDBlock.from_block_lines(text.splitlines(keepends=True))
        block.append_record("key", "value")
        output = io.StringIO()
        block.write(cast(Any, output))
        assert output.getvalue() == (
            "literal title\n" + text.split("\n", 1)[1] + "> <key>\nvalue\n\n$$$$\n"
        )
        assert (
            next(SDBlock.from_lines(output.getvalue().splitlines(keepends=True)))
            == block
        )


def test_external_iteration_and_write_errors_propagate() -> None:
    """Do not wrap external producer or sink failures as chemistry errors."""

    def failing_lines() -> Iterator[str]:
        yield "title"
        raise OSError("producer failed")

    with pytest.raises(OSError, match="^producer failed$"):
        SDBlock.from_block_lines(failing_lines())
    with pytest.raises(OSError, match="^producer failed$"):
        next(SDBlock.from_lines(failing_lines()))

    class FailingWriter:
        def write(self, text: str) -> int:
            raise OSError("sink failed")

    with pytest.raises(OSError, match="^sink failed$"):
        SDBlock("t", deque(), deque()).write(cast(Any, FailingWriter()))


def test_parse_sdf_reads_eagerly_accepts_path_types_and_normalizes_file_newlines(
    tmp_path: Path,
) -> None:
    """The path wrapper reads now, but block parsing remains a lazy generator.

    Args:
        tmp_path: Project-local pytest fixture directory.
    """
    path = tmp_path / "literal.sdf"
    original = V2000 + "$$$$\n"
    paths: tuple[Pathy, ...] = (path, str(path), os.fsencode(path))
    for path_value in paths:
        path.write_bytes(original.replace("\n", "\r\n").encode())
        blocks = parse_sdf(path_value)
        assert isinstance(blocks, GeneratorType)
        path.write_text("changed\nM  END\n$$$$\n")
        block = next(blocks)
        assert block.title == "literal title"
        assert block.mdl == deque(V2000.splitlines(keepends=True)[1:])
        assert list(blocks) == []
    missing = tmp_path / "absent.sdf"
    with pytest.raises(FileNotFoundError) as caught:
        parse_sdf(missing)
    assert caught.value.errno == 2
    assert caught.value.filename == str(missing)
    with pytest.raises(TypeError, match="expected str, bytes or os.PathLike object"):
        parse_sdf(cast(Any, io.StringIO(original)))


def test_external_iterables_are_accepted_without_sequence_access() -> None:
    """Factories require iteration, not indexing, length, or repeatable traversal."""

    class Lines:
        def __iter__(self) -> Iterator[str]:
            yield from V2000.splitlines()
            yield "$$$$"

    input_lines: Iterable[str] = Lines()
    block = next(SDBlock.from_lines(input_lines))
    assert (block.num_atoms(), block.num_bonds()) == (2, 1)


def test_title_unicode_whitespace_and_unchanged_identity() -> None:
    """Keep Python whitespace, unchanged string objects and lone surrogates."""
    assert SDBlock.from_block_lines(["\u001c\u00a0title\u2003\u001f"]).title == "title"
    title = "".join(["plain", " title"])
    assert SDBlock.from_block_lines([title]).title is title
    assert SDBlock.from_block_lines(["\ud800 title "]).title == "\ud800 title"


def test_mdl_marker_callback_runs_after_yield() -> None:
    """Do not inspect the terminator until the caller resumes the generator."""
    calls: list[str] = []

    class Marker(str):
        def startswith(
            self,
            prefix: str | tuple[str, ...],
            start: SupportsIndex | None = 0,
            end: SupportsIndex | None = None,
        ) -> bool:
            calls.append(str(self))
            return super().startswith(prefix, start, end)

    lines = iter(SDBlock.parse_mdl([Marker("M  END")]))
    assert next(lines) == "M  END"
    assert calls == []
    assert list(lines) == []
    assert calls == ["M  END"]


def test_marker_stopiteration_retains_generator_error_cause() -> None:
    """Keep Python's generator error and original cause for marker failures."""

    class Marker(str):
        def startswith(
            self,
            prefix: str | tuple[str, ...],
            start: SupportsIndex | None = 0,
            end: SupportsIndex | None = None,
        ) -> bool:
            raise StopIteration("marker failed")

    with pytest.raises(
        RuntimeError, match="^generator raised StopIteration$"
    ) as caught:
        next(iter(SDBlock.parse_metadata([Marker("value")])))
    cause = caught.value.__cause__
    assert isinstance(cause, StopIteration)
    assert cause.args == ("marker failed",)
    assert cause is caught.value.__context__


@pytest.mark.parametrize(
    ("method", "values"),
    [
        ("parse_mdl", ["raw"]),
        ("parse_metadata", ["raw"]),
        ("from_lines", ["title", "$$$$"]),
    ],
)
def test_generator_source_cycles_are_collectable(
    method: str, values: list[str]
) -> None:
    """Collect retained sources when a public generator is abandoned.

    Args:
        method: public SDBlock generator method
        values: input lines supplied to the generator
    """

    class Source:
        def __init__(self) -> None:
            self.values = iter(values)
            self.generator: object | None = None

        def __iter__(self) -> Source:
            return self

        def __next__(self) -> str:
            return next(self.values)

    source = Source()
    generator = iter(getattr(SDBlock, method)(source))
    source.generator = generator
    source_ref = weakref.ref(source)
    next(generator)
    del source, generator
    gc.collect()
    assert source_ref() is None


def test_factory_keeps_input_block_alive_through_yield() -> None:
    """Retain the input block even when a subclass constructs an unrelated block."""
    references: list[weakref.ReferenceType[str]] = []

    class Line(str):
        pass

    def source() -> Iterator[str]:
        for value in ("title", "$$$$"):
            line = Line(value)
            references.append(weakref.ref(line))
            yield line

    class Factory(SDBlock):
        @classmethod
        def from_block_lines(cls, lines: Iterable[str]) -> SDBlock:
            return SDBlock("unrelated", deque(), deque())

    blocks = Factory.from_lines(source())
    assert next(blocks).title == "unrelated"
    assert references[0]() is not None
    assert list(blocks) == []
    gc.collect()
    assert references[0]() is None
