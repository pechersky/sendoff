"""Exercise the private Rust metadata and writing algorithms directly."""

from __future__ import annotations

import gc
import importlib
import inspect
import weakref
from collections import deque
from typing import Any, Iterator, SupportsIndex
from unittest.mock import MagicMock

import pytest

from sendoff.sdblock import SDBlock

native = importlib.import_module("sendoff.native")


def test_native_sddata_private_signatures() -> None:
    """Expose the approved parameter names on the real native callables."""
    signatures = {
        "records_iter": ("block",),
        "write": ("block", "outh", "with_newlines"),
        "append_record": ("block", "record_name", "value"),
    }
    for name, parameters in signatures.items():
        assert tuple(inspect.signature(getattr(native, name)).parameters) == parameters


def test_native_records_are_lazy_and_preserve_grouping_and_header_rules() -> None:
    """Keep lazy metadata lookup, grouping, value joining, and header parsing."""
    events: list[str] = []

    class Block:
        current: deque[str]

        @property
        def metadata(self) -> deque[str]:
            events.append("metadata")
            return self.current

        @metadata.setter
        def metadata(self, value: deque[str]) -> None:
            self.current = value

    block = Block()
    iterator = getattr(native, "records_iter")(block)
    block.metadata = deque(
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
    )
    assert events == []
    assert list(iterator) == [
        (" key ", "first\nsecond"),
        ("key", ""),
        ("key", "last\n> <not-a-new-record>\nstill value"),
    ]
    assert events == ["metadata"]


def test_native_records_keep_live_deque_errors_and_lazy_parse_failures() -> None:
    """Preserve deferred malformed-header and live-deque failures."""

    class Block:
        def __init__(self, metadata: deque[Any]) -> None:
            self.metadata = metadata

    malformed = getattr(native, "records_iter")(
        Block(deque(["> missing brackets", "value"]))
    )
    with pytest.raises(IndexError, match="^list index out of range$"):
        next(malformed)
    with pytest.raises(StopIteration):
        next(malformed)

    metadata = deque(["> <one>", "value", "", "> <two>", "value", ""])
    records = getattr(native, "records_iter")(Block(metadata))
    assert next(records) == ("one", "value")
    metadata.append("changed")
    with pytest.raises(RuntimeError, match="^deque mutated during iteration$"):
        next(records)

    class StopOnStrip:
        def strip(self) -> str:
            raise StopIteration("strip failed")

    records = getattr(native, "records_iter")(Block(deque([StopOnStrip()])))
    with pytest.raises(
        RuntimeError, match="^generator raised StopIteration$"
    ) as caught:
        next(records)
    assert isinstance(caught.value.__cause__, StopIteration)
    assert caught.value.__cause__.args == ("strip failed",)
    assert caught.value.__cause__ is caught.value.__context__


def test_native_records_retain_block_through_final_yield() -> None:
    """Keep the source alive until the iterator is resumed to exhaustion."""

    class Block:
        def __init__(self) -> None:
            self.metadata = deque(["> <key>", "value"])

    block = Block()
    block_ref = weakref.ref(block)
    records = getattr(native, "records_iter")(block)
    del block

    assert next(records) == ("key", "value")
    assert block_ref() is not None
    with pytest.raises(StopIteration):
        next(records)
    assert block_ref() is None


def test_native_records_cleanup_rejects_reentrant_finalizer_safely() -> None:
    """Clear references on errors without panicking on finalizer reentry."""
    records: Any = None
    reentrant_errors: list[str] = []

    class StopOnStrip:
        def strip(self) -> str:
            raise ValueError("strip failed")

    class Block:
        def __init__(self) -> None:
            self.metadata = deque([StopOnStrip()])

        def __del__(self) -> None:
            try:
                next(records)
            except ValueError as error:
                reentrant_errors.append(str(error))

    block = Block()
    records = getattr(native, "records_iter")(block)
    del block
    with pytest.raises(ValueError, match="^strip failed$"):
        next(records)
    assert reentrant_errors == ["generator already executing"]
    with pytest.raises(StopIteration):
        next(records)


@pytest.mark.parametrize(
    ("case", "expected"),
    [
        ("different-key", [("one", "v"), ("two", "w")]),
        ("same-key", [("one", "v")]),
        ("initial-key-stop", []),
    ],
)
def test_native_records_resume_after_inner_key_stop(
    case: str, expected: list[tuple[str, str]]
) -> None:
    """Match groupby when bool raises StopIteration inside a record chunk.

    Args:
        case: Which grouping transition to exercise.
        expected: Literal records produced by the Python grouping behavior.
    """

    class StopBool(str):
        def strip(self, chars: str | None = None) -> str:
            return self

        def __bool__(self) -> bool:
            raise StopIteration("key stopped")

    metadata_by_case = {
        "different-key": [
            "> <one>",
            "v",
            StopBool("stop"),
            "",
            "> <two>",
            "w",
        ],
        "same-key": ["> <one>", "v", StopBool("stop"), "> <two>", "w"],
        "initial-key-stop": [StopBool("stop"), "> <one>", "w"],
    }
    records = getattr(native, "records_iter")(
        SDBlock("title", deque(), deque(metadata_by_case[case]))
    )
    assert list(records) == expected


def test_native_records_parse_exact_text_with_python_whitespace() -> None:
    """Trim and split normal Unicode metadata with Python-compatible whitespace."""
    block = SDBlock(
        "title",
        deque(),
        deque(
            [
                "\u001c",
                "\u00a0> prefix < name\u00a0> suffix\u2003",
                "\u2003first\u001c",
                "\u00a0second\u00a0",
                "\u001f",
            ]
        ),
    )
    assert list(getattr(native, "records_iter")(block)) == [
        (" name\u00a0", "first\nsecond")
    ]

    untouched_value = "".join(["plain", " value"])
    record = next(
        getattr(native, "records_iter")(
            SDBlock("title", deque(), deque(["> <key>", untouched_value]))
        )
    )
    assert record == ("key", untouched_value)
    assert record[1] is untouched_value


def test_native_records_only_parse_headers_at_group_heads() -> None:
    """Treat header-looking values as values without calling their protocols."""

    class Value(str):
        def strip(self, chars: str | None = None) -> str:
            return self

        def startswith(
            self,
            prefix: str | tuple[str, ...],
            start: SupportsIndex | None = None,
            end: SupportsIndex | None = None,
        ) -> bool:
            raise AssertionError("value startswith was called")

    block = SDBlock(
        "",
        deque(),
        deque(["> <key>\n", Value("> malformed value\n"), "\n"]),
    )
    assert list(block.records()) == [("key", "> malformed value\n")]

    block = SDBlock(
        "",
        deque(),
        deque(["> <key>", Value("value"), ""]),
    )
    assert list(block.records()) == [("key", "value")]


def test_native_records_keep_lone_surrogate_text() -> None:
    """Retain Python's handling for strings Rust cannot encode as UTF-8."""
    surrogate = chr(0xD800)
    record = next(
        getattr(native, "records_iter")(
            SDBlock("title", deque(), deque([f"> <{surrogate}>", surrogate]))
        )
    )
    assert record == ("\ud800", "\ud800")
    assert record[1] is surrogate


def test_native_records_preserve_falsey_string_groups() -> None:
    """Match Python grouping for falsey headers and both key transitions."""

    class FalseText(str):
        def strip(self, chars: str | None = None) -> str:
            return self

        def __bool__(self) -> bool:
            return False

    block = SDBlock(
        "title",
        deque(),
        deque(
            [
                FalseText("> <first>"),
                FalseText("one"),
                FalseText("two"),
                "> <normal>",
                "value",
                FalseText("> <last>"),
                FalseText("tail"),
            ]
        ),
    )
    expected = [("first", "one\ntwo"), ("normal", "value"), ("last", "tail")]
    assert list(block.records()) == expected
    assert list(getattr(native, "records_iter")(block)) == expected


def test_native_records_trace_the_block_and_live_metadata_iterator() -> None:
    """Trace retained block and metadata iterator references through GC."""

    class Block:
        def __init__(self) -> None:
            self.metadata: Metadata
            self.records: object | None = None

    class Metadata:
        def __init__(self, owner: Block) -> None:
            self.owner = owner
            self.values = iter(["> <key>"])

        def __iter__(self) -> Metadata:
            return self

        def __next__(self) -> str:
            return next(self.values)

    block = Block()
    metadata = Metadata(block)
    block.metadata = metadata
    records = getattr(native, "records_iter")(block)
    block.records = records
    block_ref = weakref.ref(block)
    assert next(records) == ("key", "")
    del block, metadata, records
    gc.collect()
    assert block_ref() is None


def test_native_write_preserves_call_order_granularity_and_live_truthiness() -> None:
    """Preserve print ordering, write granularity, and per-line truthiness."""
    events: list[str] = []
    raw = MagicMock(name="raw")
    raw.endswith.return_value = False
    already = MagicMock(name="already")
    already.endswith.return_value = True
    key = MagicMock(name="key")
    key.endswith.return_value = False
    value = MagicMock(name="value")
    value.endswith.return_value = False

    def mdl_source() -> Iterator[Any]:
        events.append("mdl:iter")
        yield raw
        yield already
        events.append("mdl:end")

    def metadata_source() -> Iterator[Any]:
        events.append("metadata:iter")
        yield key
        yield value

    block = MagicMock()
    block.title = "title"
    block.mdl = ["old mdl"]
    block.metadata = ["old metadata"]
    writer = MagicMock()

    def write(text: Any) -> int:
        if writer.write.call_count == 2:
            events.append("title:printed")
            block.mdl = mdl_source()
            block.metadata = metadata_source()
        return 1

    writer.write.side_effect = write
    with_newlines = MagicMock()
    with_newlines.__bool__.side_effect = [True, False, True, False]
    getattr(native, "write")(block, writer, with_newlines)
    assert [call.args[0] for call in writer.write.call_args_list] == [
        "title",
        "\n",
        raw,
        "\n",
        already,
        key,
        "\n",
        value,
        "$$$$\n",
    ]
    assert events.index("title:printed") < events.index("mdl:iter")
    assert events.index("mdl:end") < events.index("metadata:iter")
    assert with_newlines.__bool__.call_count == 4
    raw.endswith.assert_called_once_with("\n")
    already.endswith.assert_not_called()
    key.endswith.assert_called_once_with("\n")
    value.endswith.assert_not_called()


def test_native_write_propagates_partial_writer_failures() -> None:
    """Keep partial output and stop before writing later lines after an error."""
    events: list[str] = []

    class Writer:
        def write(self, text: str) -> int:
            events.append(text)
            if text == "raw":
                raise OSError("sink failed")
            return 1

    class Block:
        @property
        def title(self) -> str:
            return "title"

        @property
        def mdl(self) -> list[str]:
            events.append("mdl")
            return ["raw"]

        @property
        def metadata(self) -> list[str]:
            events.append("metadata")
            return ["metadata"]

    with pytest.raises(OSError, match="^sink failed$"):
        getattr(native, "write")(Block(), Writer(), True)
    assert events == ["title", "\n", "mdl", "metadata", "raw"]


def test_native_write_handles_plain_string_newlines() -> None:
    """Use native text checks while retaining print and per-line writes."""
    writes: list[str] = []

    class Writer:
        def write(self, text: str) -> int:
            writes.append(text)
            return len(text)

    block = SDBlock(
        "title",
        deque(["raw", "already\n"]),
        deque(["metadata"]),
    )
    getattr(native, "write")(block, Writer(), True)
    assert writes == [
        "title",
        "\n",
        "raw",
        "\n",
        "already\n",
        "metadata",
        "\n",
        "$$$$\n",
    ]


def test_native_append_keeps_formatting_and_three_append_side_effects() -> None:
    """Preserve Python formatting and the three separate append calls."""
    events: list[tuple[int, Any]] = []

    class Sink:
        def __init__(self, index: int) -> None:
            self.index = index

        def append(self, value: Any) -> None:
            events.append((self.index, value))

    class Block:
        index = 0

        @property
        def metadata(self) -> Sink:
            self.index += 1
            return Sink(self.index)

    class Name:
        def __format__(self, spec: str) -> str:
            assert spec == ""
            return "formatted"

    getattr(native, "append_record")(Block(), Name(), "value\ninner")
    assert events == [
        (1, "> <formatted>\n"),
        (2, "value\ninner\n"),
        (3, "\n"),
    ]

    metadata: deque[str] = deque()

    class StoredBlock:
        metadata: deque[str]

    stored = StoredBlock()
    stored.metadata = metadata
    with pytest.raises(TypeError):
        getattr(native, "append_record")(stored, "name", object())
    assert metadata == deque(["> <name>\n"])

    stored.metadata = deque()
    getattr(native, "append_record")(stored, "naïve", "first\nsecond")
    assert stored.metadata == deque(["> <naïve>\n", "first\nsecond\n", "\n"])

    surrogate = chr(0xD800)
    stored.metadata = deque()
    getattr(native, "append_record")(stored, surrogate, surrogate)
    assert stored.metadata == deque(["> <\ud800>\n", "\ud800\n", "\n"])
