"""Exercise the private Rust framing algorithms directly."""

from __future__ import annotations

import gc
import importlib
import inspect
import weakref
from collections import deque
from typing import Iterable, Iterator

import pytest

from sendoff.sdblock import SDBlock

native = importlib.import_module("sendoff.native")


def test_native_framing_private_signatures() -> None:
    """Expose the approved parameter names on the real native callables."""
    signatures = {
        "mdl_iter": ("lines",),
        "metadata_iter": ("lines",),
        "blocks_iter": ("cls", "lines"),
    }
    for name, parameters in signatures.items():
        assert tuple(inspect.signature(getattr(native, name)).parameters) == parameters


def test_native_section_iterators_are_lazy_and_keep_delimiter_timing() -> None:
    """Delay input iteration and preserve deferred marker checks."""
    events: list[str] = []

    class Source:
        def __init__(self, lines: list[str]) -> None:
            self.lines = iter(lines)

        def __iter__(self) -> Source:
            events.append("iter")
            return self

        def __next__(self) -> str:
            line = next(self.lines)
            events.append(line)
            return line

    mdl = getattr(native, "mdl_iter")(Source(["raw", "M  END suffix", "tail"]))
    assert events == []
    assert iter(mdl) is mdl
    assert events == []
    assert next(mdl) == "raw"
    assert events == ["iter", "raw"]
    assert next(mdl) == "M  END suffix"
    assert events == ["iter", "raw", "M  END suffix"]
    with pytest.raises(StopIteration):
        next(mdl)
    assert events == ["iter", "raw", "M  END suffix"]

    metadata_source = Source([" value ", "$$$$ suffix", "tail"])
    metadata = getattr(native, "metadata_iter")(metadata_source)
    assert events == ["iter", "raw", "M  END suffix"]
    assert list(metadata) == [" value "]
    assert next(metadata_source) == "tail"
    assert list(getattr(native, "mdl_iter")([" M  END", "tail"])) == [
        " M  END",
        "tail",
    ]
    assert list(
        getattr(native, "metadata_iter")([" $$$$ leading", "$$$$ suffix", "tail"])
    ) == [" $$$$ leading"]


def test_native_block_title_uses_python_unicode_whitespace() -> None:
    """Strip exact text titles with Python's complete whitespace definition."""

    class Parser(SDBlock):
        @classmethod
        def parse_mdl(cls, lines: Iterable[str]) -> list[str]:
            return []

        @classmethod
        def parse_metadata(cls, lines: Iterable[str]) -> list[str]:
            return []

    block = Parser.from_block_lines(["\u001c\u00a0title\u2003\u001f"])
    assert type(block) is SDBlock
    assert block.title == "title"

    untouched_title = "".join(["plain", " title"])
    block = Parser.from_block_lines([untouched_title])
    assert block.title is untouched_title

    surrogate_title = f"{chr(0xD800)} title "
    block = Parser.from_block_lines([surrogate_title])
    assert block.title == "\ud800 title"


def test_native_mdl_yields_before_marker_callback_and_preserves_pep479() -> None:
    """Match generator StopIteration errors and delimiter-check timing."""

    class Marker:
        def __init__(self) -> None:
            self.calls = 0

        def startswith(self, prefix: str) -> bool:
            self.calls += 1
            return prefix == "M  END"

    marker = Marker()
    iterator = getattr(native, "mdl_iter")([marker])
    assert next(iterator) is marker
    assert marker.calls == 0
    with pytest.raises(StopIteration):
        next(iterator)
    assert marker.calls == 1

    class StopOnIter:
        def __iter__(self) -> Iterator[str]:
            raise StopIteration("iter failed")
            yield ""

    with pytest.raises(
        RuntimeError, match="^generator raised StopIteration$"
    ) as caught:
        next(getattr(native, "mdl_iter")(StopOnIter()))
    assert isinstance(caught.value.__cause__, StopIteration)
    assert caught.value.__cause__.args == ("iter failed",)
    assert caught.value.__cause__ is caught.value.__context__

    class StopOnStartswith:
        def startswith(self, prefix: str) -> bool:
            raise StopIteration("marker failed")

    with pytest.raises(
        RuntimeError, match="^generator raised StopIteration$"
    ) as caught:
        next(getattr(native, "metadata_iter")([StopOnStartswith()]))
    assert isinstance(caught.value.__cause__, StopIteration)
    assert caught.value.__cause__.args == ("marker failed",)
    assert caught.value.__cause__ is caught.value.__context__


def test_native_block_factory_uses_hooks_and_one_shared_iterator() -> None:
    """Use one shared source iterator for both dynamic class hooks."""
    events: list[str] = []

    class Hooks(SDBlock):
        @classmethod
        def parse_mdl(cls, lines: Iterable[str]) -> Iterator[str]:
            for line in lines:
                events.append(f"mdl:{line}")
                yield line
                if line.startswith("M  END"):
                    return

        @classmethod
        def parse_metadata(cls, lines: Iterable[str]) -> Iterator[str]:
            for line in lines:
                events.append(f"metadata:{line}")
                if line.startswith("$$$$"):
                    return
                yield line

    lines = iter(["  title  ", "raw", "M  END", "> <key>", "$$$$", "tail"])
    block = Hooks.from_block_lines(lines)
    assert type(block) is SDBlock
    assert block.title == "title"
    assert block.mdl == deque(["raw", "M  END"])
    assert block.metadata == deque(["> <key>"])
    assert events == ["mdl:raw", "mdl:M  END", "metadata:> <key>", "metadata:$$$$"]
    assert next(lines) == "tail"

    missing_end = Hooks.from_block_lines(["t", "raw", "> <key>", "v", "$$$$"])
    assert missing_end.mdl == deque(["raw", "> <key>", "v", "$$$$"])
    assert missing_end.metadata == deque()

    with pytest.raises(StopIteration) as caught:
        Hooks.from_block_lines([])
    assert caught.value.args == ()


def test_native_blocks_iterator_is_lazy_and_dispatches_factory() -> None:
    """Yield each completed block before reading more input."""
    consumed: list[str] = []
    calls: list[list[str]] = []

    def source() -> Iterator[str]:
        for line in ("title", "M  END", "$$$$ suffix", "tail"):
            consumed.append(line)
            yield line

    class Custom:
        @classmethod
        def from_block_lines(cls, lines: deque[str]) -> tuple[str, list[str]]:
            materialized = list(lines)
            calls.append(materialized)
            return ("custom", materialized)

    iterator = getattr(native, "blocks_iter")(Custom, source())
    assert consumed == []
    assert next(iterator) == ("custom", ["title", "M  END", "$$$$ suffix"])
    assert calls == [["title", "M  END", "$$$$ suffix"]]
    assert consumed == ["title", "M  END", "$$$$ suffix"]
    assert list(iterator) == []
    assert consumed == ["title", "M  END", "$$$$ suffix", "tail"]

    class Stopped:
        @classmethod
        def from_block_lines(cls, lines: deque[str]) -> None:
            raise StopIteration("factory stopped")

    with pytest.raises(
        RuntimeError, match="^generator raised StopIteration$"
    ) as caught:
        next(getattr(native, "blocks_iter")(Stopped, ["t", "$$$$"]))
    cause = caught.value.__cause__
    assert isinstance(cause, StopIteration)
    assert cause is caught.value.__context__
    assert cause.args == ("factory stopped",)


@pytest.mark.parametrize(
    ("factory_name", "values"),
    [
        ("mdl_iter", ["raw"]),
        ("metadata_iter", ["raw"]),
        ("blocks_iter", ["title", "$$$$"]),
    ],
)
def test_native_framing_iterators_trace_retained_python_references(
    factory_name: str, values: list[str]
) -> None:
    """Trace lazy iterator references through cyclic garbage collection.

    Args:
        factory_name: The private iterator factory to exercise.
        values: Lines consumed by the selected iterator.
    """

    class Source:
        def __init__(self, source_values: list[str]) -> None:
            self.values = iter(source_values)
            self.iterator: object | None = None

        def __iter__(self) -> Source:
            return self

        def __next__(self) -> str:
            return next(self.values)

    class Factory:
        @classmethod
        def from_block_lines(cls, lines: deque[str]) -> tuple[str, list[str]]:
            return "block", list(lines)

    source = Source(values)
    if factory_name == "blocks_iter":
        iterator = getattr(native, factory_name)(Factory, source)
    else:
        iterator = getattr(native, factory_name)(source)
    source.iterator = iterator
    source_ref = weakref.ref(source)
    next(iterator)
    del source, iterator
    gc.collect()
    assert source_ref() is None
