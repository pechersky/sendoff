"""Exercise the private Rust framing algorithms directly."""

from __future__ import annotations

import gc
import importlib
import inspect
import io
import os
import weakref
from collections import deque
from pathlib import Path
from typing import Iterator

import pytest

_native = importlib.import_module("sendoff._native")


def test_native_framing_private_signatures() -> None:
    """Expose the approved parameter names on the real native callables."""
    signatures = {
        "_mdl_iter": ("lines",),
        "_metadata_iter": ("lines",),
        "_from_block_lines": ("cls", "block_type", "lines"),
        "_blocks_iter": ("cls", "lines"),
        "_read_sdf_lines": ("sdfpath",),
    }
    for name, parameters in signatures.items():
        assert tuple(inspect.signature(getattr(_native, name)).parameters) == parameters


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

    mdl = getattr(_native, "_mdl_iter")(Source(["raw", "M  END suffix", "tail"]))
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
    metadata = getattr(_native, "_metadata_iter")(metadata_source)
    assert events == ["iter", "raw", "M  END suffix"]
    assert list(metadata) == [" value "]
    assert next(metadata_source) == "tail"
    assert list(getattr(_native, "_mdl_iter")([" M  END", "tail"])) == [
        " M  END",
        "tail",
    ]
    assert list(
        getattr(_native, "_metadata_iter")([" $$$$ leading", "$$$$ suffix", "tail"])
    ) == [" $$$$ leading"]


def test_native_mdl_yields_before_marker_callback_and_preserves_pep479() -> None:
    """Match generator StopIteration errors and delimiter-check timing."""

    class Marker:
        def __init__(self) -> None:
            self.calls = 0

        def startswith(self, prefix: str) -> bool:
            self.calls += 1
            return prefix == "M  END"

    marker = Marker()
    iterator = getattr(_native, "_mdl_iter")([marker])
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
        next(getattr(_native, "_mdl_iter")(StopOnIter()))
    assert isinstance(caught.value.__cause__, StopIteration)
    assert caught.value.__cause__.args == ("iter failed",)
    assert caught.value.__cause__ is caught.value.__context__

    class StopOnStartswith:
        def startswith(self, prefix: str) -> bool:
            raise StopIteration("marker failed")

    with pytest.raises(
        RuntimeError, match="^generator raised StopIteration$"
    ) as caught:
        next(getattr(_native, "_metadata_iter")([StopOnStartswith()]))
    assert isinstance(caught.value.__cause__, StopIteration)
    assert caught.value.__cause__.args == ("marker failed",)
    assert caught.value.__cause__ is caught.value.__context__


def test_native_block_factory_uses_hooks_and_one_shared_iterator() -> None:
    """Use one shared source iterator for both dynamic class hooks."""
    events: list[str] = []

    class BaseBlock:
        def __init__(self, title: str, mdl: deque[str], metadata: deque[str]) -> None:
            self.title = title
            self.mdl = mdl
            self.metadata = metadata

    class Hooks:
        @classmethod
        def parse_mdl(cls, lines: Iterator[str]) -> Iterator[str]:
            for line in lines:
                events.append(f"mdl:{line}")
                yield line
                if line.startswith("M  END"):
                    return

        @classmethod
        def parse_metadata(cls, lines: Iterator[str]) -> Iterator[str]:
            for line in lines:
                events.append(f"metadata:{line}")
                if line.startswith("$$$$"):
                    return
                yield line

    lines = iter(["  title  ", "raw", "M  END", "> <key>", "$$$$", "tail"])
    block = getattr(_native, "_from_block_lines")(Hooks, BaseBlock, lines)
    assert type(block) is BaseBlock
    assert block.title == "title"
    assert block.mdl == deque(["raw", "M  END"])
    assert block.metadata == deque(["> <key>"])
    assert events == ["mdl:raw", "mdl:M  END", "metadata:> <key>", "metadata:$$$$"]
    assert next(lines) == "tail"

    missing_end = getattr(_native, "_from_block_lines")(
        Hooks, BaseBlock, ["t", "raw", "> <key>", "v", "$$$$"]
    )
    assert missing_end.mdl == deque(["raw", "> <key>", "v", "$$$$"])
    assert missing_end.metadata == deque()

    with pytest.raises(StopIteration) as caught:
        getattr(_native, "_from_block_lines")(Hooks, BaseBlock, [])
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

    iterator = getattr(_native, "_blocks_iter")(Custom, source())
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
        next(getattr(_native, "_blocks_iter")(Stopped, ["t", "$$$$"]))
    cause = caught.value.__cause__
    assert isinstance(cause, StopIteration)
    assert cause is caught.value.__context__
    assert cause.args == ("factory stopped",)


@pytest.mark.parametrize(
    ("factory_name", "values"),
    [
        ("_mdl_iter", ["raw"]),
        ("_metadata_iter", ["raw"]),
        ("_blocks_iter", ["title", "$$$$"]),
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
    if factory_name == "_blocks_iter":
        iterator = getattr(_native, factory_name)(Factory, source)
    else:
        iterator = getattr(_native, factory_name)(source)
    source.iterator = iterator
    source_ref = weakref.ref(source)
    next(iterator)
    del source, iterator
    gc.collect()
    assert source_ref() is None


def test_native_read_sdf_lines_uses_default_open_and_rejects_handles() -> None:
    """Delegate path handling and text decoding to builtins.open."""
    read_lines = getattr(_native, "_read_sdf_lines")
    expected = read_lines("LICENSE")
    assert read_lines(Path("LICENSE")) == expected
    assert read_lines(os.fsencode("LICENSE")) == expected
    with pytest.raises(TypeError, match="expected str, bytes or os.PathLike object"):
        read_lines(io.StringIO("contents"))
