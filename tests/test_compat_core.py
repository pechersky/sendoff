"""Guard the public Python contract while separating the Rust core."""

from collections import deque
from typing import Iterator, SupportsIndex

import pytest

from sendoff.ctable import CTable
from sendoff.sdblock import SDBlock
from tests.compat_literals import V3000


@pytest.mark.parametrize("value", ["> malformed value", "> <another header>"])
def test_record_value_is_not_a_group_header(value: str) -> None:
    """Treat header-looking text inside the value group as ordinary data.

    Args:
        value: text that must remain a value rather than start a record
    """
    block = SDBlock("", deque(), deque(["> <key>", value, ""]))
    assert list(block.records()) == [("key", value)]


def test_record_values_do_not_call_header_protocol() -> None:
    """Only the first line in a group is inspected as a record header."""

    class Value(str):
        def strip(self, chars: str | None = None) -> str:
            return self

        def startswith(
            self,
            prefix: str | tuple[str, ...],
            start: SupportsIndex | None = 0,
            end: SupportsIndex | None = None,
        ) -> bool:
            raise AssertionError("a value is not a group header")

    block = SDBlock("", deque(), deque(["> <key>", Value("value"), ""]))
    assert list(block.records()) == [("key", "value")]


def test_metadata_checks_each_marker_once() -> None:
    """Do not repeat the previous line's marker protocol on the next resume."""
    calls: list[str] = []
    iterator: Iterator[str]

    class Line(str):
        def startswith(
            self,
            prefix: str | tuple[str, ...],
            start: SupportsIndex | None = 0,
            end: SupportsIndex | None = None,
        ) -> bool:
            calls.append(str(self))
            if self == "value":
                with pytest.raises(ValueError, match="^generator already executing$"):
                    next(iterator)
            return super().startswith(prefix, start, end)

    iterator = iter(SDBlock.parse_metadata([Line("value"), Line("$$$$")]))
    assert list(iterator) == ["value"]
    assert calls == ["value", "$$$$"]


@pytest.mark.parametrize("position", [0, 1, 2])
def test_record_join_keeps_all_values_with_surrogates(position: int) -> None:
    """Retain preceding and following values when Python must join the text.

    Args:
        position: Location of the non-UTF-8 Python string in the value group.
    """
    values = ["first", "second"]
    values.insert(position, "\ud800")
    block = SDBlock("", deque(), deque(["> <key>", *values, ""]))
    assert list(block.records()) == [("key", "\n".join(values))]


def test_renumber_preserves_unchanged_line_objects() -> None:
    """Keep existing line objects outside replaced counts/atom/bond sections."""
    table = CTable(V3000.splitlines())
    unchanged = table.lines[0], table.lines[1], table.lines[2], table.lines[-1]
    table.renumber_indices()
    restored = table.lines[0], table.lines[1], table.lines[2], table.lines[-1]
    assert all(actual is original for actual, original in zip(restored, unchanged))
