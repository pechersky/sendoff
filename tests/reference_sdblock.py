"""Test-only reference parser: main's pure-Python SDBlock, generic over str/bytes.

This mirrors ``sendoff/sdblock.py`` from the 0.1.x pure-Python release and is
the oracle for the differential tests. It encodes the accepted 1.0 behavior
changes that an unterminated final block is still yielded and that a block
without ``M  END`` does not keep its ``$$$$`` line as MDL. Divergences that are
accepted but not modeled here are kept out of the generated inputs instead.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from itertools import groupby
from typing import Generic, Iterable, Iterator, TypeVar

AnyStr = TypeVar("AnyStr", str, bytes)


def _lit(sample: AnyStr, text: str) -> AnyStr:
    """Return ``text`` as the same type as ``sample``.

    Args:
        sample: a str or bytes value whose type to match
        text: the ASCII literal to convert

    Returns:
        ``text`` itself for str samples, otherwise its encoding.
    """
    if isinstance(sample, str):
        return text
    return text.encode()


@dataclass
class RefBlock(Generic[AnyStr]):
    """A molecule block parsed by the reference implementation."""

    title: AnyStr
    mdl: list[AnyStr] = field(default_factory=list)
    metadata: list[AnyStr] = field(default_factory=list)

    @classmethod
    def from_block_lines(cls, lines: list[AnyStr]) -> RefBlock[AnyStr]:
        """Split block lines into a stripped title, MDL, and metadata.

        Args:
            lines: one block's lines, including any ``$$$$`` delimiter

        Returns:
            The parsed reference block.
        """
        iterlines = iter(lines)
        first = next(iterlines)
        mend, delimiter = _lit(first, "M  END"), _lit(first, "$$$$")
        mdl: list[AnyStr] = []
        for line in iterlines:
            if line.startswith(delimiter):
                break  # 1.0: main appended it, so write() doubled the delimiter
            mdl.append(line)
            if line.startswith(mend):
                break
        metadata: list[AnyStr] = []
        for line in iterlines:
            if line.startswith(delimiter):
                break
            metadata.append(line)
        return cls(first.strip(), mdl, metadata)

    def records(self) -> list[tuple[AnyStr, AnyStr]]:
        """Return (name, value) records exactly as main's generator yields them.

        Returns:
            Records in order, multiline values joined with newlines.
        """
        records: list[tuple[AnyStr, AnyStr]] = []
        for _, chunk in groupby((line.strip() for line in self.metadata), bool):
            record_line = next(chunk)
            if not record_line.startswith(_lit(record_line, "> ")):
                continue
            name = (
                record_line.split(_lit(record_line, "> "), 1)[1]
                .strip()
                .rsplit(_lit(record_line, ">"), 1)[0]
                .split(_lit(record_line, "<"), 1)[1]
            )
            records.append((name, _lit(record_line, "\n").join(chunk)))
        return records

    def append_record(self, name: AnyStr, value: AnyStr) -> None:
        """Append a record as main's append_record does.

        Args:
            name: record name
            value: record value without a trailing newline
        """
        newline = _lit(name, "\n")
        self.metadata.extend(
            [_lit(name, "> <") + name + _lit(name, ">\n"), value + newline, newline]
        )

    def molblock(self) -> AnyStr:
        """Return the title and MDL lines as main's idiom joins them.

        Returns:
            The stripped title, a newline, and the MDL lines, each
            terminated by a newline as ``write`` terminates them.
        """
        newline = _lit(self.title, "\n")
        return self.title[:0].join(
            [self.title, newline]
            + [line if line.endswith(newline) else line + newline for line in self.mdl]
        )

    def serialize(self) -> AnyStr:
        """Return what main's ``write(outh)`` writes, as one value.

        Returns:
            The serialized block.
        """
        newline = _lit(self.title, "\n")
        parts: list[AnyStr] = [self.title, newline]
        for line in self.mdl + self.metadata:
            parts.append(line)
            if not line.endswith(newline):
                parts.append(newline)
        parts.append(_lit(self.title, "$$$$\n"))
        return self.title[:0].join(parts)


def ref_from_lines(lines: Iterable[AnyStr]) -> Iterator[RefBlock[AnyStr]]:
    """Group lines into blocks at each ``$$$$`` line, keeping a trailing block.

    Args:
        lines: str or bytes lines

    Yields:
        Reference blocks in input order.
    """
    block: list[AnyStr] = []
    for line in lines:
        block.append(line)
        if line.startswith(_lit(line, "$$$$")):
            yield RefBlock.from_block_lines(block)
            block = []
    if block:
        yield RefBlock.from_block_lines(block)
