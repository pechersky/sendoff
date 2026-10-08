"""Test-only reference parser: main's pure-Python SDBlock, generic over str/bytes.

This mirrors ``sendoff/sdblock.py`` from the 0.1.x pure-Python release and is
the oracle for the differential tests. It encodes the accepted native-parser
behavior changes that an unterminated final block is still yielded, that a
block without ``M  END`` does not keep its ``$$$$`` line as MDL, and that line
endings survive: a ``$$$$`` line that is exactly ``$$$$`` and a newline is kept,
and lines the block writes or ends itself end like its title line (a line
ending in a bare ``\r`` gets just ``\n``). An appended record follows a blank
line, so that it does not run on from a last record without one, and its
header copies the first named header's line around the name. A block whose
first line starts with ``$$$$`` is an empty record, as the CTfile spec reads
any such line as a terminator, rather than a ``$$$$`` title that each write
doubled. Divergences that are accepted but not modeled here are kept out of
the generated inputs instead.
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
    newline: AnyStr
    trailer: AnyStr
    mdl: list[AnyStr] = field(default_factory=list)
    metadata: list[AnyStr] = field(default_factory=list)
    titled: bool = True

    @classmethod
    def from_block_lines(cls, lines: list[AnyStr]) -> RefBlock[AnyStr]:
        """Split block lines into a stripped title, MDL, and metadata.

        Args:
            lines: one block's lines, including any ``$$$$`` delimiter

        Returns:
            The parsed reference block.
        """
        lf, crlf = _lit(lines[0], "\n"), _lit(lines[0], "\r\n")
        title = lines[0] if lines[0].endswith(lf) else lines[0] + lf
        newline = crlf if title.endswith(crlf) else lf
        lines = [
            (
                line
                if line.endswith(lf)
                else line + (lf if line.endswith(crlf[:1]) else newline)
            )
            for line in [title] + lines[1:]
        ]
        iterlines = iter(lines)
        first = next(iterlines)
        mend, delimiter = _lit(first, "M  END"), _lit(first, "$$$$")
        trailer = delimiter + newline
        if first.startswith(delimiter):
            return cls(first[:0], newline, trailer, titled=False)
        if len(lines) > 1 and lines[-1] in (delimiter + lf, delimiter + crlf):
            trailer = lines[-1]
        mdl: list[AnyStr] = []
        for line in iterlines:
            if line.startswith(delimiter):
                break  # main appended it, so write() doubled the delimiter
            mdl.append(line)
            if line.startswith(mend):
                break
        metadata: list[AnyStr] = []
        for line in iterlines:
            if line.startswith(delimiter):
                break
            metadata.append(line)
        return cls(first.strip(), newline, trailer, mdl, metadata)

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

        Unlike main, a blank line comes first if the SD data ends in a line
        that is not blank, and the header copies the first named header.

        Args:
            name: record name
            value: record value without a trailing newline
        """
        newline = self.newline
        self.titled = True
        if self.metadata and self.metadata[-1].strip():
            self.metadata.append(newline)
        before, after = self._template() or (_lit(name, "> <"), _lit(name, ">"))
        self.metadata.extend(
            [before + name + after + newline, value + newline, newline]
        )

    def _template(self) -> tuple[AnyStr, AnyStr] | None:
        """Return the first named header line before and after its name.

        Returns:
            The text through ``<`` and from the closing ``>``, without the
            line ending, or None without a named header or a closing ``>``.
        """
        for _, chunk in groupby(self.metadata, lambda line: bool(line.strip())):
            line = next(chunk)
            field = line.strip()
            if not field.startswith(_lit(line, "> ")):
                continue
            if _lit(line, "<") not in field[2:].rsplit(_lit(line, ">"), 1)[0]:
                continue  # no name
            if _lit(line, ">") not in field[2:]:
                return None
            start, end = line.find(_lit(line, "<")), line.rfind(_lit(line, ">"))
            text = line.removesuffix(_lit(line, "\n")).removesuffix(_lit(line, "\r"))
            return line[: start + 1], text[end:]
        return None

    def molblock(self) -> AnyStr:
        """Return the title and MDL lines as main's idiom joins them.

        Returns:
            The stripped title, the title line's ending, and the MDL lines,
            or nothing for an empty record.
        """
        if not self.titled:
            return self.title
        return self.title[:0].join([self.title, self.newline] + self.mdl)

    def serialize(self) -> AnyStr:
        """Return what main's ``write(outh)`` writes, as one value.

        Returns:
            The serialized block.
        """
        if not self.titled:
            return self.trailer
        return self.title[:0].join(
            [self.title, self.newline] + self.mdl + self.metadata + [self.trailer]
        )


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
