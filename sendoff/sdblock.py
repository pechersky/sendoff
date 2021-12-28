"""Separating an SDF into molecule "SDBlock"s of chemical data and metadata."""
from __future__ import annotations

import os
from collections import deque
from dataclasses import dataclass
from enum import Enum
from io import TextIOWrapper
from itertools import chain, groupby
from typing import Iterable, Iterator, Tuple, Union

Pathy = Union[str, bytes, os.PathLike]


class CTableFormat(Enum):
    """The format a CTAB comes in, either V2000 or V3000."""

    V2000 = "V2000"
    V3000 = "V3000"


class CTable:
    """Handle a Connection table (CTAB) from a single molecule block.

    Parsed in based on https://en.wikipedia.org/wiki/Chemical_table_file.
    As part of parsing, record the title, the source line, the comment line,
    and the counts line. The rest of the lines are kept but not parsed.
    The counts line is parsed to infer the format and the number of atoms.
    Note, the number of atoms in the counts line might not actually
    be the number of lines of atoms.
    """

    lines: deque[str]
    title: str
    source: str
    comment: str
    counts: str
    format: CTableFormat
    num_atoms: int

    def __init__(self, lines: Iterable[str]) -> None:
        """Parse in lines representing the CTAB, including the title line.

        Args:
            lines: an Iterable of str that comprises the CTAB
                The title line is parsed in by the SDBlock, and should
                be supplied as the first in lines.

        """
        self.lines = deque(lines)
        iterlines = iter(self.lines)
        self.title = next(iterlines).strip()
        self.source = next(iterlines)
        self.comment = next(iterlines)
        self.counts = next(iterlines).strip()
        self.format = self.parse_format(self.counts)
        if self.format is CTableFormat.V3000:
            next(iterlines)  # CTAB BEGIN line
            self.counts = next(iterlines).strip()
            self.num_atoms = self.parse_v3000_counts(self.counts)
        else:
            self.num_atoms = self.parse_v2000_counts(self.counts)
        return

    @staticmethod
    def parse_format(line: str) -> CTableFormat:
        """Parse a v2000 counts line according to get the format.

        The line could be from a V3000 block, where it is used for
        compatibility purposes. In that case, the actual
        counts line comes later.

        Can raise KeyError: when the CTableFormat could not be parsed
                from the counts line

        Args:
            line: a counts line parsed in as part of the CTAB block.
                Expected to be whitespace stripped.

        Returns:
            A CTAB format, based on the end of the counts line
        """
        sline = line.split()
        ctformat = CTableFormat[sline[-1]]
        return ctformat

    @staticmethod
    def parse_v2000_counts(line: str) -> int:
        """Parse a v2000 counts line according to get the number of atoms.

        The line could be from a V3000 block, where it is used for
        compatibility purposes. In that case, the actual
        counts line comes later, and this is not the right function to call.

        Note, V2000 maxes out at 999 atoms because the counts line
        can only handle 3 characters for the number of atoms.

        Can raise ValueError: when the number could not be parsed
                from the counts line

        Args:
            line: a counts line parsed in as part of the CTAB block.
                Expected to be whitespace stripped.

        Returns:
            The number of atoms indicated by the counts line.
                Not the actual number of atom lines further down.
        """
        num_atoms = int(line[:3])
        return num_atoms

    @staticmethod
    def parse_v3000_counts(line: str) -> int:
        """Parse a v3000 counts line according to get the number of atoms.

        Can raise ValueError: when the number could not be parsed
                from the counts line

        Args:
            line: a counts line parsed in as part of the CTAB block.
                Expected to be whitespace stripped.

        Returns:
            The number of atoms indicated by the counts line.
                Not the actual number of atom lines further down.
        """
        sline = line.split()
        num_atoms = int(sline[3])
        return num_atoms


@dataclass
class SDBlock:
    """Handle a single molecule block from an SD file."""

    title: str
    mdl: deque[str]
    metadata: deque[str]

    @classmethod
    def parse_mdl(cls, lines: Iterable[str]) -> Iterable[str]:
        """Parse in an MDL block from the ingested lines.

        Args:
            lines: an Iterable of str that comprises the MDL block

        Yields:
            str lines comprising the MDL block
        """
        for line in lines:
            yield line
            if line.startswith("M  END"):
                return

    @classmethod
    def parse_metadata(cls, lines: Iterable[str]) -> Iterable[str]:
        """Parse in the SDF metadata from the ingested lines.

        Args:
            lines: an Iterable of str that metadata, optionally
                including the $$$$ delimiter

        Yields:
            str lines comprising the metadata, excluding the $$$$ delimiter
        """
        for line in lines:
            if line.startswith("$$$$"):
                return
            yield line

    @classmethod
    def from_block_lines(cls, lines: Iterable[str]) -> SDBlock:
        """Parse in an SDBlock from a sequence of lines.

        Args:
            lines: an Iterable of str that
                optionally have whitespace that is stripped

        Returns:
            An SDBlock with a parsed in title and lines
        """
        iterlines = iter(lines)
        title = next(iterlines).strip()
        mdl = deque(cls.parse_mdl(iterlines))
        metadata = deque(cls.parse_metadata(iterlines))
        return SDBlock(title, mdl, metadata)

    @classmethod
    def from_lines(cls, lines: Iterable[str]) -> Iterator[SDBlock]:
        """Parse in SDBlocks from a sequence of lines.

        Args:
            lines: an Iterable of str that optionally
                have whitespace that is stripped,
                likely separated by the $$$$ SD delimiter

        Yields:
            SDBlocks parsed in from the lines
        """
        block: deque[str] = deque()
        for line in lines:
            block.append(line)
            if line.startswith("$$$$"):
                yield cls.from_block_lines(block)
                block = deque()

    def records(self) -> Iterable[Tuple[str, str]]:
        """Generate SD metadata records one by one.

        Yields:
            Tuples of str, str of record, value.
                The value includes any newlines if it is multiline
        """
        # make data chunks by breaking on empty lines
        for _, chunk in groupby((line.strip() for line in self.metadata), bool):
            record_line: str = next(chunk)
            if not record_line.startswith("> "):
                continue
            # something of the form `> ___<___>___` where `_` is anything
            record_name = (
                record_line.split("> ", 1)[1].strip().rsplit(">")[0].split("<", 1)[1]
            )
            yield record_name, str.join("\n", chunk)

    def write(self, outh: TextIOWrapper, with_newlines: bool = True) -> None:
        """Write an SDBlock to a file-like handle.

        Args:
            outh: file-like handle, such as what is returned by open(..., "w")
            with_newlines: each line is written with a trailing "\\n".
                The newline character is appended only if it wasn't in the line

        """
        print(self.title, file=outh)
        for line in chain(self.mdl, self.metadata):
            outh.write(line)
            if with_newlines and not line.endswith("\n"):
                outh.write("\n")
        outh.write("$$$$\n")

    def append_record(self, record_name: str, value: str) -> None:
        """Append a field and value to the SD data record.

        Args:
            record_name: str for field's name
            value: str for value, with no terminating newline

        """
        self.metadata.append(f"> <{record_name}>\n")
        self.metadata.append(value + "\n")
        self.metadata.append("\n")

    def num_atoms(self) -> int:
        """Get number of atoms as indicated in the MDL block.

        The mdl lines are parsed in to generate a CTable.
        Based on the format parsed, the number of atoms is parsed
        from the counts line. Note, this is not the number of atom lines.

        Returns:
            The number of atoms, parsed in from the counts line.
        """
        ctable = CTable(chain([self.title], self.mdl))
        return ctable.num_atoms


def parse_sdf(sdfpath: Pathy) -> Iterator[SDBlock]:
    """Parse in SDBlocks from a path. Wrapper around SDBlock.from_lines.

    Args:
        sdfpath: file-like handle, such as what is returned by open(..., "r")

    Returns:
        An Iterator of SDBlocks parsed in from the lines in the file
    """
    return SDBlock.from_lines(open(sdfpath).readlines())
