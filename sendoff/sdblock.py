"""Separating an SDF into molecule "SDBlock"s of chemical data and metadata."""

from __future__ import annotations

import os
from collections import deque
from dataclasses import dataclass
from itertools import chain
from typing import Iterable, Iterator, TextIO, Tuple, Union

import sendoff.native as native
from sendoff.ctable import CTable

Pathy = Union[str, bytes, "os.PathLike[str]"]


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
        for line in native.mdl_iter(lines):
            yield line

    @classmethod
    def parse_metadata(cls, lines: Iterable[str]) -> Iterable[str]:
        """Parse in the SDF metadata from the ingested lines.

        Args:
            lines: an Iterable of str that metadata, optionally
                including the $$$$ delimiter

        Yields:
            str lines comprising the metadata, excluding the $$$$ delimiter
        """
        for line in native.metadata_iter(lines):
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
        for record in native.records_iter(self):
            yield record

    def write(self, outh: TextIO, with_newlines: bool = True) -> None:
        """Write an SDBlock to a file-like handle.

        Args:
            outh: file-like handle, such as what is returned by open(..., "w")
            with_newlines: each line is written with a trailing "\\n".
                The newline character is appended only if it wasn't in the line

        """
        native.write(self, outh, with_newlines)

    def append_record(self, record_name: str, value: str) -> None:
        """Append a field and value to the SD data record.

        Args:
            record_name: str for field's name
            value: str for value, with no terminating newline

        """
        native.append_record(self, record_name, value)

    def ctable(self) -> CTable:
        """Parse out the underlying CTable object from the SDBlock.

        Returns:
            The CTable object, titled with the SDBlock.title,
                parsed in from the SDBlock.mdl, but without validation.
        """
        return CTable(chain([self.title], self.mdl))

    def num_atoms(self) -> int:
        """Get number of atoms as indicated in the MDL block.

        The mdl lines are parsed in to generate a CTable.
        Based on the format parsed, the number of atoms is parsed
        from the counts line. Note, this is not the number of atom lines.

        Returns:
            The number of atoms, parsed in from the counts line.
        """
        return self.ctable().num_atoms

    def num_bonds(self) -> int:
        """Get number of bonds as indicated in the MDL block.

        The mdl lines are parsed in to generate a CTable.
        Based on the format parsed, the number of bonds is parsed
        from the counts line. Note, this is not the number of bond lines.

        Returns:
            The number of bonds, parsed in from the counts line.
        """
        return self.ctable().num_bonds

    def renumber_indices(self) -> None:
        """Renumber atom and bond indices to be 1-indexed and in order.

        Utilizes the underlying CTable's methods, and overwrites
        this object's mdl.
        """
        ctable = self.ctable()
        ctable.renumber_indices()
        self.mdl = ctable.lines
        self.mdl.popleft()


def parse_sdf(sdfpath: Pathy) -> Iterator[SDBlock]:
    """Parse in SDBlocks from a path. Wrapper around SDBlock.from_lines.

    Args:
        sdfpath: file-like handle, such as what is returned by open(..., "r")

    Returns:
        An Iterator of SDBlocks parsed in from the lines in the file
    """
    return SDBlock.from_lines(open(sdfpath).readlines())
