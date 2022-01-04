"""Reading and validating connection tables from SDFs."""
from __future__ import annotations

from collections import deque
from enum import Enum
from typing import Iterable, Tuple


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
    num_bonds: int

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
            self.num_atoms, self.num_bonds = self.parse_v3000_counts(self.counts)
        else:
            self.num_atoms, self.num_bonds = self.parse_v2000_counts(self.counts)
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
    def parse_v2000_counts(line: str) -> Tuple[int, int]:
        """Parse a v2000 counts line according to get the number of atoms and bonds.

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
            A tuple:
                (the number of atoms indicated by the counts line.
                    Not the actual number of atom lines further down.,
                the number of bonds indicated by the counts line.
                    Not the actual number of bonds lines further down.)
        """
        num_atoms = int(line[:3])
        num_bonds = int(line[3:6])
        return num_atoms, num_bonds

    @staticmethod
    def parse_v3000_counts(line: str) -> Tuple[int, int]:
        """Parse a v3000 counts line according to get the number of atoms and bonds.

        Can raise ValueError: when the number could not be parsed
                from the counts line

        Args:
            line: a counts line parsed in as part of the CTAB block.
                Expected to be whitespace stripped.

        Returns:
            A tuple:
                (the number of atoms indicated by the counts line.
                    Not the actual number of atom lines further down.,
                the number of bonds indicated by the counts line.
                    Not the actual number of bonds lines further down.)
        """
        sline = line.split()
        num_atoms = int(sline[3])
        num_bonds = int(sline[4])
        return num_atoms, num_bonds
