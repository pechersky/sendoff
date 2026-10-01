"""Reading and validating connection tables from SDFs."""

from __future__ import annotations

import sys
from collections import deque
from enum import Enum
from typing import Iterable, Tuple

import sendoff.native as native


class CTableFormat(Enum):
    """The format a CTAB comes in, either V2000 or V3000."""

    V2000 = "V2000"
    V3000 = "V3000"


class IndicesMismatchError(Exception):
    """When then number of atom or bond lines does not match the counts line."""


class IndicesOutOfOrderError(Exception):
    """When then atom or bond lines are not in increasing index order."""


class IndicesDuplicateError(Exception):
    """When then atom or bond lines have duplicate indices."""


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
        native.ctable_init(self, lines, CTableFormat.V3000)

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
        return native.parse_format(line, CTableFormat)

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
        return native.parse_v2000_counts(line)

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
        return native.parse_v3000_counts(line)

    def atomlines(self) -> Iterable[str]:
        """Get atom lines in the atom table, assumed to be after 7 lines.

        Returns:
            A single-use iterable of the atom lines
        """
        return native.atomlines(self)

    def bondlines(self) -> Iterable[str]:
        """Get bond lines in the bond table.

        We need to get the atomlines first.
        TODO: Make the usage friendlier to iteration, so that
        other methods don't end up calling bondlines twice.

        Returns:
            A single-use iterable of the bond lines
        """
        return native.bondlines(self)

    def valid_atom_indices(self, strict: bool = False) -> bool:
        """Validate that the atom lines match the counts line.

        If strict, make sure they are 1-indexed and in order.
        This can break if the V3000 block has "-" terminated lines,
            which means that the next line is a continuation of the previous.

        Args:
            strict: the indices start with 1, and increment by one.

        Raises:
            IndicesDuplicateError: if an atom line met has an index seen before
            IndicesMismatchError: if number of atom lines does not match count line
            IndicesOutOfOrderError: if strict, and indices are not in 1-indexed order
            NotImplementedError: if trying to validate a V2000 format table

        Returns:
            If all the checks pass, return True.
        """
        return native.valid_atom_indices(
            self, strict, CTableFormat.V3000, sys.modules[__name__]
        )

    def valid_bond_indices(self, strict: bool = False) -> bool:
        """Validate that the bond lines match the counts line.

        If strict, make sure they are 1-indexed and in order.
        This can break if the V3000 block has "-" terminated lines,
            which means that the next line is a continuation of the previous.

        Args:
            strict: the indices start with 1, and increment by one.

        Raises:
            IndicesDuplicateError: if an bond line met has an index seen before
            IndicesMismatchError: if number of bond lines does not match count line
            IndicesOutOfOrderError: if strict, and indices are not in 1-indexed order
            NotImplementedError: if trying to validate a V2000 format table

        Returns:
            If all the checks pass, return True.
        """
        return native.valid_bond_indices(
            self, strict, CTableFormat.V3000, sys.modules[__name__]
        )

    def renumber_indices(self) -> None:
        """Renumber the indices in the block, including changing counts line.

        Iterating through the atom lines, replace the indexes into a
        1-indexed order, keeping track of what maps to what.
        Then, renumber the bond lines, making sure to map the new atom indices
        properly.
        At the end, regenerate the counts line to match the number of atom and bond
        lines.
        This will run regardless of whether the indices are already valid.
        It cannot fix duplicate atom indices properly if the duplicate one
        occurs in the bond lines, and will raise an error.
        The lines are then in-place replaced within the CTable object.
        Only implemented for V3000 tables.

        This can break if the V3000 block has "-" terminated lines,
            which means that the next line is a continuation of the previous.

        Raises:
            IndicesDuplicateError: if there was a duplicate atom index, and it was used
                somewhere in a bond line. Raised, because it is not clear which of the
                original indices to use in the remapping.
            NotImplementedError: if trying to renumber in a V2000 format table
        """
        native.renumber_ctable(self, CTableFormat.V3000, IndicesDuplicateError)
