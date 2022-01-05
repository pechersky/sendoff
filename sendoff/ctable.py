"""Reading and validating connection tables from SDFs."""
from __future__ import annotations

import itertools as itt
from collections import defaultdict, deque
from enum import Enum
from typing import Iterable, Tuple


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
        self.lines = deque(lines)
        iterlines = iter(self.lines)
        self.title = next(iterlines).strip()
        self.source = next(iterlines)
        self.comment = next(iterlines)
        self.counts = next(iterlines)
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
        sline = line.strip().split()
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

    def atomlines(self) -> Iterable[str]:
        """Get atom lines in the atom table, assumed to be after 7 lines.

        Returns:
            A single-use iterable of the atom lines
        """
        atomlines = itt.takewhile(
            lambda x: not str.startswith(x, "M  V30 END ATOM"),
            # title source comment compat begin counts begin
            itt.islice(self.lines, 7, None),
        )
        return atomlines

    def bondlines(self) -> Iterable[str]:
        """Get bond lines in the bond table.

        We need to get the atomlines first.
        TODO: Make the usage friendlier to iteration, so that
        other methods don't end up calling bondlines twice.

        Returns:
            A single-use iterable of the bond lines
        """
        bondlines = itt.takewhile(
            lambda x: not str.startswith(x, "M  V30 END BOND"),
            # islice(..., 1, None) means to drop one
            itt.islice(
                itt.dropwhile(
                    lambda x: not str.startswith(x, "M  V30 BEGIN BOND"), self.lines
                ),
                1,
                None,
            ),
        )
        return bondlines

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
        if self.format is not CTableFormat.V3000:
            raise NotImplementedError
        seen_indices: set[int] = set()
        for line_ix, line in enumerate(self.atomlines()):
            sline = line.split()
            atom_ix = int(sline[2])
            if strict and line_ix + 1 != atom_ix:
                raise IndicesOutOfOrderError("atoms")
            if atom_ix in seen_indices:
                raise IndicesDuplicateError("atoms")
            seen_indices.add(atom_ix)
        if len(seen_indices) < self.num_atoms:
            raise IndicesMismatchError("fewer atom lines than count line")
        if len(seen_indices) > self.num_atoms:
            raise IndicesMismatchError("more atom lines than count line")
        return True

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
        if self.format is not CTableFormat.V3000:
            raise NotImplementedError
        seen_indices: set[int] = set()
        for line_ix, line in enumerate(self.bondlines()):
            sline = line.split()
            bond_ix = int(sline[2])
            if strict and line_ix + 1 != bond_ix:
                raise IndicesOutOfOrderError("bonds")
            if bond_ix in seen_indices:
                raise IndicesDuplicateError("bonds")
            seen_indices.add(bond_ix)
        if len(seen_indices) < self.num_bonds:
            raise IndicesMismatchError("fewer bond lines than count line")
        if len(seen_indices) > self.num_bonds:
            raise IndicesMismatchError("more bond lines than count line")
        return True

    def renumber_indices(self) -> None:  # noqa: max-complexity: 13
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
        if self.format is not CTableFormat.V3000:
            raise NotImplementedError
        new_atomlines: deque[str] = deque()
        new_bondlines: deque[str] = deque()
        # old_ix: [new_ix, new_ix2, ...]
        atom_index_mapping: defaultdict[int, list[int]] = defaultdict(list)
        prefix = "M  V30 "
        len_prefix = len(prefix)
        for line_ix, atomline in enumerate(self.atomlines()):
            sline = atomline.split()
            atom_ix = int(sline[2])
            new_ix = line_ix + 1
            atom_index_mapping[atom_ix].append(new_ix)
            # can't use str.join because we want to preserve all whitespace properly
            # but we do assume that `M  V30 ` is correctly done by the spec
            skip_chars = len_prefix + len(sline[2])
            nline = f"{prefix}{new_ix}" + atomline[skip_chars:]
            new_atomlines.append(nline)
        for line_ix, bondline in enumerate(self.bondlines()):
            # here, we take less care to retain the whitespace and just
            # reconstruct the line with single whitespace
            sline = bondline.split()
            new_ix = line_ix + 1
            old_from_ix = int(sline[4])
            old_to_ix = int(sline[5])
            if 1 < len(atom_index_mapping[old_from_ix]) or 1 < len(
                atom_index_mapping[old_to_ix]
            ):
                raise IndicesDuplicateError("atom index mapping in bond")
            new_from_ix = atom_index_mapping[old_from_ix][0]
            new_to_ix = atom_index_mapping[old_to_ix][0]
            trailing = "\n" if bondline[-1] == "\n" else ""
            # not doing anything to the bond order
            nline = f"{prefix}{new_ix} {sline[3]} {new_from_ix} {new_to_ix}{trailing}"
            new_bondlines.append(nline)
        trailing = "\n" if self.counts[-1] == "\n" else ""
        scounts = self.counts.split()
        # suffix could be empty, but if not, prepend with space to help with constuction
        # but for now, assume that count lines are well-formed and aren't missing specs
        suffix = " " + str.join(" ", scounts[5:])
        new_counts = f"{prefix}COUNTS {len(new_atomlines)} {len(new_bondlines)}{suffix}"
        new_lines: deque[str] = deque()
        appending = True
        for line in self.lines:
            if appending:
                new_lines.append(line)
            if line.startswith("M  V30 COUNTS"):
                # unappend, place ours
                new_lines.pop()
                new_lines.append(new_counts)
            if line.startswith("M  V30 BEGIN ATOM"):
                # stop appending the old atom lines, iterate through them
                appending = False
            if line.startswith("M  V30 END ATOM"):
                appending = True
                new_lines.extend(new_atomlines)
                new_lines.append(line)
            if line.startswith("M  V30 BEGIN BOND"):
                # stop appending the old bond lines, iterate through them
                appending = False
            if line.startswith("M  V30 END BOND"):
                appending = True
                new_lines.extend(new_bondlines)
                new_lines.append(line)
        self.lines = new_lines
        return
