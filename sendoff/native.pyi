from collections.abc import Sequence
from typing import Any, Final, final

__version__: Final[str]

@final
class CTable:
    """
    Handle a Connection table (CTAB) from a single molecule block.

    Parsed in based on https://en.wikipedia.org/wiki/Chemical_table_file.
    As part of parsing, record the title, the source line, the comment line,
    and the counts line. The rest of the lines are kept but not parsed.
    The counts line is parsed to infer the format and the number of atoms.
    Note, the number of atoms in the counts line might not actually
    be the number of lines of atoms.
    """

    def __new__(cls, /, lines: Any) -> CTable:
        """
        Parse in lines representing the CTAB, including the title line.

        Args:
            lines: an Iterable of str that comprises the CTAB
                The title line is parsed in by the SDBlock, and should
                be supplied as the first in lines.
        """
    @property
    def lines(self, /) -> list[str]:
        """
        Return a copy of the original connection-table lines.
        """
    @property
    def num_atoms(self, /) -> int:
        """
        Return the atom count from the counts line.
        """
    @property
    def num_bonds(self, /) -> int:
        """
        Return the bond count from the counts line.
        """
    def renumber_indices(self, /) -> None:
        """
        Renumber the indices in the block, including changing counts line.

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
    def valid_atom_indices(self, /, strict: "bool" = False) -> bool:
        """
        Validate that the atom lines match the counts line.

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
    def valid_bond_indices(self, /, strict: "bool" = False) -> bool:
        """
        Validate that the bond lines match the counts line.

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

@final
class SDBlock:
    """
    Handle a single molecule block from an SD file.
    """

    def append_record(self, /, record_name: str, value: str) -> None:
        """
        Append a field and value to the SD data record.

        Args:
            record_name: field name
            value: field value, with no terminating newline
        """
    def ctable(self, /) -> CTable:
        """
        Parse out the underlying CTable object from the SDBlock.
        """
    @staticmethod
    def from_lines(lines: Any) -> Any:
        """
        Parse in SDBlocks from a sequence of lines.

        Args:
            lines: an Iterable of str that optionally have whitespace that is stripped,
                likely separated by the $$$$ SD delimiter

        Returns:
            An iterator of SDBlocks parsed in from the lines
        """
    def num_atoms(self, /) -> int:
        """
        Get the number of atoms indicated in the MDL block.
        """
    def num_bonds(self, /) -> int:
        """
        Get the number of bonds indicated in the MDL block.
        """
    def records(self, /) -> list[tuple[str, str]]:
        """
        Generate SD metadata records one by one.

        Returns:
            Tuples of record name and value. Multiline values include newlines.
        """
    def renumber_indices(self, /) -> None:
        """
        Renumber atom and bond indices to be 1-indexed and in order.

        Utilizes the underlying CTable's renumbering and replaces this object's MDL lines.
        """
    @property
    def title(self, /) -> str:
        """
        Return the parsed molecule title.
        """
    def write(self, /, outh: Any, with_newlines: bool = True) -> None:
        """
        Write an SDBlock to a file-like handle.

        Args:
            outh: file-like handle, such as what is returned by open(..., "w")
            with_newlines: each line is written with a trailing "\\n".
                The newline character is appended only if it was not in the line.
        """

def append_record(record_name: str, value: str) -> tuple[str, str, str]: ...
def parse_mdl(lines: Any) -> list[str]: ...
def parse_metadata(lines: Any) -> list[str]: ...
def records(metadata: Sequence[str]) -> list[tuple[str, str]]: ...
def serialize_sdfblock(
    title: str,
    mdl_lines: Sequence[str],
    metadata_lines: Sequence[str],
    with_newlines: bool,
) -> str: ...
