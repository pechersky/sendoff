import os
from typing import Final, Iterable, Iterator, Protocol, TypeVar, final

_T = TypeVar("_T")

DEBUG_BUILD: Final[bool]
"""True when compiled without optimizations; never benchmark such a build."""

__version__: Final[str]

class SupportsTextWrite(Protocol):
    """A text file-like handle accepting str payloads."""

    def write(self, data: str, /) -> object: ...

class SupportsRead(Protocol):
    """A text or binary file-like handle read in chunks."""

    def read(self, size: int, /) -> str | bytes: ...

class SupportsBytesWrite(Protocol):
    """A binary file-like handle accepting bytes payloads."""

    def write(self, data: bytes, /) -> object: ...

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

    def __new__(cls, /, lines: Iterable[str]) -> CTable:
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

    def append_record(self, /, record_name: str | bytes, value: str | bytes) -> None:
        """
        Append a field and value to the SD data record.

        Args:
            record_name: field name of the same type as the input lines
            value: field value of the same type, without a terminating newline
        """

    def ctable(self, /) -> CTable:
        """
        Parse out the underlying CTable object from the SDBlock.

        Returns:
            The CTable titled with this SDBlock's title, parsed from its MDL
            lines without validation.

        Raises:
            NotImplementedError: for byte-backed SDBlocks
        """

    @staticmethod
    def from_lines(
        lines: Iterable[str] | Iterable[bytes] | SupportsRead,
    ) -> _SDBlockIterator:
        """
        Parse in SDBlocks from a sequence of lines or a file object.

        File objects (anything with ``read1`` or ``read``) are read in 1 MiB
        chunks, so the handle's position runs ahead of the blocks produced.
        Line iterables are consumed only up to each block's $$$$ line.

        Args:
            lines: homogeneous str or bytes lines separated by the $$$$
                delimiter, or a text or binary file object

        Returns:
            A lazy iterator of SDBlocks whose values have the same type as the input

        The iterator does not close the caller-owned lines iterable.
        """

    def num_atoms(self, /) -> int:
        """
        Get the number of atoms as indicated in the MDL block.

        The count is read from the counts line and may differ from the number
        of atom lines.

        Returns:
            The number of atoms parsed from the counts line.
        """

    def num_bonds(self, /) -> int:
        """
        Get the number of bonds as indicated in the MDL block.

        The count is read from the counts line and may differ from the number
        of bond lines.

        Returns:
            The number of bonds parsed from the counts line.
        """

    def records(self, /) -> list[tuple[str, str]] | list[tuple[bytes, bytes]]:
        """
        Generate SD metadata records.

        Returns:
            Tuples of input-typed name/value pairs; multiline values include newlines.
        """

    def iter_records(self, /) -> Iterator[tuple[str, str] | tuple[bytes, bytes]]:
        """
        Lazily iterate over SD metadata records in input order.

        Returned names and values are owned and have the same type as the block.
        """

    def get(
        self, /, name: str | bytes, default: _T | None = None
    ) -> str | bytes | _T | None:
        """
        Return the value of the last record named `name`, or `default`.

        Equivalent to `dict(block.records()).get(name, default)`, without
        building the other values.

        Args:
            name: record name in the block's type
            default: returned when no record has that name

        Returns:
            The record value in the block's type, or `default`

        Raises:
            TypeError: if `name` is not the block's type
        """

    def to_dict(self, /) -> dict[str, str] | dict[bytes, bytes]:
        """
        Return the records as a dict, later duplicates winning.

        Equivalent to `dict(block.records())`.
        """

    def renumber_indices(self, /) -> None:
        """
        Renumber atom and bond indices to be 1-indexed and in order.

        Utilizes the underlying CTable's renumbering and overwrites this object's MDL lines.

        Raises:
            NotImplementedError: for byte-backed SDBlocks
        """

    @property
    def title(self, /) -> str | bytes:
        """
        Return the molecule title in the input type.

        Titles, record names, and record values are stripped as str.strip()
        does for str blocks and as bytes.strip() does for bytes blocks.
        """

    @property
    def molblock(self, /) -> str | bytes:
        """
        Return the title line and CTAB through `M  END` as one value.

        Suitable for RDKit's `Chem.MolFromMolBlock` or OpenEye's
        `OEReadMolFromBytes`.
        """

    def write(self, /, outh: SupportsTextWrite | SupportsBytesWrite) -> None:
        """
        Write an SDBlock, ending with "$$$$", to a file-like handle.

        Every line is terminated with "\\n" and the block is written in a
        single call.

        Args:
            outh: text or binary file-like handle matching the input type
        """

@final
class _SDBlockIterator(Iterator[SDBlock]):
    """
    A lazy iterator of SDBlocks.
    """

    def __next__(self, /) -> SDBlock: ...
    def close(self, /) -> None:
        """
        Stop iterating and release the input source, closing files that
        parse_sdf opened.
        """

def parse_sdf(
    path: str | bytes | os.PathLike[str] | os.PathLike[bytes], /
) -> _SDBlockIterator:
    """
    Lazily parse SDBlocks from an SDF path.

    The file is opened on the first next() and read in chunks without
    holding the GIL. str paths give str values decoded as UTF-8 on access,
    with universal newlines; bytes paths give bytes split on b"\\n" only.
    Use a bytes path for feeds that are not valid UTF-8.

    Args:
        path: a str, bytes, or os.PathLike path

    Returns:
        An iterator of SDBlocks; close() closes the file early
    """
