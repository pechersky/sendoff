import os
from typing import (
    Any,
    AnyStr,
    Final,
    Generic,
    Iterable,
    Iterator,
    Mapping,
    Protocol,
    TypeVar,
    final,
    overload,
)

_T = TypeVar("_T")
_AnyStr_co = TypeVar("_AnyStr_co", str, bytes, covariant=True)
_AnyStr_contra = TypeVar("_AnyStr_contra", str, bytes, contravariant=True)

DEBUG_BUILD: Final[bool]
"""
True when compiled without optimizations; never benchmark such a build.
"""

__version__: Final[str]

__all__ = [
    "CTable",
    "SDBlock",
    "_SDBlockIterator",
    "_SDRecordIterator",
    "parse_sdf",
    "write_sdf",
    "__version__",
    "DEBUG_BUILD",
]

class _SupportsRead(Protocol[_AnyStr_co]):
    """A text or binary file-like handle read in chunks."""

    def read(self, size: int, /) -> _AnyStr_co: ...

class _SupportsWrite(Protocol[_AnyStr_contra]):
    """A text or binary file-like handle accepting str or bytes payloads."""

    def write(self, data: _AnyStr_contra, /) -> object: ...

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

    def fix_counts(self, /) -> None:
        """
        Set the counts line to the number of atom and bond lines.

        Updates num_atoms and num_bonds to match, and keeps the rest of the
        counts line, such as the Sgroup count, as written. Use it after
        adding or removing atom or bond lines.

        Raises:
            NotImplementedError: for a V2000 format table, "-" continued atom
                or bond lines, or CTABs nested in RGROUP or TEMPLATE blocks
            ValueError: if an ATOM or BOND section has no END line, or the root
                CTAB has more than one of either
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
        Renumber the atom and bond indices in the block.

        Iterating through the atom lines, replace the indexes into a
        1-indexed order, keeping track of what maps to what.
        Then, renumber the bond lines, making sure to map the new atom indices
        properly, and keeping the rest of each line, such as CFG, as written.
        The counts line, num_atoms, and num_bonds stay as they are; call
        fix_counts after adding or removing atom or bond lines.
        This will run regardless of whether the indices are already valid.
        It cannot fix duplicate atom indices properly if the duplicate one
        occurs in the bond lines, and will raise an error.
        The lines are then in-place replaced within the CTable object.
        Only implemented for V3000 tables.

        Indices that other data refers to (atom ATTCHORD, bond ENDPTS,
        SGROUP, COLLECTION, OBJ3D, and LINKNODE lines) are not remapped, so
        renumbering refuses to change them. It also refuses "-" continued
        atom or bond lines and CTABs nested in RGROUP or TEMPLATE blocks.

        Raises:
            IndicesDuplicateError: if there was a duplicate atom index, and it was used
                somewhere in a bond line. Raised, because it is not clear which of the
                original indices to use in the remapping.
            NotImplementedError: if trying to renumber in a V2000 format table,
                or a table that renumbering would corrupt, as above
            ValueError: if an ATOM or BOND section has no END line, or the root
                CTAB has more than one of either
        """

    def valid_atom_indices(self, /, strict: bool = False) -> bool:
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
            ValueError: if the ATOM section has no END line, or the root
                CTAB has more than one ATOM section

        Returns:
            If all the checks pass, return True.
        """

    def valid_bond_indices(self, /, strict: bool = False) -> bool:
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
            ValueError: if the BOND section has no END line, or the root
                CTAB has more than one BOND section

        Returns:
            If all the checks pass, return True.
        """

@final
class SDBlock(Generic[AnyStr]):
    """
    Handle a single molecule block from an SD file.

    `SDBlock[str]` blocks come from text input and `SDBlock[bytes]` blocks
    from binary input; names, values, and lines have the block's type.
    """

    def __bool__(self, /) -> bool:
        """
        A block is true even without records, as before it had a length.
        """

    def __contains__(self, name: AnyStr, /) -> bool:
        """
        Whether any record is named `name`: `name in block.names()`.

        Like `get`, every record header is checked, so a malformed one raises
        even after a match.

        Raises:
            TypeError: if `name` is not the block's type
            ValueError: if a record header has no `<name>`
        """

    def __copy__(self, /) -> SDBlock[AnyStr]:
        """
        Return an independent copy; a block holds no other objects, so
        `copy.copy` and `copy.deepcopy` are the same.
        """

    def __deepcopy__(self, memo: dict[int, Any], /) -> SDBlock[AnyStr]: ...
    def __eq__(self, value: object, /) -> bool: ...
    def __len__(self, /) -> int:
        """
        The number of records, duplicates included: `len(block.records())`.
        """

    def __ne__(self, value: object, /) -> bool: ...
    def __new__(
        cls,
        /,
        sdf: AnyStr,
        records: (
            Mapping[AnyStr, AnyStr] | Iterable[tuple[AnyStr, AnyStr]] | None
        ) = None,
    ) -> SDBlock[AnyStr]:
        """
        Parse one block from SDF str or bytes, such as a molblock.

        The block ends at an optional `$$$$` line, after which only
        whitespace may follow. str input has universal newlines.

        Args:
            sdf: a molblock or one SDF record, as str or bytes
            records: SD data to append after any in `sdf`, as a mapping or
                an iterable of (name, value) tuples of the same type

        Raises:
            ValueError: if `sdf` holds more than one block
            TypeError: if `sdf` is not str or bytes, or a record is not its type
        """

    def __reduce__(self, /) -> tuple[Any, ...]:
        """
        Pickle the block's exact state rather than its SDF text, which
        does not always parse back to the same block.
        """

    def __repr__(self, /) -> str: ...
    @classmethod
    def _from_state(
        cls, /, data: bytes, title_end: int, mdl_end: int, text: bool
    ) -> SDBlock[Any]:
        """
        Rebuild a block pickled by `__reduce__`.

        Raises:
            ValueError: if the state cannot come from a block
        """

    def append_record(self, /, record_name: AnyStr, value: AnyStr) -> None:
        """
        Append a field and value of the same type as the input block.

        Args:
            record_name: field name
            value: field value, with no terminating newline

        Raises:
            TypeError: if `record_name` or `value` is not the block's type
        """

    def ctable(self: SDBlock[str], /) -> CTable:
        """
        Parse out the underlying CTable object from a text SDBlock.

        Returns:
            The CTable titled with this block's title, parsed from its MDL
            lines without validation

        Raises:
            NotImplementedError: for byte-backed SDBlocks
        """

    def fix_counts(self: SDBlock[str], /) -> None:
        """
        Set the counts line to the number of atom and bond lines.

        Uses CTable.fix_counts and replaces this block's MDL lines.

        Raises:
            NotImplementedError: for byte-backed SDBlocks
        """

    @overload
    @staticmethod
    def from_lines(
        lines: Iterable[str] | _SupportsRead[str],
    ) -> _SDBlockIterator[str]: ...
    @overload
    @staticmethod
    def from_lines(
        lines: Iterable[bytes] | _SupportsRead[bytes],
    ) -> _SDBlockIterator[bytes]: ...
    @overload
    @staticmethod
    def from_lines(
        lines: (
            Iterable[str] | Iterable[bytes] | _SupportsRead[str] | _SupportsRead[bytes]
        ),
    ) -> _SDBlockIterator[Any]:
        """
        Parse SDBlocks from homogeneous str or bytes lines or a file object.

        File objects are read in 128 KiB chunks with `read1`, when available,
        or `read`, so the handle's position runs ahead of the blocks
        produced. Line iterables are consumed only up to each block's `$$$$`
        line, and the iterator does not close them.

        Args:
            lines: homogeneous str or bytes lines separated by the `$$$$`
                delimiter, or a text or binary file object

        Returns:
            An iterator that produces owned SDBlocks as input is consumed

        Raises:
            TypeError: for a str or bytes argument; use `from_sdf` for those
        """

    @overload
    @staticmethod
    def from_sdf(sdf: str) -> _SDBlockIterator[str]: ...
    @overload
    @staticmethod
    def from_sdf(sdf: bytes) -> _SDBlockIterator[bytes]: ...
    @overload
    @staticmethod
    def from_sdf(sdf: str | bytes) -> _SDBlockIterator[Any]:
        """
        Lazily parse SDBlocks from SDF str or bytes held in memory.

        Equivalent to `from_lines(io.StringIO(sdf))` or
        `from_lines(io.BytesIO(sdf))`, without splitting lines: str input
        has universal newlines, and bytes are split on b"\\n" only.

        Args:
            sdf: SDF content as str or bytes

        Returns:
            An iterator of SDBlocks of the input type

        Raises:
            TypeError: if `sdf` is not str or bytes
        """

    @overload
    def get(self, /, name: AnyStr, default: None = None) -> AnyStr | None: ...
    @overload
    def get(self, /, name: AnyStr, default: _T) -> AnyStr | _T:
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

    def iter_records(self, /) -> _SDRecordIterator[AnyStr]:
        """
        Lazily iterate over SD data records in input order.

        Names and values are owned and have the block's type.
        """

    @property
    def mdl(self, /) -> AnyStr:
        """
        Return the CTAB after the title, through `M  END`, as one value.

        Every line ends in a newline, so `needle in block.mdl` searches the
        CTAB and `block.mdl.splitlines()` gives its lines. `block.molblock`
        is `block.title`, the title line's ending, then `block.mdl`.
        """

    @property
    def molblock(self, /) -> AnyStr:
        """
        Return the title line and CTAB through `M  END` as one value.

        Suitable for RDKit's `Chem.MolFromMolBlock` or OpenEye's
        `OEReadMolFromBytes`.
        """

    @molblock.setter
    def molblock(self, /, molblock: AnyStr) -> None:
        """
        Replace the title line and CTAB, keeping the SD data.

        For example, `block.molblock = Chem.MolToMolBlock(mol)`. Text
        molblocks have CRLF and CR line endings translated to LF.

        Raises:
            ValueError: if the molblock has no `M  END` line, has more than
                whitespace after it, or has a line starting with `$$$$`
            TypeError: if the molblock is not the block's type
        """

    def names(self, /) -> list[AnyStr]:
        """
        Return the record names in input order, duplicates included.

        Equivalent to `[name for name, _ in block.records()]`, without
        building the values. Blocks need not share names, so collect them
        across a file with `Counter(n for b in blocks for n in b.names())`.
        """

    def num_atoms(self, /) -> int:
        """
        Get the number of atoms declared by the counts line.

        The count may differ from the number of atom lines.
        """

    def num_bonds(self, /) -> int:
        """
        Get the number of bonds declared by the counts line.

        The count may differ from the number of bond lines.
        """

    def records(self, /) -> list[tuple[AnyStr, AnyStr]]:
        """
        Return all SD data records.

        Returns:
            Tuples of record name and value in the input type; multiline
            values include their inner newlines
        """

    def remove_record(self, /, name: AnyStr) -> int:
        """
        Remove every record named `name`.

        Args:
            name: record name in the block's type

        Returns:
            How many records were removed

        Raises:
            TypeError: if `name` is not the block's type
        """

    def renumber_indices(self: SDBlock[str], /) -> None:
        """
        Renumber atom and bond indices to be 1-indexed and in order.

        Uses the CTable renumbering and replaces this block's MDL lines. The
        counts line stays as it is; see fix_counts.

        Raises:
            NotImplementedError: for byte-backed SDBlocks
        """

    def replace_records(
        self, /, records: Mapping[AnyStr, AnyStr] | Iterable[tuple[AnyStr, AnyStr]]
    ) -> None:
        """
        Replace all SD data with the given records.

        Args:
            records: a mapping, or an iterable of (name, value) tuples, in
                the block's type

        Raises:
            TypeError: if a name or value is not the block's type
        """

    def set_record(self, /, name: AnyStr, value: AnyStr) -> None:
        """
        Set a record, keeping the position of the first with that name.

        Later records with the same name are removed; if there is none, the
        record is appended.

        Args:
            name: record name in the block's type
            value: record value in the block's type

        Raises:
            TypeError: if `name` or `value` is not the block's type
        """

    @property
    def title(self, /) -> AnyStr:
        """
        Return the molecule title in the input type.

        Titles, record names, and record values are stripped as `str.strip()`
        does for str blocks and as `bytes.strip()` does for bytes blocks.
        """

    @title.setter
    def title(self, /, title: AnyStr) -> None:
        """
        Replace the title line.

        Raises:
            ValueError: if the title contains a line break
            TypeError: if the title is not the block's type
        """

    def to_dict(self, /) -> dict[AnyStr, AnyStr]:
        """
        Return the records as a dict, later duplicates winning.

        Equivalent to `dict(block.records())`.
        """

    def to_sdf(self, /) -> AnyStr:
        """
        Return the block as SDF text ending with $$$$, in the input type.

        The inverse of `SDBlock(sdf)`, and what `write` writes.
        """

    def write(self, /, outh: _SupportsWrite[AnyStr]) -> None:
        """
        Write the block, ending with $$$$, to a file-like handle.

        Lines end as they were read, so blocks from a bytes path keep "\\r\\n";
        str paths read universal newlines, so their lines end in "\\n". Lines
        the block adds or ends (a stripped title, new records, a missing $$$$,
        a last line without "\\n") end like the title line. The block is
        written in one `write` call.

        Args:
            outh: a text or binary file-like handle matching the input type
        """

@final
class _SDBlockIterator(Iterator[SDBlock[AnyStr]]):
    """
    Lazily parse SDBlocks from an input source.
    """

    def __iter__(self, /) -> _SDBlockIterator[AnyStr]: ...
    def __next__(self, /) -> SDBlock[AnyStr]: ...
    def close(self, /) -> None:
        """
        Stop iterating and release the input source, closing files that
        `parse_sdf` opened.
        """

@final
class _SDRecordIterator(Iterator[tuple[AnyStr, AnyStr]]):
    """
    Lazily iterate over one block's SD data records.
    """

    def __iter__(self, /) -> _SDRecordIterator[AnyStr]: ...
    def __next__(self, /) -> tuple[AnyStr, AnyStr]: ...

@overload
def parse_sdf(path: str | os.PathLike[str]) -> _SDBlockIterator[str]: ...
@overload
def parse_sdf(path: bytes | os.PathLike[bytes]) -> _SDBlockIterator[bytes]: ...
@overload
def parse_sdf(
    path: str | bytes | os.PathLike[str] | os.PathLike[bytes],
) -> _SDBlockIterator[Any]:
    """
    Lazily parse SDBlocks from an SDF path.

    The file is opened on the first `next()`. Regular files on local Linux
    filesystems are memory-mapped and split a batch of blocks at a time
    without holding the GIL; pipes, network filesystems, and other files are
    read in chunks. As with any mapped file, truncating the file while
    iterating can kill the process with SIGBUS, so write results to a new
    path. str paths give str values decoded as UTF-8 on access, with
    universal newlines; bytes paths give bytes split on b"\\n" only. Use a
    bytes path for feeds that are not valid UTF-8: from a str path, their
    text raises UnicodeDecodeError with a note naming the block and record.

    Args:
        path: a str, bytes, or os.PathLike path

    Returns:
        An iterator of SDBlocks; `close()` releases the file early
    """

@overload
def write_sdf(
    dest: str | bytes | os.PathLike[str] | os.PathLike[bytes],
    blocks: Iterable[SDBlock[Any]],
) -> int: ...
@overload
def write_sdf(dest: _SupportsWrite[AnyStr], blocks: Iterable[SDBlock[AnyStr]]) -> int:
    """
    Write SDBlocks to a path or a file object, and return how many.

    Blocks are written as `SDBlock.write` would, but 64-256 KiB at a time:
    to a path without holding the GIL, or to a file object with one
    `write` call per chunk, as str for text blocks or bytes for byte-backed
    ones. Text blocks must be valid UTF-8.

    Args:
        dest: a str, bytes, or os.PathLike path to create or truncate, or a
            text or binary file object
        blocks: an iterable of SDBlocks, all str or all bytes for a file
            object

    Returns:
        The number of blocks written
    """
