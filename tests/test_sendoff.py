# mypy: allow-untyped-decorators
"""Test sendoff parsing of metadata on generated SDFs."""
import io
import os
from itertools import chain
from pathlib import Path
from typing import Generator, Iterable, Iterator, cast

import pytest

from sendoff.ctable import (
    IndicesDuplicateError,
    IndicesMismatchError,
    IndicesOutOfOrderError,
)
from sendoff.sdblock import Pathy, SDBlock, parse_sdf


def test_single_mol_read(single_mol_sdf: Pathy) -> None:
    """An sdf of a single molecule contains one molecule.

    Args:
        single_mol_sdf: pytest fixture of a Path to the sdf
    """
    mols = list(parse_sdf(single_mol_sdf))
    assert len(mols) == 1


def test_double_mol_read(double_mol_sdf: Pathy) -> None:
    """An sdf of two molecule contains two molecules.

    Args:
        double_mol_sdf: pytest fixture of a Path to the sdf
    """
    mols = list(parse_sdf(double_mol_sdf))
    assert len(mols) == 2


@pytest.mark.parametrize("trailing_newline", ["", "\n"])
def test_unterminated_final_record_is_returned(
    tmp_path: Path, trailing_newline: str
) -> None:
    """Return a final block without a delimiter, with or without a newline.

    Args:
        tmp_path: Temporary directory for the SDF fixture.
        trailing_newline: Newline appended after the unterminated value.
    """
    sdf = tmp_path / "unterminated.sdf"
    sdf.write_text(
        "Closed\n\nM  END\n$$$$\n"
        "Title\n\nM  END\n> <Record>\nValue" + trailing_newline
    )

    blocks = list(parse_sdf(sdf))

    assert [block.title for block in blocks] == ["Closed", "Title"]
    assert blocks[-1].records() == [("Record", "Value")]


def test_empty_input_has_no_blocks() -> None:
    """Do not synthesize a block for empty input."""
    assert list(SDBlock.from_lines([])) == []


def test_bytes_lines_preserve_invalid_bytes() -> None:
    """Byte-backed blocks retain non-UTF-8 title and metadata octets."""
    raw = [
        b"\xff title \n",
        b"\n",
        b"M  END\n",
        b"> <field\xfe>\n",
        b"value\xfd\n",
        b"\n",
        b"$$$$\n",
        b"second\n",
        b"\n",
        b"M  END\n",
    ]

    blocks = SDBlock.from_lines(raw)
    block = next(blocks)
    tail = next(blocks)

    assert block.title == b"\xff title"
    assert block.records() == [(b"field\xfe", b"value\xfd")]
    assert tail.title == b"second"
    output = io.BytesIO()
    block.write(output)
    assert b"\xff" in output.getvalue()
    assert b"field\xfe" in output.getvalue()
    assert b"value\xfd" in output.getvalue()
    reparsed = next(SDBlock.from_lines(output.getvalue().splitlines(keepends=True)))
    assert reparsed.title == block.title
    assert reparsed.records() == block.records()

    with pytest.raises(NotImplementedError, match="byte-backed"):
        block.ctable()


def test_bytes_path_uses_binary_input(tmp_path: Path) -> None:
    """A bytes filesystem path yields bytes and preserves invalid UTF-8.

    Args:
        tmp_path: Temporary directory for the SDF fixture.
    """
    path = tmp_path / "bytes.sdf"
    path.write_bytes(b"\xff title\n\nM  END\n> <field>\nvalue\xfe\n\n$$$$\n")

    block = next(parse_sdf(os.fsencode(path)))

    assert block.title == b"\xff title"
    assert block.records() == [(b"field", b"value\xfe")]


def test_text_path_preserves_text_values(tmp_path: Path) -> None:
    """A str filesystem path continues to use text decoding.

    Args:
        tmp_path: Temporary directory for the SDF fixture.
    """
    path = tmp_path / "text.sdf"
    path.write_text("title é\n\nM  END\n> <field>\nvalue é\n\n$$$$\n")

    block = next(parse_sdf(path))

    assert block.title == "title é"
    assert block.records() == [("field", "value é")]


def test_from_lines_rejects_mixed_text_and_bytes() -> None:
    """A single line source cannot mix text and bytes."""
    mixed_lines = cast(Iterable[str], [b"title\n", "M  END\n"])
    with pytest.raises(TypeError, match="mixed str and bytes"):
        list(SDBlock.from_lines(mixed_lines))


def test_from_lines_yields_incrementally_and_retains_blocks() -> None:
    """Only consume through the next delimiter and retain already parsed data."""
    consumed: list[str] = []
    first_record = [
        "First\n",
        "\n",
        "M  END\n",
        "> <Record>\n",
        "one\n",
        "\n",
        "$$$$\n",
    ]
    second_record = ["Second\n", "\n", "M  END\n", "$$$$\n"]

    def lines() -> Iterator[str]:
        for line in first_record + second_record:
            consumed.append(line)
            yield line

    blocks = SDBlock.from_lines(lines())

    assert consumed == []
    assert iter(blocks) is blocks
    first = next(blocks)
    assert consumed == first_record

    second = next(blocks)

    assert first.title == "First"
    assert first.records() == [("Record", "one")]
    assert second.title == "Second"
    assert consumed == first_record + second_record


def test_from_lines_surfaces_source_errors_on_advance() -> None:
    """Yield completed blocks before surfacing a later source error."""

    def lines() -> Iterator[str]:
        yield "First\n"
        yield "\n"
        yield "M  END\n"
        yield "$$$$\n"
        raise OSError("read failed")

    blocks = SDBlock.from_lines(lines())

    assert next(blocks).title == "First"
    with pytest.raises(OSError, match="read failed"):
        next(blocks)
    with pytest.raises(StopIteration):
        next(blocks)


def test_from_lines_does_not_close_caller_owned_iterable() -> None:
    """Abandoning the parser leaves the caller's source open."""
    closed: list[bool] = []

    def lines() -> Generator[str, None, None]:
        try:
            yield "First\n"
            yield "\n"
            yield "M  END\n"
            yield "$$$$\n"
        finally:
            closed.append(True)

    source = lines()
    blocks = SDBlock.from_lines(source)
    assert next(blocks).title == "First"
    del blocks

    assert closed == []
    source.close()
    assert closed == [True]


def test_parse_sdf_closes_owned_file_when_iterator_closes(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Closing the path parser exits its file context.

    Args:
        monkeypatch: Fixture used to replace the module-level open.
    """
    first_record = "First\n\nM  END\n$$$$\n"
    stream = io.StringIO(first_record + "Second\n\nM  END\n$$$$\n")

    def fake_open(path: str | bytes, mode: str) -> io.StringIO:
        assert path == "unused.sdf"
        assert mode == "r"
        return stream

    monkeypatch.setattr("builtins.open", fake_open)
    blocks = parse_sdf("unused.sdf")

    assert next(blocks).title == "First"
    assert stream.tell() == len(first_record)
    assert not stream.closed
    blocks.close()
    assert stream.closed


def test_bytes_block_requires_binary_append_values() -> None:
    """Appending to a bytes block cannot introduce implicitly encoded text."""
    block = next(SDBlock.from_lines([b"title\n", b"\n", b"M  END\n"]))

    with pytest.raises(TypeError, match="must be appended as bytes"):
        block.append_record("field", b"value")

    block.append_record(b"field", b"value")
    assert block.records() == [(b"field", b"value")]


def test_titled_mol_read(single_titled_mol_sdf: Pathy) -> None:
    """An sdf of a single titled molecule contains it with the expected title.

    Args:
        single_titled_mol_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock = next(parse_sdf(single_titled_mol_sdf))
    assert mol.title == "Title"


@pytest.mark.xfail  # sendoff does not support SDF titles that start with $$$$
def test_first_delimiter_titled_mol_read(
    double_first_delimiter_titled_mol_sdf: Pathy,
) -> None:
    """SDF with 1st molecule delimiter titled, containing with expected titles.

    Args:
        double_first_delimiter_titled_mol_sdf: pytest fixture of a Path to sdf
    """
    supp = parse_sdf(double_first_delimiter_titled_mol_sdf)
    mol: SDBlock = next(supp)
    other_mol: SDBlock = next(supp)
    assert mol.title == "$$$$"
    assert other_mol.title == "Title"


def test_second_delimiter_titled_mol_read(
    double_second_delimiter_titled_mol_sdf: Pathy,
) -> None:
    """SDF with 2nd molecule delimiter titled, containing with expected titles.

    Args:
        double_second_delimiter_titled_mol_sdf: pytest fixture of a Path to sdf
    """
    supp = parse_sdf(double_second_delimiter_titled_mol_sdf)
    other_mol: SDBlock = next(supp)
    mol: SDBlock = next(supp)
    assert other_mol.title == "Title"
    assert mol.title == "$$$$"


def test_single_record_mol_read(single_record_mol_sdf: Pathy) -> None:
    """An sdf of a delimiter record valued molecule contains the expected data.

    Args:
        single_record_mol_sdf: pytest fixture of a Path to the sdf
    """
    mol = next(parse_sdf(single_record_mol_sdf))
    records = list(mol.records())
    assert ("Record", "Value") in records


def test_iter_records_preserves_order_duplicates_and_multiline_values() -> None:
    """Iterate records on demand without changing list-returning records()."""
    block = next(
        SDBlock.from_lines(
            [
                "Title\n",
                "\n",
                "M  END\n",
                "> <Record>\n",
                "first\n",
                "second\n",
                "\n",
                "> <Record>\n",
                "third\n",
                "\n",
                "$$$$\n",
            ]
        )
    )
    iterator = block.iter_records()
    assert iter(iterator) is iterator
    first = next(iterator)
    second = next(iterator)

    assert first == ("Record", "first\nsecond")
    assert second == ("Record", "third")
    assert first == ("Record", "first\nsecond")
    records = block.records()
    assert isinstance(records, list)
    assert records == [first, second]


def test_iter_records_defers_malformed_field_error() -> None:
    """A later malformed record does not prevent reading earlier fields."""
    block = next(
        SDBlock.from_lines(
            [
                "Title\n",
                "\n",
                "M  END\n",
                "> <Good>\n",
                "value\n",
                "\n",
                "> missing name\n",
                "other\n",
                "\n",
                "$$$$\n",
            ]
        )
    )
    iterator = block.iter_records()

    assert next(iterator) == ("Good", "value")
    with pytest.raises(ValueError, match="missing a name"):
        next(iterator)
    with pytest.raises(StopIteration):
        next(iterator)
    with pytest.raises(ValueError, match="missing a name"):
        block.records()


def test_iter_records_preserves_bytes_and_owns_its_block() -> None:
    """Byte record iteration remains lazy and retains source data safely."""
    block = next(
        SDBlock.from_lines(
            [
                b"Title\n",
                b"\n",
                b"M  END\n",
                b"> <Record>\n",
                b"first\xfe\n",
                b"second\n",
                b"\n",
                b"> <Record>\n",
                b"third\n",
                b"\n",
                b"$$$$\n",
            ]
        )
    )
    iterator = block.iter_records()
    del block

    record = next(iterator)

    assert record == (b"Record", b"first\xfe\nsecond")
    assert list(iterator) == [(b"Record", b"third")]


def test_string_and_bytes_whitespace_follow_their_input_rules() -> None:
    """Text uses Unicode trim while bytes preserve non-ASCII whitespace."""
    text_block = next(
        SDBlock.from_lines(
            [
                " \u2003Title\u2003 \n",
                "\n",
                "M  END\n",
                "> <Record>\u2003\n",
                "\u2003Text\u2003\n",
                "\n",
                "$$$$\n",
            ]
        )
    )
    byte_block = next(
        SDBlock.from_lines(
            [
                b"\xc2\xa0Title\xc2\xa0\n",
                b"\n",
                b"M  END\n",
                b"> <Record>\n",
                b"\xc2\xa0Bytes\xc2\xa0\n",
                b"\n",
                b"$$$$\n",
            ]
        )
    )

    assert text_block.title == "Title"
    assert text_block.records() == [("Record", "Text")]
    assert byte_block.title == b"\xc2\xa0Title\xc2\xa0"
    assert byte_block.records() == [(b"Record", b"\xc2\xa0Bytes\xc2\xa0")]


@pytest.mark.xfail  # sendoff does not support metadata vals starting with $$$$
def test_delimiter_record_mol_read(
    single_delimiter_record_mol_sdf: Pathy,
) -> None:
    """An sdf of a delimiter record valued molecule contains the expected data.

    Args:
        single_delimiter_record_mol_sdf: pytest fixture of a Path to the sdf
    """
    mol = next(parse_sdf(single_delimiter_record_mol_sdf))
    records = list(mol.records())
    assert ("Record", "$$$$") in records


def test_multiline_record_mol_read(
    single_multiline_record_mol_sdf: Pathy,
) -> None:
    """An sdf of a multiline record valued molecule contains the expected data.

    Args:
        single_multiline_record_mol_sdf: pytest fixture of a Path to the sdf
    """
    mol = next(parse_sdf(single_multiline_record_mol_sdf))
    records = list(mol.records())
    assert ("Record", "0\n1\n2\n3") in records


def test_empty_string_record_mol_read(
    single_empty_string_record_mol_sdf: Pathy,
) -> None:
    """An sdf of an empty string record valued molecule contains the expected data.

    Args:
        single_empty_string_record_mol_sdf: pytest fixture of a Path to the sdf
    """
    mol = next(parse_sdf(single_empty_string_record_mol_sdf))
    records = list(mol.records())
    assert ("Record", "") in records


def test_multiline_record_name_mol_read(
    single_multiline_record_name_mol_sdf: Pathy,
) -> None:
    """An sdf of a multiline record name molecule cannot contain the expected data.

    Args:
        single_multiline_record_name_mol_sdf: fixture of a Path to the sdf
    """
    mol = next(parse_sdf(single_multiline_record_name_mol_sdf))
    records = list(mol.records())
    assert "Rec\nord" not in map(str, chain(zip(*records)))


def test_single_mol_newline_write(single_mol_sdf: Pathy) -> None:
    """An sdf block gets written with as many lines as the input.

    Args:
        single_mol_sdf: pytest fixture of a Path to the sdf
    """
    mols = list(parse_sdf(single_mol_sdf))
    buffer = io.StringIO()
    for mol in mols:
        mol.write(buffer)
    assert len(buffer.getvalue().splitlines()) == len(
        open(single_mol_sdf).read().splitlines()
    )


def test_single_mol_newline_write_splitlines(single_mol_sdf: Pathy) -> None:
    """An sdf block gets written with as many lines as a splitlines input.

    Args:
        single_mol_sdf: pytest fixture of a Path to the sdf
    """
    mols = list(SDBlock.from_lines(open(single_mol_sdf).read().splitlines()))
    buffer = io.StringIO()
    for mol in mols:
        mol.write(buffer)
    assert len(buffer.getvalue().splitlines()) == len(
        open(single_mol_sdf).read().splitlines()
    )


def test_single_mol_newline_write_splitlines_no_trailing(single_mol_sdf: Pathy) -> None:
    """An sdf block written with no trailing newlines has fewer lines than input.

    Args:
        single_mol_sdf: pytest fixture of a Path to the sdf
    """
    mols = list(SDBlock.from_lines(open(single_mol_sdf).read().splitlines()))
    buffer = io.StringIO()
    for mol in mols:
        mol.write(buffer, with_newlines=False)
    assert len(buffer.getvalue().splitlines()) < len(
        open(single_mol_sdf).read().splitlines()
    )


def test_right_angle_bracket_record_name_mol_read(
    single_right_angle_bracket_record_name_mol_sdf: Pathy,
) -> None:
    """An sdf with a record name with a ">" is present in the read data as expected.

    Args:
        single_right_angle_bracket_record_name_mol_sdf: fixture of a Path to the sdf
    """
    mol: SDBlock = next(parse_sdf(single_right_angle_bracket_record_name_mol_sdf))
    records = list(mol.records())
    assert records == [("Rec>ord", "Value")]


def test_0_atoms_mol_num_atoms(single_0_atoms_mol_sdf: Pathy) -> None:
    """An sdf block with 0 atoms written parses as 0 atoms.

    Args:
        single_0_atoms_mol_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock = next(parse_sdf(single_0_atoms_mol_sdf))
    assert mol.num_atoms() == 0


def test_star_atom_mol_num_atoms(single_star_atom_mol_sdf: Pathy) -> None:
    """An sdf block with a single star atom written parses as 1 atom.

    Args:
        single_star_atom_mol_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock = next(parse_sdf(single_star_atom_mol_sdf))
    assert mol.num_atoms() == 1


def test_999_atoms_mol_num_atoms(single_999_atoms_mol_sdf: Pathy) -> None:
    """An sdf block with 999 atoms written parses as 999 atoms.

    Args:
        single_999_atoms_mol_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock = next(parse_sdf(single_999_atoms_mol_sdf))
    assert mol.num_atoms() == 999


def test_1001_atoms_mol_num_atoms(single_1001_atoms_mol_sdf: Pathy) -> None:
    """An sdf block with 1001 atoms written parses as 1001 atoms.

    Args:
        single_1001_atoms_mol_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock = next(parse_sdf(single_1001_atoms_mol_sdf))
    assert mol.num_atoms() == 1001


def test_0_atoms_mol_num_bonds(single_0_atoms_mol_sdf: Pathy) -> None:
    """An sdf block with 0 atoms written parses as 0 bonds.

    Args:
        single_0_atoms_mol_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock = next(parse_sdf(single_0_atoms_mol_sdf))
    assert mol.num_bonds() == 0


def test_star_atom_mol_num_bonds(single_star_atom_mol_sdf: Pathy) -> None:
    """An sdf block with a single star atom written parses as 0 bonds.

    Args:
        single_star_atom_mol_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock = next(parse_sdf(single_star_atom_mol_sdf))
    assert mol.num_bonds() == 0


def test_999_atoms_mol_num_bonds(single_999_atoms_mol_sdf: Pathy) -> None:
    """An sdf block with 999 atoms written parses as 998 bonds.

    Args:
        single_999_atoms_mol_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock = next(parse_sdf(single_999_atoms_mol_sdf))
    assert mol.num_bonds() == 998


def test_1001_atoms_mol_num_bonds(single_1001_atoms_mol_sdf: Pathy) -> None:
    """An sdf block with 1001 atoms written parses as 1000 bonds.

    Args:
        single_1001_atoms_mol_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock = next(parse_sdf(single_1001_atoms_mol_sdf))
    assert mol.num_bonds() == 1000


def test_single_mol_v2000_sdf_valid_atom_indices_not_implemented(
    single_mol_v2000_sdf: Pathy,
) -> None:
    """A v2000 sdf block cannot yet be validated for atom indices.

    Args:
        single_mol_v2000_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock = next(parse_sdf(single_mol_v2000_sdf))
    ctable = mol.ctable()
    with pytest.raises(NotImplementedError):
        assert ctable.valid_atom_indices()


def test_single_mol_v2000_sdf_valid_bond_indices_not_implemented(
    single_mol_v2000_sdf: Pathy,
) -> None:
    """A v2000 sdf block cannot yet be validated for bond indices.

    Args:
        single_mol_v2000_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock = next(parse_sdf(single_mol_v2000_sdf))
    ctable = mol.ctable()
    with pytest.raises(NotImplementedError):
        assert ctable.valid_bond_indices()


def test_single_large_atom_index_valid_atom_indices(
    single_large_atom_index_v3000_sdf: Pathy,
) -> None:
    """An sdf block with large atom indices has valid indices.

    That means that the number of atoms matches the counts line,
    and each index is unique and integer.
    Does not check that the indices are 1-indexed and in order.

    Args:
        single_large_atom_index_v3000_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock = next(parse_sdf(single_large_atom_index_v3000_sdf))
    ctable = mol.ctable()
    assert ctable.valid_atom_indices()


def test_single_large_atom_index_strict_invalid_atom_indices(
    single_large_atom_index_v3000_sdf: Pathy,
) -> None:
    """An sdf block with large atom indices has invalid indices, when strict.

    That means that the number of atoms matches the counts line,
    and each index is unique and integer.
    However, the indices are not 1-indexed and in order.

    Args:
        single_large_atom_index_v3000_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock = next(parse_sdf(single_large_atom_index_v3000_sdf))
    ctable = mol.ctable()
    with pytest.raises(IndicesOutOfOrderError, match="atoms"):
        ctable.valid_atom_indices(strict=True)


def test_single_large_bond_index_valid_bond_indices(
    single_large_bond_index_v3000_sdf: Pathy,
) -> None:
    """An sdf block with large bond indices has valid indices.

    That means that the number of bonds matches the counts line,
    and each index is unique and integer.
    Does not check that the indices are 1-indexed and in order.

    Args:
        single_large_bond_index_v3000_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock = next(parse_sdf(single_large_bond_index_v3000_sdf))
    ctable = mol.ctable()
    assert ctable.valid_bond_indices()


def test_single_large_bond_index_strict_invalid_bond_indices(
    single_large_bond_index_v3000_sdf: Pathy,
) -> None:
    """An sdf block with large bond indices has invalid indices, when strict.

    That means that the number of bonds matches the counts line,
    and each index is unique and integer.
    However, the indices are not 1-indexed and in order.

    Args:
        single_large_bond_index_v3000_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock = next(parse_sdf(single_large_bond_index_v3000_sdf))
    ctable = mol.ctable()
    with pytest.raises(IndicesOutOfOrderError, match="bonds"):
        ctable.valid_bond_indices(strict=True)


def test_single_missing_atom_line_invalid_atom_indices(
    single_missing_atom_line_v3000_sdf: Pathy,
) -> None:
    """An sdf block with a missing atom line has invalid indices.

    That means that the number of atoms is fewer than the counts line.

    Args:
        single_missing_atom_line_v3000_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock = next(parse_sdf(single_missing_atom_line_v3000_sdf))
    ctable = mol.ctable()
    with pytest.raises(IndicesMismatchError, match="fewer atom lines than count line"):
        ctable.valid_atom_indices()


def test_single_extra_atom_line_invalid_atom_indices(
    single_extra_atom_line_v3000_sdf: Pathy,
) -> None:
    """An sdf block with a extra atom line has invalid indices.

    That means that the number of atoms is more than the counts line.

    Args:
        single_extra_atom_line_v3000_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock = next(parse_sdf(single_extra_atom_line_v3000_sdf))
    ctable = mol.ctable()
    with pytest.raises(IndicesMismatchError, match="more atom lines than count line"):
        ctable.valid_atom_indices()


def test_single_missing_bond_line_invalid_atom_indices(
    single_missing_bond_line_v3000_sdf: Pathy,
) -> None:
    """An sdf block with a missing bond line has invalid indices.

    That means that the number of bonds is fewer than the counts line.

    Args:
        single_missing_bond_line_v3000_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock = next(parse_sdf(single_missing_bond_line_v3000_sdf))
    ctable = mol.ctable()
    with pytest.raises(IndicesMismatchError, match="fewer bond lines than count line"):
        ctable.valid_bond_indices()


def test_single_extra_bond_line_invalid_atom_indices(
    single_extra_bond_line_v3000_sdf: Pathy,
) -> None:
    """An sdf block with a extra bond line has invalid indices.

    That means that the number of bonds is more than the counts line.

    Args:
        single_extra_bond_line_v3000_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock = next(parse_sdf(single_extra_bond_line_v3000_sdf))
    ctable = mol.ctable()
    with pytest.raises(IndicesMismatchError, match="more bond lines than count line"):
        ctable.valid_bond_indices()


def test_single_duplicate_atom_index_invalid_atom_indices(
    single_duplicate_atom_index_v3000_sdf: Pathy,
) -> None:
    """An sdf block with a duplicate atom index has invalid indices.

    Args:
        single_duplicate_atom_index_v3000_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock = next(parse_sdf(single_duplicate_atom_index_v3000_sdf))
    ctable = mol.ctable()
    with pytest.raises(IndicesDuplicateError, match="atoms"):
        ctable.valid_atom_indices()


def test_single_duplicate_bond_index_invalid_bond_indices(
    single_duplicate_bond_index_v3000_sdf: Pathy,
) -> None:
    """An sdf block with a duplicate bond index has invalid indices.

    Args:
        single_duplicate_bond_index_v3000_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock = next(parse_sdf(single_duplicate_bond_index_v3000_sdf))
    ctable = mol.ctable()
    with pytest.raises(IndicesDuplicateError, match="bonds"):
        ctable.valid_bond_indices()


def test_single_shuffled_atom_lines_valid_atom_indices(
    single_shuffled_atom_lines_v3000_sdf: Pathy,
) -> None:
    """An sdf block with shuffled atom lines has valid indices.

    That means that the number of atoms matches the counts line,
    and each index is unique and integer.
    Does not check that the indices are 1-indexed and in order.

    Args:
        single_shuffled_atom_lines_v3000_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock = next(parse_sdf(single_shuffled_atom_lines_v3000_sdf))
    ctable = mol.ctable()
    assert ctable.valid_atom_indices()


def test_single_shuffled_atom_lines_strict_invalid_atom_indices(
    single_shuffled_atom_lines_v3000_sdf: Pathy,
) -> None:
    """An sdf block with shuffled atom lines has invalid indices, when strict.

    That means that the number of atoms matches the counts line,
    and each index is unique and integer.
    However, the indices are not 1-indexed and in order.

    Args:
        single_shuffled_atom_lines_v3000_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock = next(parse_sdf(single_shuffled_atom_lines_v3000_sdf))
    ctable = mol.ctable()
    with pytest.raises(IndicesOutOfOrderError, match="atoms"):
        ctable.valid_atom_indices(strict=True)


def test_single_shuffled_bond_lines_valid_bond_indices(
    single_shuffled_bond_lines_v3000_sdf: Pathy,
) -> None:
    """An sdf block with shuffled bond lines has valid indices.

    That means that the number of bonds matches the counts line,
    and each index is unique and integer.
    Does not check that the indices are 1-indexed and in order.

    Args:
        single_shuffled_bond_lines_v3000_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock = next(parse_sdf(single_shuffled_bond_lines_v3000_sdf))
    ctable = mol.ctable()
    assert ctable.valid_bond_indices()


def test_single_shuffled_bond_lines_strict_invalid_bond_indices(
    single_shuffled_bond_lines_v3000_sdf: Pathy,
) -> None:
    """An sdf block with shuffled bond lines has invalid indices, when strict.

    That means that the number of bonds matches the counts line,
    and each index is unique and integer.
    However, the indices are not 1-indexed and in order.

    Args:
        single_shuffled_bond_lines_v3000_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock = next(parse_sdf(single_shuffled_bond_lines_v3000_sdf))
    ctable = mol.ctable()
    with pytest.raises(IndicesOutOfOrderError, match="bonds"):
        ctable.valid_bond_indices(strict=True)


def test_single_mol_v2000_sdf_renumber_indices_not_implemented(
    single_mol_v2000_sdf: Pathy,
) -> None:
    """A v2000 sdf block cannot yet have indices renumbered.

    Args:
        single_mol_v2000_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock = next(parse_sdf(single_mol_v2000_sdf))
    with pytest.raises(NotImplementedError):
        mol.renumber_indices()


def test_single_large_atom_index_strict_valid_atom_indices_renumbered(
    single_large_atom_index_v3000_sdf: Pathy,
) -> None:
    """An sdf block with large atom indices has strict valid indices, after renumbering.

    That means that the number of atoms matches the counts line,
    and each index is unique and integer.
    The indices are 1-indexed and in order after renumbering.

    Args:
        single_large_atom_index_v3000_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock = next(parse_sdf(single_large_atom_index_v3000_sdf))
    mol.renumber_indices()
    ctable = mol.ctable()
    assert ctable.valid_atom_indices(strict=True)


def test_single_large_bond_index_strict_valid_bond_indices_renumbered(
    single_large_bond_index_v3000_sdf: Pathy,
) -> None:
    """An sdf block with large bond indices has strict valid indices, after renumbering.

    That means that the number of atoms matches the counts line,
    and each index is unique and integer.
    The indices are 1-indexed and in order after renumbering.

    Args:
        single_large_bond_index_v3000_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock = next(parse_sdf(single_large_bond_index_v3000_sdf))
    mol.renumber_indices()
    ctable = mol.ctable()
    assert ctable.valid_bond_indices(strict=True)


def test_single_missing_atom_line_valid_atom_indices_after_renumbering(
    single_missing_atom_line_v3000_sdf: Pathy,
) -> None:
    """An sdf block with a missing atom line has stric valid indices, after renumbering.

    Args:
        single_missing_atom_line_v3000_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock = next(parse_sdf(single_missing_atom_line_v3000_sdf))
    mol.renumber_indices()
    ctable = mol.ctable()
    assert ctable.valid_atom_indices(strict=True)


def test_single_extra_atom_line_valid_atom_indices_after_renumbering(
    single_extra_atom_line_v3000_sdf: Pathy,
) -> None:
    """An sdf block with a extra atom line has strict valid indices, after renumbering.

    Args:
        single_extra_atom_line_v3000_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock = next(parse_sdf(single_extra_atom_line_v3000_sdf))
    mol.renumber_indices()
    ctable = mol.ctable()
    assert ctable.valid_atom_indices(strict=True)


def test_single_missing_bond_line_valid_atom_indices_after_renumbering(
    single_missing_bond_line_v3000_sdf: Pathy,
) -> None:
    """An sdf block with a missing bond line has strict valid indices after renumbering.

    Args:
        single_missing_bond_line_v3000_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock = next(parse_sdf(single_missing_bond_line_v3000_sdf))
    mol.renumber_indices()
    ctable = mol.ctable()
    assert ctable.valid_bond_indices(strict=True)


def test_single_extra_bond_line_valid_atom_indices_after_renumbering(
    single_extra_bond_line_v3000_sdf: Pathy,
) -> None:
    """An sdf block with a extra bond line has strict valid indices, after renumbering.

    Args:
        single_extra_bond_line_v3000_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock = next(parse_sdf(single_extra_bond_line_v3000_sdf))
    mol.renumber_indices()
    ctable = mol.ctable()
    assert ctable.valid_bond_indices(strict=True)


def test_single_duplicate_atom_index_invalid_atom_indices_even_renumbering(
    single_duplicate_atom_index_v3000_sdf: Pathy,
) -> None:
    """An sdf block with a duplicate atom index in bonds is invalid in renumbering.

    Args:
        single_duplicate_atom_index_v3000_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock = next(parse_sdf(single_duplicate_atom_index_v3000_sdf))
    with pytest.raises(IndicesDuplicateError, match="atom index mapping in bond"):
        mol.renumber_indices()


def test_single_duplicate_bond_index_valid_bond_indices_after_renumbering(
    single_duplicate_bond_index_v3000_sdf: Pathy,
) -> None:
    """An sdf block w/ duplicate bond index has strict valid indices after renumbering.

    Args:
        single_duplicate_bond_index_v3000_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock = next(parse_sdf(single_duplicate_bond_index_v3000_sdf))
    mol.renumber_indices()
    ctable = mol.ctable()
    assert ctable.valid_bond_indices(strict=True)


def test_single_shuffled_atom_lines_strict_valid_atom_indices_after_renumbering(
    single_shuffled_atom_lines_v3000_sdf: Pathy,
) -> None:
    """An sdf block w/ shuffled atom lines has strict valid indices, after renumbering.

    That means that the number of atoms matches the counts line,
    and each index is unique and integer.
    After renumbering, the indices are 1-indexed and in order.

    Args:
        single_shuffled_atom_lines_v3000_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock = next(parse_sdf(single_shuffled_atom_lines_v3000_sdf))
    mol.renumber_indices()
    ctable = mol.ctable()
    assert ctable.valid_atom_indices(strict=True)


def test_single_shuffled_bond_lines_strict_valid_bond_indices_after_renumbering(
    single_shuffled_bond_lines_v3000_sdf: Pathy,
) -> None:
    """An sdf block w/ shuffled bond lines has strict valid indices, after renumbering.

    That means that the number of bonds matches the counts line,
    and each index is unique and integer.
    After renumbering, the indices are 1-indexed and in order.

    Args:
        single_shuffled_bond_lines_v3000_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock = next(parse_sdf(single_shuffled_bond_lines_v3000_sdf))
    mol.renumber_indices()
    ctable = mol.ctable()
    assert ctable.valid_bond_indices(strict=True)


def test_single_mol_v3000_sdf_idempotent_renumber(
    single_mol_v3000_sdf: Pathy,
) -> None:
    """A V3000 sdf block generated by RDKit, renumbered, is the same as the input.

    Args:
        single_mol_v3000_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock = next(parse_sdf(single_mol_v3000_sdf))
    mol.renumber_indices()
    buffer = io.StringIO()
    mol.write(buffer)
    assert buffer.getvalue() == open(single_mol_v3000_sdf).read()


def test_single_mol_v3000_sdf_idempotent_renumber_with_splitlines(
    single_mol_v3000_sdf: Pathy,
) -> None:
    """A V3000 sdf block generated by RDKit, renumbered, is the same as the input.

    Args:
        single_mol_v3000_sdf: pytest fixture of a Path to the sdf
    """
    mol: SDBlock = next(
        SDBlock.from_lines(open(single_mol_v3000_sdf).read().splitlines())
    )
    mol.renumber_indices()
    buffer = io.StringIO()
    mol.write(buffer)
    assert buffer.getvalue() == open(single_mol_v3000_sdf).read()


V3000_HEADER = "t\n  sendoff\n\n  0  0  0     0  0            999 V3000\n"


def v3000_block(
    atom: str = "M  V30 1 C 0 0 0 0",
    bond: str = "M  V30 1 1 1 2",
    counts: str = "M  V30 COUNTS 2 1 0 0 0",
) -> SDBlock:
    """Build a two-atom V3000 block with one replaceable line of each kind.

    Args:
        atom: the first atom line
        bond: the only bond line
        counts: the V3000 counts line

    Returns:
        The parsed block.
    """
    text = (
        V3000_HEADER
        + f"M  V30 BEGIN CTAB\n{counts}\nM  V30 BEGIN ATOM\n{atom}\n"
        + "M  V30 2 O 0 0 0 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\n"
        + f"{bond}\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n$$$$\n"
    )
    return next(SDBlock.from_lines(text.splitlines()))


@pytest.mark.parametrize(
    "text, error",
    [
        ("t\nM  END\n$$$$\n", ValueError),
        (V3000_HEADER + "M  END\n$$$$\n", ValueError),
        ("t\n\n\n   \nM  END\n$$$$\n", IndexError),
        ("t\n\n\n  0  0  0     0  0            999 V9\nM  END\n$$$$\n", KeyError),
        ("t\n\n\n  x  0  0     0  0            999 V2000\nM  END\n$$$$\n", ValueError),
    ],
)
def test_malformed_ctable_header_raises_exception(
    text: str, error: type[Exception]
) -> None:
    """Malformed CTable headers raise ordinary exceptions, not Rust panics.

    Args:
        text: SDF text with a malformed CTable
        error: the expected exception type
    """
    block = next(SDBlock.from_lines(text.splitlines()))
    with pytest.raises(error):
        block.num_atoms()


def test_malformed_v3000_counts_raises_index_error() -> None:
    """A V3000 counts line without atom and bond fields raises IndexError."""
    with pytest.raises(IndexError):
        v3000_block(counts="M  V30 COUNTS 1").num_atoms()


@pytest.mark.parametrize(
    "atom, bond, error",
    [
        ("M  V30", "M  V30 1 1 1 2", IndexError),
        ("M  V30 x C 0 0 0 0", "M  V30 1 1 1 2", ValueError),
        ("M  V30 1 C 0 0 0 0", "M  V30", IndexError),
        ("M  V30 1 C 0 0 0 0", "M  V30 x 1 1 2", ValueError),
    ],
)
def test_malformed_index_lines_raise_exception(
    atom: str, bond: str, error: type[Exception]
) -> None:
    """Atom and bond lines without integer indices raise ordinary exceptions.

    Args:
        atom: the first atom line
        bond: the only bond line
        error: the expected exception type
    """
    ctable = v3000_block(atom=atom, bond=bond).ctable()
    with pytest.raises(error):
        ctable.valid_atom_indices()
        ctable.valid_bond_indices()


@pytest.mark.parametrize(
    "atom, bond, error",
    [
        ("M  V30", "M  V30 1 1 1 2", IndexError),
        ("M  V30 x C 0 0 0 0", "M  V30 1 1 1 2", ValueError),
        ("M  V30 1 C 0 0 0 0", "M  V30 1 1 1", IndexError),
        ("M  V30 1 C 0 0 0 0", "M  V30 1 1 1 x", ValueError),
        ("M  V30 1 C 0 0 0 0", "M  V30 1 1 1 9", IndexError),
    ],
)
def test_malformed_renumber_input_raises_exception(
    atom: str, bond: str, error: type[Exception]
) -> None:
    """Renumbering malformed or dangling lines raises and leaves the block intact.

    Args:
        atom: the first atom line
        bond: the only bond line
        error: the expected exception type
    """
    block = v3000_block(atom=atom, bond=bond)
    before = block.ctable().lines
    with pytest.raises(error):
        block.renumber_indices()
    assert block.ctable().lines == before


def test_renumber_slices_atom_lines_by_character() -> None:
    """Renumbering skips characters like Python slicing, even mid-UTF-8."""
    block = v3000_block(atom="a b 1 \u4e2d")
    block.renumber_indices()
    assert "M  V30 1" in block.ctable().lines
