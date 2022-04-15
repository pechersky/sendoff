# mypy: allow-untyped-decorators
"""Test sendoff parsing of metadata on generated SDFs."""
import io
from itertools import chain

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
    """An sdf block written with no trailing newlines has fewer lines than splitlines input.

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
    """An sdf block with a missing bond line has strict valid indices, after renumbering.

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
    """An sdf block with a duplicate bond index has strict valid indices after renumbering.

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
    """An sdf block with shuffled atom lines has strict valid indices, after renumbering.

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
    """An sdf block with shuffled bond lines has strict valid indices, after renumbering.

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
