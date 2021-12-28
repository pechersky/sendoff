# mypy: allow-untyped-decorators
"""Test sendoff parsing of metadata on generated SDFs."""
import io
from itertools import chain

import pytest

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


def test_num_atoms_in_double_mol_sdf(double_mol_sdf: Pathy) -> None:
    """The position in the SDF file is how many atoms are in the record.

    Args:
        double_mol_sdf: pytest fixture of a Path to the sdf
    """
    for idx, mol in enumerate(parse_sdf(double_mol_sdf)):
        assert mol.num_atoms() == idx + 1
    assert idx == 1


def test_num_atoms_in_double_mol_v3000(double_mol_v3000: Pathy) -> None:
    """The position in the SDF file is how many atoms are in the record.

    Args:
        double_mol_v3000: pytest fixture of a Path to the sdf
    """
    for idx, mol in enumerate(parse_sdf(double_mol_v3000)):
        assert mol.num_atoms() == idx + 1
    assert idx == 1


def test_num_atoms_999_atom_mol(test_999_atom_mol: Pathy) -> None:
    """All molecules in the file contain 999 atoms. And both V2000 and V3000.

    Args:
        test_999_atom_mol: pytest fixture of a Path to the sdf
    """
    data = open(test_999_atom_mol).read()
    assert "V2000" in data
    assert "V3000" in data
    for idx, mol in enumerate(parse_sdf(test_999_atom_mol)):
        assert mol.num_atoms() == 999
    assert idx == 1


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
