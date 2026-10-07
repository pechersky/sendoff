"""Test rdkit parsing of metadata on generated SDFs."""

import os
from pathlib import Path

import pytest
from rdkit import Chem

from sendoff.sdblock import Pathy, parse_sdf
from tests.fixtures import INDEX_DEFECTS


def test_single_mol_read(single_mol_sdf: Pathy) -> None:
    """An sdf of a single molecule contains one molecule.

    Args:
        single_mol_sdf: pytest fixture of a Path to the sdf
    """
    mols = list(Chem.SDMolSupplier(str(single_mol_sdf)))
    assert len(mols) == 1


def test_double_mol_read(double_mol_sdf: Pathy) -> None:
    """An sdf of two molecule contains two molecules.

    Args:
        double_mol_sdf: pytest fixture of a Path to the sdf
    """
    mols = list(Chem.SDMolSupplier(str(double_mol_sdf)))
    assert len(mols) == 2


def test_titled_mol_read(single_titled_mol_sdf: Pathy) -> None:
    """An sdf of a single titled molecule contains it with the expected title.

    Args:
        single_titled_mol_sdf: pytest fixture of a Path to the sdf
    """
    mol = next(Chem.SDMolSupplier(str(single_titled_mol_sdf)))
    assert mol.GetProp("_Name") == "Title"


def test_delimiter_titled_mol_read(
    single_delimiter_titled_mol_sdf: Pathy,
) -> None:
    """An sdf of a delimiter titled molecule contains it with the expected title.

    Args:
        single_delimiter_titled_mol_sdf: pytest fixture of a Path to the sdf
    """
    mol = next(Chem.SDMolSupplier(str(single_delimiter_titled_mol_sdf)))
    assert mol.GetProp("_Name") == "$$$$"


def test_first_delimiter_titled_mol_read(
    double_first_delimiter_titled_mol_sdf: Pathy,
) -> None:
    """SDF with 1st molecule delimiter titled, containing with expected titles.

    Args:
        double_first_delimiter_titled_mol_sdf: pytest fixture of a Path to sdf
    """
    supp = Chem.SDMolSupplier(str(double_first_delimiter_titled_mol_sdf))
    mol = next(supp)
    other_mol = next(supp)
    assert mol.GetProp("_Name") == "$$$$"
    assert other_mol.GetProp("_Name") == "Title"


def test_second_delimiter_titled_mol_read(
    double_second_delimiter_titled_mol_sdf: Pathy,
) -> None:
    """SDF with 2nd molecule delimiter titled, containing with expected titles.

    Args:
        double_second_delimiter_titled_mol_sdf: pytest fixture of a Path to sdf
    """
    supp = Chem.SDMolSupplier(str(double_second_delimiter_titled_mol_sdf))
    other_mol = next(supp)
    mol = next(supp)
    assert other_mol.GetProp("_Name") == "Title"
    assert mol.GetProp("_Name") == "$$$$"


def test_single_record_mol_read(single_record_mol_sdf: Pathy) -> None:
    """An sdf of a delimiter record valued molecule contains the expected data.

    Args:
        single_record_mol_sdf: pytest fixture of a Path to the sdf
    """
    mol = next(Chem.SDMolSupplier(str(single_record_mol_sdf)))
    assert mol.GetProp("Record") == "Value"


def test_delimiter_record_mol_read(
    single_delimiter_record_mol_sdf: Pathy,
) -> None:
    """An sdf of a delimiter record valued molecule contains the expected data.

    Args:
        single_delimiter_record_mol_sdf: pytest fixture of a Path to the sdf
    """
    mol = next(Chem.SDMolSupplier(str(single_delimiter_record_mol_sdf)))
    assert mol.GetProp("Record") == "$$$$"


def test_multiline_record_mol_read(
    single_multiline_record_mol_sdf: Pathy,
) -> None:
    """An sdf of a multiline record valued molecule contains the expected data.

    Args:
        single_multiline_record_mol_sdf: pytest fixture of a Path to the sdf
    """
    mol = next(Chem.SDMolSupplier(str(single_multiline_record_mol_sdf)))
    assert mol.GetProp("Record") == "0\n1\n2\n3"


def test_empty_string_record_mol_read(
    single_empty_string_record_mol_sdf: Pathy,
) -> None:
    """An sdf of an empty string record valued molecule contains the expected data.

    Args:
        single_empty_string_record_mol_sdf: pytest fixture of a Path to the sdf
    """
    mol = next(Chem.SDMolSupplier(str(single_empty_string_record_mol_sdf)))
    assert mol.GetProp("Record") == ""


def test_multiline_record_name_mol_read(
    single_multiline_record_name_mol_sdf: Pathy,
) -> None:
    """An sdf of a multiline record name molecule cannot contain the expected data.

    Args:
        single_multiline_record_name_mol_sdf: fixture of a Path to the sdf
    """
    mol = next(Chem.SDMolSupplier(str(single_multiline_record_name_mol_sdf)))
    assert not mol.HasProp("Rec\nord")


def test_right_angle_bracket_record_name_mol_read_field_name_changed(
    single_right_angle_bracket_record_name_mol_sdf: Pathy,
) -> None:
    """An sdf with a record name with a ">" is not present in the read data.

    Args:
        single_right_angle_bracket_record_name_mol_sdf: fixture of a Path to the sdf
    """
    mol = next(Chem.SDMolSupplier(str(single_right_angle_bracket_record_name_mol_sdf)))
    assert not mol.HasProp("Rec>ord")


def test_right_angle_bracket_record_name_mol_read_data_at_modified_field(
    single_right_angle_bracket_record_name_mol_sdf: Pathy,
) -> None:
    """An sdf with a record name with a ">" contains data at a modified field name.

    Args:
        single_right_angle_bracket_record_name_mol_sdf: fixture of a Path to the sdf
    """
    mol = next(Chem.SDMolSupplier(str(single_right_angle_bracket_record_name_mol_sdf)))
    assert mol.GetProp("Rec") == "Value"


def test_0_atoms_mol_num_atoms(single_0_atoms_mol_sdf: Pathy) -> None:
    """An sdf block with 0 atoms written parses as 0 atoms.

    Args:
        single_0_atoms_mol_sdf: pytest fixture of a Path to the sdf
    """
    mol = next(Chem.SDMolSupplier(str(single_0_atoms_mol_sdf)))
    assert mol.GetNumAtoms() == 0


def test_star_atom_mol_num_atoms(single_star_atom_mol_sdf: Pathy) -> None:
    """An sdf block with a single star atom written parses as 1 atom.

    Args:
        single_star_atom_mol_sdf: pytest fixture of a Path to the sdf
    """
    mol = next(Chem.SDMolSupplier(str(single_star_atom_mol_sdf)))
    assert mol.GetNumAtoms() == 1


def test_999_atoms_mol_num_atoms(single_999_atoms_mol_sdf: Pathy) -> None:
    """An sdf block with 999 atoms written parses as 999 atoms.

    Args:
        single_999_atoms_mol_sdf: pytest fixture of a Path to the sdf
    """
    mol = next(Chem.SDMolSupplier(str(single_999_atoms_mol_sdf)))
    assert mol.GetNumAtoms() == 999


def test_1001_atoms_mol_num_atoms(single_1001_atoms_mol_sdf: Pathy) -> None:
    """An sdf block with 1001 atoms written parses as 1001 atoms.

    Args:
        single_1001_atoms_mol_sdf: pytest fixture of a Path to the sdf
    """
    mol = next(Chem.SDMolSupplier(str(single_1001_atoms_mol_sdf)))
    assert mol.GetNumAtoms() == 1001


def test_0_atoms_mol_num_bonds(single_0_atoms_mol_sdf: Pathy) -> None:
    """An sdf block with 0 atoms written parses as 0 bonds.

    Args:
        single_0_atoms_mol_sdf: pytest fixture of a Path to the sdf
    """
    mol = next(Chem.SDMolSupplier(str(single_0_atoms_mol_sdf)))
    assert mol.GetNumBonds() == 0


def test_star_atom_mol_num_bonds(single_star_atom_mol_sdf: Pathy) -> None:
    """An sdf block with a single star atom written parses as 0 bonds.

    Args:
        single_star_atom_mol_sdf: pytest fixture of a Path to the sdf
    """
    mol = next(Chem.SDMolSupplier(str(single_star_atom_mol_sdf)))
    assert mol.GetNumBonds() == 0


def test_999_atoms_mol_num_bonds(single_999_atoms_mol_sdf: Pathy) -> None:
    """An sdf block with 999 atoms written parses as 998 bonds.

    Args:
        single_999_atoms_mol_sdf: pytest fixture of a Path to the sdf
    """
    mol = next(Chem.SDMolSupplier(str(single_999_atoms_mol_sdf)))
    assert mol.GetNumBonds() == 998


def test_1001_atoms_mol_num_bonds(single_1001_atoms_mol_sdf: Pathy) -> None:
    """An sdf block with 1001 atoms written parses as 1000 bonds.

    Args:
        single_1001_atoms_mol_sdf: pytest fixture of a Path to the sdf
    """
    mol = next(Chem.SDMolSupplier(str(single_1001_atoms_mol_sdf)))
    assert mol.GetNumBonds() == 1000


@pytest.mark.parametrize(
    "defect, readable",
    [
        ("large_atom_index", True),
        ("large_bond_index", True),
        ("missing_atom_line", False),
        ("extra_atom_line", False),
        ("missing_bond_line", False),
        ("extra_bond_line", False),
        ("duplicate_atom_index", True),
        ("duplicate_bond_index", True),
        ("shuffled_atom_lines", True),
        ("shuffled_bond_lines", True),
    ],
)
def test_index_defects_mol_read(defect: str, readable: bool) -> None:
    """Large, duplicate or unordered indices read, but wrong counts do not.

    Args:
        defect: the name of the defect in INDEX_DEFECTS
        readable: whether RDKit reads the molecule
    """
    supplier = Chem.SDMolSupplier()
    supplier.SetData(INDEX_DEFECTS[defect].sdf())
    mol = next(supplier)
    if readable:
        assert (mol.GetNumAtoms(), mol.GetNumBonds()) == (9, 9)
    else:
        assert mol is None


def test_molblock_round_trips_through_rdkit(double_mol_sdf: Pathy) -> None:
    """Each block's molblock parses to the molecule RDKit reads from the file.

    Args:
        double_mol_sdf: pytest fixture of a Path to the sdf
    """
    expected = [
        Chem.MolToSmiles(mol) for mol in Chem.SDMolSupplier(str(double_mol_sdf))
    ]
    for path in (str(double_mol_sdf), os.fsencode(double_mol_sdf)):
        smiles = [
            Chem.MolToSmiles(Chem.MolFromMolBlock(block.molblock))
            for block in parse_sdf(path)
        ]
        assert smiles == expected


def test_molblock_setter_swaps_structures_from_rdkit(
    double_mol_sdf: Pathy, tmp_path: Path
) -> None:
    """Assigning an RDKit molblock keeps records, and RDKit reads the result.

    Args:
        double_mol_sdf: pytest fixture of a Path to the sdf
        tmp_path: pytest fixture for writing files to a temp directory
    """
    outpath = tmp_path / "swapped.sdf"
    with open(outpath, "w") as handle:
        for block in parse_sdf(str(double_mol_sdf)):
            mol = Chem.MolFromSmiles("CCO")
            mol.SetProp("_Name", f"{block.title!s}-swapped")
            block.molblock = Chem.MolToMolBlock(mol)
            block.set_record("Original", block.title)
            block.write(handle)
    swapped = list(Chem.SDMolSupplier(str(outpath)))
    assert [Chem.MolToSmiles(mol) for mol in swapped] == ["CCO", "CCO"]
    assert all(mol.GetProp("_Name").endswith("-swapped") for mol in swapped)
    assert all(mol.GetProp("Original") == mol.GetProp("_Name") for mol in swapped)
