"""Test rdkit parsing of metadata on generated SDFs."""
from rdkit import Chem

from sendoff.sdblock import Pathy


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


def test_single_large_atom_index_valid_atom_indices(
    single_large_atom_index_v3000_sdf: Pathy,
) -> None:
    """An sdf block with large atom indices has valid indices.

    That means that the number of atoms matches the counts line,
    and each index is unique and integer. That means RDKit does not complain.
    Does not check that the indices are 1-indexed and in order.

    Args:
        single_large_atom_index_v3000_sdf: pytest fixture of a Path to the sdf
    """
    mol = next(Chem.SDMolSupplier(str(single_large_atom_index_v3000_sdf)))
    assert mol.GetNumAtoms() == 9


def test_single_large_bond_index_valid_bond_indices(
    single_large_bond_index_v3000_sdf: Pathy,
) -> None:
    """An sdf block with large bond indices has valid indices.

    That means that the number of bonds matches the counts line,
    and each index is unique and integer. That means RDKit does not complain.
    Does not check that the indices are 1-indexed and in order.

    Args:
        single_large_bond_index_v3000_sdf: pytest fixture of a Path to the sdf
    """
    mol = next(Chem.SDMolSupplier(str(single_large_bond_index_v3000_sdf)))
    assert mol.GetNumBonds() == 9


def test_single_missing_atom_line_invalid_atom_indices(
    single_missing_atom_line_v3000_sdf: Pathy,
) -> None:
    """An sdf block with a missing atom line has invalid indices.

    That means that the number of atoms is fewer than the counts line,
    and RDKit won't read the molecule.

    Args:
        single_missing_atom_line_v3000_sdf: pytest fixture of a Path to the sdf
    """
    mol = next(Chem.SDMolSupplier(str(single_missing_atom_line_v3000_sdf)))
    assert mol is None


def test_single_extra_atom_line_invalid_atom_indices(
    single_extra_atom_line_v3000_sdf: Pathy,
) -> None:
    """An sdf block with a extra atom line has invalid indices.

    That means that the number of atoms is more than the counts line,
    and RDKit won't read the molecule.

    Args:
        single_extra_atom_line_v3000_sdf: pytest fixture of a Path to the sdf
    """
    mol = next(Chem.SDMolSupplier(str(single_extra_atom_line_v3000_sdf)))
    assert mol is None


def test_single_missing_bond_line_invalid_atom_indices(
    single_missing_bond_line_v3000_sdf: Pathy,
) -> None:
    """An sdf block with a missing bond line has invalid indices.

    That means that the number of bonds is fewer than the counts line,
    and RDKit won't read the molecule.

    Args:
        single_missing_bond_line_v3000_sdf: pytest fixture of a Path to the sdf
    """
    mol = next(Chem.SDMolSupplier(str(single_missing_bond_line_v3000_sdf)))
    assert mol is None


def test_single_extra_bond_line_invalid_atom_indices(
    single_extra_bond_line_v3000_sdf: Pathy,
) -> None:
    """An sdf block with a extra bond line has invalid indices.

    That means that the number of bonds is more than the counts line,
    and RDKit won't read the molecule.

    Args:
        single_extra_bond_line_v3000_sdf: pytest fixture of a Path to the sdf
    """
    mol = next(Chem.SDMolSupplier(str(single_extra_bond_line_v3000_sdf)))
    assert mol is None


def test_single_duplicate_atom_index_valid_atom_indices(
    single_duplicate_atom_index_v3000_sdf: Pathy,
) -> None:
    """An sdf block with a duplicate atom index has valid indices.

    That means that the number of atoms is weird,
    yet RDKit will read the molecule.

    Args:
        single_duplicate_atom_index_v3000_sdf: pytest fixture of a Path to the sdf
    """
    mol = next(Chem.SDMolSupplier(str(single_duplicate_atom_index_v3000_sdf)))
    assert mol.GetNumAtoms() == 9


def test_single_duplicate_bond_index_valid_bond_indices(
    single_duplicate_bond_index_v3000_sdf: Pathy,
) -> None:
    """An sdf block with a duplicate bond index has valid indices.

    That means that the number of bonds is weird,
    yet RDKit will read the molecule.

    Args:
        single_duplicate_bond_index_v3000_sdf: pytest fixture of a Path to the sdf
    """
    mol = next(Chem.SDMolSupplier(str(single_duplicate_bond_index_v3000_sdf)))
    assert mol.GetNumBonds() == 9


def test_single_shuffled_atom_lines_valid_atom_indices(
    single_shuffled_atom_lines_v3000_sdf: Pathy,
) -> None:
    """An sdf block with shuffled atom lines has valid indices.

    That means that the number of atoms is weird,
    yet RDKit will read the molecule.

    Args:
        single_shuffled_atom_lines_v3000_sdf: pytest fixture of a Path to the sdf
    """
    mol = next(Chem.SDMolSupplier(str(single_shuffled_atom_lines_v3000_sdf)))
    assert mol.GetNumAtoms() == 9


def test_single_shuffled_bond_lines_valid_bond_indices(
    single_shuffled_bond_lines_v3000_sdf: Pathy,
) -> None:
    """An sdf block with shuffled bond lines has valid indices.

    That means that the number of bonds is weird,
    yet RDKit will read the molecule.

    Args:
        single_shuffled_bond_lines_v3000_sdf: pytest fixture of a Path to the sdf
    """
    mol = next(Chem.SDMolSupplier(str(single_shuffled_bond_lines_v3000_sdf)))
    assert mol.GetNumBonds() == 9
