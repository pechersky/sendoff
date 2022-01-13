# mypy: allow-untyped-decorators
"""Provide fixtures of SDFs with well- or misbehaved metadata."""
from pathlib import Path
from typing import TYPE_CHECKING

import pytest
from rdkit import Chem

from sendoff.sdblock import parse_sdf

if TYPE_CHECKING:
    from _pytest.fixtures import FixtureRequest as __FixtureRequest

    # This is to provide a proper type; we are subclassing Any, which causes mypy warn
    class FixtureRequest(__FixtureRequest):  # type: ignore
        """FixtureRequest wrapper class to provide param attribute."""

        param: str

else:
    from _pytest.fixtures import FixtureRequest


@pytest.fixture(params=["rdkitV2000", "rdkitV3000"])
def single_mol_sdf(request: FixtureRequest, tmp_path: Path) -> Path:
    """Write a single molecule with no metadata to an sdf.

    Args:
        request: pytest fixture configuration handling param passing
        tmp_path: pytest fixture for writing files to a temp directory

    Returns:
        Path to the sdf
    """
    mol = Chem.MolFromSmiles("C")
    outpath = tmp_path / "input.sdf"
    with Chem.SDWriter(str(outpath)) as sdw:
        sdw.SetForceV3000(request.param == "rdkitV3000")
        sdw.write(mol)
    return outpath


@pytest.fixture(params=["rdkitV2000", "rdkitV3000"])
def double_mol_sdf(request: FixtureRequest, tmp_path: Path) -> Path:
    """Write two molecules with no metadata to an sdf.

    Args:
        request: pytest fixture configuration handling param passing
        tmp_path: pytest fixture for writing files to a temp directory

    Returns:
        Path to the sdf
    """
    methane = Chem.MolFromSmiles("C")
    ethane = Chem.MolFromSmiles("CC")
    outpath = tmp_path / "input.sdf"
    with Chem.SDWriter(str(outpath)) as sdw:
        sdw.SetForceV3000(request.param == "rdkitV3000")
        sdw.write(methane)
        sdw.write(ethane)
    return outpath


@pytest.fixture(params=["rdkitV2000", "rdkitV3000"])
def single_0_atoms_mol_sdf(request: FixtureRequest, tmp_path: Path) -> Path:
    """Write a single molecule with 0 atoms and no metadata to an sdf.

    This relies on the fact that RDKit does not write hydrogens by default.

    Args:
        request: pytest fixture configuration handling param passing
        tmp_path: pytest fixture for writing files to a temp directory

    Returns:
        Path to the sdf
    """
    mol = Chem.MolFromSmiles("C" * 0)
    outpath = tmp_path / "input.sdf"
    with Chem.SDWriter(str(outpath)) as sdw:
        sdw.SetForceV3000(request.param == "rdkitV3000")
        sdw.write(mol)
    return outpath


@pytest.fixture(params=["rdkitV2000", "rdkitV3000"])
def single_star_atom_mol_sdf(request: FixtureRequest, tmp_path: Path) -> Path:
    """Write a single molecule with a single star atom and no metadata to an sdf.

    This relies on the fact that RDKit does not write hydrogens by default.

    Args:
        request: pytest fixture configuration handling param passing
        tmp_path: pytest fixture for writing files to a temp directory

    Returns:
        Path to the sdf
    """
    mol = Chem.MolFromSmiles("*")
    outpath = tmp_path / "input.sdf"
    with Chem.SDWriter(str(outpath)) as sdw:
        sdw.SetForceV3000(request.param == "rdkitV3000")
        sdw.write(mol)
    return outpath


@pytest.fixture(params=["rdkitV2000", "rdkitV3000"])
def single_999_atoms_mol_sdf(request: FixtureRequest, tmp_path: Path) -> Path:
    """Write a single molecule with 999 atoms and no metadata to an sdf.

    This relies on the fact that RDKit does not write hydrogens by default.

    Args:
        request: pytest fixture configuration handling param passing
        tmp_path: pytest fixture for writing files to a temp directory

    Returns:
        Path to the sdf
    """
    mol = Chem.MolFromSmiles("C" * 999)
    outpath = tmp_path / "input.sdf"
    with Chem.SDWriter(str(outpath)) as sdw:
        sdw.SetForceV3000(request.param == "rdkitV3000")
        sdw.write(mol)
    return outpath


@pytest.fixture
def single_1001_atoms_mol_sdf(tmp_path: Path) -> Path:
    """Write a single molecule with 1001 atoms and no metadata to an sdf.

    This relies on the fact that RDKit does not write hydrogens by default,
    and that RDKit falls back to V3000 for molecules of more than 999 atoms.

    Args:
        tmp_path: pytest fixture for writing files to a temp directory

    Returns:
        Path to the sdf
    """
    mol = Chem.MolFromSmiles("C" * 1001)
    outpath = tmp_path / "input.sdf"
    with Chem.SDWriter(str(outpath)) as sdw:
        sdw.write(mol)
    return outpath


@pytest.fixture(params=["rdkitV2000", "rdkitV3000"])
def single_titled_mol_sdf(request: FixtureRequest, tmp_path: Path) -> Path:
    """Write a single molecule with only title metadata to an sdf.

    Args:
        request: pytest fixture configuration handling param passing
        tmp_path: pytest fixture for writing files to a temp directory

    Returns:
        Path to the sdf
    """
    mol = Chem.MolFromSmiles("C")
    mol.SetProp("_Name", "Title")
    outpath = tmp_path / "input.sdf"
    with Chem.SDWriter(str(outpath)) as sdw:
        sdw.SetForceV3000(request.param == "rdkitV3000")
        sdw.write(mol)
    return outpath


@pytest.fixture(params=["rdkitV2000", "rdkitV3000"])
def single_delimiter_titled_mol_sdf(request: FixtureRequest, tmp_path: Path) -> Path:
    """Write a single molecule with delimiter title metadata to an sdf.

    Args:
        request: pytest fixture configuration handling param passing
        tmp_path: pytest fixture for writing files to a temp directory

    Returns:
        Path to the sdf
    """
    mol = Chem.MolFromSmiles("C")
    mol.SetProp("_Name", "$$$$")
    outpath = tmp_path / "input.sdf"
    with Chem.SDWriter(str(outpath)) as sdw:
        sdw.SetForceV3000(request.param == "rdkitV3000")
        sdw.write(mol)
    return outpath


@pytest.fixture(params=["rdkitV2000", "rdkitV3000"])
def double_first_delimiter_titled_mol_sdf(
    request: FixtureRequest, tmp_path: Path
) -> Path:
    """Write two molecules with the first as delimiter title metadata to an sdf.

    Args:
        request: pytest fixture configuration handling param passing
        tmp_path: pytest fixture for writing files to a temp directory

    Returns:
        Path to the sdf
    """
    mol = Chem.MolFromSmiles("C")
    mol.SetProp("_Name", "$$$$")
    other_mol = Chem.MolFromSmiles("CC")
    other_mol.SetProp("_Name", "Title")
    outpath = tmp_path / "input.sdf"
    with Chem.SDWriter(str(outpath)) as sdw:
        sdw.SetForceV3000(request.param == "rdkitV3000")
        sdw.write(mol)
        sdw.write(other_mol)
    return outpath


@pytest.fixture(params=["rdkitV2000", "rdkitV3000"])
def double_second_delimiter_titled_mol_sdf(
    request: FixtureRequest, tmp_path: Path
) -> Path:
    """Write two molecules with the second as delimiter title metadata to an sdf.

    Args:
        request: pytest fixture configuration handling param passing
        tmp_path: pytest fixture for writing files to a temp directory

    Returns:
        Path to the sdf
    """
    other_mol = Chem.MolFromSmiles("CC")
    other_mol.SetProp("_Name", "Title")
    mol = Chem.MolFromSmiles("C")
    mol.SetProp("_Name", "$$$$")
    outpath = tmp_path / "input.sdf"
    with Chem.SDWriter(str(outpath)) as sdw:
        sdw.SetForceV3000(request.param == "rdkitV3000")
        sdw.write(other_mol)
        sdw.write(mol)
    return outpath


@pytest.fixture(params=["rdkitV2000", "rdkitV3000", "sendoff"])
def single_record_mol_sdf(request: FixtureRequest, tmp_path: Path) -> Path:
    """Write a single molecule with only record value metadata to an sdf.

    Args:
        request: pytest fixture configuration handling param passing
        tmp_path: pytest fixture for writing files to a temp directory

    Returns:
        Path to the sdf
    """
    record, value = ("Record", "Value")
    mol = Chem.MolFromSmiles("C")
    outpath = tmp_path / "input.sdf"
    with Chem.SDWriter(str(outpath)) as sdw:
        sdw.write(mol)
    if request.param.startswith("rdkit"):
        mols = list(Chem.SDMolSupplier(str(outpath)))
        with Chem.SDWriter(str(outpath)) as outh:
            sdw.SetForceV3000(request.param == "rdkitV3000")
            for mol in mols:
                mol.SetProp(record, value)
                outh.write(mol)
    elif request.param == "sendoff":
        mols = list(parse_sdf(outpath))
        with open(outpath, "w") as outh:
            for mol in mols:
                mol.append_record(record, value)
                mol.write(outh)
    return outpath


@pytest.fixture(params=["rdkitV2000", "rdkitV3000", "sendoff"])
def single_delimiter_record_mol_sdf(request: FixtureRequest, tmp_path: Path) -> Path:
    """Write a single molecule with delimiter record value metadata to an sdf.

    Args:
        request: pytest fixture configuration handling param passing
        tmp_path: pytest fixture for writing files to a temp directory

    Returns:
        Path to the sdf
    """
    record, value = ("Record", "$$$$")
    mol = Chem.MolFromSmiles("C")
    outpath = tmp_path / "input.sdf"
    with Chem.SDWriter(str(outpath)) as sdw:
        sdw.write(mol)
    if request.param.startswith("rdkit"):
        mols = list(Chem.SDMolSupplier(str(outpath)))
        with Chem.SDWriter(str(outpath)) as outh:
            sdw.SetForceV3000(request.param == "rdkitV3000")
            for mol in mols:
                mol.SetProp(record, value)
                outh.write(mol)
    elif request.param == "sendoff":
        mols = list(parse_sdf(outpath))
        with open(outpath, "w") as outh:
            for mol in mols:
                mol.append_record(record, value)
                mol.write(outh)
    return outpath


@pytest.fixture(params=["rdkitV2000", "rdkitV3000", "sendoff"])
def single_multiline_record_mol_sdf(request: FixtureRequest, tmp_path: Path) -> Path:
    """Write a single molecule with a multiline record value metadata to an sdf.

    Args:
        request: pytest fixture configuration handling param passing
        tmp_path: pytest fixture for writing files to a temp directory

    Returns:
        Path to the sdf
    """
    record, value = ("Record", str.join("\n", map(str, [0, 1, 2, 3])))
    mol = Chem.MolFromSmiles("C")
    outpath = tmp_path / "input.sdf"
    with Chem.SDWriter(str(outpath)) as sdw:
        sdw.write(mol)
    if request.param.startswith("rdkit"):
        mols = list(Chem.SDMolSupplier(str(outpath)))
        with Chem.SDWriter(str(outpath)) as outh:
            sdw.SetForceV3000(request.param == "rdkitV3000")
            for mol in mols:
                mol.SetProp(record, value)
                outh.write(mol)
    elif request.param == "sendoff":
        mols = list(parse_sdf(outpath))
        with open(outpath, "w") as outh:
            for mol in mols:
                mol.append_record(record, value)
                mol.write(outh)
    return outpath


@pytest.fixture(params=["rdkitV2000", "rdkitV3000", "sendoff"])
def single_empty_string_record_mol_sdf(request: FixtureRequest, tmp_path: Path) -> Path:
    """Write a single molecule with an empty string record value to an sdf.

    Args:
        request: pytest fixture configuration handling param passing
        tmp_path: pytest fixture for writing files to a temp directory

    Returns:
        Path to the sdf
    """
    record, value = ("Record", "")
    mol = Chem.MolFromSmiles("C")
    outpath = tmp_path / "input.sdf"
    with Chem.SDWriter(str(outpath)) as sdw:
        sdw.write(mol)
    if request.param.startswith("rdkit"):
        mols = list(Chem.SDMolSupplier(str(outpath)))
        with Chem.SDWriter(str(outpath)) as outh:
            sdw.SetForceV3000(request.param == "rdkitV3000")
            for mol in mols:
                mol.SetProp(record, value)
                outh.write(mol)
    elif request.param == "sendoff":
        mols = list(parse_sdf(outpath))
        with open(outpath, "w") as outh:
            for mol in mols:
                mol.append_record(record, value)
                mol.write(outh)
    return outpath


@pytest.fixture(params=["rdkitV2000", "rdkitV3000", "sendoff"])
def single_multiline_record_name_mol_sdf(
    request: FixtureRequest, tmp_path: Path
) -> Path:
    """Write a single molecule with a multiline string record name to an sdf.

    Args:
        request: pytest fixture configuration handling param passing
        tmp_path: pytest fixture for writing files to a temp directory

    Returns:
        Path to the sdf
    """
    record, value = ("Rec\nord", "Value")
    mol = Chem.MolFromSmiles("C")
    outpath = tmp_path / "input.sdf"
    with Chem.SDWriter(str(outpath)) as sdw:
        sdw.write(mol)
    if request.param.startswith("rdkit"):
        mols = list(Chem.SDMolSupplier(str(outpath)))
        with Chem.SDWriter(str(outpath)) as outh:
            sdw.SetForceV3000(request.param == "rdkitV3000")
            for mol in mols:
                mol.SetProp(record, value)
                outh.write(mol)
    elif request.param == "sendoff":
        mols = list(parse_sdf(outpath))
        with open(outpath, "w") as outh:
            for mol in mols:
                mol.append_record(record, value)
                mol.write(outh)
    return outpath


@pytest.fixture(params=["rdkitV2000", "rdkitV3000", "sendoff"])
def single_right_angle_bracket_record_name_mol_sdf(
    request: FixtureRequest, tmp_path: Path
) -> Path:
    """Write a single molecule with record name with a right angle bracket to an sdf.

    Args:
        request: pytest fixture configuration handling param passing
        tmp_path: pytest fixture for writing files to a temp directory

    Returns:
        Path to the sdf
    """
    record, value = ("Rec>ord", "Value")
    mol = Chem.MolFromSmiles("C")
    outpath = tmp_path / "input.sdf"
    with Chem.SDWriter(str(outpath)) as sdw:
        sdw.write(mol)
    if request.param.startswith("rdkit"):
        mols = list(Chem.SDMolSupplier(str(outpath)))
        with Chem.SDWriter(str(outpath)) as outh:
            sdw.SetForceV3000(request.param == "rdkitV3000")
            for mol in mols:
                mol.SetProp(record, value)
                outh.write(mol)
    elif request.param == "sendoff":
        mols = list(parse_sdf(outpath))
        with open(outpath, "w") as outh:
            for mol in mols:
                mol.append_record(record, value)
                mol.write(outh)
    return outpath


@pytest.fixture()
def single_large_atom_index_v3000_sdf(tmp_path: Path) -> Path:
    """Write a single molecule to a v3000 sdf with large atom indices.

    Args:
        tmp_path: pytest fixture for writing files to a temp directory

    Returns:
        Path to the sdf
    """
    sdf_text = """


  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 9 9 0 0 0
M  V30 BEGIN ATOM
M  V30 2101315825 C 87.71 -95.64 0 0
M  V30 2101315840 C 87.71 -81.29 0 0
M  V30 2101315860 C 100.18 -74.09 0 0
M  V30 2101315861 C 100.18 -59.69 0 0
M  V30 2101315862 C 87.71 -52.49 0 0
M  V30 2101315863 C 75.24 -59.69 0 0
M  V30 2101315864 C 75.24 -74.09 0 0
M  V30 2101315865 C 87.71 -38.09 0 0
M  V30 2101315866 O 100.18 -30.89 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2101315825 2101315840
M  V30 2 1 2101315840 2101315860
M  V30 3 2 2101315860 2101315861
M  V30 4 1 2101315861 2101315862
M  V30 5 2 2101315862 2101315863
M  V30 6 1 2101315863 2101315864
M  V30 7 2 2101315864 2101315840
M  V30 8 1 2101315862 2101315865
M  V30 9 2 2101315865 2101315866
M  V30 END BOND
M  V30 END CTAB
M  END
$$$$

"""
    outpath = tmp_path / "input.sdf"
    with open(outpath, "w") as outh:
        outh.write(sdf_text)
    return outpath


@pytest.fixture()
def single_large_bond_index_v3000_sdf(tmp_path: Path) -> Path:
    """Write a single molecule to a v3000 sdf with large bond indices.

    Args:
        tmp_path: pytest fixture for writing files to a temp directory

    Returns:
        Path to the sdf
    """
    sdf_text = """


  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 9 9 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 87.71 -95.64 0 0
M  V30 2 C 87.71 -81.29 0 0
M  V30 3 C 100.18 -74.09 0 0
M  V30 4 C 100.18 -59.69 0 0
M  V30 5 C 87.71 -52.49 0 0
M  V30 6 C 75.24 -59.69 0 0
M  V30 7 C 75.24 -74.09 0 0
M  V30 8 C 87.71 -38.09 0 0
M  V30 9 O 100.18 -30.89 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 2101315886 1 1 2
M  V30 2101315911 1 2 3
M  V30 2101315912 2 3 4
M  V30 2101315913 1 4 5
M  V30 2101315914 2 5 6
M  V30 2101315915 1 6 7
M  V30 2101315916 2 7 2
M  V30 2101315917 1 5 8
M  V30 2101315918 2 8 9
M  V30 END BOND
M  V30 END CTAB
M  END
$$$$
"""
    outpath = tmp_path / "input.sdf"
    with open(outpath, "w") as outh:
        outh.write(sdf_text)
    return outpath


@pytest.fixture()
def single_missing_atom_line_v3000_sdf(tmp_path: Path) -> Path:
    """Write a single molecule to a v3000 sdf with a missing atom line.

    Args:
        tmp_path: pytest fixture for writing files to a temp directory

    Returns:
        Path to the sdf
    """
    sdf_text = """


  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 10 9 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 87.71 -95.64 0 0
M  V30 2 C 87.71 -81.29 0 0
M  V30 3 C 100.18 -74.09 0 0
M  V30 4 C 100.18 -59.69 0 0
M  V30 5 C 87.71 -52.49 0 0
M  V30 6 C 75.24 -59.69 0 0
M  V30 7 C 75.24 -74.09 0 0
M  V30 8 C 87.71 -38.09 0 0
M  V30 9 O 100.18 -30.89 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 2 3 4
M  V30 4 1 4 5
M  V30 5 2 5 6
M  V30 6 1 6 7
M  V30 7 2 7 2
M  V30 8 1 5 8
M  V30 9 2 8 9
M  V30 END BOND
M  V30 END CTAB
M  END
$$$$
"""
    outpath = tmp_path / "input.sdf"
    with open(outpath, "w") as outh:
        outh.write(sdf_text)
    return outpath


@pytest.fixture()
def single_extra_atom_line_v3000_sdf(tmp_path: Path) -> Path:
    """Write a single molecule to a v3000 sdf with an extra atom line.

    Args:
        tmp_path: pytest fixture for writing files to a temp directory

    Returns:
        Path to the sdf
    """
    sdf_text = """


  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 8 9 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 87.71 -95.64 0 0
M  V30 2 C 87.71 -81.29 0 0
M  V30 3 C 100.18 -74.09 0 0
M  V30 4 C 100.18 -59.69 0 0
M  V30 5 C 87.71 -52.49 0 0
M  V30 6 C 75.24 -59.69 0 0
M  V30 7 C 75.24 -74.09 0 0
M  V30 8 C 87.71 -38.09 0 0
M  V30 9 O 100.18 -30.89 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 2 3 4
M  V30 4 1 4 5
M  V30 5 2 5 6
M  V30 6 1 6 7
M  V30 7 2 7 2
M  V30 8 1 5 8
M  V30 9 2 8 9
M  V30 END BOND
M  V30 END CTAB
M  END
$$$$
"""
    outpath = tmp_path / "input.sdf"
    with open(outpath, "w") as outh:
        outh.write(sdf_text)
    return outpath


@pytest.fixture()
def single_missing_bond_line_v3000_sdf(tmp_path: Path) -> Path:
    """Write a single molecule to a v3000 sdf with a missing bond line.

    Args:
        tmp_path: pytest fixture for writing files to a temp directory

    Returns:
        Path to the sdf
    """
    sdf_text = """


  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 9 10 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 87.71 -95.64 0 0
M  V30 2 C 87.71 -81.29 0 0
M  V30 3 C 100.18 -74.09 0 0
M  V30 4 C 100.18 -59.69 0 0
M  V30 5 C 87.71 -52.49 0 0
M  V30 6 C 75.24 -59.69 0 0
M  V30 7 C 75.24 -74.09 0 0
M  V30 8 C 87.71 -38.09 0 0
M  V30 9 O 100.18 -30.89 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 2 3 4
M  V30 4 1 4 5
M  V30 5 2 5 6
M  V30 6 1 6 7
M  V30 7 2 7 2
M  V30 8 1 5 8
M  V30 9 2 8 9
M  V30 END BOND
M  V30 END CTAB
M  END
$$$$
"""
    outpath = tmp_path / "input.sdf"
    with open(outpath, "w") as outh:
        outh.write(sdf_text)
    return outpath


@pytest.fixture()
def single_extra_bond_line_v3000_sdf(tmp_path: Path) -> Path:
    """Write a single molecule to a v3000 sdf with an extra bond line.

    Args:
        tmp_path: pytest fixture for writing files to a temp directory

    Returns:
        Path to the sdf
    """
    sdf_text = """


  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 9 8 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 87.71 -95.64 0 0
M  V30 2 C 87.71 -81.29 0 0
M  V30 3 C 100.18 -74.09 0 0
M  V30 4 C 100.18 -59.69 0 0
M  V30 5 C 87.71 -52.49 0 0
M  V30 6 C 75.24 -59.69 0 0
M  V30 7 C 75.24 -74.09 0 0
M  V30 8 C 87.71 -38.09 0 0
M  V30 9 O 100.18 -30.89 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 2 3 4
M  V30 4 1 4 5
M  V30 5 2 5 6
M  V30 6 1 6 7
M  V30 7 2 7 2
M  V30 8 1 5 8
M  V30 9 2 8 9
M  V30 END BOND
M  V30 END CTAB
M  END
$$$$
"""
    outpath = tmp_path / "input.sdf"
    with open(outpath, "w") as outh:
        outh.write(sdf_text)
    return outpath


@pytest.fixture()
def single_duplicate_atom_index_v3000_sdf(tmp_path: Path) -> Path:
    """Write a single molecule to a v3000 sdf with a duplicate atom index.

    Args:
        tmp_path: pytest fixture for writing files to a temp directory

    Returns:
        Path to the sdf
    """
    sdf_text = """


  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 9 9 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 87.71 -95.64 0 0
M  V30 2 C 87.71 -81.29 0 0
M  V30 3 C 100.18 -74.09 0 0
M  V30 4 C 100.18 -59.69 0 0
M  V30 5 C 87.71 -52.49 0 0
M  V30 6 C 75.24 -59.69 0 0
M  V30 7 C 75.24 -74.09 0 0
M  V30 1 C 87.71 -38.09 0 0
M  V30 9 O 100.18 -30.89 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 2 3 4
M  V30 4 1 4 5
M  V30 5 2 5 6
M  V30 6 1 6 7
M  V30 7 2 7 2
M  V30 8 1 5 1
M  V30 9 2 1 9
M  V30 END BOND
M  V30 END CTAB
M  END
$$$$
    """
    outpath = tmp_path / "input.sdf"
    with open(outpath, "w") as outh:
        outh.write(sdf_text)
    return outpath


@pytest.fixture()
def single_duplicate_bond_index_v3000_sdf(tmp_path: Path) -> Path:
    """Write a single molecule to a v3000 sdf with a duplicate bond index.

    Args:
        tmp_path: pytest fixture for writing files to a temp directory

    Returns:
        Path to the sdf
    """
    sdf_text = """


  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 9 9 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 87.71 -95.64 0 0
M  V30 2 C 87.71 -81.29 0 0
M  V30 3 C 100.18 -74.09 0 0
M  V30 4 C 100.18 -59.69 0 0
M  V30 5 C 87.71 -52.49 0 0
M  V30 6 C 75.24 -59.69 0 0
M  V30 7 C 75.24 -74.09 0 0
M  V30 8 C 87.71 -38.09 0 0
M  V30 9 O 100.18 -30.89 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 2 3 4
M  V30 4 1 4 5
M  V30 5 2 5 6
M  V30 6 1 6 7
M  V30 7 2 7 2
M  V30 1 1 5 8
M  V30 9 2 8 9
M  V30 END BOND
M  V30 END CTAB
M  END
$$$$
"""
    outpath = tmp_path / "input.sdf"
    with open(outpath, "w") as outh:
        outh.write(sdf_text)
    return outpath


@pytest.fixture()
def single_shuffled_atom_lines_v3000_sdf(tmp_path: Path) -> Path:
    """Write a single molecule to a v3000 sdf with atom lines not in index order.

    Args:
        tmp_path: pytest fixture for writing files to a temp directory

    Returns:
        Path to the sdf
    """
    sdf_text = """


  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 9 9 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 87.71 -95.64 0 0
M  V30 2 C 87.71 -81.29 0 0
M  V30 3 C 100.18 -74.09 0 0
M  V30 4 C 100.18 -59.69 0 0
M  V30 5 C 87.71 -52.49 0 0
M  V30 6 C 75.24 -59.69 0 0
M  V30 8 C 87.71 -38.09 0 0
M  V30 7 C 75.24 -74.09 0 0
M  V30 9 O 100.18 -30.89 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 2 3 4
M  V30 4 1 4 5
M  V30 5 2 5 6
M  V30 6 1 6 7
M  V30 7 2 7 2
M  V30 8 1 5 8
M  V30 9 2 8 9
M  V30 END BOND
M  V30 END CTAB
M  END
$$$$
"""
    outpath = tmp_path / "input.sdf"
    with open(outpath, "w") as outh:
        outh.write(sdf_text)
    return outpath


@pytest.fixture()
def single_shuffled_bond_lines_v3000_sdf(tmp_path: Path) -> Path:
    """Write a single molecule to a v3000 sdf with bond lines not in index order.

    Args:
        tmp_path: pytest fixture for writing files to a temp directory

    Returns:
        Path to the sdf
    """
    sdf_text = """


  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 9 9 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 87.71 -95.64 0 0
M  V30 2 C 87.71 -81.29 0 0
M  V30 3 C 100.18 -74.09 0 0
M  V30 4 C 100.18 -59.69 0 0
M  V30 5 C 87.71 -52.49 0 0
M  V30 6 C 75.24 -59.69 0 0
M  V30 7 C 75.24 -74.09 0 0
M  V30 8 C 87.71 -38.09 0 0
M  V30 9 O 100.18 -30.89 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 2 3 4
M  V30 4 1 4 5
M  V30 5 2 5 6
M  V30 6 1 6 7
M  V30 8 1 5 8
M  V30 7 2 7 2
M  V30 9 2 8 9
M  V30 END BOND
M  V30 END CTAB
M  END
$$$$
"""
    outpath = tmp_path / "input.sdf"
    with open(outpath, "w") as outh:
        outh.write(sdf_text)
    return outpath


@pytest.fixture()
def single_mol_v2000_sdf(tmp_path: Path) -> Path:
    """Write a single molecule with no metadata to an sdf, forced to v2000.

    Args:
        tmp_path: pytest fixture for writing files to a temp directory

    Returns:
        Path to the sdf
    """
    mol = Chem.MolFromSmiles("C")
    outpath = tmp_path / "input.sdf"
    with Chem.SDWriter(str(outpath)) as sdw:
        sdw.SetForceV3000(False)
        sdw.write(mol)
    return outpath


@pytest.fixture()
def single_mol_v3000_sdf(tmp_path: Path) -> Path:
    """Write a single molecule with no metadata to an sdf, forced to v3000.

    Args:
        tmp_path: pytest fixture for writing files to a temp directory

    Returns:
        Path to the sdf
    """
    mol = Chem.MolFromSmiles("C")
    outpath = tmp_path / "input.sdf"
    with Chem.SDWriter(str(outpath)) as sdw:
        sdw.SetForceV3000(True)
        sdw.write(mol)
    return outpath
