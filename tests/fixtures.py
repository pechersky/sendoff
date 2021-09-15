"""Provide fixtures of SDFs with well- or misbehaved metadata."""
import pytest
from rdkit import Chem


@pytest.fixture
def methane():
    """Provide a molecule with no metadata.

    Returns:
        Methane, as an rdkit.Chem.rdchem.Mol
    """
    return Chem.MolFromSmiles("C")


@pytest.fixture
def ethane():
    """Provide a molecule with no metadata.

    Returns:
        Ethane, as an rdkit.Chem.rdchem.Mol
    """
    return Chem.MolFromSmiles("CC")


@pytest.fixture
def single_mol_sdf(methane, tmp_path):
    """Write a single molecule with no metadata to an sdf.

    Args:
        methane: pytest fixture providing a single molecule
        tmp_path: pytest fixture for writing files to a temp directory

    Returns:
        Path to the sdf
    """
    outpath = tmp_path / "input.sdf"
    with Chem.SDWriter(str(outpath)) as sdw:
        sdw.write(methane)
    return outpath


@pytest.fixture
def double_mol_sdf(methane, ethane, tmp_path):
    """Write two molecules with no metadata to an sdf.

    Args:
        methane: pytest fixture providing a single molecule
        ethane: pytest fixture providing a single molecule
        tmp_path: pytest fixture for writing files to a temp directory

    Returns:
        Path to the sdf
    """
    outpath = tmp_path / "input.sdf"
    with Chem.SDWriter(str(outpath)) as sdw:
        sdw.write(methane)
        sdw.write(ethane)
    return outpath


@pytest.fixture
def single_titled_mol_sdf(methane, tmp_path):
    """Write a single molecule with only title metadata to an sdf.

    Args:
        methane: pytest fixture providing a single molecule
        tmp_path: pytest fixture for writing files to a temp directory

    Returns:
        Path to the sdf
    """
    mol = methane
    mol.SetProp("_Name", "Title")
    outpath = tmp_path / "input.sdf"
    with Chem.SDWriter(str(outpath)) as sdw:
        sdw.write(mol)
    return outpath


@pytest.fixture
def single_delimiter_titled_mol_sdf(methane, tmp_path):
    """Write a single molecule with delimiter title metadata to an sdf.

    Args:
        methane: pytest fixture providing a single molecule
        tmp_path: pytest fixture for writing files to a temp directory

    Returns:
        Path to the sdf
    """
    mol = methane
    mol.SetProp("_Name", "$$$$")
    outpath = tmp_path / "input.sdf"
    with Chem.SDWriter(str(outpath)) as sdw:
        sdw.write(mol)
    return outpath


@pytest.fixture
def single_record_mol_sdf(methane, tmp_path):
    """Write a single molecule with only record value metadata to an sdf.

    Args:
        methane: pytest fixture providing a single molecule
        tmp_path: pytest fixture for writing files to a temp directory

    Returns:
        Path to the sdf
    """
    mol = methane
    mol.SetProp("Record", "Value")
    outpath = tmp_path / "input.sdf"
    with Chem.SDWriter(str(outpath)) as sdw:
        sdw.write(mol)
    return outpath


@pytest.fixture
def single_delimiter_record_mol_sdf(methane, tmp_path):
    """Write a single molecule with delimiter record value metadata to an sdf.

    Args:
        methane: pytest fixture providing a single molecule
        tmp_path: pytest fixture for writing files to a temp directory

    Returns:
        Path to the sdf
    """
    mol = methane
    mol.SetProp("Record", "$$$$")
    outpath = tmp_path / "input.sdf"
    with Chem.SDWriter(str(outpath)) as sdw:
        sdw.write(mol)
    return outpath


@pytest.fixture
def single_multiline_record_mol_sdf(methane, tmp_path):
    """Write a single molecule with a multiline record value metadata to an sdf.

    Args:
        methane: pytest fixture providing a single molecule
        tmp_path: pytest fixture for writing files to a temp directory

    Returns:
        Path to the sdf
    """
    mol = methane
    value = str.join("\n", map(str, [0, 1, 2, 3]))
    mol.SetProp("Record", value)
    outpath = tmp_path / "input.sdf"
    with Chem.SDWriter(str(outpath)) as sdw:
        sdw.write(mol)
    return outpath


@pytest.fixture
def single_empty_string_record_mol_sdf(methane, tmp_path):
    """Write a single molecule with an empty string record value to an sdf.

    Args:
        methane: pytest fixture providing a single molecule
        tmp_path: pytest fixture for writing files to a temp directory

    Returns:
        Path to the sdf
    """
    mol = methane
    mol.SetProp("Record", "")
    outpath = tmp_path / "input.sdf"
    with Chem.SDWriter(str(outpath)) as sdw:
        sdw.write(mol)
    return outpath
