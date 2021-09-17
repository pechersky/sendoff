"""Provide fixtures of SDFs with well- or misbehaved metadata."""
from pathlib import Path

import pytest
from rdkit import Chem


@pytest.fixture
def single_mol_sdf(tmp_path: Path) -> Path:
    """Write a single molecule with no metadata to an sdf.

    Args:
        tmp_path: pytest fixture for writing files to a temp directory

    Returns:
        Path to the sdf
    """
    mol = Chem.MolFromSmiles("C")
    outpath = tmp_path / "input.sdf"
    with Chem.SDWriter(str(outpath)) as sdw:
        sdw.write(mol)
    return outpath


@pytest.fixture
def double_mol_sdf(tmp_path: Path) -> Path:
    """Write two molecules with no metadata to an sdf.

    Args:
        tmp_path: pytest fixture for writing files to a temp directory

    Returns:
        Path to the sdf
    """
    methane = Chem.MolFromSmiles("C")
    ethane = Chem.MolFromSmiles("CC")
    outpath = tmp_path / "input.sdf"
    with Chem.SDWriter(str(outpath)) as sdw:
        sdw.write(methane)
        sdw.write(ethane)
    return outpath


@pytest.fixture
def single_titled_mol_sdf(tmp_path: Path) -> Path:
    """Write a single molecule with only title metadata to an sdf.

    Args:
        tmp_path: pytest fixture for writing files to a temp directory

    Returns:
        Path to the sdf
    """
    mol = Chem.MolFromSmiles("C")
    mol.SetProp("_Name", "Title")
    outpath = tmp_path / "input.sdf"
    with Chem.SDWriter(str(outpath)) as sdw:
        sdw.write(mol)
    return outpath


@pytest.fixture
def single_delimiter_titled_mol_sdf(tmp_path: Path) -> Path:
    """Write a single molecule with delimiter title metadata to an sdf.

    Args:
        tmp_path: pytest fixture for writing files to a temp directory

    Returns:
        Path to the sdf
    """
    mol = Chem.MolFromSmiles("C")
    mol.SetProp("_Name", "$$$$")
    outpath = tmp_path / "input.sdf"
    with Chem.SDWriter(str(outpath)) as sdw:
        sdw.write(mol)
    return outpath


@pytest.fixture
def double_first_delimiter_titled_mol_sdf(tmp_path: Path) -> Path:
    """Write two molecules with the first as delimiter title metadata to an sdf.

    Args:
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
        sdw.write(mol)
        sdw.write(other_mol)
    return outpath


@pytest.fixture
def double_second_delimiter_titled_mol_sdf(tmp_path: Path) -> Path:
    """Write two molecules with the second as delimiter title metadata to an sdf.

    Args:
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
        sdw.write(other_mol)
        sdw.write(mol)
    return outpath


@pytest.fixture
def single_record_mol_sdf(tmp_path: Path) -> Path:
    """Write a single molecule with only record value metadata to an sdf.

    Args:
        tmp_path: pytest fixture for writing files to a temp directory

    Returns:
        Path to the sdf
    """
    mol = Chem.MolFromSmiles("C")
    mol.SetProp("Record", "Value")
    outpath = tmp_path / "input.sdf"
    with Chem.SDWriter(str(outpath)) as sdw:
        sdw.write(mol)
    return outpath


@pytest.fixture
def single_delimiter_record_mol_sdf(tmp_path: Path) -> Path:
    """Write a single molecule with delimiter record value metadata to an sdf.

    Args:
        tmp_path: pytest fixture for writing files to a temp directory

    Returns:
        Path to the sdf
    """
    mol = Chem.MolFromSmiles("C")
    mol.SetProp("Record", "$$$$")
    outpath = tmp_path / "input.sdf"
    with Chem.SDWriter(str(outpath)) as sdw:
        sdw.write(mol)
    return outpath


@pytest.fixture
def single_multiline_record_mol_sdf(tmp_path: Path) -> Path:
    """Write a single molecule with a multiline record value metadata to an sdf.

    Args:
        tmp_path: pytest fixture for writing files to a temp directory

    Returns:
        Path to the sdf
    """
    mol = Chem.MolFromSmiles("C")
    value = str.join("\n", map(str, [0, 1, 2, 3]))
    mol.SetProp("Record", value)
    outpath = tmp_path / "input.sdf"
    with Chem.SDWriter(str(outpath)) as sdw:
        sdw.write(mol)
    return outpath


@pytest.fixture
def single_empty_string_record_mol_sdf(tmp_path: Path) -> Path:
    """Write a single molecule with an empty string record value to an sdf.

    Args:
        tmp_path: pytest fixture for writing files to a temp directory

    Returns:
        Path to the sdf
    """
    mol = Chem.MolFromSmiles("C")
    mol.SetProp("Record", "")
    outpath = tmp_path / "input.sdf"
    with Chem.SDWriter(str(outpath)) as sdw:
        sdw.write(mol)
    return outpath
