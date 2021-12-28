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
