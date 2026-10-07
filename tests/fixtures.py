# mypy: allow-untyped-decorators
"""Provide fixtures of SDFs with well- or misbehaved metadata."""
from collections.abc import Sequence
from pathlib import Path
from typing import TYPE_CHECKING, NamedTuple

import pytest
from rdkit import Chem

from sendoff.sdblock import parse_sdf

if TYPE_CHECKING:
    from _pytest.fixtures import FixtureRequest as __FixtureRequest

    # This is to provide a proper type; we are subclassing Any, which causes mypy warn
    class FixtureRequest(__FixtureRequest):
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


V3000_ATOMS = [
    "C 87.71 -95.64 0 0",
    "C 87.71 -81.29 0 0",
    "C 100.18 -74.09 0 0",
    "C 100.18 -59.69 0 0",
    "C 87.71 -52.49 0 0",
    "C 75.24 -59.69 0 0",
    "C 75.24 -74.09 0 0",
    "C 87.71 -38.09 0 0",
    "O 100.18 -30.89 0 0",
]


# Each bond's type and the positions, from 1, of its two atoms in V3000_ATOMS.
V3000_BONDS = [
    (1, 1, 2),
    (1, 2, 3),
    (2, 3, 4),
    (1, 4, 5),
    (2, 5, 6),
    (1, 6, 7),
    (2, 7, 2),
    (1, 5, 8),
    (2, 8, 9),
]


class V3000Mol(NamedTuple):
    """A 9-atom V3000 molecule, whose counts and indices tests can break.

    Attributes:
        counts: the atom and bond counts on the COUNTS line
        atom_ids: the index of each atom in V3000_ATOMS
        bond_ids: the index of each bond in V3000_BONDS
        atom_order: the positions in V3000_ATOMS, from 0, in line order
        bond_order: the positions in V3000_BONDS, from 0, in line order
    """

    counts: tuple[int, int] = (9, 9)
    atom_ids: Sequence[int] = range(1, 10)
    bond_ids: Sequence[int] = range(1, 10)
    atom_order: Sequence[int] = range(9)
    bond_order: Sequence[int] = range(9)

    def sdf(self) -> str:
        """Write the molecule as an untitled sdf block without metadata.

        Returns:
            The sdf text.
        """
        atoms = [f"M  V30 {self.atom_ids[i]} {V3000_ATOMS[i]}" for i in self.atom_order]
        bonds = []
        for i in self.bond_order:
            kind, first, second = V3000_BONDS[i]
            first_id, second_id = self.atom_ids[first - 1], self.atom_ids[second - 1]
            bonds.append(f"M  V30 {self.bond_ids[i]} {kind} {first_id} {second_id}")
        lines = [
            "",
            "",
            "",
            "  0  0  0     0  0            999 V3000",
            "M  V30 BEGIN CTAB",
            "M  V30 COUNTS {} {} 0 0 0".format(*self.counts),
            "M  V30 BEGIN ATOM",
            *atoms,
            "M  V30 END ATOM",
            "M  V30 BEGIN BOND",
            *bonds,
            "M  V30 END BOND",
            "M  V30 END CTAB",
            "M  END",
            "$$$$",
        ]
        return str.join("", (line + "\n" for line in lines))


SEVENTH_AND_EIGHTH_SWAPPED = (0, 1, 2, 3, 4, 5, 7, 6, 8)


INDEX_DEFECTS = {
    "large_atom_index": V3000Mol(
        atom_ids=(2101315825, 2101315840, *range(2101315860, 2101315867))
    ),
    "large_bond_index": V3000Mol(bond_ids=(2101315886, *range(2101315911, 2101315919))),
    "missing_atom_line": V3000Mol(counts=(10, 9)),
    "extra_atom_line": V3000Mol(counts=(8, 9)),
    "missing_bond_line": V3000Mol(counts=(9, 10)),
    "extra_bond_line": V3000Mol(counts=(9, 8)),
    "duplicate_atom_index": V3000Mol(atom_ids=(1, 2, 3, 4, 5, 6, 7, 1, 9)),
    "duplicate_bond_index": V3000Mol(bond_ids=(1, 2, 3, 4, 5, 6, 7, 1, 9)),
    "shuffled_atom_lines": V3000Mol(atom_order=SEVENTH_AND_EIGHTH_SWAPPED),
    "shuffled_bond_lines": V3000Mol(bond_order=SEVENTH_AND_EIGHTH_SWAPPED),
}


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
