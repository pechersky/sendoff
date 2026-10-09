# mypy: allow-untyped-decorators
"""Test CTAB operations: counts, index checks, renumbering, and CTable."""
import io
import typing
from pathlib import Path

import pytest

from sendoff.ctable import (
    CTable,
    CTableError,
    IndicesDuplicateError,
    IndicesMismatchError,
    IndicesOutOfOrderError,
    MalformedCTableError,
)
from sendoff.sdblock import SDBlock, parse_sdf, write_sdf
from tests.fixtures import INDEX_DEFECTS, CountedMol


def test_counts_mol_read(counted_mol_sdf: CountedMol) -> None:
    """An sdf block parses with its molecule's atom and bond counts.

    Args:
        counted_mol_sdf: pytest fixture of a Path to the sdf, and its counts
    """
    mol: SDBlock[str] = next(parse_sdf(counted_mol_sdf.path))
    assert (mol.num_atoms, mol.num_bonds) == (
        counted_mol_sdf.atoms,
        counted_mol_sdf.bonds,
    )


def defect_block(defect: str) -> SDBlock[str]:
    """Parse the 9-atom V3000 molecule with one of INDEX_DEFECTS.

    Args:
        defect: the name of the defect

    Returns:
        The parsed block.
    """
    return next(SDBlock.from_sdf(INDEX_DEFECTS[defect].sdf()))


IndexCheck = tuple[type[CTableError], str] | None


@pytest.mark.parametrize(
    "defect, kind, loose, strict",
    [
        (
            "large_atom_index",
            "atom",
            None,
            (
                IndicesOutOfOrderError,
                "CTable line 8: atom index 2101315825 is out of order",
            ),
        ),
        (
            "large_bond_index",
            "bond",
            None,
            (
                IndicesOutOfOrderError,
                "CTable line 19: bond index 2101315886 is out of order",
            ),
        ),
        (
            "missing_atom_line",
            "atom",
            (
                IndicesMismatchError,
                "fewer atom lines (9) than the counts line declares (10)",
            ),
            (
                IndicesMismatchError,
                "fewer atom lines (9) than the counts line declares (10)",
            ),
        ),
        (
            "extra_atom_line",
            "atom",
            (
                IndicesMismatchError,
                "more atom lines (9) than the counts line declares (8)",
            ),
            (
                IndicesMismatchError,
                "more atom lines (9) than the counts line declares (8)",
            ),
        ),
        (
            "missing_bond_line",
            "bond",
            (
                IndicesMismatchError,
                "fewer bond lines (9) than the counts line declares (10)",
            ),
            (
                IndicesMismatchError,
                "fewer bond lines (9) than the counts line declares (10)",
            ),
        ),
        (
            "extra_bond_line",
            "bond",
            (
                IndicesMismatchError,
                "more bond lines (9) than the counts line declares (8)",
            ),
            (
                IndicesMismatchError,
                "more bond lines (9) than the counts line declares (8)",
            ),
        ),
        # Strict checking finds the repeated index out of order first.
        (
            "duplicate_atom_index",
            "atom",
            (IndicesDuplicateError, "CTable line 15: atom index 1 is a duplicate"),
            (IndicesOutOfOrderError, "CTable line 15: atom index 1 is out of order"),
        ),
        (
            "duplicate_bond_index",
            "bond",
            (IndicesDuplicateError, "CTable line 26: bond index 1 is a duplicate"),
            (IndicesOutOfOrderError, "CTable line 26: bond index 1 is out of order"),
        ),
        (
            "shuffled_atom_lines",
            "atom",
            None,
            (IndicesOutOfOrderError, "CTable line 14: atom index 8 is out of order"),
        ),
        (
            "shuffled_bond_lines",
            "bond",
            None,
            (IndicesOutOfOrderError, "CTable line 25: bond index 8 is out of order"),
        ),
    ],
)
def test_index_checks_find_each_defect(
    defect: str, kind: str, loose: IndexCheck, strict: IndexCheck
) -> None:
    """Checking the defective section raises, and the other section is valid.

    Without strict, indices need only match the counts line and be unique;
    strict also needs them to count up from 1 in line order.

    Args:
        defect: the name of the defect in INDEX_DEFECTS
        kind: the section with the defect, "atom" or "bond"
        loose: the expected error type and message without strict, if any
        strict: the expected error type and message with strict, if any
    """
    ctable = defect_block(defect).ctable()
    for is_strict, expected in [(False, loose), (True, strict)]:
        check = getattr(ctable, f"valid_{kind}_indices")
        if expected is None:
            assert check(strict=is_strict)
            continue
        error, message = expected
        with pytest.raises(error) as raised:
            check(strict=is_strict)
        assert type(raised.value) is error
        assert str(raised.value) == message
    other = "bond" if kind == "atom" else "atom"
    assert getattr(ctable, f"valid_{other}_indices")(strict=True)


@pytest.mark.parametrize(
    "defect, kind, mismatch",
    [
        ("large_atom_index", "atom", None),
        ("large_bond_index", "bond", None),
        (
            "missing_atom_line",
            "atom",
            "fewer atom lines (9) than the counts line declares (10)",
        ),
        (
            "extra_atom_line",
            "atom",
            "more atom lines (9) than the counts line declares (8)",
        ),
        (
            "missing_bond_line",
            "bond",
            "fewer bond lines (9) than the counts line declares (10)",
        ),
        (
            "extra_bond_line",
            "bond",
            "more bond lines (9) than the counts line declares (8)",
        ),
        ("duplicate_bond_index", "bond", None),
        ("shuffled_atom_lines", "atom", None),
        ("shuffled_bond_lines", "bond", None),
    ],
)
def test_renumbering_mends_index_defects(
    defect: str, kind: str, mismatch: str | None
) -> None:
    """Renumbered indices are strictly valid, once the counts line matches.

    Renumbering leaves the counts line; fix_counts makes it match.

    Args:
        defect: the name of the defect in INDEX_DEFECTS
        kind: the section with the defect, "atom" or "bond"
        mismatch: the counts mismatch left after renumbering, if any
    """
    block = defect_block(defect)
    block.renumber_indices()
    if mismatch is not None:
        with pytest.raises(IndicesMismatchError) as raised:
            getattr(block.ctable(), f"valid_{kind}_indices")(strict=True)
        assert str(raised.value) == mismatch
        block.fix_counts()
    ctable = block.ctable()
    assert ctable.valid_atom_indices(strict=True)
    assert ctable.valid_bond_indices(strict=True)


def test_renumbering_refuses_bonds_to_a_duplicate_atom_index() -> None:
    """Renumbering cannot tell which atom with a duplicate index a bond names."""
    block = defect_block("duplicate_atom_index")
    with pytest.raises(IndicesDuplicateError) as raised:
        block.renumber_indices()
    assert str(raised.value) == (
        "CTable line 19: bond atom 1 is the index of several atoms"
    )


def test_single_mol_v3000_sdf_idempotent_renumber(single_mol_v3000_sdf: Path) -> None:
    """A V3000 sdf block generated by RDKit, renumbered, is the same as the input.

    Whether parsed from the file or from its lines without endings.

    Args:
        single_mol_v3000_sdf: pytest fixture of a Path to the sdf
    """
    text = single_mol_v3000_sdf.read_text()
    for mol in (
        next(parse_sdf(single_mol_v3000_sdf)),
        next(SDBlock.from_lines(text.splitlines())),
    ):
        mol.renumber_indices()
        buffer = io.StringIO()
        mol.write(buffer)
        assert buffer.getvalue() == text


V3000_HEADER = "t\n  sendoff\n\n  0  0  0     0  0            999 V3000\n"


def v3000_block(
    atom: str = "M  V30 1 C 0 0 0 0",
    bond: str = "M  V30 1 1 1 2",
    counts: str = "M  V30 COUNTS 2 1 0 0 0",
    extra: str = "",
    after: str = "",
) -> SDBlock[str]:
    """Build a two-atom V3000 block with one replaceable line of each kind.

    Args:
        atom: the first atom line
        bond: the only bond line
        counts: the V3000 counts line
        extra: lines to add after the bond section, inside the CTAB
        after: lines to add after the CTAB

    Returns:
        The parsed block.
    """
    text = (
        V3000_HEADER
        + f"M  V30 BEGIN CTAB\n{counts}\nM  V30 BEGIN ATOM\n{atom}\n"
        + "M  V30 2 O 0 0 0 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\n"
        + f"{bond}\nM  V30 END BOND\n{extra}M  V30 END CTAB\n{after}M  END\n$$$$\n"
    )
    return next(SDBlock.from_lines(text.splitlines()))


CTABLE_OPERATIONS = [
    "ctable",
    "num_atoms",
    "num_bonds",
    "valid_atom_indices",
    "valid_bond_indices",
    "renumber_indices",
    "fix_counts",
]


def run_ctable_operation(block: SDBlock[typing.Any], operation: str) -> object:
    """Run an SDBlock CTable operation, or a validation on the block's CTable.

    Args:
        block: the block to read or change
        operation: one of CTABLE_OPERATIONS

    Returns:
        What the operation returns.
    """
    if operation.startswith("valid_"):
        return getattr(block.ctable(), operation)()
    if operation.startswith("num_"):
        return getattr(block, operation)
    return getattr(block, operation)()


V3000_CTAB = "M  V30 BEGIN CTAB\n"


@pytest.mark.parametrize("binary", [False, True], ids=["str", "bytes"])
@pytest.mark.parametrize("operation", CTABLE_OPERATIONS)
@pytest.mark.parametrize(
    "header, message",
    [
        ("t\n", "CTable line 3: missing; the table has only 2 lines"),
        (V3000_HEADER, "CTable line 6: missing; the table has only 5 lines"),
        ("t\n\n\n   \n", "CTable line 4: no format field"),
        (
            "t\n\n\n  0  0  0     0  0            999 V9\n",
            'CTable line 4: format "V9" is not V2000 or V3000',
        ),
        (
            "t\n\n\n  x  0  0     0  0            999 V2000\n",
            'CTable line 4: atom count "  x" is not an integer',
        ),
        (
            "t\n\n\n  0  y  0     0  0            999 V2000\n",
            'CTable line 4: bond count "  y" is not an integer',
        ),
        (
            V3000_HEADER + V3000_CTAB + "M  V30 COUNTS\n",
            "CTable line 6: no atom count field",
        ),
        (
            V3000_HEADER + V3000_CTAB + "M  V30 COUNTS 1\n",
            "CTable line 6: no bond count field",
        ),
        (
            V3000_HEADER + V3000_CTAB + "M  V30 COUNTS x 0 0 0 0\n",
            'CTable line 6: atom count "x" is not an integer',
        ),
        (
            V3000_HEADER + V3000_CTAB + "M  V30 COUNTS 0 y 0 0 0\n",
            'CTable line 6: bond count "y" is not an integer',
        ),
    ],
)
def test_malformed_headers_fail_every_operation(
    header: str, message: str, operation: str, binary: bool
) -> None:
    """Every CTable operation reports a malformed header line by number.

    Args:
        header: the CTable lines before M  END
        message: the expected error message
        operation: one of CTABLE_OPERATIONS
        binary: whether the block is a bytes block
    """
    text = header + "M  END\n$$$$\n"
    block: SDBlock[typing.Any] = next(
        SDBlock.from_sdf(text.encode() if binary else text)
    )
    with pytest.raises(MalformedCTableError) as raised:
        run_ctable_operation(block, operation)
    assert str(raised.value) == message
    assert block.to_sdf() == (text.encode() if binary else text)


RENUMBERING = {"renumber_indices"}


ATOM_INDICES = {"valid_atom_indices", "renumber_indices"}


BOND_INDICES = {"valid_bond_indices", "renumber_indices"}


CTableErrors = dict[str, tuple[type[Exception], str]]


def malformed(operations: set[str], message: str) -> CTableErrors:
    """Expect each of `operations` to raise MalformedCTableError(message).

    Args:
        operations: the operations that read the malformed line
        message: the expected error message

    Returns:
        The expected error of each failing operation.
    """
    return {operation: (MalformedCTableError, message) for operation in operations}


def ambiguous_bond(bond: str) -> tuple[dict[str, str], CTableErrors]:
    """Expect errors for a bond to atom 2 when two atoms have index 2.

    Args:
        bond: the bond line, which names atom 2

    Returns:
        v3000_block arguments and the expected error of each failing operation.
    """
    return {"atom": "M  V30 2 C 0 0 0 0", "bond": bond}, {
        "valid_atom_indices": (
            IndicesDuplicateError,
            "CTable line 9: atom index 2 is a duplicate",
        ),
        "renumber_indices": (
            IndicesDuplicateError,
            "CTable line 12: bond atom 2 is the index of several atoms",
        ),
    }


@pytest.mark.parametrize("operation", CTABLE_OPERATIONS)
@pytest.mark.parametrize(
    "lines, errors",
    [
        (
            {"atom": "M  V30"},
            malformed(ATOM_INDICES, "CTable line 8: no atom index field"),
        ),
        (
            {"atom": "M  V30 x C 0 0 0 0"},
            malformed(ATOM_INDICES, 'CTable line 8: atom index "x" is not an integer'),
        ),
        (
            {"bond": "M  V30"},
            malformed(BOND_INDICES, "CTable line 12: no bond index field"),
        ),
        (
            {"bond": "M  V30 x 1 1 2"},
            malformed(BOND_INDICES, 'CTable line 12: bond index "x" is not an integer'),
        ),
        (
            {"bond": "M  V30 1"},
            malformed(RENUMBERING, "CTable line 12: no bond type field"),
        ),
        (
            {"bond": "M  V30 1 1"},
            malformed(RENUMBERING, "CTable line 12: no first bond atom field"),
        ),
        (
            {"bond": "M  V30 1 1 1"},
            malformed(RENUMBERING, "CTable line 12: no second bond atom field"),
        ),
        (
            {"bond": "M  V30 1 1 x 2"},
            malformed(
                RENUMBERING, 'CTable line 12: first bond atom "x" is not an integer'
            ),
        ),
        (
            {"bond": "M  V30 1 1 1 x"},
            malformed(
                RENUMBERING, 'CTable line 12: second bond atom "x" is not an integer'
            ),
        ),
        (
            {"bond": "M  V30 1 1 9 2"},
            malformed(RENUMBERING, "CTable line 12: bond atom 9 is not an atom index"),
        ),
        (
            {"bond": "M  V30 1 1 1 9"},
            malformed(RENUMBERING, "CTable line 12: bond atom 9 is not an atom index"),
        ),
        # A continued line is one logical line, which only renumbering refuses.
        (
            {"atom": "M  V30 1 C 0 0 0 0 -\nM  V30 CHG=1"},
            {
                "renumber_indices": (
                    NotImplementedError,
                    "renumbering does not support continued (-) atom or bond "
                    "lines yet",
                )
            },
        ),
        # Continuing a section's last line is malformed, even for renumbering.
        (
            {"bond": "M  V30 1 1 1 2 -"},
            malformed(
                BOND_INDICES | {"fix_counts"},
                'CTable line 12: continues ("-") onto the "M  V30 END BOND" line',
            ),
        ),
        # A duplicated atom index outranks a missing one in the same bond.
        ambiguous_bond("M  V30 1 1 9 2"),
        ambiguous_bond("M  V30 1 1 2 9"),
    ],
)
def test_malformed_lines_fail_the_operations_that_read_them(
    lines: dict[str, str], errors: CTableErrors, operation: str
) -> None:
    """Operations that read a bad atom or bond line report it; the rest work.

    Failing operations leave the block as it was.

    Args:
        lines: v3000_block arguments that replace an atom or bond line
        errors: the expected error type and message of each failing operation
        operation: one of CTABLE_OPERATIONS
    """
    block = v3000_block(**lines)
    text = block.to_sdf()
    if operation not in errors:
        run_ctable_operation(block, operation)
        return
    error, message = errors[operation]
    with pytest.raises(error) as raised:
        run_ctable_operation(block, operation)
    assert type(raised.value) is error
    assert str(raised.value) == message
    assert block.to_sdf() == text


@pytest.mark.parametrize(
    "operation, message",
    [
        ("valid_atom_indices", "validating atom indices"),
        ("valid_bond_indices", "validating bond indices"),
        ("renumber_indices", "renumbering"),
        ("fix_counts", "fixing counts"),
    ],
)
def test_v3000_operations_name_themselves_on_v2000_tables(
    operation: str, message: str
) -> None:
    """V3000-only operations raise NotImplementedError, not a CTableError.

    Args:
        operation: one of CTABLE_OPERATIONS
        message: the operation's name in the error message
    """
    text = "t\n\n\n  0  0  0     0  0            999 V2000\nM  END\n$$$$\n"
    with pytest.raises(NotImplementedError) as raised:
        run_ctable_operation(next(SDBlock.from_sdf(text)), operation)
    assert str(raised.value) == f"{message} is only implemented for V3000 tables"


def test_ctable_errors_are_value_errors() -> None:
    """Every error about a table's contents is a CTableError and a ValueError."""
    for error in (
        MalformedCTableError,
        IndicesDuplicateError,
        IndicesMismatchError,
        IndicesOutOfOrderError,
    ):
        assert issubclass(error, CTableError)
    assert issubclass(CTableError, ValueError)


def test_renumber_keeps_atom_lines_after_the_index_field() -> None:
    """Renumbering keeps whatever follows an atom's index field, even non-ASCII."""
    block = v3000_block(atom="a b 1 \u4e2d")
    block.renumber_indices()
    assert "M  V30 1 \u4e2d" in block.ctable().lines


def test_renumber_reads_index_fields_after_extra_spaces() -> None:
    """Atom and bond lines with extra spaces renumber consistently."""
    block = v3000_block(atom="M  V30  10 C 0 0 0 0", bond="M  V30   7 1  10  2")
    block.renumber_indices()
    lines = block.ctable().lines
    assert "M  V30 1 C 0 0 0 0" in lines
    assert "M  V30 1 1 1 2" in lines
    assert block.ctable().valid_atom_indices(strict=True)


def test_renumber_leaves_counts() -> None:
    """Renumbering leaves a counts line that disagrees with the table."""
    block = v3000_block(
        "M  V30 3 C 0 0 0 0", "M  V30 1 1 3 2", counts="M  V30 COUNTS 5 3 0 0 0"
    )
    ctable = block.ctable()
    ctable.renumber_indices()
    assert "M  V30 1 C 0 0 0 0" in ctable.lines
    assert "M  V30 COUNTS 5 3 0 0 0" in ctable.lines
    assert (ctable.num_atoms, ctable.num_bonds) == (5, 3)
    block.renumber_indices()
    assert (block.num_atoms, block.num_bonds) == (5, 3)


@pytest.mark.parametrize("ending", ["\n", "\r\n", ""])
def test_fix_counts_sets_counts_to_table_lines(ending: str) -> None:
    """Fixing counts rewrites only the two counts, keeping the rest of the line.

    Args:
        ending: what ends each line given to the CTable
    """
    counts = "M  V30 COUNTS  5 3  1 0 1 REGNO=7"
    block = v3000_block(counts=counts)
    lines = [line + ending for line in block.to_sdf().splitlines()[:-1]]
    ctable = CTable(lines)
    assert (ctable.num_atoms, ctable.num_bonds) == (5, 3)
    ctable.fix_counts()
    fixed = "M  V30 COUNTS  2 1  1 0 1 REGNO=7" + ending
    assert ctable.lines == [
        fixed if line == counts + ending else line for line in lines
    ]
    assert (ctable.num_atoms, ctable.num_bonds) == (2, 1)
    reparsed = CTable(ctable.lines)
    assert (reparsed.num_atoms, reparsed.num_bonds) == (2, 1)
    block.fix_counts()
    assert (block.num_atoms, block.num_bonds) == (2, 1)
    assert block.to_sdf() == v3000_block(counts=fixed.rstrip()).to_sdf()


def test_fix_counts_reads_section_markers() -> None:
    """Fixing counts counts lines between markers, and absent sections as empty."""
    ctab = "M  V30 BEGIN CTAB\nM  V30 COUNTS 3 1 0 0 0\nM  V30 END CTAB\nM  END\n"
    block = next(SDBlock.from_sdf(V3000_HEADER + ctab + "$$$$\n"))
    block.fix_counts()
    assert (block.num_atoms, block.num_bonds) == (0, 0)


@pytest.mark.parametrize("section", ["ATOM", "BOND"])
@pytest.mark.parametrize(
    "operation", ["valid_indices", "renumber_indices", "fix_counts"]
)
def test_sections_need_one_begin_and_one_end_line(operation: str, section: str) -> None:
    """An unterminated or repeated section is malformed and changes nothing.

    Args:
        operation: the method to call; valid_indices validates the section
        section: the section to break
    """
    text = v3000_block().to_sdf()
    begin, end = f"M  V30 BEGIN {section}", f"M  V30 END {section}"
    # The repeated begin line follows the section's end line.
    begin_line, repeat_line = (7, 11) if section == "ATOM" else (11, 14)
    for broken, message in [
        (
            text.replace(end + "\n", ""),
            f'CTable line {begin_line}: "{begin}" has no "{end}" line',
        ),
        (
            text.replace(end + "\n", f"{end}\n{begin}\n{end}\n"),
            f'CTable line {repeat_line}: more than one "{begin}" line in the root CTAB',
        ),
    ]:
        block = next(SDBlock.from_sdf(broken))
        if operation == "valid_indices":
            call = getattr(block.ctable(), f"valid_{section.lower()}_indices")
        else:
            call = getattr(block, operation)
        with pytest.raises(MalformedCTableError) as raised:
            call()
        assert str(raised.value) == message
        assert block.to_sdf() == broken


def test_sections_are_the_root_ctabs() -> None:
    """Validation reads the root CTAB's sections, not those of nested CTABs."""
    root = (
        "M  V30 BEGIN CTAB\nM  V30 COUNTS 1 0 0 0 0\nM  V30 BEGIN ATOM\n"
        + "M  V30 1 R# 0 0 0 0 RGROUPS=(1 1)\nM  V30 END ATOM\nM  V30 END CTAB\n"
    )
    rgroup = (
        "M  V30 BEGIN RGROUP 1\nM  V30 BEGIN CTAB\nM  V30 COUNTS 2 1 0 0 0\n"
        + "M  V30 BEGIN ATOM\nM  V30 1 C 0 0 0 0\nM  V30 2 C 0 0 0 0\n"
        + "M  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 1 1 1 2\nM  V30 END BOND\n"
        + "M  V30 END CTAB\nM  V30 END RGROUP\n"
    )
    text = V3000_HEADER + root + rgroup + "M  END\n$$$$\n"
    ctable = next(SDBlock.from_sdf(text)).ctable()
    assert ctable.valid_atom_indices(strict=True)
    assert ctable.valid_bond_indices(strict=True)


def test_bytes_blocks_edit_their_ctable_as_latin_1() -> None:
    """Bytes blocks read and edit their CTable as Latin-1 text, byte for byte."""
    sdf = b"\xff" + v3000_block(
        "M  V30 3 C 0 0 0 0", "M  V30 1 1 3 2", counts="M  V30 COUNTS 5 3 0 0 0"
    ).to_sdf().encode().replace(b" 2 O 0 0 0 0", b" 2 O 0 0 0 0 X=\xe9")
    raw = next(SDBlock.from_sdf(sdf))
    text = next(SDBlock.from_sdf(sdf.decode("latin-1")))
    lines = raw.ctable().lines
    assert lines[0] == "\xfft"
    assert lines == text.ctable().lines
    for block in (raw, text):
        block.renumber_indices()
        block.fix_counts()
    assert raw.to_sdf() == text.to_sdf().encode("latin-1")
    assert b"M  V30 1 C 0 0 0 0\n" in raw.to_sdf()
    assert b"M  V30 2 O 0 0 0 0 X=\xe9\n" in raw.to_sdf()
    assert (raw.num_atoms, raw.num_bonds) == (2, 1)
    assert raw.title == b"\xfft"


@pytest.mark.parametrize("errors", ["replace", "backslashreplace", "surrogateescape"])
def test_lenient_blocks_edit_their_ctable_byte_for_byte(
    tmp_path: Path,
    errors: typing.Literal["replace", "backslashreplace", "surrogateescape"],
) -> None:
    """Lenient str blocks show bad CTAB bytes decoded, and edit them as read.

    Args:
        tmp_path: Temporary directory for the SDF fixture.
        errors: the handler the block decodes with
    """
    sdf = b"\xff" + v3000_block(
        "M  V30 3 C 0 0 0 0 X=?", "M  V30 1 1 3 2", counts="M  V30 COUNTS 5 3 0 0 0"
    ).to_sdf().encode().replace(b"X=?", b"X=\xe9").replace(
        b" 2 O 0 0 0 0", " 2 O 0 0 0 0 X=\xe9".encode()
    )
    path = tmp_path / "bad.sdf"
    path.write_bytes(sdf)
    (block,) = parse_sdf(path, errors=errors)
    lines = block.ctable().lines
    shown = "\ufffd" if errors == "surrogateescape" else b"\xff".decode(errors=errors)
    assert lines[0] == shown + "t"
    assert "M  V30 2 O 0 0 0 0 X=\xe9" in lines
    raw = next(SDBlock.from_sdf(sdf))
    for edited in (block, raw):
        edited.renumber_indices()
        edited.fix_counts()
    assert b"M  V30 1 C 0 0 0 0 X=\xe9\n" in raw.to_sdf()
    write_sdf(path, [block])
    assert path.read_bytes() == raw.to_sdf()


@pytest.mark.parametrize("ending", ["\n", "\r\n", ""])
def test_renumber_keeps_line_endings(ending: str) -> None:
    """Renumbering an already numbered table keeps its lines, endings included.

    Args:
        ending: what ends each line given to the CTable
    """
    lines = [line + ending for line in v3000_block().to_sdf().splitlines()[:-1]]
    ctable = CTable(lines)
    ctable.renumber_indices()
    assert ctable.lines == lines


def test_renumber_keeps_bond_properties() -> None:
    """Renumbering remaps bond atoms and keeps the properties after them."""
    block = v3000_block(
        atom="M  V30 3 C 0 0 0 0 CHG=1", bond="M  V30 4 1 3 2 CFG=1  STBOX=1"
    )
    block.renumber_indices()
    lines = block.ctable().lines
    assert "M  V30 1 C 0 0 0 0 CHG=1" in lines
    assert "M  V30 1 1 1 2 CFG=1  STBOX=1" in lines


@pytest.mark.parametrize(
    "atom, bond, extra, construct",
    [
        ("M  V30 {} C 0 0 0 0 ATTCHORD=(2 2 1)", "M  V30 1 1 {} 2", "", "ATTCHORD"),
        (
            "M  V30 {} C 0 0 0 0",
            "M  V30 1 9 {} 2 ENDPTS=(1 2) ATTACH=ALL",
            "",
            "ENDPTS",
        ),
        (
            "M  V30 {} C 0 0 0 0",
            "M  V30 1 1 {} 2",
            "M  V30 BEGIN SGROUP\nM  V30 1 SUP 0 ATOMS=(1 {})\nM  V30 END SGROUP\n",
            "SGROUP",
        ),
        (
            "M  V30 1 C 0 0 0 0",
            "M  V30 {} 1 1 2",
            "M  V30 begin sgroup\nM  V30 1 SUP 0 XBONDS=(1 {})\nM  V30 end sgroup\n",
            "SGROUP",
        ),
        (
            "M  V30 {} C 0 0 0 0",
            "M  V30 1 1 {} 2",
            "M  V30 BEGIN COLLECTION\nM  V30 MDLV30/STEABS ATOMS=(1 {})\n"
            + "M  V30 END COLLECTION\n",
            "COLLECTION",
        ),
        (
            "M  V30 {} C 0 0 0 0",
            "M  V30 1 1 {} 2",
            "M  V30 LINKNODE 1 3 1 {} 2\n",
            "LINKNODE",
        ),
    ],
)
def test_renumber_refuses_to_change_referenced_indices(
    atom: str, bond: str, extra: str, construct: str
) -> None:
    """Renumbering refuses to change indices that it would not remap.

    Args:
        atom: the first atom line, formatted with its index
        bond: the only bond line, formatted with an index to change
        extra: lines referring to that index, formatted with it
        construct: what the error names
    """
    block = v3000_block(atom.format(3), bond.format(3), extra=extra.format(3))
    before = block.to_sdf()
    with pytest.raises(NotImplementedError, match=construct):
        block.renumber_indices()
    assert block.to_sdf() == before
    numbered = v3000_block(atom.format(1), bond.format(1), extra=extra.format(1))
    before = numbered.to_sdf()
    numbered.renumber_indices()
    assert numbered.to_sdf() == before


@pytest.mark.parametrize(
    "operation, name",
    [("renumber_indices", "renumbering"), ("fix_counts", "fixing counts")],
)
def test_refuses_nested_ctabs(operation: str, name: str) -> None:
    """Renumbering and fixing counts refuse CTABs nested in RGROUP blocks.

    Args:
        operation: the SDBlock method to call
        name: the operation's name in the error message
    """
    block = v3000_block(
        after="M  V30 BEGIN RGROUP 1\nM  V30 BEGIN CTAB\nM  V30 COUNTS 1 0 0 0 0\n"
        + "M  V30 BEGIN ATOM\nM  V30 1 C 0 0 0 0\nM  V30 END ATOM\n"
        + "M  V30 END CTAB\nM  V30 END RGROUP\n"
    )
    before = block.to_sdf()
    with pytest.raises(NotImplementedError) as raised:
        getattr(block, operation)()
    assert str(raised.value) == (
        f"{name} does not support CTABs nested in RGROUP or TEMPLATE blocks yet"
    )
    assert block.to_sdf() == before


@pytest.mark.parametrize("ending", ["\n", "\r\n", ""])
def test_continued_lines_are_one_logical_line(ending: str) -> None:
    """Validation and fix_counts join a "-" ended line and the next line.

    The "-" and the next line's "M  V30 " marker, if it has one, are dropped.

    Args:
        ending: what ends each line given to the CTable
    """
    block = v3000_block(
        atom="M  V30 -\nM  V30 1 C 0 0 0 -\nM  V30 0 CHG=1",
        bond="M  V30 1 1 -\n1 2",
        counts="M  V30 COUNTS 9 9 0 0 0",
    )
    lines = [line + ending for line in block.to_sdf().splitlines()[:-1]]
    ctable = CTable(lines)
    ctable.fix_counts()
    assert (ctable.num_atoms, ctable.num_bonds) == (2, 1)
    assert ctable.valid_atom_indices(strict=True)
    assert ctable.valid_bond_indices(strict=True)
    block.fix_counts()
    assert (block.num_atoms, block.num_bonds) == (2, 1)


def test_continued_lines_are_numbered_by_their_first_line() -> None:
    """An error in a logical line names the line it starts on."""
    block = v3000_block(
        atom="M  V30 1 C 0 0 -\nM  V30 0 0\nM  V30 -\nM  V30 1 N 0 0 0 0"
    )
    with pytest.raises(IndicesDuplicateError) as raised:
        block.ctable().valid_atom_indices()
    assert str(raised.value) == "CTable line 10: atom index 1 is a duplicate"


def test_ctable_lines_are_a_fresh_list() -> None:
    """Each CTable.lines access is a new list that does not alias the table."""
    table = v3000_block(atom="a b 1 C").ctable()
    lines = table.lines
    assert lines == table.lines and lines is not table.lines
    lines.clear()
    assert table.lines and table.lines[0] == "t"
