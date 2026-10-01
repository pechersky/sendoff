"""Freeze raw connection-table behavior without using RDKit as an API oracle."""

from __future__ import annotations

import io
import itertools
from collections import deque
from typing import Any, Callable, Iterator, cast

import pytest

from sendoff.ctable import (
    CTable,
    CTableFormat,
    IndicesDuplicateError,
    IndicesMismatchError,
    IndicesOutOfOrderError,
)
from sendoff.sdblock import SDBlock
from tests.compat_literals import V2000, V3000


def test_constructor_eager_copy_and_raw_header_attributes() -> None:
    """Construction exhausts external iterables, copies lines, and only strips title."""
    for text, expected_format in (
        (V2000, CTableFormat.V2000),
        (V3000, CTableFormat.V3000),
    ):
        raw = deque(text.splitlines(keepends=True))
        source = iter(raw)
        table = CTable(source)
        assert list(source) == []
        assert table.lines == raw
        assert table.lines is not raw
        assert table.title == "literal title"
        assert table.source == "  literal source  \n"
        assert table.comment == " literal comment \n"
        assert table.format is expected_format
        assert (table.num_atoms, table.num_bonds) == (2, 1)
        assert table.counts == (
            raw[3] if expected_format is CTableFormat.V2000 else raw[5].strip()
        )
        raw.clear()
        assert table.lines
    consumed: list[str] = []

    def short_input() -> Iterator[str]:
        for line in ("only title", "source"):
            consumed.append(line)
            yield line

    with pytest.raises(StopIteration) as caught:
        CTable(short_input())
    assert caught.value.args == ()
    assert consumed == ["only title", "source"]


def test_invalid_header_still_exhausts_external_input_before_parsing() -> None:
    """Even a format failure occurs only after materializing the complete iterable."""
    external = iter(["t", "source", "comment", "V4000", "unparsed body", "tail"])
    with pytest.raises(KeyError) as caught:
        CTable(external)
    assert caught.value.args == ("V4000",)
    assert list(external) == []


@pytest.mark.parametrize("size", [0, 1, 2, 3])
def test_short_header_raises_plain_stopiteration(size: int) -> None:
    """Short CTable headers do not acquire new domain-specific errors.

    Args:
        size: Number of available header lines.
    """
    with pytest.raises(StopIteration) as caught:
        CTable(V2000.splitlines()[:size])
    assert str(caught.value) == ""
    assert caught.value.args == ()


@pytest.mark.parametrize(
    ("parser", "line", "error_type", "message"),
    [
        (CTable.parse_format, "", IndexError, "list index out of range"),
        (CTable.parse_format, "  0  0 V4000", KeyError, "'V4000'"),
        (
            CTable.parse_v2000_counts,
            "abc  2 V2000",
            ValueError,
            "invalid literal for int() with base 10: 'abc'",
        ),
        (
            CTable.parse_v2000_counts,
            "  2",
            ValueError,
            "invalid literal for int() with base 10: ''",
        ),
        (
            CTable.parse_v3000_counts,
            "M  V30 COUNTS 2",
            IndexError,
            "list index out of range",
        ),
        (
            CTable.parse_v3000_counts,
            "M  V30 COUNTS nope 1",
            ValueError,
            "invalid literal for int() with base 10: 'nope'",
        ),
    ],
)
def test_counts_parser_error_types_and_exact_messages(
    parser: Callable[[str], Any], line: str, error_type: type[Exception], message: str
) -> None:
    """Freeze malformed count/format diagnostics at their existing boundaries.

    Args:
        parser: Exposed static parser under test.
        line: Malformed raw count line.
        error_type: Current exception class.
        message: Current exact message.
    """
    with pytest.raises(error_type) as caught:
        parser(line)
    assert str(caught.value) == message
    assert caught.value.args == ("V4000" if error_type is KeyError else message,)


def test_counts_parsing_is_positional_and_permissive() -> None:
    """Numeric parsers do not validate chemistry, markers, or nonnegative counts."""
    assert CTable.parse_format("\t arbitrary prefix V3000 \n") is CTableFormat.V3000
    assert CTable.parse_v2000_counts("999001ignored") == (999, 1)
    assert CTable.parse_v2000_counts(" -1 -2 V2000") == (-1, -2)
    assert CTable.parse_v3000_counts("wrong words here -1 +2 ignored") == (-1, 2)
    table = CTable(
        ["t", "s", "c", "  0  0 V3000", "not a begin marker", "wrong words here 0 0"]
    )
    assert table.counts == "wrong words here 0 0"
    assert table.valid_atom_indices()
    assert table.valid_bond_indices()
    with pytest.raises(StopIteration):
        CTable(["t", "s", "c", "  0  0 V3000", "begin"])


def test_v3000_raw_iterators_are_single_use_repeatable_and_marker_based() -> None:
    """Raw iterators are takewhile objects, independent of declared counts."""
    table = CTable(V3000.splitlines(keepends=True))
    atoms = iter(table.atomlines())
    bonds = iter(table.bondlines())
    assert type(atoms) is itertools.takewhile
    assert type(bonds) is itertools.takewhile
    assert iter(atoms) is atoms
    assert list(atoms) == ["M  V30 1 C 0 0 0 0 CHG=1\n", "M  V30 2 O 1 0 0 0\n"]
    assert list(atoms) == []
    assert list(bonds) == ["M  V30 1 1 1 2 CFG=1\n"]
    assert list(bonds) == []
    assert list(table.bondlines()) == ["M  V30 1 1 1 2 CFG=1\n"]
    table.num_atoms = 500
    assert len(list(table.atomlines())) == 2
    del table.lines[9]
    assert list(table.atomlines())[-1] == "M  END\n"


def test_raw_iterator_creation_captures_deque_iteration_before_consumption() -> None:
    """Mutating table.lines invalidates already created iterators, even before next."""
    for method_name in ("atomlines", "bondlines"):
        table = CTable(V3000.splitlines())
        lines = iter(getattr(table, method_name)())
        table.lines.append("changed")
        with pytest.raises(RuntimeError, match="^deque mutated during iteration$"):
            next(lines)
    table = CTable(V3000.splitlines())
    old_atoms = iter(table.atomlines())
    table.lines = deque(V2000.splitlines())
    assert list(old_atoms) == ["M  V30 1 C 0 0 0 0 CHG=1", "M  V30 2 O 1 0 0 0"]
    assert list(table.atomlines()) == ["M  END"]


def test_v2000_raw_iteration_known_failure_and_unimplemented_operations() -> None:
    """Freeze sendoff-712.2 separately from the intentionally unsupported validators."""
    table = CTable(V2000.splitlines(keepends=True))
    assert list(table.atomlines()) == ["M  END\n"]
    assert list(table.bondlines()) == []
    single_atom = CTable(["t", "s", "c", "  1  0 V2000", "one opaque atom", "M  END"])
    assert list(single_atom.atomlines()) == []
    operations: tuple[Callable[[], object], ...] = (
        table.valid_atom_indices,
        table.valid_bond_indices,
        table.renumber_indices,
    )
    for operation in operations:
        with pytest.raises(NotImplementedError) as caught:
            operation()
        assert caught.value.args == ()
        assert str(caught.value) == ""


@pytest.mark.parametrize(
    ("section", "change", "strict", "error_type", "message"),
    [
        (
            "atom",
            "fewer",
            False,
            IndicesMismatchError,
            "fewer atom lines than count line",
        ),
        (
            "atom",
            "more",
            False,
            IndicesMismatchError,
            "more atom lines than count line",
        ),
        (
            "bond",
            "fewer",
            False,
            IndicesMismatchError,
            "fewer bond lines than count line",
        ),
        (
            "bond",
            "more",
            False,
            IndicesMismatchError,
            "more bond lines than count line",
        ),
        ("atom", "duplicate", False, IndicesDuplicateError, "atoms"),
        ("bond", "duplicate", False, IndicesDuplicateError, "bonds"),
        ("atom", "duplicate", True, IndicesOutOfOrderError, "atoms"),
        ("bond", "duplicate", True, IndicesOutOfOrderError, "bonds"),
    ],
)
def test_validation_diagnostics_and_strict_error_precedence(
    section: str, change: str, strict: bool, error_type: type[Exception], message: str
) -> None:
    """Freeze literal failures, including strict ordering before duplicate detection.

    Args:
        section: Atom or bond section to alter.
        change: Focused mismatch or duplication.
        strict: Whether ordering is checked.
        error_type: Existing exception class.
        message: Existing exact error message.
    """
    table = CTable(V3000.splitlines(keepends=True))
    position = 8 if section == "atom" else 11
    if change == "fewer":
        del table.lines[position]
    elif change == "more":
        new_line = "M  V30 3 C 0 0 0 0\n" if section == "atom" else "M  V30 2 1 1 2\n"
        table.lines.insert(position + 1, new_line)
    elif section == "atom":
        table.lines[position] = "M  V30 1 O 1 0 0 0\n"
    else:
        table.lines.insert(position + 1, table.lines[position])
    validator = (
        table.valid_atom_indices if section == "atom" else table.valid_bond_indices
    )
    with pytest.raises(error_type) as caught:
        validator(strict=strict)
    assert caught.value.args == (message,)
    assert str(caught.value) == message


def test_validation_checks_indices_not_atom_chemistry_or_bond_endpoints() -> None:
    """Malformed chemistry is accepted until an index token itself cannot be parsed."""
    table = CTable(V3000.splitlines())
    table.lines[7] = "M  V30 -7 imaginary chemistry"
    table.lines[8] = "M  V30 0 also opaque"
    table.lines[11] = "M  V30 -5 impossible bond endpoints"
    assert table.valid_atom_indices()
    assert table.valid_bond_indices()
    with pytest.raises(IndicesOutOfOrderError, match="^atoms$"):
        table.valid_atom_indices(strict=True)
    for validator, position in (
        (table.valid_atom_indices, 7),
        (table.valid_bond_indices, 11),
    ):
        for raw, error_type, message in (
            ("short", IndexError, "list index out of range"),
            (
                "M  V30 nope C",
                ValueError,
                "invalid literal for int() with base 10: 'nope'",
            ),
        ):
            table.lines[position] = raw
            with pytest.raises(error_type) as caught:
                validator()
            assert str(caught.value) == message


def test_continuation_lines_keep_the_documented_validation_limitation() -> None:
    """Atom continuations are treated as extra index-bearing lines, not joined."""
    table = CTable(V3000.splitlines())
    table.lines[7] = "M  V30 1 C 0 0 0 0 -"
    table.lines.insert(8, "M  V30 CHG=1")
    with pytest.raises(ValueError) as caught:
        table.valid_atom_indices()
    assert str(caught.value) == "invalid literal for int() with base 10: 'CHG=1'"
    original = table.lines
    with pytest.raises(ValueError, match="CHG=1"):
        table.renumber_indices()
    assert table.lines is original


def test_ctable_does_not_wrap_external_producer_errors() -> None:
    """Eager materialization propagates producer failures before header parsing."""

    def failing_lines() -> Iterator[str]:
        yield "title"
        raise OSError("producer failed")

    with pytest.raises(OSError, match="^producer failed$"):
        CTable(failing_lines())


def test_direct_line_edits_leave_cached_attributes_stale_until_reparsing() -> None:
    """Parsed table headers/counts are snapshots; raw iterators use live lines."""
    table = CTable(V3000.splitlines(keepends=True))
    table.lines[0] = "new title\n"
    table.lines[5] = "M  V30 COUNTS 9 8 0 0 0\n"
    table.lines[7] = "M  V30 7 N 0 0 0 0\n"
    assert table.title == "literal title"
    assert (table.num_atoms, table.num_bonds) == (2, 1)
    assert next(iter(table.atomlines())) == "M  V30 7 N 0 0 0 0\n"
    reparsed = CTable(table.lines)
    assert reparsed.title == "new title"
    assert (reparsed.num_atoms, reparsed.num_bonds) == (9, 8)
    table.num_atoms = 9
    with pytest.raises(
        IndicesMismatchError, match="^fewer atom lines than count line$"
    ):
        table.valid_atom_indices()
    table.lines = deque(V2000.splitlines())
    assert table.format is CTableFormat.V3000
    assert CTable(table.lines).format is CTableFormat.V2000


def test_renumber_exact_lines_property_loss_and_stale_counts() -> None:
    """Freeze sendoff-712.3/.4, raw counts newline loss, and deque replacement."""
    text = (
        V3000.replace("COUNTS 2 1", "COUNTS 5 4")
        .replace("M  V30 1 C", "M  V30 7  C")
        .replace("M  V30 2 O", "M  V30 12 O")
        .replace("M  V30 1 1 1 2 CFG=1", "M  V30 9   2 7 12 CFG=1")
    )
    table = CTable(text.splitlines(keepends=True))
    old_lines = table.lines
    old_counts = table.counts
    assert cast(Callable[[], object], table.renumber_indices)() is None
    expected = deque(V3000.splitlines(keepends=True))
    expected[5] = "M  V30 COUNTS 2 1 0 0 0"
    expected[7] = "M  V30 1  C 0 0 0 0 CHG=1\n"
    expected[11] = "M  V30 1 2 1 2\n"
    assert table.lines == expected
    assert table.lines is not old_lines
    assert old_lines == deque(text.splitlines(keepends=True))
    assert table.counts == old_counts == "M  V30 COUNTS 5 4 0 0 0"
    assert (table.num_atoms, table.num_bonds) == (5, 4)
    with pytest.raises(
        IndicesMismatchError, match="^fewer atom lines than count line$"
    ):
        table.valid_atom_indices(strict=True)
    reparsed = CTable(table.lines)
    assert (reparsed.num_atoms, reparsed.num_bonds) == (2, 1)
    assert reparsed.valid_atom_indices(strict=True)
    assert reparsed.valid_bond_indices(strict=True)


def test_block_renumber_replaces_mdl_preserves_metadata_and_reparses_counts() -> None:
    """Blocks differ from direct table use by reconstructing count snapshots."""
    block = SDBlock.from_block_lines(
        V3000.replace("COUNTS 2 1", "COUNTS 5 4").splitlines(keepends=True)
    )
    old_mdl = block.mdl
    old_metadata = block.metadata
    old_title = block.title
    assert cast(Callable[[], object], block.renumber_indices)() is None
    assert block.mdl is not old_mdl
    assert block.metadata is old_metadata
    assert block.title == old_title
    assert block.mdl[0] == "  literal source  \n"
    assert (block.num_atoms(), block.num_bonds()) == (2, 1)
    assert old_mdl[4] == "M  V30 COUNTS 5 4 0 0 0\n"


@pytest.mark.parametrize(
    ("replacement", "error_type", "message"),
    [
        ("M  V30 1 O 1 0 0 0", IndicesDuplicateError, "atom index mapping in bond"),
        ("M  V30 7 O 1 0 0 0", IndexError, "list index out of range"),
    ],
)
def test_renumber_failures_do_not_replace_original_lines(
    replacement: str, error_type: type[Exception], message: str
) -> None:
    """Duplicate mapped atoms and missing endpoints fail without committing new lines.

    Args:
        replacement: Atom line causing ambiguous or absent bond endpoint mapping.
        error_type: Existing exception class.
        message: Existing exact diagnostic.
    """
    table = CTable(V3000.splitlines())
    table.lines[8] = replacement
    original = table.lines
    before = deque(original)
    with pytest.raises(error_type) as caught:
        table.renumber_indices()
    assert str(caught.value) == message
    assert caught.value.args == (message,)
    assert table.lines is original
    assert table.lines == before


def test_unused_duplicate_atom_indices_can_be_renumbered() -> None:
    """Duplicates without bond references are repairable, unlike mapped duplicates."""
    table = CTable(V3000.splitlines())
    table.lines[8] = "M  V30 1 O 1 0 0 0"
    del table.lines[11]
    assert cast(Callable[[], object], table.renumber_indices)() is None
    reparsed = CTable(table.lines)
    assert reparsed.valid_atom_indices(strict=True)
    assert reparsed.valid_bond_indices(strict=True)
    assert (reparsed.num_atoms, reparsed.num_bonds) == (2, 0)


def test_renumber_retains_opaque_sections_and_handles_absent_count_suffix() -> None:
    """Opaque sections survive; a short counts suffix becomes one trailing space."""
    table = CTable(V3000.replace("COUNTS 2 1 0 0 0", "COUNTS 2 1").splitlines())
    extension = ["M  V30 BEGIN SGROUP", "opaque extension", "M  V30 END SGROUP"]
    table.lines.insert(13, extension[0])
    table.lines.insert(14, extension[1])
    table.lines.insert(15, extension[2])
    table.renumber_indices()
    assert table.lines[5] == "M  V30 COUNTS 2 1 "
    assert list(table.lines)[13:16] == extension
    assert table.lines[-1] == "M  END"


def test_renumber_raw_write_loses_counts_newline() -> None:
    """Keep the separately deferred raw-output defect, even on valid indices."""
    block = SDBlock.from_block_lines(V3000.splitlines(keepends=True))
    block.renumber_indices()
    output = io.StringIO()
    block.write(cast(Any, output), with_newlines=False)
    assert output.getvalue() == (
        "literal title\n"
        "  literal source  \n"
        " literal comment \n"
        "  0  0  0  0  0  0  0  0  0  0999 V3000\n"
        "M  V30 BEGIN CTAB\n"
        "M  V30 COUNTS 2 1 0 0 0M  V30 BEGIN ATOM\n"
        "M  V30 1 C 0 0 0 0 CHG=1\n"
        "M  V30 2 O 1 0 0 0\n"
        "M  V30 END ATOM\n"
        "M  V30 BEGIN BOND\n"
        "M  V30 1 1 1 2\n"
        "M  V30 END BOND\n"
        "M  V30 END CTAB\n"
        "M  END\n"
        "$$$$\n"
    )
