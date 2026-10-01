"""Exercise public CTable index operations against independent cases."""

from __future__ import annotations

import itertools
import sys
from collections import deque
from typing import Any, Callable, Iterable, Iterator, SupportsIndex, cast

import pytest

from sendoff.ctable import CTable
from tests.compat_literals import V2000, V3000
import sendoff.ctable as ctable_module


def outcome(table: object, operation: str, strict: object = False) -> object:
    """Capture the operation result or exact legacy exception.

    Args:
        table: Original or native-backed Python table.
        operation: Validation or renumbering method.
        strict: Uncoerced validation operand.

    Returns:
        Result or exception type, arguments and message.
    """
    try:
        method = getattr(table, operation)
        return method() if operation == "renumber_indices" else method(strict)
    except Exception as error:
        return type(error), error.args, str(error)


def assert_public_behavior(
    lines: Iterable[str], operation: str, strict: object = False
) -> None:
    """Exercise results, raw deque replacement and snapshot state.

    Args:
        lines: Raw connection table lines.
        operation: Operation to compare.
        strict: Uncoerced validation operand.
    """
    raw = deque(lines)
    table = CTable(raw)
    original = table.lines
    result = outcome(table, operation, strict)
    assert type(table.lines) is deque
    assert original == raw
    assert (table.lines is original) == (
        operation != "renumber_indices" or result is not None
    )
    assert type(table.atomlines()) is itertools.takewhile
    assert type(table.bondlines()) is itertools.takewhile


@pytest.mark.parametrize(
    "operation",
    ["valid_atom_indices", "valid_bond_indices", "renumber_indices"],
)
@pytest.mark.parametrize("strict", [False, True])
@pytest.mark.parametrize(
    "change",
    [
        "normal",
        "v2000",
        "duplicate_atoms",
        "duplicate_bonds",
        "missing_atoms",
        "missing_bonds",
        "extra_atoms",
        "extra_bonds",
        "missing_endpoint",
        "duplicate_before_missing_endpoint",
        "bad_atom",
        "bad_bond",
        "short_atom",
        "short_bond",
        "continuation",
        "stale_counts",
        "empty_suffix",
        "unterminated_atoms",
        "unterminated_bonds",
        "opaque",
    ],
)
def test_literal_parity(change: str, strict: bool, operation: str) -> None:
    """Freeze results, diagnostics, check order and all legacy output defects.

    Args:
        change: Focused raw-line mutation.
        strict: Require one-based ordering.
        operation: Validation or renumbering.
    """
    lines = deque((V2000 if change == "v2000" else V3000).splitlines(keepends=True))
    replacements = {
        "duplicate_atoms": (8, lines[7]),
        "missing_endpoint": (11, "M  V30 1 1 99 2\n"),
        "bad_atom": (7, "M  V30 nope nope nope nope"),
        "bad_bond": (11, "M  V30 nope nope nope nope"),
        "short_atom": (7, "short"),
        "short_bond": (11, "short"),
        "stale_counts": (5, "M  V30 COUNTS 5 4 0 0 0\n"),
        "empty_suffix": (5, "M  V30 COUNTS 2 1\n"),
    }
    insertions = {
        "duplicate_bonds": (12, lines[11] if change != "v2000" else ""),
        "extra_atoms": (9, "M  V30 3 N 0 0 0 0\n"),
        "extra_bonds": (12, "M  V30 2 1 1 2 CFG=9\n"),
        "continuation": (8, "M  V30 CHG=1"),
        "opaque": (13, "opaque SGROUP\n"),
    }
    deletions = {"missing_atoms": 8, "missing_bonds": 11}
    truncations = {"unterminated_atoms": 9, "unterminated_bonds": 12}
    if change in replacements:
        position, replacement = replacements[change]
        lines[position] = replacement
    elif change in insertions:
        position, replacement = insertions[change]
        lines.insert(position, replacement)
    elif change in deletions:
        del lines[deletions[change]]
    elif change in truncations:
        lines = deque(list(lines)[: truncations[change]])
    elif change == "duplicate_before_missing_endpoint":
        lines[8] = lines[7]
        lines[11] = "M  V30 1 1 99 1\n"
    assert_public_behavior(lines, operation, strict)


@pytest.mark.parametrize(
    "token", ["-7", "+0007", "7_0", "٧", "７", str(2**257), "-" + str(2**257)]
)
@pytest.mark.parametrize("spacing", [" ", "\t", "\u001c", "\u0085", "\u2003"])
def test_integer_and_unicode_parity(token: str, spacing: str) -> None:
    """Use exact integers and Python whitespace, including big mapped endpoints.

    Args:
        token: Valid arbitrary-precision Python integer spelling.
        spacing: Python-recognized token separator.
    """
    text = V3000.replace("M  V30 1 C", f"M  V30 {token} C").replace(
        "M  V30 1 1 1 2 CFG=1", f"M  V30 1 1 {token} 2 CFG=1"
    )
    lines = text.splitlines(keepends=True)
    lines[7] = lines[7].replace(" ", spacing)
    for operation in ("valid_atom_indices", "valid_bond_indices", "renumber_indices"):
        assert_public_behavior(lines, operation)


def test_digit_limit_and_surrogates() -> None:
    """Keep native int diagnostics and allow lone surrogates in opaque strings."""
    lines = V3000.splitlines()
    lines[7] += "\ud800"
    assert_public_behavior(lines, "renumber_indices")
    token = "9" * (sys.get_int_max_str_digits() + 1)
    lines[7] = f"M  V30 {token} C"
    assert_public_behavior(lines, "valid_atom_indices")
    assert_public_behavior(lines, "renumber_indices")


@pytest.mark.parametrize("counts", ["", "opaque", "\ud800", "M V30 COUNTS 8 9"])
def test_cached_counts_failures_and_snapshots(counts: str) -> None:
    """Use cached counts independently from mutable raw counts lines.

    Args:
        counts: Direct replacement of the parsed count snapshot.
    """
    table = CTable(V3000.splitlines())
    original = table.lines
    table.counts = counts
    result = outcome(table, "renumber_indices")
    assert (table.lines is original) == (result is not None)
    assert table.counts == counts
    assert (table.num_atoms, table.num_bonds) == (2, 1)


def test_renumber_preserves_line_identity_and_call_order() -> None:
    """Keep untouched Python lines and consume each source section once."""
    events: list[str] = []

    class Lines(deque[str]):
        def __iter__(self) -> Iterator[str]:
            events.append("lines")
            return super().__iter__()

    class Receiver(CTable):
        def atomlines(self) -> Iterable[str]:
            events.append("atoms")
            return super().atomlines()

        def bondlines(self) -> Iterable[str]:
            events.append("bonds")
            return super().bondlines()

    table = Receiver(V3000.splitlines())
    table.lines = Lines(table.lines)
    title, source, comment = table.lines[0], table.lines[1], table.lines[2]
    table.renumber_indices()
    assert table.lines[0] is title
    assert table.lines[1] is source
    assert table.lines[2] is comment
    assert events == ["atoms", "lines", "bonds", "lines", "lines"]


def test_empty_tables_do_not_evaluate_strict() -> None:
    """Defer strict truthiness until an actual index has been parsed."""

    class Strict:
        def __bool__(self) -> bool:
            raise AssertionError("strict evaluated")

    lines = V3000.splitlines()
    del lines[11]
    del lines[7:9]
    lines[5] = "M V30 COUNTS 0 0"
    for operation in ("valid_atom_indices", "valid_bond_indices"):
        assert_public_behavior(lines, operation, Strict())


def test_strict_callback_preserves_live_deque_mutation_error() -> None:
    """Do not snapshot raw lines or hide iterator errors from callback mutations."""
    table = CTable(V3000.splitlines())
    original = table.lines

    class Strict:
        def __bool__(self) -> bool:
            table.lines.append("external mutation")
            return False

    result = outcome(table, "valid_atom_indices", Strict())
    assert table.lines is original
    assert result == (
        RuntimeError,
        ("deque mutated during iteration",),
        "deque mutated during iteration",
    )


@pytest.mark.parametrize("operation", ["valid_atom_indices", "renumber_indices"])
def test_string_subclass_protocol_order(operation: str) -> None:
    """Dispatch overrides rather than treating str subclasses as normal strings.

    Args:
        operation: Validation or renumbering.
    """
    events: list[object] = []

    class Tokens(list[str]):
        def __getitem__(self, key: Any) -> Any:
            events.append(("token", key))
            return super().__getitem__(key)

        def __del__(self) -> None:
            events.append("release tokens")

    class Line(str):
        def split(
            self, sep: str | None = None, maxsplit: SupportsIndex = -1
        ) -> list[str]:
            """Record dynamic splitting.

            Args:
                sep: Token separator.
                maxsplit: Maximum splits.

            Returns:
                Observable token sequence.
            """
            events.append(("split", str(self)))
            return Tokens(super().split(sep, maxsplit))

        def __getitem__(self, key: Any) -> Any:
            events.append(("line", key))
            return super().__getitem__(key)

        def startswith(self, prefix: Any, *args: Any) -> bool:
            """Record dynamic prefix checks.

            Args:
                prefix: Prefix to check.
                *args: Optional string bounds.

            Returns:
                Whether the prefix matches.
            """
            events.append(("prefix", prefix))
            return super().startswith(prefix, *args)

    table = CTable(V3000.splitlines())
    table.lines = deque(Line(line) for line in table.lines)
    table.counts = Line(table.counts)
    events.clear()
    result = outcome(table, operation)
    assert result is (None if operation == "renumber_indices" else True)
    assert events


@pytest.mark.parametrize("operation", ["valid_atom_indices", "valid_bond_indices"])
def test_truthiness_and_count_comparison_protocols(operation: str) -> None:
    """Check truthiness per line and compare fresh cached attributes in order.

    Args:
        operation: Atom or bond validation.
    """
    events: list[object] = []

    class Strict:
        def __bool__(self) -> bool:
            events.append("strict")
            return False

    class Count:
        def __gt__(self, size: int) -> bool:
            events.append(("less", size))
            return False

        def __lt__(self, size: int) -> bool:
            events.append(("greater", size))
            return False

    table = CTable(V3000.splitlines())
    setattr(
        table,
        "num_atoms" if operation == "valid_atom_indices" else "num_bonds",
        Count(),
    )
    assert outcome(table, operation, Strict()) is True
    assert events.count("strict") == (2 if operation == "valid_atom_indices" else 1)


def test_uncoerced_integer_and_format_protocols() -> None:
    """Honor custom split sequences, int operands, slices and bond formatting."""
    events: list[object] = []

    class Token:
        def __int__(self) -> int:
            events.append("int")
            return 2**50000

        def __len__(self) -> int:
            events.append("len")
            return 1

    class Order:
        def __format__(self, spec: str) -> str:
            events.append(("format", spec))
            return "opaque"

    class Line(str):
        def split(
            self, sep: str | None = None, maxsplit: SupportsIndex = -1
        ) -> list[Any]:
            """Supply non-string integer and formatting operands.

            Args:
                sep: Token separator.
                maxsplit: Maximum splits.

            Returns:
                Protocol-bearing tokens.
            """
            events.append(("split", str(self)))
            tokens: list[Any] = super().split(sep, maxsplit)
            if " C " in self:
                tokens[2] = Token()
            elif "CFG" in self:
                tokens[3] = Order()
                tokens[4] = Token()
            return tokens

    table = CTable(V3000.splitlines())
    table.lines[7] = Line(table.lines[7])
    table.lines[11] = Line(table.lines[11])
    assert outcome(table, "renumber_indices") is None
    assert events


@pytest.mark.parametrize("stage", ["strict", "split", "producer", "prefix", "setter"])
def test_original_failure_instance_and_atomicity(stage: str) -> None:
    """Return original callback errors without replacing the live raw deque.

    Args:
        stage: Callback failing within validation or renumbering.
    """
    failure = OSError(stage)

    class Strict:
        def __bool__(self) -> bool:
            raise failure

    class Line(str):
        def split(
            self, sep: str | None = None, maxsplit: SupportsIndex = -1
        ) -> list[str]:
            """Inject the selected splitting failure.

            Args:
                sep: Token separator.
                maxsplit: Maximum splits.

            Returns:
                Tokens if splitting succeeds.

            Raises:
                failure: The selected original callback exception.
            """
            if stage == "split":
                raise failure
            return super().split(sep, maxsplit)

        def startswith(self, prefix: Any, *args: Any) -> bool:
            """Inject the selected prefix failure.

            Args:
                prefix: Prefix to check.
                *args: Optional string bounds.

            Returns:
                Whether the prefix matches.

            Raises:
                failure: The selected original callback exception.
            """
            if stage == "prefix":
                raise failure
            return super().startswith(prefix, *args)

    class Receiver(CTable):
        def atomlines(self) -> Iterable[str]:
            """Return a failing external producer when requested.

            Returns:
                Atom lines.

            Raises:
                failure: The selected original callback exception.
            """
            if stage == "producer":
                raise failure
            return super().atomlines()

        def __setattr__(self, name: str, value: object) -> None:
            if stage == "setter" and name == "lines" and hasattr(self, "lines"):
                raise failure
            super().__setattr__(name, value)

    table = Receiver(V3000.splitlines())
    table.lines[7] = Line(table.lines[7])
    original = table.lines
    before = deque(original)
    with pytest.raises(OSError) as caught:
        if stage == "strict":
            validate = cast(Callable[[CTable, object], bool], CTable.valid_atom_indices)
            validate(table, Strict())
        else:
            table.renumber_indices()
    assert caught.value is failure
    assert table.lines is original
    assert table.lines == before


@pytest.mark.parametrize("name", ["valid_atom_indices", "valid_bond_indices"])
@pytest.mark.parametrize(
    "error_name",
    ["IndicesMismatchError", "IndicesOutOfOrderError", "IndicesDuplicateError"],
)
def test_live_module_class_lookup(
    name: str, error_name: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Resolve the current module class after index-processing callbacks.

    Args:
        name: Native validator to exercise.
        error_name: Class replaced during the algorithm.
        monkeypatch: Restore the existing public exception classes afterward.
    """

    class Changed(Exception):
        """A dynamically replaced error class."""

    class Strict:
        def __bool__(self) -> bool:
            monkeypatch.setattr(ctable_module, error_name, Changed)
            return error_name == "IndicesOutOfOrderError"

    table = CTable(V3000.splitlines())
    atoms = name == "valid_atom_indices"
    position = 7 if atoms else 11
    if error_name == "IndicesMismatchError":
        setattr(table, "num_atoms" if atoms else "num_bonds", 99)
    elif error_name == "IndicesOutOfOrderError":
        table.lines[position] = "M V30 9"
    else:
        table.lines.insert(position + 1, table.lines[position])
    with pytest.raises(Changed):
        getattr(table, name)(Strict())
