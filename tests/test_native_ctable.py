"""Exercise real native CTable exports against the frozen Python protocols."""

from __future__ import annotations

import builtins
import inspect
import itertools
import pickle
from collections import deque
from decimal import Decimal
from importlib import import_module
from typing import Any, Callable, Iterable, Iterator, SupportsIndex, Tuple

import pytest

from sendoff.ctable import CTable, CTableFormat, IndicesMismatchError
from tests.compat_literals import V2000, V3000

native = import_module("sendoff.native")
init: Callable[..., None] = getattr(native, "ctable_init")
parse_format: Callable[..., CTableFormat] = getattr(native, "parse_format")
v2000_counts: Callable[..., Tuple[int, int]] = getattr(native, "parse_v2000_counts")
v3000_counts: Callable[..., Tuple[int, int]] = getattr(native, "parse_v3000_counts")
atomlines: Callable[..., Iterator[str]] = getattr(native, "atomlines")
bondlines: Callable[..., Iterator[str]] = getattr(native, "bondlines")


class NativeTable(CTable):
    """Provide test-only thin adapters without editing the public facade."""

    def __init__(self, lines: Iterable[str]) -> None:
        init(self, lines, CTableFormat.V3000)

    @staticmethod
    def parse_format(line: str) -> CTableFormat:
        """Delegate format parsing to Rust.

        Args:
            line: counts line

        Returns:
            Parsed format.
        """
        return parse_format(line, CTableFormat)

    @staticmethod
    def parse_v2000_counts(line: str) -> Tuple[int, int]:
        """Delegate V2000 counts to Rust.

        Args:
            line: counts line

        Returns:
            Atom and bond counts.
        """
        return v2000_counts(line)

    @staticmethod
    def parse_v3000_counts(line: str) -> Tuple[int, int]:
        """Delegate V3000 counts to Rust.

        Args:
            line: counts line

        Returns:
            Atom and bond counts.
        """
        return v3000_counts(line)

    def atomlines(self) -> Iterable[str]:
        """Return the native atom-line iterator.

        Returns:
            Atom lines.
        """
        return atomlines(self)

    def bondlines(self) -> Iterable[str]:
        """Return the native bond-line iterator.

        Returns:
            Bond lines.
        """
        return bondlines(self)


def test_plain_counts_and_indices_do_not_call_python_int(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Require ordinary counts and index processing to remain in Rust.

    Args:
        monkeypatch: temporary Python integer-constructor guard
    """
    table = CTable(V3000.splitlines())
    table.num_atoms = len(list(table.atomlines()))
    table.num_bonds = len(list(table.bondlines()))

    def fail(value: object) -> int:
        raise AssertionError("Python int called on normal input")

    with monkeypatch.context() as boundary:
        boundary.setattr(builtins, "int", fail)
        counts = (
            CTable.parse_v2000_counts("002001 V2000"),
            CTable.parse_v3000_counts("M V30 COUNTS 2 1 0 0 0"),
        )
        table.valid_atom_indices()
        table.valid_bond_indices()
        table.renumber_indices()
    assert counts == ((2, 1), (2, 1))
    assert type(table.num_atoms) is type(table.num_bonds) is int


@pytest.mark.parametrize(
    "line,v3000",
    [
        ("002001 V2000", False),
        ("1_22_3", False),
        (" \u0661  \uff12 ", False),
        ("\x1c1  2 ", False),
        ("\ud80001 02", False),
        ("", False),
        ("001", False),
        ("M V30 COUNTS 2 1", True),
        ("M\x1cV30\x1dCOUNTS\x1e2\x1f1", True),
        ("M V30 COUNTS 1_2 \u0663", True),
        (f"M V30 COUNTS {2**140} {-2**160}", True),
        ("M V30 COUNTS + 1", True),
        ("M V30 COUNTS invalid", True),
        ("M V30 COUNTS", True),
        ("M V30 COUNTS 1", True),
        ("M V30 COUNTS \ud800 1", True),
    ],
)
def test_rust_counts_keep_python_integer_and_slice_semantics(
    line: str, v3000: bool
) -> None:
    """Compare against independent Python expressions, including failures.

    Args:
        line: valid, malformed or Unicode counts line
        v3000: whitespace-delimited rather than fixed-width counts
    """
    reference: Callable[[], Tuple[int, int]] = (
        (lambda: (int(line.split()[3]), int(line.split()[4])))
        if v3000
        else (lambda: (int(line[:3]), int(line[3:6])))
    )
    parser = v3000_counts if v3000 else v2000_counts
    implementations: tuple[Callable[[], Tuple[int, int]], ...] = (
        reference,
        lambda: parser(line),
    )
    outcomes: list[object] = []
    for implementation in implementations:
        try:
            outcomes.append(implementation())
        except (ValueError, IndexError) as error:
            outcomes.append((type(error), error.args))
    assert outcomes[0] == outcomes[1]


def test_unchanged_strip_preserves_raw_string_identity() -> None:
    """Keep Python str.strip identity while still trimming changed fields."""
    lines = V3000.splitlines()
    lines[0] = "".join(["unchanged ", "title"])
    table = CTable(lines)
    assert table.title is table.lines[0] is lines[0]
    assert table.counts is lines[5]

    lines[0] = "\x1c" + lines[0] + "\x1d"
    lines[5] = "\x1e" + lines[5] + "\x1f\n"
    trimmed = CTable(lines)
    assert trimmed.title == table.title
    assert trimmed.title is not lines[0]
    assert trimmed.counts == table.counts
    assert trimmed.counts is not lines[5]


def test_exports_and_keyword_arguments() -> None:
    """Expose required positional-or-keyword operands, not coercing signatures."""
    expected = {
        "ctable_init": ("table", "lines", "v3000"),
        "parse_format": ("line", "formats"),
        "parse_v2000_counts": ("line",),
        "parse_v3000_counts": ("line",),
        "atomlines": ("table",),
        "bondlines": ("table",),
    }
    for name, operands in expected.items():
        function = getattr(native, name)
        assert inspect.isbuiltin(function)
        parameters = inspect.signature(function).parameters
        assert tuple(parameters) == operands
        assert all(
            parameter.kind is inspect.Parameter.POSITIONAL_OR_KEYWORD
            and parameter.default is inspect.Parameter.empty
            for parameter in parameters.values()
        )
    table = CTable.__new__(CTable)
    assert (
        getattr(native, "ctable_init")(
            table=table, lines=V3000.splitlines(), v3000=CTableFormat.V3000
        )
        is None
    )
    assert parse_format(line="V2000", formats=CTableFormat) is CTableFormat.V2000
    assert v2000_counts(line="002001") == (2, 1)
    assert v3000_counts(line="M V30 COUNTS 2 1") == (2, 1)
    assert type(atomlines(table=table)) is itertools.takewhile
    assert type(bondlines(table=table)) is itertools.takewhile


@pytest.mark.parametrize("text", [V2000, V3000])
def test_eager_copy_headers_snapshots_and_public_pickle(text: str) -> None:
    """Copy before parsing and leave original public objects and snapshots intact.

    Args:
        text: The existing literal connection table.
    """
    raw = deque(text.splitlines(keepends=True))
    producer = iter(raw)
    table = CTable.__new__(CTable)
    init(table, producer, CTableFormat.V3000)
    assert list(producer) == []
    assert vars(table) == vars(CTable(raw))
    assert table.lines is not raw
    assert table.format is (CTableFormat.V2000 if text is V2000 else CTableFormat.V3000)
    assert type(table.num_atoms) is type(table.num_bonds) is int
    old_fields = {key: value for key, value in vars(table).items() if key != "lines"}
    table.lines[0] = "new title\n"
    table.lines[3 if table.format is CTableFormat.V2000 else 5] = (
        "009008 V2000\n" if table.format is CTableFormat.V2000 else "M V30 COUNTS 9 8\n"
    )
    table.lines[7] = "replacement body\n"
    assert {key: value for key, value in vars(table).items() if key != "lines"} == (
        old_fields
    )
    reparsed = NativeTable(table.lines)
    assert (reparsed.title, reparsed.num_atoms, reparsed.num_bonds) == (
        "new title",
        9,
        8,
    )
    restored = pickle.loads(pickle.dumps(table))
    assert type(restored) is CTable
    assert restored.format is table.format
    assert vars(restored) == vars(table)
    error = IndicesMismatchError("unchanged")
    assert type(pickle.loads(pickle.dumps(error))) is IndicesMismatchError
    raw.clear()
    assert table.lines


@pytest.mark.parametrize(
    ("lines", "error_type"),
    [
        *[(V2000.splitlines()[:size], StopIteration) for size in range(4)],
        *[(V3000.splitlines()[:size], StopIteration) for size in (4, 5)],
        (["t", "s", "c", ""], IndexError),
        (["t", "s", "c", "V4000", "tail"], KeyError),
        (["t", "s", "c", "abc  2 V2000"], ValueError),
        (["t", "s", "c", "V3000", "ignored", "M V30 COUNTS 1"], IndexError),
        (["t", "s", "c", "V3000", "ignored", "M V30 COUNTS nope 1"], ValueError),
        ([1, "s", "c", "002001 V2000"], AttributeError),
    ],
)
def test_malformed_constructor_partial_state(
    lines: list[object], error_type: type[Exception]
) -> None:
    """Keep exact diagnostics, assignment order and eager external exhaustion.

    Args:
        lines: Malformed or short constructor input.
        error_type: The legacy builtin error class.
    """
    tables = [CTable.__new__(CTable), CTable.__new__(CTable)]
    errors = []
    initializers: tuple[Callable[..., None], ...] = (
        CTable.__init__,
        NativeTable.__init__,
    )
    for table, initializer in zip(tables, initializers):
        producer = iter(lines)
        with pytest.raises(error_type) as caught:
            initializer(table, producer)
        errors.append(caught.value)
        assert list(producer) == []
    assert type(errors[0]) is type(errors[1])
    assert errors[0].args == errors[1].args
    assert str(errors[0]) == str(errors[1])
    assert list(vars(tables[0])) == list(vars(tables[1]))
    assert vars(tables[0]) == vars(tables[1])
    if error_type is StopIteration:
        assert errors[1].args == ()
        assert errors[1].__cause__ is None


@pytest.mark.parametrize(
    "field",
    [
        "lines",
        "title",
        "source",
        "comment",
        "counts",
        "format",
        "num_atoms",
        "num_bonds",
    ],
)
def test_receiver_assignment_failure_and_attribute_order(field: str) -> None:
    """Preserve dynamic lookup order, setter errors and partial assignments.

    Args:
        field: The constructor assignment that raises the original exception.
    """
    events: list[tuple[str, str]] = []
    failure = OSError("setter failed")

    class Receiver(CTable):
        def __getattribute__(self, name: str) -> Any:
            events.append(("get", name))
            return super().__getattribute__(name)

        def __setattr__(self, name: str, value: object) -> None:
            events.append(("set", name))
            if name == field:
                raise failure
            super().__setattr__(name, value)

    initializers: tuple[Callable[..., None], ...] = (
        CTable.__init__,
        NativeTable.__init__,
    )
    states, lookups = [], []
    for initializer in initializers:
        events.clear()
        table = Receiver.__new__(Receiver)
        with pytest.raises(OSError) as caught:
            initializer(table, V3000.splitlines())
        assert caught.value is failure
        lookups.append(events.copy())
        states.append(vars(table))
    assert lookups[0] == lookups[1]
    assert states[0] == states[1]
    assert list(states[0]) == list(states[1])


def test_producer_errors_and_consumption_precede_all_assignments() -> None:
    """Propagate the identical producer error without replacing existing state."""
    table = CTable(V2000.splitlines())
    old = vars(table).copy()
    failure = OSError("producer failed")

    def producer() -> Iterator[str]:
        yield "title"
        raise failure

    with pytest.raises(OSError) as caught:
        init(table, producer(), CTableFormat.V3000)
    assert caught.value is failure
    assert vars(table) == old
    assert table.lines is old["lines"]
    consumed: list[str] = []

    class Receiver(CTable):
        @staticmethod
        def parse_format(line: str) -> CTableFormat:
            assert consumed == ["t", "s", "c", "V2000", "tail"]
            raise failure

    def complete() -> Iterator[str]:
        for line in ["t", "s", "c", "V2000", "tail"]:
            consumed.append(line)
            yield line

    receiver = Receiver.__new__(Receiver)
    with pytest.raises(OSError) as caught:
        init(receiver, complete(), CTableFormat.V3000)
    assert caught.value is failure
    assert list(vars(receiver)) == ["lines", "title", "source", "comment", "counts"]


@pytest.mark.parametrize("format_value", [CTableFormat.V3000, "V3000", None])
def test_dynamic_parser_identity_dispatch_and_arbitrary_unpack(
    format_value: object,
) -> None:
    """Respect subclass results and complete unpacking before numeric assignments.

    Args:
        format_value: A real enum or a nonidentical result from a subclass hook.
    """
    events: list[object] = []
    atoms, bonds = object(), object()

    class Pair(tuple[Any, ...]):
        def __del__(self) -> None:
            events.append("pair released")

        def __iter__(self) -> Iterator[object]:
            events.append("iter")
            for value in (atoms, bonds):
                events.append(value)
                yield value
            events.append("exhausted")

    class Receiver(CTable):
        @staticmethod
        def parse_format(line: str) -> Any:
            events.append(("format", line))
            return format_value

        @staticmethod
        def parse_v2000_counts(line: str) -> Any:
            events.append(("v2000", line))
            return Pair((99, 88))

        @staticmethod
        def parse_v3000_counts(line: str) -> Any:
            events.append(("v3000", line))
            return Pair((99, 88))

        def __setattr__(self, name: str, value: object) -> None:
            events.append(("set", name))
            super().__setattr__(name, value)

    receiver = Receiver.__new__(Receiver)
    init(receiver, V3000.splitlines(keepends=True), CTableFormat.V3000)
    counts = (
        V3000.splitlines()[5]
        if format_value is CTableFormat.V3000
        else V3000.splitlines(keepends=True)[3]
    )
    assert (
        "v3000" if format_value is CTableFormat.V3000 else "v2000",
        counts,
    ) in events
    assert events[-7:] == [
        "iter",
        atoms,
        bonds,
        "exhausted",
        "pair released",
        ("set", "num_atoms"),
        ("set", "num_bonds"),
    ]
    assert receiver.num_atoms is atoms
    assert receiver.num_bonds is bonds


@pytest.mark.parametrize(
    ("result", "error_type"),
    [
        (None, TypeError),
        (7, TypeError),
        (Decimal("1"), TypeError),
        ({}, ValueError),
        ([], ValueError),
        ([1], ValueError),
        ([1, 2, 3, 4], ValueError),
        ((1, 2, 3), ValueError),
        ({"a": 1, "b": 2, "c": 3}, ValueError),
        ("abc", ValueError),
        (type("NotIterable", (), {})(), TypeError),
        (type("Disabled", (), {"__iter__": None})(), TypeError),
        (type("é" * 110, (), {})(), TypeError),
        (type("x" + "é" * 110, (), {})(), TypeError),
    ],
)
def test_count_hook_unpack_diagnostics(
    result: object, error_type: type[Exception]
) -> None:
    """Do not coerce subclass count results or assign a partial unpacked pair.

    Args:
        result: An arbitrary count hook result that cannot unpack into two values.
        error_type: The precise builtin unpack failure.
    """

    class Receiver(CTable):
        @staticmethod
        def parse_v2000_counts(line: str) -> Any:
            return result

    tables = [Receiver.__new__(Receiver), Receiver.__new__(Receiver)]
    errors = []
    initializers: tuple[Callable[..., None], ...] = (
        CTable.__init__,
        NativeTable.__init__,
    )
    for table, initializer in zip(tables, initializers):
        with pytest.raises(error_type) as caught:
            initializer(table, V2000.splitlines())
        errors.append(caught.value)
    assert type(errors[0]) is type(errors[1])
    assert errors[0].args == errors[1].args
    assert str(errors[0]) == str(errors[1])
    assert vars(tables[0]) == vars(tables[1])
    assert not hasattr(tables[1], "num_atoms")


@pytest.mark.parametrize("step", [0, 1, 2, 3])
@pytest.mark.parametrize("error_type", [OSError, TypeError, StopIteration])
def test_count_hook_iteration_errors(step: int, error_type: type[Exception]) -> None:
    """Match iterator creation, two-item lookahead, and original callback errors.

    Args:
        step: The iteration call that fails, where zero means iterator creation.
        error_type: An external exception, including unpack exhaustion.
    """
    failure = error_type("callback failed")
    events: list[int] = []

    class Pair:
        def __iter__(self) -> Iterator[object]:
            events.append(0)
            if step == 0:
                raise failure
            return self

        def __next__(self) -> object:
            events.append(len(events))
            if len(events) - 1 == step:
                raise failure
            return object()

    class Receiver(CTable):
        @staticmethod
        def parse_v2000_counts(line: str) -> Any:
            return Pair()

    initializers: tuple[Callable[..., None], ...] = (
        CTable.__init__,
        NativeTable.__init__,
    )
    for initializer in initializers:
        events.clear()
        receiver = Receiver.__new__(Receiver)
        if error_type is StopIteration and step == 3:
            initializer(receiver, V2000.splitlines())
            assert hasattr(receiver, "num_bonds")
        else:
            expected = error_type
            if error_type is StopIteration and step:
                expected = ValueError
            with pytest.raises(expected) as caught:
                initializer(receiver, V2000.splitlines())
            if expected is error_type:
                assert caught.value is failure
            else:
                assert str(caught.value) == (
                    f"not enough values to unpack (expected 2, got {step - 1})"
                )
            assert not hasattr(receiver, "num_atoms")
        assert events == list(range(step + 1))


@pytest.mark.parametrize(
    ("parser", "legacy", "line", "expected"),
    [
        (parse_format, CTable.parse_format, "\tignored V3000 \n", CTableFormat.V3000),
        (v2000_counts, CTable.parse_v2000_counts, "999001ignored", (999, 1)),
        (v2000_counts, CTable.parse_v2000_counts, " -1 -2 ignored", (-1, -2)),
        (v2000_counts, CTable.parse_v2000_counts, " ٢٣ ４５", (23, 45)),
        (v3000_counts, CTable.parse_v3000_counts, "wrong words here -1 +2", (-1, 2)),
        (v3000_counts, CTable.parse_v3000_counts, "a b c ٢٣ ４５", (23, 45)),
        (
            v3000_counts,
            CTable.parse_v3000_counts,
            f"a b c {10**100} 1_000",
            (10**100, 1000),
        ),
    ],
)
def test_positional_permissive_unicode_and_unbounded_counts(
    parser: Callable[..., object],
    legacy: Callable[[str], object],
    line: str,
    expected: object,
) -> None:
    """Delegate exact slicing, Unicode digits and integer semantics to Python.

    Args:
        parser: The real native export.
        legacy: The frozen parser.
        line: Permissive count or format input.
        expected: The original Python value or enum.
    """
    value = parser(line, CTableFormat) if parser is parse_format else parser(line)
    assert value == legacy(line) == expected
    if parser is parse_format:
        assert value is expected
    else:
        assert type(value) is tuple
        assert all(type(number) is int for number in value)


@pytest.mark.parametrize(
    ("parser", "legacy", "line", "error_type"),
    [
        (parse_format, CTable.parse_format, "", IndexError),
        (parse_format, CTable.parse_format, "V4000", KeyError),
        (parse_format, CTable.parse_format, "\ud800", KeyError),
        (v2000_counts, CTable.parse_v2000_counts, "abc  2", ValueError),
        (v2000_counts, CTable.parse_v2000_counts, "  2", ValueError),
        (v2000_counts, CTable.parse_v2000_counts, "\ud800  1", ValueError),
        (v3000_counts, CTable.parse_v3000_counts, "M V30 COUNTS 2", IndexError),
        (v3000_counts, CTable.parse_v3000_counts, "M V30 COUNTS nope 1", ValueError),
        (v3000_counts, CTable.parse_v3000_counts, "a b c \ud800 1", ValueError),
        (
            v3000_counts,
            CTable.parse_v3000_counts,
            "a b c " + "1" * 5000 + " 1",
            ValueError,
        ),
        (v2000_counts, CTable.parse_v2000_counts, b"002001", None),
        (v3000_counts, CTable.parse_v3000_counts, b"a b c 2 1", None),
        (parse_format, CTable.parse_format, 1, AttributeError),
        (v2000_counts, CTable.parse_v2000_counts, 1, TypeError),
        (v3000_counts, CTable.parse_v3000_counts, 1, AttributeError),
    ],
)
def test_parser_protocols_and_exact_malformed_diagnostics(
    parser: Callable[..., object],
    legacy: Callable[..., object],
    line: object,
    error_type: type[Exception] | None,
) -> None:
    """Retain permissive Python object operands and builtin error text.

    Args:
        parser: The native parser.
        legacy: The original parser.
        line: Malformed or non-string input.
        error_type: The builtin error, or None for accepted bytes.
    """
    args = (line, CTableFormat) if parser is parse_format else (line,)
    if error_type is None:
        assert parser(*args) == legacy(line)
    else:
        with pytest.raises(error_type) as reference:
            legacy(line)
        with pytest.raises(error_type) as actual:
            parser(*args)
        assert type(actual.value) is type(reference.value)
        assert actual.value.args == reference.value.args
        assert str(actual.value) == str(reference.value)


def test_dynamic_primitive_hooks_and_exact_slice_objects() -> None:
    """Call subclass strip/split/index hooks while bypassing only startswith hooks."""
    events: list[object] = []

    class Text(str):
        def strip(self, chars: str | None = None) -> Text:
            events.append("strip")
            return self

        def split(self, sep: str | None = None, maxsplit: SupportsIndex = -1) -> Any:
            events.append("split")
            return Tokens()

        def __getitem__(self, key: Any) -> Any:
            events.append(key)
            return " ٢" if key.start is None else "-3"

    class Tokens:
        def __getitem__(self, key: int) -> str:
            events.append(key)
            return {-1: "V3000", 3: "٢", 4: "-3"}[key]

    assert parse_format(Text("ignored"), CTableFormat) is CTableFormat.V3000
    assert events == ["strip", "split", -1]
    events.clear()
    assert v2000_counts(Text("ignored")) == (2, -3)
    assert events == [slice(None, 3, None), slice(3, 6, None)]
    events.clear()
    assert v3000_counts(Text("ignored")) == (2, -3)
    assert events == ["split", 3, 4]
    sentinel = object()

    class Formats:
        def __getitem__(self, name: str) -> object:
            assert name == "V3000"
            return sentinel

    assert parse_format(Text("ignored"), Formats()) is sentinel


@pytest.mark.parametrize("method", ["atomlines", "bondlines"])
@pytest.mark.parametrize("text", [V2000, V3000])
def test_real_single_use_raw_iterators_and_preserved_v2000_defect(
    method: str, text: str
) -> None:
    """Keep exact itertools identity, independent calls, ignored declared counts.

    Args:
        method: The raw section to traverse.
        text: The existing V2000 or V3000 fixture.
    """
    table = NativeTable(text.splitlines(keepends=True))
    reference = CTable(table.lines)
    expected = list(getattr(reference, method)())
    table.num_atoms = table.num_bonds = -999
    raw = getattr(table, method)()
    assert type(raw) is itertools.takewhile
    assert iter(raw) is raw
    assert list(raw) == expected
    assert list(raw) == []
    assert list(getattr(table, method)()) == expected
    assert table.lines == reference.lines
    if text is V2000:
        assert expected == (["M  END\n"] if method == "atomlines" else [])


@pytest.mark.parametrize("method", ["atomlines", "bondlines"])
@pytest.mark.parametrize("started", [False, True])
def test_call_time_capture_mutation_and_reassignment(
    method: str, started: bool
) -> None:
    """Capture old deques immediately and preserve their real mutation detection.

    Args:
        method: The section iterator factory.
        started: Whether to consume a value before structural mutation.
    """
    table = NativeTable(V3000.splitlines())
    raw = iter(getattr(table, method)())
    if started:
        next(raw)
    table.lines.append("mutation")
    with pytest.raises(RuntimeError) as caught:
        next(raw)
    assert caught.value.args == ("deque mutated during iteration",)
    table = NativeTable(V3000.splitlines())
    old = iter(getattr(table, method)())
    expected = list(getattr(CTable(table.lines), method)())
    table.lines = deque(V2000.splitlines())
    assert list(old) == expected
    assert list(getattr(table, method)()) == (
        ["M  END"] if method == "atomlines" else []
    )


@pytest.mark.parametrize("method", ["atomlines", "bondlines"])
def test_raw_prefixes_spill_marker_consumption_and_protocol_timing(method: str) -> None:
    """Avoid snapshots, consume stopping markers and call the str descriptor.

    Args:
        method: The raw iterator factory.
    """
    function = atomlines if method == "atomlines" else bondlines

    class Text(str):
        def startswith(self, prefix: Any, *args: Any) -> bool:
            raise AssertionError("instance startswith must not be called")

    lines = V3000.splitlines()
    lines[7] = Text("opaque \ud800 atom")
    lines[9] = Text("M  V30 END ATOM suffix")
    lines[10] = Text("M  V30 BEGIN BOND suffix")
    lines[12] = Text("M  V30 END BOND suffix")
    table = CTable.__new__(CTable)
    source = iter(lines)
    calls: list[str] = []

    class Lines:
        def __iter__(self) -> Iterator[str]:
            calls.append("iter")
            return source

    setattr(table, "lines", Lines())
    raw = function(table)
    assert calls == ["iter"]
    assert list(raw) == (lines[7:9] if method == "atomlines" else lines[11:12])
    assert next(source) == (lines[10] if method == "atomlines" else lines[13])
    assert list(raw) == []
    table.lines = deque(line for line in lines if "END ATOM" not in line)
    assert list(atomlines(table))[-1] == "M  END"
    table.lines = deque(line for line in lines if "BEGIN BOND" not in line)
    assert list(bondlines(table)) == []
    table.lines = deque(lines)
    table.lines[9] = " M  V30 END ATOM"
    table.lines[12] = " M  V30 END BOND"
    assert list(function(table))[-1] == "M  END"
    setattr(table, "lines", deque([*lines[:7], 1, *lines[8:]]))
    actual = function(table)
    reference = getattr(CTable, method)(table)
    with pytest.raises(TypeError) as expected:
        next(iter(reference))
    with pytest.raises(TypeError) as caught:
        next(actual)
    assert caught.value.args == expected.value.args
    assert str(caught.value) == str(expected.value)
