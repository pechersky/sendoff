"""Characterize owned public Python surfaces, including accidental behavior."""

from __future__ import annotations

import inspect
import pickle
from collections import deque
from dataclasses import asdict, fields, is_dataclass, replace
from enum import Enum
from typing import Any, Iterable, cast, get_args

import pytest

import sendoff
import sendoff.ctable as ctable_module
import sendoff.sdblock as sdblock_module
from sendoff.ctable import (
    CTable,
    CTableFormat,
    IndicesDuplicateError,
    IndicesMismatchError,
    IndicesOutOfOrderError,
)
from sendoff.sdblock import Pathy, SDBlock, parse_sdf
from tests.compat_literals import V2000


def test_owned_module_paths_and_typing() -> None:
    """Preserve owned imports without inventing top-level class reexports."""
    assert isinstance(sendoff.__version__, str)
    assert not hasattr(sendoff, "SDBlock")
    assert not hasattr(sendoff, "CTable")
    assert SDBlock.__module__ == "sendoff.sdblock"
    assert parse_sdf.__module__ == "sendoff.sdblock"
    assert ctable_module.CTable is getattr(sdblock_module, "CTable") is CTable
    for exported in (
        CTable,
        CTableFormat,
        IndicesDuplicateError,
        IndicesMismatchError,
        IndicesOutOfOrderError,
    ):
        assert exported.__module__ == "sendoff.ctable"
        assert getattr(ctable_module, exported.__name__) is exported
    assert sdblock_module.Pathy is Pathy
    assert parse_sdf.__annotations__ == {
        "sdfpath": "Pathy",
        "return": "Iterator[SDBlock]",
    }
    assert get_args(Pathy)[:2] == (str, bytes)
    assert get_args(Pathy)[2].__forward_arg__ == "os.PathLike[str]"
    assert SDBlock.__annotations__ == {
        "title": "str",
        "mdl": "deque[str]",
        "metadata": "deque[str]",
    }
    assert CTable.__annotations__ == {
        "lines": "deque[str]",
        "title": "str",
        "source": "str",
        "comment": "str",
        "counts": "str",
        "format": "CTableFormat",
        "num_atoms": "int",
        "num_bonds": "int",
    }


def test_constructor_and_method_call_signatures() -> None:
    """Keep required arguments, keyword names, and boolean defaults callable."""
    assert tuple(inspect.signature(SDBlock).parameters) == ("title", "mdl", "metadata")
    assert tuple(inspect.signature(CTable).parameters) == ("lines",)
    assert tuple(inspect.signature(parse_sdf).parameters) == ("sdfpath",)
    for factory in (SDBlock.from_lines, SDBlock.from_block_lines):
        assert tuple(inspect.signature(factory).parameters) == ("lines",)
    signatures: tuple[tuple[Any, tuple[str, ...]], ...] = (
        (SDBlock.parse_mdl, ("lines",)),
        (SDBlock.parse_metadata, ("lines",)),
        (SDBlock.records, ("self",)),
        (SDBlock.write, ("self", "outh", "with_newlines")),
        (SDBlock.append_record, ("self", "record_name", "value")),
        (SDBlock.ctable, ("self",)),
        (SDBlock.num_atoms, ("self",)),
        (SDBlock.num_bonds, ("self",)),
        (SDBlock.renumber_indices, ("self",)),
        (CTable.parse_format, ("line",)),
        (CTable.parse_v2000_counts, ("line",)),
        (CTable.parse_v3000_counts, ("line",)),
        (CTable.atomlines, ("self",)),
        (CTable.bondlines, ("self",)),
        (CTable.valid_atom_indices, ("self", "strict")),
        (CTable.valid_bond_indices, ("self", "strict")),
        (CTable.renumber_indices, ("self",)),
    )
    for method, parameters in signatures:
        signature = inspect.signature(method)
        assert tuple(signature.parameters) == parameters
        assert all(
            parameter.kind is inspect.Parameter.POSITIONAL_OR_KEYWORD
            for parameter in signature.parameters.values()
        )
    assert inspect.signature(SDBlock.write).parameters["with_newlines"].default is True
    for validator in (CTable.valid_atom_indices, CTable.valid_bond_indices):
        assert inspect.signature(validator).parameters["strict"].default is False
    assert SDBlock(title="t", mdl=deque(), metadata=deque()).title == "t"
    assert CTable(lines=V2000.splitlines()).num_atoms == 2
    with pytest.raises(TypeError, match="missing 3 required positional arguments"):
        cast(Any, SDBlock)()
    with pytest.raises(TypeError, match="unexpected keyword argument 'other'"):
        cast(Any, SDBlock)("t", deque(), deque(), other=True)


def test_dataclass_fields_repr_equality_and_mutation() -> None:
    """Freeze dataclass value behavior and reference-preserving construction."""
    mdl = deque(["M  END\n"])
    metadata = deque(["> <key>\n", "value\n", "\n"])
    block = SDBlock(" raw title ", mdl, metadata)
    assert is_dataclass(cast(object, SDBlock))
    assert [field.name for field in fields(block)] == ["title", "mdl", "metadata"]
    if hasattr(SDBlock, "__match_args__"):
        assert getattr(SDBlock, "__match_args__") == ("title", "mdl", "metadata")
    assert block.title == " raw title "
    assert block.mdl is mdl
    assert block.metadata is metadata
    assert repr(block) == (
        "SDBlock(title=' raw title ', mdl=deque(['M  END\\n']), "
        "metadata=deque(['> <key>\\n', 'value\\n', '\\n']))"
    )
    equal = SDBlock(" raw title ", deque(mdl), deque(metadata))
    assert block == equal
    with pytest.raises(TypeError, match="unhashable type: 'SDBlock'"):
        hash(block)
    with pytest.raises(TypeError, match="not supported"):
        cast(Any, block) < equal
    mdl.append("opaque")
    assert block != equal
    block.title = "changed"
    block.metadata = deque()
    assert block.title == "changed"
    assert metadata == deque(["> <key>\n", "value\n", "\n"])
    setattr(block, "application_state", 7)
    assert getattr(block, "application_state") == 7


def test_dataclass_helpers_and_no_runtime_type_validation() -> None:
    """Keep replace aliasing, asdict copies, and unchecked annotated fields."""
    block = SDBlock("t", deque(["raw"]), deque(["metadata"]))
    copied = asdict(block)
    assert copied == {
        "title": "t",
        "mdl": deque(["raw"]),
        "metadata": deque(["metadata"]),
    }
    assert copied["mdl"] is not block.mdl
    assert copied["metadata"] is not block.metadata
    changed = replace(block, title="other")
    assert changed.mdl is block.mdl
    assert changed.metadata is block.metadata
    unchecked = SDBlock(cast(str, 7), cast(deque[str], ["raw"]), deque())
    assert cast(Any, unchecked.title) == 7
    assert isinstance(unchecked.mdl, list)


def test_owned_objects_pickle_with_their_python_module_paths() -> None:
    """Base block/table/enum/error objects retain ordinary Python pickle round trips."""
    block = SDBlock.from_block_lines(V2000.splitlines(keepends=True))
    restored_block = pickle.loads(pickle.dumps(block))
    assert type(restored_block) is SDBlock
    assert restored_block == block
    assert restored_block.mdl is not block.mdl
    assert restored_block.metadata is not block.metadata
    table = block.ctable()
    table.num_atoms = 99
    restored_table = pickle.loads(pickle.dumps(table))
    assert type(restored_table) is CTable
    assert restored_table.lines == table.lines
    assert restored_table.num_atoms == 99
    assert CTable(restored_table.lines).num_atoms == 2
    assert pickle.loads(pickle.dumps(CTableFormat.V2000)) is CTableFormat.V2000
    error = pickle.loads(pickle.dumps(IndicesDuplicateError("atoms")))
    assert type(error) is IndicesDuplicateError
    assert error.args == ("atoms",)


def test_subclass_factories_and_exact_class_equality() -> None:
    """Class parse hooks dispatch, but inherited factories return base blocks."""

    class CustomBlock(SDBlock):
        @classmethod
        def parse_mdl(cls, lines: Iterable[str]) -> Iterable[str]:
            yield from super().parse_mdl(lines)
            yield "custom mdl"

        @classmethod
        def parse_metadata(cls, lines: Iterable[str]) -> Iterable[str]:
            yield from super().parse_metadata(lines)
            yield "custom metadata"

    custom = CustomBlock("t", deque(), deque())
    base = SDBlock("t", deque(), deque())
    assert custom != base
    parsed = CustomBlock.from_block_lines(["t", "M  END"])
    assert type(parsed) is SDBlock
    assert parsed.mdl == deque(["M  END", "custom mdl"])
    assert parsed.metadata == deque(["custom metadata"])
    parsed_stream = next(CustomBlock.from_lines(["t", "M  END", "$$$$"]))
    assert type(parsed_stream) is SDBlock
    assert parsed_stream.mdl == parsed.mdl
    assert parsed_stream.metadata == parsed.metadata

    class CustomFactory(SDBlock):
        @classmethod
        def from_block_lines(cls, lines: Iterable[str]) -> SDBlock:
            block = super().from_block_lines(lines)
            return cls(block.title, block.mdl, block.metadata)

    assert isinstance(next(CustomFactory.from_lines(["t", "$$$$"])), CustomFactory)


def test_ctable_subclass_hooks_and_identity_behavior() -> None:
    """Tables support subclass parser overrides, not dataclass equality."""

    class CustomTable(CTable):
        @staticmethod
        def parse_v2000_counts(line: str) -> tuple[int, int]:
            return 8, 9

    table = CustomTable(V2000.splitlines())
    assert (table.num_atoms, table.num_bonds) == (8, 9)
    assert not is_dataclass(table)
    assert table != CustomTable(V2000.splitlines())
    assert table == table
    assert isinstance(hash(table), int)
    setattr(table, "application_state", 7)
    assert getattr(table, "application_state") == 7
    table.source = "changed"
    assert table.source == "changed"
    assert table.lines[1] == "  literal source  "


def test_enum_members_names_values_and_errors() -> None:
    """Formats remain ordinary two-member enums rather than str enums."""
    assert list(CTableFormat) == [CTableFormat.V2000, CTableFormat.V3000]
    assert list(CTableFormat.__members__) == ["V2000", "V3000"]
    for member in CTableFormat:
        assert isinstance(member, Enum)
        assert not isinstance(member, str)
        assert member.name == member.value
        assert CTableFormat[member.name] is member
        assert CTableFormat(member.value) is member
        assert str(member) == f"CTableFormat.{member.name}"
        assert repr(member) == f"<CTableFormat.{member.name}: '{member.value}'>"
        assert member != cast(object, member.value)
    with pytest.raises(KeyError) as keyed:
        CTableFormat["unknown"]
    assert keyed.value.args == ("unknown",)
    assert str(keyed.value) == "'unknown'"
    with pytest.raises(ValueError, match="'unknown' is not a valid CTableFormat"):
        CTableFormat("unknown")


def test_custom_exceptions_are_distinct_plain_exceptions() -> None:
    """Keep catch relationships, arbitrary args, and exact custom message text."""
    for error_type in (
        IndicesDuplicateError,
        IndicesMismatchError,
        IndicesOutOfOrderError,
    ):
        assert error_type.__bases__ == (Exception,)
        error = error_type("atoms")
        assert error.args == ("atoms",)
        assert str(error) == "atoms"
        assert not isinstance(error, ValueError)
        assert error_type().args == ()
        assert error_type("a", "b").args == ("a", "b")
