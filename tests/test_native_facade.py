"""Require public operations to enter the native implementation."""

from __future__ import annotations

from io import BytesIO, TextIOWrapper
from typing import Callable

import pytest

import sendoff.native as native
from sendoff.ctable import CTable
from sendoff.sdblock import SDBlock, parse_sdf
from tests.compat_literals import V3000


@pytest.mark.parametrize(
    "name,operation",
    [
        ("mdl_iter", lambda table, block: list(SDBlock.parse_mdl([]))),
        ("metadata_iter", lambda table, block: list(SDBlock.parse_metadata([]))),
        ("from_block_lines", lambda table, block: SDBlock.from_block_lines([])),
        ("blocks_iter", lambda table, block: list(SDBlock.from_lines([]))),
        ("read_sdf_lines", lambda table, block: parse_sdf("unused")),
        ("records_iter", lambda table, block: list(block.records())),
        ("write", lambda table, block: block.write(TextIOWrapper(BytesIO()))),
        ("append_record", lambda table, block: block.append_record("key", "value")),
        ("ctable_init", lambda table, block: CTable([])),
        ("parse_format", lambda table, block: CTable.parse_format("")),
        ("parse_v2000_counts", lambda table, block: CTable.parse_v2000_counts("")),
        ("parse_v3000_counts", lambda table, block: CTable.parse_v3000_counts("")),
        ("atomlines", lambda table, block: list(table.atomlines())),
        ("bondlines", lambda table, block: list(table.bondlines())),
        ("valid_atom_indices", lambda table, block: table.valid_atom_indices()),
        ("valid_bond_indices", lambda table, block: table.valid_bond_indices()),
        ("renumber_ctable", lambda table, block: table.renumber_indices()),
    ],
)
def test_public_operations_require_native(
    monkeypatch: pytest.MonkeyPatch,
    name: str,
    operation: Callable[[CTable, SDBlock], object],
) -> None:
    """Propagate native failures rather than running a Python fallback.

    Args:
        monkeypatch: native function replacement
        name: native entry point
        operation: public operation invoking the entry point
    """
    table = CTable(V3000.splitlines())
    block = SDBlock.from_block_lines(V3000.splitlines())

    def fail(*args: object, **kwargs: object) -> None:
        raise RuntimeError("native entry point")

    monkeypatch.setattr(native, name, fail)
    with pytest.raises(RuntimeError, match="^native entry point$"):
        operation(table, block)
