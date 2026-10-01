"""Check the mandatory mixed-package bridge without implementing leaf behavior."""

from __future__ import annotations

import inspect
import subprocess
import sys
from importlib.machinery import ExtensionFileLoader
from importlib.metadata import version
from pathlib import Path

import pytest

import sendoff
import sendoff._native as _native

NATIVE_SIGNATURES = (
    ("_mdl_iter", ("lines",)),
    ("_metadata_iter", ("lines",)),
    ("_from_block_lines", ("cls", "block_type", "lines")),
    ("_blocks_iter", ("cls", "lines")),
    ("_read_sdf_lines", ("sdfpath",)),
    ("_records_iter", ("block",)),
    ("_write", ("block", "outh", "with_newlines")),
    ("_append_record", ("block", "record_name", "value")),
    ("_ctable_init", ("table", "lines", "v3000")),
    ("_parse_format", ("line", "formats")),
    ("_parse_v2000_counts", ("line",)),
    ("_parse_v3000_counts", ("line",)),
    ("_atomlines", ("table",)),
    ("_bondlines", ("table",)),
    ("_valid_atom_indices", ("table", "strict", "v3000", "errors")),
    ("_valid_bond_indices", ("table", "strict", "v3000", "errors")),
    ("_renumber_ctable", ("table", "v3000", "duplicate_error")),
)


def test_mandatory_binary_and_distribution_version() -> None:
    """Load an actual extension and coordinate Python, Cargo and package versions."""
    assert isinstance(_native.__loader__, ExtensionFileLoader)
    assert _native.__name__ == "sendoff._native"
    assert _native.__file__ is not None
    assert Path(_native.__file__).name == "_native.abi3.so"
    assert sendoff.__version__ == _native.__version__ == version("sendoff")
    assert (Path(sendoff.__file__).parent / "_native.pyi").is_file()
    assert (Path(sendoff.__file__).parent / "py.typed").is_file()


@pytest.mark.parametrize("name,parameters", NATIVE_SIGNATURES)
def test_private_bridge_call_signatures(name: str, parameters: tuple[str, ...]) -> None:
    """Keep leaf exports callable with the approved uncoerced argument names.

    Args:
        name: the private extension function name
        parameters: its approved required positional-or-keyword arguments
    """
    function = getattr(_native, name)
    assert inspect.isbuiltin(function)
    signature = inspect.signature(function)
    assert tuple(signature.parameters) == parameters
    assert all(
        argument.kind is inspect.Parameter.POSITIONAL_OR_KEYWORD
        and argument.default is inspect.Parameter.empty
        for argument in signature.parameters.values()
    )


@pytest.mark.parametrize(
    "module_name", ("sendoff", "sendoff.sdblock", "sendoff.ctable")
)
def test_missing_native_module_fails_import(module_name: str) -> None:
    """Fail public imports normally when the mandatory extension is unavailable.

    Args:
        module_name: the public import path to attempt in a fresh interpreter
    """
    package_parent = str(Path(sendoff.__file__).resolve().parent.parent)
    script = f"""
import importlib
import importlib.abc
import sys

class MissingNative(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname == "sendoff._native":
            raise ModuleNotFoundError(
                "blocked mandatory sendoff._native", name=fullname
            )

sys.path.insert(0, {package_parent!r})
sys.meta_path.insert(0, MissingNative())
importlib.import_module({module_name!r})
print("unexpected import success")
"""
    result = subprocess.run(
        [sys.executable, "-I", "-c", script],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode != 0
    assert "ModuleNotFoundError: blocked mandatory sendoff._native" in result.stderr
    assert "unexpected import success" not in result.stdout
