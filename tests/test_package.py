"""Check the mandatory mixed-package bridge without implementing leaf behavior."""

from __future__ import annotations

import shutil
import subprocess
import sys
from importlib.machinery import ExtensionFileLoader
from importlib.metadata import version
from pathlib import Path

import pytest

import sendoff
import sendoff._native as _native


def test_mandatory_binary_and_distribution_version() -> None:
    """Load an actual extension and coordinate Python, Cargo and package versions."""
    assert isinstance(_native.__loader__, ExtensionFileLoader)
    assert _native.__name__ == "sendoff._native"
    assert _native.__file__ is not None
    assert Path(_native.__file__).name == "_native.abi3.so"
    assert sendoff.__version__ == _native.__version__ == version("sendoff")
    assert (Path(sendoff.__file__).parent / "_native.pyi").is_file()


@pytest.mark.parametrize(
    "module_name", ("sendoff", "sendoff.sdblock", "sendoff.ctable")
)
def test_missing_native_module_fails_import(module_name: str, tmp_path: Path) -> None:
    """Fail public imports normally when the mandatory extension is unavailable.

    Args:
        module_name: the public import path to attempt in a fresh interpreter
        tmp_path: the isolated package directory without the native binary
    """
    shutil.copytree(
        Path(sendoff.__file__).parent,
        tmp_path / "sendoff",
        ignore=shutil.ignore_patterns("*.so", "__pycache__"),
    )
    script = f"import sys; sys.path.insert(0, {str(tmp_path)!r}); import {module_name}"
    result = subprocess.run(
        [sys.executable, "-I", "-c", script],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode != 0
    assert "ModuleNotFoundError: No module named 'sendoff._native'" in result.stderr
