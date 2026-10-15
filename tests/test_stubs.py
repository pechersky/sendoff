"""Check the hand-written native stub against the compiled extension."""

import subprocess
import sys
from importlib import resources
from pathlib import Path

import pytest

# Generic[AnyStr] supplies __class_getitem__ in the stub. PyO3 lists every
# export, private iterators included, in the runtime __all__; sendoff's own
# __all__ is the public one.
ALLOWLIST = """\
sendoff.native.SDBlock.__class_getitem__
sendoff.native.__all__
"""
# Python 3.11 has the buffer slots but not their PEP 688 methods.
if sys.version_info < (3, 12):
    ALLOWLIST += """\
sendoff.native.SDBlock.__buffer__
sendoff.native.SDBlock.__release_buffer__
"""


def test_native_stub_matches_extension(tmp_path: Path) -> None:
    """Stubtest finds no drift between native.pyi and the built module.

    Args:
        tmp_path: directory for the stubtest allowlist
    """
    pytest.importorskip("mypy.stubtest")
    allowlist = tmp_path / "allowlist.txt"
    allowlist.write_text(ALLOWLIST, encoding="utf-8")
    completed = subprocess.run(
        [
            sys.executable,
            "-m",
            "mypy.stubtest",
            "sendoff.native",
            "--allowlist",
            str(allowlist),
        ],
        capture_output=True,
        check=False,
        cwd=Path(__file__).parents[1],
        encoding="utf-8",
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr


def test_package_is_marked_typed() -> None:
    """The PEP 561 marker ships, so downstream type checkers read the stub."""
    assert (resources.files("sendoff") / "py.typed").is_file()


def test_package_exports_the_api() -> None:
    """The package exports what its submodules define, as the same objects."""
    import sendoff
    from sendoff import ctable, sdblock

    for name in sendoff.__all__:
        module = sdblock if hasattr(sdblock, name) else ctable
        assert getattr(sendoff, name) is getattr(module, name)
