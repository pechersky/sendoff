"""Check the hand-written native stub against the compiled extension."""

import subprocess
import sys
from pathlib import Path

import pytest

# Generic[AnyStr] supplies __class_getitem__ in the stub, and iter_records is
# typed as Iterator[tuple] rather than by its private iterator class.
ALLOWLIST = """\
sendoff.native.SDBlock.__class_getitem__
sendoff.native._SDRecordIterator
"""


def test_native_stub_matches_extension(tmp_path: Path) -> None:
    """Stubtest finds no drift between native.pyi and the built module.

    Args:
        tmp_path: directory for the stubtest allowlist
    """
    pytest.importorskip("mypy.stubtest")
    allowlist = tmp_path / "allowlist.txt"
    allowlist.write_text(ALLOWLIST)
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
        text=True,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr
