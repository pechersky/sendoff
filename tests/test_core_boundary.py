"""Build the Rust core independently of the Python extension."""

import subprocess
from pathlib import Path


def test_core_builds_without_python(tmp_path: Path) -> None:
    """Keep the implementation independent of PyO3 and Python dependencies.

    Args:
        tmp_path: directory for the compiled Rust library
    """
    core = Path(__file__).resolve().parents[1] / "rust" / "src" / "core" / "mod.rs"
    subprocess.run(
        [
            "rustc",
            "--edition=2024",
            "--crate-name",
            "sendoff_core",
            "--crate-type",
            "lib",
            str(core),
            "--out-dir",
            str(tmp_path),
        ],
        check=True,
    )
