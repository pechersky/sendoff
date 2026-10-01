"""Inspect and smoke-test an installed sendoff wheel without source imports."""

from __future__ import annotations

import os
import subprocess
import sys
import tempfile
from pathlib import Path
from zipfile import ZipFile


def inspect_wheel() -> None:
    """Check the built wheel's abi3 tags and extension metadata.

    Raises:
        AssertionError: If wheel metadata is incorrect.
    """
    wheels = sorted(Path("wheelhouse").glob("sendoff-*.whl"))
    assert len(wheels) == 1, wheels
    wheel = wheels[0]
    _, interpreter, abi, platforms = wheel.stem.rsplit("-", 3)
    expected_platforms = set(os.environ["EXPECTED_PLATFORM_TAGS"].split("."))
    assert (interpreter, abi) == ("cp311", "abi3"), wheel.name
    assert set(platforms.split(".")) == expected_platforms, wheel.name

    with ZipFile(wheel) as archive:
        names = archive.namelist()
        extensions = [name for name in names if name.endswith(".so")]
        metadata_files = [name for name in names if name.endswith(".dist-info/WHEEL")]
        assert extensions == ["sendoff/native.abi3.so"], extensions
        assert len(metadata_files) == 1, metadata_files
        metadata = archive.read(metadata_files[0]).decode()
        assert "Root-Is-Purelib: false" in metadata.splitlines()
        tags = {
            line.removeprefix("Tag: ")
            for line in metadata.splitlines()
            if line.startswith("Tag: ")
        }
        assert tags == {f"cp311-abi3-{tag}" for tag in expected_platforms}, tags

    # ponytail: importing on the matching runner proves binary architecture.
    print("PASS", wheel.name, sorted(expected_platforms))


def smoke_installed_wheel() -> None:
    """Exercise installed imports and small public SDBlock/CTable operations.

    Raises:
        AssertionError: If imports, class identities, or operations differ.
    """
    import platform
    import sysconfig
    from importlib.machinery import EXTENSION_SUFFIXES, ExtensionFileLoader
    from importlib.metadata import version

    import sendoff
    import sendoff.ctable as ctable_module
    import sendoff.native as native
    import sendoff.sdblock as sdblock_module
    from sendoff.ctable import CTable
    from sendoff.sdblock import SDBlock

    assert sys.implementation.name == "cpython"
    assert sys.version_info >= (3, 11)
    expected_minor = os.environ["EXPECTED_PYTHON_MINOR"]
    if expected_minor != "latest":
        assert sys.version_info[:2] == (3, int(expected_minor))
    else:
        assert sys.version_info[:2] != (3, 11)
    assert not sysconfig.get_config_var("Py_GIL_DISABLED")
    assert platform.machine() == os.environ["EXPECTED_MACHINE"]
    assert isinstance(native.__loader__, ExtensionFileLoader)
    assert native.__file__ is not None
    assert any(native.__file__.endswith(suffix) for suffix in EXTENSION_SUFFIXES)
    assert sendoff.__file__ is not None
    site_packages = Path(sysconfig.get_paths()["purelib"]).resolve()
    assert Path(sendoff.__file__).resolve().is_relative_to(site_packages)
    assert Path(native.__file__).resolve().is_relative_to(site_packages)
    assert sendoff.__version__ == native.__version__ == version("sendoff")
    assert ctable_module.CTable is getattr(sdblock_module, "CTable") is CTable
    assert sdblock_module.SDBlock is SDBlock

    block = SDBlock.from_block_lines(
        [
            "wheel smoke",
            "  sendoff",
            "literal SDF",
            "  1  0  0  0  0  0            999 V2000",
            "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0",
            "M  END",
        ]
    )
    table = block.ctable()
    assert type(block) is SDBlock
    assert type(table) is CTable
    assert (block.num_atoms(), table.num_atoms, table.num_bonds) == (1, 1, 0)
    block.append_record("CI", "installed")
    assert list(block.records()) == [("CI", "installed")]
    try:
        CTable.parse_format("UNKNOWN")
    except KeyError as error:
        assert error.args == ("UNKNOWN",)
    else:
        raise AssertionError("invalid format error was hidden")
    print("PASS", sys.version.split()[0], platform.machine(), native.__file__)


def parity_installed_wheel() -> None:
    """Run the RDKit-independent Python corpus against the installed wheel."""
    smoke_installed_wheel()
    tests = Path(__file__).resolve().parent
    files = sorted(tests.glob("test_native_*.py"))
    files += sorted(tests.glob("test_compat_*.py"))
    files.append(tests / "test_package.py")
    with tempfile.TemporaryDirectory(prefix="sendoff-parity-") as directory:
        subprocess.run(
            [
                sys.executable,
                "-I",
                "-m",
                "pytest",
                "--import-mode=importlib",
                "--noconftest",
                "-o",
                "addopts=",
                "-q",
                *(str(path) for path in files),
            ],
            cwd=directory,
            check=True,
        )


def main() -> None:
    """Run the requested wheel check mode.

    Raises:
        SystemExit: If the command line does not select a known check mode.
    """
    if len(sys.argv) != 2:
        raise SystemExit("usage: ci_wheel_smoke.py inspect|smoke|parity")
    if sys.argv[1] == "inspect":
        inspect_wheel()
    elif sys.argv[1] == "smoke":
        smoke_installed_wheel()
    elif sys.argv[1] == "parity":
        parity_installed_wheel()
    else:
        raise SystemExit(f"unknown mode: {sys.argv[1]}")


if __name__ == "__main__":
    main()
