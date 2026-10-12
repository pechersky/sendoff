# mypy: allow-untyped-decorators
"""Test the python -m sendoff.cli.<command> shell helpers."""
import importlib
import io
import subprocess
import sys
from pathlib import Path

import pytest

MOL = b"\n  sendoff\n\n  0  0  0  0  0  0  0  0  0  0999 V2000\nM  END\n"
GBK = "\N{GREEK SMALL LETTER BETA}-Bisabolene".encode("gbk")
SDF = (
    b"first" + MOL + b"> <ID>\nA-1\n\n> <Name>\n" + GBK + b"\n\n$$$$\n"
    b"second" + MOL + b"> <Note>\nline one\nline two\n\n$$$$\n"
    b"third" + MOL + b">  <ID>  (3)\nA-3\n\n$$$$\n"
)


@pytest.fixture
def sdf(tmp_path: Path) -> Path:
    """Write a ragged SD file with non-UTF-8 bytes.

    Args:
        tmp_path: pytest fixture of a temporary directory

    Returns:
        The file's path
    """
    path = tmp_path / "in.sdf"
    path.write_bytes(SDF)
    return path


def cli(
    capsysbinary: pytest.CaptureFixture[bytes],
    monkeypatch: pytest.MonkeyPatch,
    command: str,
    *argv: str | Path,
    stdin: bytes = b"",
) -> tuple[object, bytes, bytes]:
    """Run a command in-process.

    Args:
        capsysbinary: pytest fixture capturing binary output
        monkeypatch: pytest fixture to replace stdin
        command: module name under sendoff.cli
        argv: arguments
        stdin: bytes to read from stdin

    Returns:
        Exit status, stdout, and stderr
    """
    monkeypatch.setattr(sys, "stdin", io.TextIOWrapper(io.BytesIO(stdin)))
    main = importlib.import_module(f"sendoff.cli.{command}").main
    try:
        status = main([str(arg) for arg in argv])
    except SystemExit as exc:
        status = exc.code
    out, err = capsysbinary.readouterr()
    return status, out, err


def test_count(
    capsysbinary: pytest.CaptureFixture[bytes],
    monkeypatch: pytest.MonkeyPatch,
    sdf: Path,
) -> None:
    """Count prints one number, or COUNT<TAB>FILE lines and a total.

    Args:
        capsysbinary: pytest fixture capturing binary output
        monkeypatch: pytest fixture to replace stdin
        sdf: pytest fixture of a Path to the sdf
    """
    assert cli(capsysbinary, monkeypatch, "count", sdf) == (0, b"3\n", b"")
    assert cli(capsysbinary, monkeypatch, "count", stdin=SDF) == (0, b"3\n", b"")
    assert cli(
        capsysbinary, monkeypatch, "count", sdf, "-", stdin=SDF[: SDF.index(b"third")]
    ) == (
        0,
        b"3\t%s\n2\t-\n5\ttotal\n" % bytes(sdf),
        b"",
    )


def test_missing_file(
    capsysbinary: pytest.CaptureFixture[bytes],
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """A missing file fails with status 1 and no traceback.

    Args:
        capsysbinary: pytest fixture capturing binary output
        monkeypatch: pytest fixture to replace stdin
        tmp_path: pytest fixture of a temporary directory
    """
    status, out, err = cli(capsysbinary, monkeypatch, "count", tmp_path / "nope")
    assert (status, out) == (1, b"")
    assert err.startswith(b"python -m sendoff.cli.count: [Errno 2]")


def test_module_dies_quietly_on_closed_pipe(sdf: Path) -> None:
    """The module runs with python -m and stops silently when stdout closes.

    Args:
        sdf: pytest fixture of a Path to the sdf
    """
    count = [sys.executable, "-m", "sendoff.cli.count", str(sdf)]
    assert subprocess.run(count, capture_output=True, check=True).stdout == b"3\n"
    many = [sys.executable, "-m", "sendoff.cli.count", *[str(sdf)] * 5000]
    with subprocess.Popen(many, stdout=subprocess.PIPE, stderr=subprocess.PIPE) as p:
        assert p.stdout is not None and p.stderr is not None
        p.stdout.read(1)
        p.stdout.close()
        assert p.stderr.read() == b""
    assert p.returncode != 0


def test_get(
    capsysbinary: pytest.CaptureFixture[bytes],
    monkeypatch: pytest.MonkeyPatch,
    sdf: Path,
) -> None:
    """Get prints one raw line per block, blank or --missing where absent.

    Args:
        capsysbinary: pytest fixture capturing binary output
        monkeypatch: pytest fixture to replace stdin
        sdf: pytest fixture of a Path to the sdf
    """
    assert cli(capsysbinary, monkeypatch, "get", "ID", sdf) == (0, b"A-1\n\nA-3\n", b"")
    assert cli(
        capsysbinary, monkeypatch, "get", "Name", "--missing", "NA", stdin=SDF
    ) == (
        0,
        b"%s\nNA\nNA\n" % GBK,
        b"",
    )
    assert cli(capsysbinary, monkeypatch, "get", "-z", "Note", sdf) == (
        0,
        b"\0line one\nline two\0\0",
        b"",
    )
    assert cli(capsysbinary, monkeypatch, "get", "-t", sdf) == (
        0,
        b"first\nsecond\nthird\n",
        b"",
    )
    assert cli(capsysbinary, monkeypatch, "get", "-t", stdin=SDF)[1] == (
        b"first\nsecond\nthird\n"
    )


@pytest.mark.parametrize("argv", [(), ("-t", "ID", "in.sdf")])
def test_get_usage(
    capsysbinary: pytest.CaptureFixture[bytes],
    monkeypatch: pytest.MonkeyPatch,
    argv: tuple[str, ...],
) -> None:
    """Get takes NAME or -t, not both.

    Args:
        capsysbinary: pytest fixture capturing binary output
        monkeypatch: pytest fixture to replace stdin
        argv: bad arguments
    """
    status, out, err = cli(capsysbinary, monkeypatch, "get", *argv)
    assert (status, out) == (2, b"")
    assert b"usage: python -m sendoff.cli.get" in err


def test_malformed_block(
    capsysbinary: pytest.CaptureFixture[bytes],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A malformed block fails with its number after the earlier output.

    Args:
        capsysbinary: pytest fixture capturing binary output
        monkeypatch: pytest fixture to replace stdin
    """
    bad = SDF.replace(b"> <Note>", b"> Note")
    assert cli(capsysbinary, monkeypatch, "get", "ID", stdin=bad) == (
        1,
        b"A-1\n",
        b"python -m sendoff.cli.get: block 2: "
        b"SD data field header is missing a name\n",
    )
