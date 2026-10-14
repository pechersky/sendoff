# mypy: allow-untyped-decorators
"""Test the python -m sendoff.cli.<command> shell helpers."""
import errno
import importlib
import io
import os
import signal
import subprocess
import sys
from pathlib import Path

import pytest

from sendoff.cli._common import STDIO, Blocks, failures, parser, run

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
    monkeypatch.setattr(
        sys, "stdin", io.TextIOWrapper(io.BytesIO(stdin), encoding="utf-8")
    )
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
    assert err.startswith(
        (
            b"python -m sendoff.cli.count: [Errno 2]",
            b"python -m sendoff.cli.count: [WinError 2]",
        )
    )


# Runs the module as Windows would, where writes to a closed pipe raise.
WITHOUT_SIGPIPE = (
    "import runpy, signal; vars(signal).pop('SIGPIPE', None); "
    "runpy.run_module('sendoff.cli.count', run_name='__main__')"
)


@pytest.mark.parametrize(
    "module",
    [["-m", "sendoff.cli.count"], ["-c", WITHOUT_SIGPIPE]],
    ids=["sigpipe", "no-sigpipe"],
)
def test_module_dies_quietly_on_closed_pipe(sdf: Path, module: list[str]) -> None:
    """The module runs with python -m and stops silently when stdout closes.

    Args:
        sdf: pytest fixture of a Path to the sdf
        module: how python runs the count module
    """
    count = [sys.executable, *module, str(sdf)]
    assert subprocess.run(count, capture_output=True, check=True).stdout == b"3\n"
    read_fd, write_fd = os.pipe()
    os.close(read_fd)
    try:
        closed = subprocess.run(count, stdout=write_fd, stderr=subprocess.PIPE)
    finally:
        os.close(write_fd)
    assert closed.stderr == b""
    assert closed.returncode != 0


def test_run_without_sigpipe(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """A closed pipe passes failures() and exits 1; other OSErrors raise.

    Args:
        monkeypatch: pytest fixture to remove SIGPIPE and replace stdout
        tmp_path: pytest fixture of a temporary directory
    """

    def closed(_: object) -> int:
        with failures(parser("count", ""), Blocks(STDIO)):
            raise BrokenPipeError(errno.EPIPE, "Broken pipe")

    def full(_: object) -> int:
        raise OSError(errno.ENOSPC, "No space left on device")

    monkeypatch.delattr(signal, "SIGPIPE", raising=False)
    # run() points stdout's descriptor at os.devnull.
    with open(tmp_path / "stdout", "w", encoding="utf-8") as stdout:
        monkeypatch.setattr(sys, "stdout", stdout)
        with pytest.raises(SystemExit) as exited:
            run(closed)
        assert exited.value.code == 1
        with pytest.raises(OSError, match="No space"):
            run(full)


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
    assert cli(capsysbinary, monkeypatch, "get", "ID", "-z", sdf) == (
        0,
        b"A-1\0\0A-3\0",
        b"",
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


def test_set(
    capsysbinary: pytest.CaptureFixture[bytes],
    monkeypatch: pytest.MonkeyPatch,
    sdf: Path,
) -> None:
    """Set replaces the first record in place or appends one.

    Args:
        capsysbinary: pytest fixture capturing binary output
        monkeypatch: pytest fixture to replace stdin
        sdf: pytest fixture of a Path to the sdf
    """
    replaced = (
        SDF.replace(b"A-1", b"X")
        .replace(b"A-3", b"X")
        .replace(b"two\n\n$$$$", b"two\n\n> <ID>\nX\n\n$$$$")
    )
    assert cli(capsysbinary, monkeypatch, "set", "ID", "--value", "X", sdf) == (
        0,
        replaced,
        b"",
    )
    appended = SDF.replace(b"\n\n$$$$", b"\n\n>  <ID>\nX\n\n$$$$")
    assert cli(
        capsysbinary,
        monkeypatch,
        "set",
        "ID",
        "--value",
        "X",
        "--append",
        "--header",
        "padded",
        stdin=SDF,
    ) == (0, appended, b"")


@pytest.mark.parametrize("null", [False, True])
def test_set_from(
    capsysbinary: pytest.CaptureFixture[bytes],
    monkeypatch: pytest.MonkeyPatch,
    sdf: Path,
    tmp_path: Path,
    null: bool,
) -> None:
    """Set --from takes a value per block; empty ones are set or skipped.

    Args:
        capsysbinary: pytest fixture capturing binary output
        monkeypatch: pytest fixture to replace stdin
        sdf: pytest fixture of a Path to the sdf
        tmp_path: pytest fixture of a temporary directory
        null: whether values end with NUL
    """
    end = b"\0" if null else b"\n"
    values = tmp_path / "values"
    values.write_bytes(b"x1%s%sx3%s" % (end, end, end))
    z = ["-z"] if null else []
    expected = SDF.replace(b"A-1", b"x1").replace(b"A-3", b"x3")
    assert cli(capsysbinary, monkeypatch, "set", "ID", *z, "--from", values, sdf) == (
        0,
        expected.replace(b"two\n\n$$$$", b"two\n\n> <ID>\n\n\n$$$$"),
        b"",
    )
    assert cli(
        capsysbinary,
        monkeypatch,
        "set",
        "ID",
        *z,
        "--skip-empty",
        "--from",
        "-",
        sdf,
        stdin=values.read_bytes(),
    ) == (0, expected, b"")


@pytest.mark.parametrize(
    "values, error",
    [
        (b"a\nb\n", b"block 3: stdin has only 2 values"),
        (b"a\nb\nc\nd", b"block 3: no more blocks, but stdin has 4 values"),
    ],
)
def test_set_from_mismatch(
    capsysbinary: pytest.CaptureFixture[bytes],
    monkeypatch: pytest.MonkeyPatch,
    sdf: Path,
    values: bytes,
    error: bytes,
) -> None:
    """Set --from fails unless there is one value per block.

    Args:
        capsysbinary: pytest fixture capturing binary output
        monkeypatch: pytest fixture to replace stdin
        sdf: pytest fixture of a Path to the sdf
        values: --from input
        error: expected message after the program name
    """
    status, _, err = cli(
        capsysbinary, monkeypatch, "set", "ID", "--from", "-", sdf, stdin=values
    )
    assert (status, err) == (1, b"python -m sendoff.cli.set: %s\n" % error)


def test_set_from_stdin_needs_file(
    capsysbinary: pytest.CaptureFixture[bytes],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Set cannot read both values and blocks from stdin.

    Args:
        capsysbinary: pytest fixture capturing binary output
        monkeypatch: pytest fixture to replace stdin
    """
    status, out, err = cli(
        capsysbinary, monkeypatch, "set", "ID", "--from", "-", stdin=SDF
    )
    assert (status, out) == (2, b"")
    assert err.endswith(b"error: --from - needs a FILE\n")


def test_title(
    capsysbinary: pytest.CaptureFixture[bytes],
    monkeypatch: pytest.MonkeyPatch,
    sdf: Path,
) -> None:
    """Title sets titles from --value, --field, or get -t output.

    Args:
        capsysbinary: pytest fixture capturing binary output
        monkeypatch: pytest fixture to replace stdin
        sdf: pytest fixture of a Path to the sdf
    """
    first, second, third = (b"first", b"second", b"third")

    def titled(*titles: bytes) -> bytes:
        return (
            SDF.replace(first, titles[0], 1)
            .replace(second, titles[1], 1)
            .replace(third, titles[2], 1)
        )

    assert cli(capsysbinary, monkeypatch, "title", "--value", "T", sdf) == (
        0,
        titled(b"T", b"T", b"T"),
        b"",
    )
    assert cli(capsysbinary, monkeypatch, "title", "--field", "ID", stdin=SDF) == (
        0,
        titled(b"A-1", b"", b"A-3"),
        b"",
    )
    assert cli(
        capsysbinary, monkeypatch, "title", "--field", "ID", "--skip-empty", sdf
    ) == (0, titled(b"A-1", second, b"A-3"), b"")
    _, titles, _ = cli(capsysbinary, monkeypatch, "get", "-t", sdf)
    assert cli(
        capsysbinary, monkeypatch, "title", "--from", "-", sdf, stdin=titles
    ) == (0, SDF, b"")


def test_title_line_break(
    capsysbinary: pytest.CaptureFixture[bytes],
    monkeypatch: pytest.MonkeyPatch,
    sdf: Path,
) -> None:
    """A title with a line break fails with the block number.

    Args:
        capsysbinary: pytest fixture capturing binary output
        monkeypatch: pytest fixture to replace stdin
        sdf: pytest fixture of a Path to the sdf
    """
    status, _, err = cli(capsysbinary, monkeypatch, "title", "--field", "Note", sdf)
    assert status == 1
    assert err.startswith(b"python -m sendoff.cli.title: block 2: ")


def test_filter(
    capsysbinary: pytest.CaptureFixture[bytes],
    monkeypatch: pytest.MonkeyPatch,
    sdf: Path,
) -> None:
    """Filter writes the matching blocks byte for byte.

    Args:
        capsysbinary: pytest fixture capturing binary output
        monkeypatch: pytest fixture to replace stdin
        sdf: pytest fixture of a Path to the sdf
    """
    first, second, third = (block + b"$$$$\n" for block in SDF.split(b"$$$$\n")[:3])
    assert cli(capsysbinary, monkeypatch, "filter", "--has", "ID", sdf) == (
        0,
        first + third,
        b"",
    )
    assert cli(capsysbinary, monkeypatch, "filter", "-v", "--has", "ID", sdf) == (
        0,
        second,
        b"",
    )
    assert cli(capsysbinary, monkeypatch, "filter", "--eq", "ID", "A-3", stdin=SDF) == (
        0,
        third,
        b"",
    )
    if sys.getfilesystemencodeerrors() == "surrogateescape":
        name = os.fsdecode(GBK)
        assert cli(capsysbinary, monkeypatch, "filter", "--eq", "Name", name, sdf) == (
            0,
            first,
            b"",
        )
    crlf = SDF.replace(b"\n", b"\r\n")
    _, out, _ = cli(
        capsysbinary, monkeypatch, "filter", "-v", "--eq", "ID", "A-1", stdin=crlf
    )
    assert out == (second + third).replace(b"\n", b"\r\n")


def test_filter_unique(
    capsysbinary: pytest.CaptureFixture[bytes],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Filter --unique keeps the first block per value and those without one.

    Args:
        capsysbinary: pytest fixture capturing binary output
        monkeypatch: pytest fixture to replace stdin
    """
    first, second, third = (block + b"$$$$\n" for block in SDF.split(b"$$$$\n")[:3])
    repeated = first + second + first + third + second + third
    assert cli(
        capsysbinary, monkeypatch, "filter", "--unique", "ID", stdin=repeated
    ) == (
        0,
        first + second + third + second,
        b"",
    )
    assert cli(
        capsysbinary, monkeypatch, "filter", "-v", "--unique", "ID", stdin=repeated
    ) == (0, first + third, b"")


COMMANDS = [
    ("count",),
    ("get", "ID"),
    ("get", "-t"),
    ("set", "ID", "--value", "X"),
    ("title", "--field", "ID"),
    ("filter", "--has", "ID"),
    ("filter", "--unique", "ID"),
]


@pytest.mark.parametrize("command", COMMANDS)
def test_input_and_output_options(
    capsysbinary: pytest.CaptureFixture[bytes],
    monkeypatch: pytest.MonkeyPatch,
    sdf: Path,
    command: tuple[str, ...],
) -> None:
    """-i reads like FILE, -o writes what stdout would, and -o FILE is no input.

    Args:
        capsysbinary: pytest fixture capturing binary output
        monkeypatch: pytest fixture to replace stdin
        sdf: pytest fixture of a Path to the sdf
        command: module name and arguments
    """
    expected = cli(capsysbinary, monkeypatch, *command, sdf)
    assert expected[0] == 0 and expected[1]
    assert cli(capsysbinary, monkeypatch, *command, "-i", sdf) == expected
    out = sdf.with_name("out")
    assert cli(capsysbinary, monkeypatch, *command, "-i", sdf, "-o", out) == (
        0,
        b"",
        b"",
    )
    assert out.read_bytes() == expected[1]
    status, _, err = cli(capsysbinary, monkeypatch, *command, sdf, "-i", sdf)
    if command[0] != "count":
        assert status == 2
        assert b"give FILE or -i, not both" in err
    status, _, err = cli(capsysbinary, monkeypatch, *command, "-i", sdf, "-o", sdf)
    assert status == 2
    assert b"is also an input" in err
    assert sdf.read_bytes() == SDF
    missing = sdf.with_name("no-dir") / "out"
    status, _, err = cli(capsysbinary, monkeypatch, *command, sdf, "-o", missing)
    assert status == 1
    assert b"[Errno 2]" in err or b"[WinError 3]" in err


def test_count_repeats_input(
    capsysbinary: pytest.CaptureFixture[bytes],
    monkeypatch: pytest.MonkeyPatch,
    sdf: Path,
) -> None:
    """Count reads every FILE and -i.

    Args:
        capsysbinary: pytest fixture capturing binary output
        monkeypatch: pytest fixture to replace stdin
        sdf: pytest fixture of a Path to the sdf
    """
    _, out, _ = cli(capsysbinary, monkeypatch, "count", sdf, "-i", sdf, "-i", sdf)
    assert out == b"3\t%s\n" % bytes(sdf) * 3 + b"9\ttotal\n"
