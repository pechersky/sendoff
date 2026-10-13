# mypy: allow-untyped-decorators
"""Test writing blocks with write_sdf and sdf_writer."""
import errno
import io
import os
import typing
from pathlib import Path
from typing import Iterator

import pytest

from sendoff.sdblock import SDBlock, parse_sdf, sdf_writer, write_sdf
from tests.fixtures import Pathy, sdf_of


def test_single_mol_newline_write(single_mol_sdf: Pathy) -> None:
    """An sdf block gets written with as many lines as the input.

    Args:
        single_mol_sdf: pytest fixture of a Path to the sdf
    """
    mols = list(parse_sdf(single_mol_sdf))
    buffer = io.StringIO()
    write_sdf(buffer, mols)
    assert len(buffer.getvalue().splitlines()) == len(
        open(single_mol_sdf).read().splitlines()
    )


def test_single_mol_newline_write_splitlines(single_mol_sdf: Pathy) -> None:
    """An sdf block gets written with as many lines as a splitlines input.

    Args:
        single_mol_sdf: pytest fixture of a Path to the sdf
    """
    mols = list(SDBlock.from_lines(open(single_mol_sdf).read().splitlines()))
    buffer = io.StringIO()
    write_sdf(buffer, mols)
    assert len(buffer.getvalue().splitlines()) == len(
        open(single_mol_sdf).read().splitlines()
    )


def test_write_sdf_writes_each_blocks_sdf(tmp_path: Path) -> None:
    """write_sdf to a path or handle writes each block's SDF, across chunks.

    Args:
        tmp_path: pytest fixture of a temporary directory
    """
    # Padded titles are stripped; 3 MiB of blocks span several chunks, and
    # non-ASCII text blocks go to a text handle on their own.
    text = str.join(
        "",
        [
            f"  t{i}  \r\nM  END\n> <A>\n{'é' * (i % 7 < 2)}{'v' * 1000}\n\n$$$$\n"
            for i in range(3000)
        ],
    )
    source = tmp_path / "input.sdf"
    source.write_text(text, newline="")
    out = tmp_path / "out.sdf"
    paths: list[Pathy] = [source, os.fsencode(source)]
    for path in paths:
        blocks = list(parse_sdf(path))
        handle = io.StringIO() if isinstance(path, Path) else io.BytesIO()
        expected = bytes.join(
            b"",
            [str(b).encode() if isinstance(b.title, str) else bytes(b) for b in blocks],
        )
        assert write_sdf(out, blocks) == 3000
        assert out.read_bytes() == expected
        assert write_sdf(handle, iter(blocks)) == 3000
        written = handle.getvalue()
        assert (written.encode() if isinstance(written, str) else written) == expected
    assert write_sdf(str(out), []) == 0
    assert out.read_bytes() == b""


def _then_raise(blocks: list[SDBlock[str]]) -> Iterator[SDBlock[str]]:
    """Produce `blocks`, then fail.

    Args:
        blocks: the blocks to produce first

    Yields:
        each of `blocks`

    Raises:
        RuntimeError: after the last block
    """
    yield from blocks
    raise RuntimeError("source failed")


@pytest.mark.parametrize("to_path", [True, False], ids=["path", "handle"])
@pytest.mark.parametrize(
    ("failing", "error", "match"),
    [
        (None, RuntimeError, "source failed"),
        ("t\nM  END\n", TypeError, "SDBlock"),
        (b"latin1", UnicodeDecodeError, "utf-8"),
        (b"t\nM  END\n", TypeError, None),
    ],
    ids=["source-error", "not-a-block", "undecodable", "bytes-block"],
)
def test_write_sdf_writes_blocks_before_a_failure(
    tmp_path: Path,
    to_path: bool,
    failing: typing.Any,
    error: type[Exception],
    match: str | None,
) -> None:
    """Blocks before a failing one are written, and the error propagates.

    Args:
        tmp_path: pytest fixture of a temporary directory
        to_path: whether to write to a path or to a file object
        failing: what follows the good blocks
        error: the exception it raises
        match: a pattern in its message
    """
    good = [SDBlock(f"{title}\nM  END\n") for title in "ab"]
    if failing == b"t\nM  END\n" and to_path:
        # Paths take str and bytes blocks alike.
        assert write_sdf(tmp_path / "mixed.sdf", [*good, SDBlock(failing)]) == 3
        return
    blocks: typing.Any = _then_raise(good)
    if failing == b"latin1":
        source = tmp_path / "latin1.sdf"
        source.write_bytes(b"t\nM  END\n> <A>\n\xe9\n\n$$$$\n")
        blocks = [*good, *parse_sdf(source)]
    elif failing is not None:
        blocks = [*good, SDBlock(failing) if isinstance(failing, bytes) else failing]
    out = tmp_path / "out.sdf"
    handle = io.StringIO()
    with pytest.raises(error, match=match) as raised:
        write_sdf(out if to_path else handle, blocks)
    if error is UnicodeDecodeError:
        assert raised.value.__notes__[0].startswith("Writing block 2 (from 0),")
    written = out.read_text() if to_path else handle.getvalue()
    assert written == "a\nM  END\n$$$$\nb\nM  END\n$$$$\n"


def test_writes_check_utf8_until_it_is_known(tmp_path: Path) -> None:
    """Path text blocks are checked on every write until one decodes.

    Args:
        tmp_path: pytest fixture of a temporary directory
    """
    source = tmp_path / "in.sdf"
    source.write_bytes(b"t\nM  END\n> <A>\n\xe9\n\n$$$$\nu\nM  END\n$$$$\n")
    bad, good = parse_sdf(source)
    out = tmp_path / "out.sdf"
    for _ in range(2):
        with pytest.raises(UnicodeDecodeError):
            write_sdf(out, [bad, bad.__copy__()])
    good["A"] = "\u00e9"
    assert sdf_of(good) == "u\nM  END\n> <A>\n\u00e9\n\n$$$$\n"
    assert write_sdf(out, [good, good]) == 2
    assert out.read_text() == sdf_of(good) * 2


def test_write_sdf_raises_os_errors(tmp_path: Path) -> None:
    """Paths raise the OSError subclass open() would, naming the path.

    Args:
        tmp_path: pytest fixture of a temporary directory
    """
    blocks = [SDBlock("t\nM  END\n")]
    locked = tmp_path / "locked"
    locked.mkdir(mode=0o500)
    dests: dict[type[OSError], Path] = {
        IsADirectoryError: tmp_path,
        FileNotFoundError: tmp_path / "missing" / "x.sdf",
    }
    if os.geteuid() != 0:
        dests[PermissionError] = locked / "x.sdf"
    for error, dest in dests.items():
        with pytest.raises(error) as raised:
            write_sdf(dest, blocks)
        assert raised.value.filename == str(dest)
    if os.path.exists("/dev/full"):
        with pytest.raises(OSError) as raised:
            write_sdf(b"/dev/full", blocks)
        assert (raised.value.errno, raised.value.filename) == (
            errno.ENOSPC,
            b"/dev/full",
        )


def test_write_sdf_checks_blocks_before_creating_the_path(tmp_path: Path) -> None:
    """A non-iterable `blocks` leaves the destination alone.

    Args:
        tmp_path: pytest fixture of a temporary directory
    """
    out = tmp_path / "out.sdf"
    out.write_text("kept")
    with pytest.raises(TypeError, match="not iterable"):
        write_sdf(out, 5)  # type: ignore[call-overload]
    assert out.read_text() == "kept"


def test_write_sdf_propagates_file_object_errors() -> None:
    """A file object's write error propagates."""
    handle = io.StringIO()
    handle.close()
    with pytest.raises(ValueError, match="closed file"):
        write_sdf(handle, [SDBlock("t\nM  END\n")])


def test_sdf_writer_writes_as_write_sdf(tmp_path: Path) -> None:
    """One pass writes to several paths and handles as write_sdf would.

    Args:
        tmp_path: pytest fixture of a temporary directory
    """
    text = str.join(
        "",
        [
            f"t{i}\nM  END\n> <A>\n{'é' * (i % 7 < 2)}{'v' * 1000}\n\n$$$$\n"
            for i in range(3000)
        ],
    )
    source = tmp_path / "input.sdf"
    source.write_text(text)
    paths: list[Pathy] = [source, os.fsencode(source)]
    for path in paths:
        blocks = list(parse_sdf(path))
        handle: typing.Any = io.StringIO() if isinstance(path, Path) else io.BytesIO()
        expected = tmp_path / "expected.sdf"
        write_sdf(expected, blocks)
        with (
            sdf_writer(tmp_path / "all.sdf") as everything,
            sdf_writer(handle) as copy,
            sdf_writer(str(tmp_path / "odd.sdf")) as odd,
        ):
            for i, block in enumerate(blocks):
                everything.write(block)
                copy.write(block)
                if i % 2:
                    odd.write(block)
        assert (tmp_path / "all.sdf").read_bytes() == expected.read_bytes()
        written = handle.getvalue()
        assert (written.encode() if isinstance(written, str) else written) == (
            expected.read_bytes()
        )
        write_sdf(expected, blocks[1::2])
        assert (tmp_path / "odd.sdf").read_bytes() == expected.read_bytes()
        assert not handle.closed


def test_sdf_writer_closes(tmp_path: Path) -> None:
    """Closing writes what is buffered, once; dropping or an error closes too.

    Args:
        tmp_path: pytest fixture of a temporary directory
    """
    block = SDBlock("t\nM  END\n")
    out = tmp_path / "out.sdf"
    writer = sdf_writer(out)
    writer.write(block)
    assert out.read_text() == ""
    writer.close()
    writer.close()
    assert out.read_text() == sdf_of(block)
    with pytest.raises(ValueError, match="closed"):
        writer.write(block)
    writer = sdf_writer(out)
    writer.write(block)
    del writer
    assert out.read_text() == sdf_of(block)
    with pytest.raises(RuntimeError), sdf_writer(out) as writer:
        for each in _then_raise([block]):
            writer.write(each)
    assert out.read_text() == sdf_of(block)
    if os.path.exists("/dev/full"):
        writer = sdf_writer("/dev/full")
        writer.write(block)
        with pytest.raises(OSError) as raised:
            writer.close()
        assert raised.value.errno == errno.ENOSPC
        writer.close()


def test_sdf_writer_rejects_a_block_and_goes_on(tmp_path: Path) -> None:
    """A rejected block raises, and the writer keeps the blocks around it.

    Args:
        tmp_path: pytest fixture of a temporary directory
    """
    source = tmp_path / "latin1.sdf"
    source.write_bytes(b"t\nM  END\n> <A>\n\xe9\n\n$$$$\n")
    (undecodable,) = parse_sdf(source)
    a, b = SDBlock("a\nM  END\n"), SDBlock("b\nM  END\n")
    handle = io.StringIO()
    with sdf_writer(tmp_path / "out.sdf") as to_path, sdf_writer(handle) as to_handle:
        for writer in (to_path, to_handle):
            writer.write(a)
            with pytest.raises(TypeError):
                writer.write("t\nM  END\n")  # type: ignore[arg-type]
        with pytest.raises(TypeError, match="mixed"):
            to_handle.write(SDBlock(b"t\nM  END\n"))  # type: ignore[arg-type]
        with pytest.raises(UnicodeDecodeError) as raised:
            to_path.write(undecodable)
        assert raised.value.__notes__[0].startswith("Writing block 1 (from 0),")
        to_path.write(b)
        to_handle.write(b)
    expected = "a\nM  END\n$$$$\nb\nM  END\n$$$$\n"
    assert (tmp_path / "out.sdf").read_text() == handle.getvalue() == expected
