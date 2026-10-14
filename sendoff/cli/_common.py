"""Plumbing shared by the ``python -m sendoff.cli.<command>`` modules."""

from __future__ import annotations

import argparse
import os
import signal
import sys
from collections.abc import Callable, Iterator, Sequence
from contextlib import contextmanager
from typing import BinaryIO, NoReturn

from sendoff import SDBlock, parse_sdf

STDIN = "-"


def parser(command: str, description: str) -> argparse.ArgumentParser:
    """Return an argument parser for a command.

    Args:
        command: module name under ``sendoff.cli``
        description: help text

    Returns:
        The parser
    """
    return argparse.ArgumentParser(
        prog=f"python -m sendoff.cli.{command}", description=description
    )


def add_file(p: argparse.ArgumentParser) -> None:
    """Add the optional input file operand, with -i and -o.

    Args:
        p: the command's parser
    """
    p.add_argument("file", nargs="?", metavar="FILE", help="SD file (default: stdin)")
    add_io(p)


def add_io(p: argparse.ArgumentParser, *, many: bool = False) -> None:
    """Add -i, an input file in place of the FILE operand, and -o.

    Args:
        p: the command's parser
        many: whether -i may repeat
    """
    p.add_argument(
        "-i",
        "--input",
        action="append" if many else "store",
        metavar="FILE",
        help="SD file, in place of FILE" + ("; repeatable" if many else ""),
    )
    p.add_argument(
        "-o", "--output", default=STDIN, metavar="FILE", help="default: stdout"
    )


def input_path(
    p: argparse.ArgumentParser, operand: str | None, option: str | None
) -> str:
    """Return the input file, from the FILE operand or -i, or stdin.

    Args:
        p: the command's parser
        operand: the FILE operand, if given
        option: -i, if given

    Returns:
        The path, or ``-`` for stdin
    """
    if operand is not None and option is not None:
        p.error("give FILE or -i, not both")
    if option is not None:
        return option
    return STDIN if operand is None else operand


def output(
    p: argparse.ArgumentParser, path: str, *inputs: str | None
) -> str | BinaryIO:
    """Return -o's path, or stdout for ``-``, refusing an input file.

    Writing over a mapped input would truncate it while it is read.

    Args:
        p: the command's parser
        path: -o
        inputs: files the command reads

    Returns:
        The path, or stdout's buffer
    """
    if path == STDIN:
        return sys.stdout.buffer
    for source in inputs:
        if source is not None and source != STDIN and same_file(source, path):
            p.error(f"-o {path} is also an input; write elsewhere, then mv")
    return path


def same_file(a: str, b: str) -> bool:
    """Return whether two paths name the same existing file.

    Args:
        a: a path
        b: another path

    Returns:
        Whether both exist and are one file
    """
    try:
        return os.path.samefile(a, b)
    except OSError:
        return False


@contextmanager
def opened(p: argparse.ArgumentParser, target: str | BinaryIO) -> Iterator[BinaryIO]:
    """Open an output from `output` for writing bytes.

    Args:
        p: the command's parser
        target: a path or a binary file

    Yields:
        A binary file, closed after if opened here
    """
    if not isinstance(target, str):
        yield target
        return
    try:
        f = open(target, "wb")  # noqa: SIM115
    except OSError as exc:
        p.exit(1, f"{p.prog}: {exc}\n")
    with f:
        yield f


class Blocks:
    """Bytes blocks of an SD file, numbering the block being worked on.

    ``borrow`` lets blocks of a mapped file read it in place, for commands
    that never write over their input.
    """

    def __init__(self, path: str, *, borrow: bool = False) -> None:
        self.path = path
        self.borrow = borrow
        self.number = 0

    def __iter__(self) -> Iterator[SDBlock[bytes]]:
        """Yield the blocks.

        Yields:
            Each block, in bytes mode
        """
        source = sys.stdin.buffer if self.path == STDIN else os.fsencode(self.path)
        blocks = parse_sdf(source, borrow=self.borrow)
        for self.number, block in enumerate(blocks, 1):
            yield block


def add_sources(
    p: argparse.ArgumentParser, what: str
) -> argparse._MutuallyExclusiveGroup:
    """Add --value and --from, with -z and --skip-empty.

    Args:
        p: the command's parser
        what: what each value becomes, for help

    Returns:
        The required group of sources, for more of them
    """
    sources = p.add_mutually_exclusive_group(required=True)
    sources.add_argument("--value", metavar="VALUE", help=f"the {what} of every block")
    sources.add_argument(
        "--from",
        dest="source",
        metavar="VALUES",
        help=f"file of one {what} per line, in block order; - for stdin",
    )
    p.add_argument("-z", "--null", action="store_true", help="VALUES end with NUL")
    p.add_argument(
        "--skip-empty",
        action="store_true",
        help=f"leave a block unchanged where its {what} is empty",
    )
    return sources


def pairs(
    p: argparse.ArgumentParser, args: argparse.Namespace, blocks: Blocks
) -> Iterator[tuple[SDBlock[bytes], bytes]]:
    """Pair each block with its value from --value or --from.

    Args:
        p: the command's parser
        args: parsed arguments from add_sources
        blocks: the input

    Returns:
        The pairs, raising ValueError if --from and the input differ in length
    """
    if args.source is None:
        value = os.fsencode(args.value)
        return ((block, value) for block in blocks)
    if args.source == STDIN == blocks.path:
        p.error("--from - needs a FILE")
    if args.source == STDIN:
        data = sys.stdin.buffer.read()
    else:
        with open(args.source, "rb") as f:
            data = f.read()
    values = data.split(b"\0" if args.null else b"\n")
    if values[-1] == b"":
        values.pop()
    return aligned(blocks, values, "stdin" if args.source == STDIN else args.source)


def aligned(
    blocks: Blocks, values: list[bytes], source: str
) -> Iterator[tuple[SDBlock[bytes], bytes]]:
    """Pair blocks with values in order.

    Args:
        blocks: the input
        values: one per block
        source: where the values came from, for errors

    Yields:
        Each block and its value

    Raises:
        ValueError: if there are fewer or more values than blocks
    """
    rest = iter(values)
    for block in blocks:
        value = next(rest, None)
        if value is None:
            raise ValueError(f"{source} has only {len(values)} values")
        yield block, value
    if blocks.number < len(values):
        raise ValueError(f"no more blocks, but {source} has {len(values)} values")


@contextmanager
def failures(p: argparse.ArgumentParser, blocks: Blocks) -> Iterator[None]:
    """Exit with status 1 and the block number on bad input.

    Args:
        p: the command's parser
        blocks: the input being read

    Yields:
        Nothing
    """
    try:
        yield
    except (OSError, ValueError) as exc:
        where = f"block {blocks.number}: " if blocks.number else ""
        p.exit(1, f"{p.prog}: {where}{exc}\n")


def run(main: Callable[[Sequence[str] | None], int]) -> NoReturn:
    """Run a command as a program, dying quietly when its reader goes away.

    Args:
        main: the command's main function
    """
    if hasattr(signal, "SIGPIPE"):
        signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    sys.exit(main(None))
