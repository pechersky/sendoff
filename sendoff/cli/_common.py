"""Plumbing shared by the ``python -m sendoff.cli.<command>`` modules."""

from __future__ import annotations

import argparse
import os
import signal
import sys
from collections.abc import Callable, Iterator, Sequence
from contextlib import contextmanager
from typing import NoReturn

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
    """Add the optional input file operand.

    Args:
        p: the command's parser
    """
    p.add_argument(
        "file",
        nargs="?",
        default=STDIN,
        metavar="FILE",
        help="SD file (default: stdin)",
    )


class Blocks:
    """Bytes blocks of an SD file, numbering the block being worked on."""

    def __init__(self, path: str) -> None:
        self.path = path
        self.number = 0

    def __iter__(self) -> Iterator[SDBlock[bytes]]:
        """Yield the blocks.

        Yields:
            Each block, in bytes mode
        """
        source = sys.stdin.buffer if self.path == STDIN else os.fsencode(self.path)
        for self.number, block in enumerate(parse_sdf(source), 1):
            yield block


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
