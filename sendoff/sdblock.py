"""Separating an SDF into molecule "SDBlock"s of chemical data and metadata."""
from __future__ import annotations

import os
from collections import deque
from dataclasses import dataclass
from io import TextIOWrapper
from itertools import groupby, tee
from typing import Any, Iterable, Iterator, Tuple, Union

Pathy = Union[str, bytes, os.PathLike[Any]]


@dataclass
class SDBlock:
    """Handle a single molecule block from an SD file."""

    title: str
    lines: Iterator[str]

    @classmethod
    def from_block_lines(cls, lines: Iterable[str]) -> SDBlock:
        """Parse in an SDBlock from a sequence of lines.

        Args:
            lines: an Iterable of str that
                optionally have whitespace that is stripped

        Returns:
            An SDBlock with a parsed in title and lines
        """
        title_gen, lines = tee(lines)
        title = next(title_gen).strip()
        return SDBlock(title, lines)

    @classmethod
    def from_lines(cls, lines: Iterable[str]) -> Iterator[SDBlock]:
        """Parse in SDBlocks from a sequence of lines.

        Args:
            lines: an Iterable of str that optionally
                have whitespace that is stripped,
                likely separated by the $$$$ SD delimiter

        Yields:
            SDBlocks parsed in from the lines
        """
        block: deque[str] = deque()
        for line in lines:
            block.append(line)
            if line.startswith("$$$$"):
                yield cls.from_block_lines(block)
                block = deque()

    def parse_mdl(self) -> Iterable[str]:
        """Parse in an MDL block from the ingested lines.

        Yields:
            str lines comprising the MDL block
        """
        for line in self.lines:
            yield line
            if line.startswith("M  END"):
                return

    def records(self) -> Iterable[Tuple[str, str]]:
        """Generate SD metadata records one by one.

        Yields:
            Tuples of str, str of record, value.
                The value includes any newlines if it is multiline
        """
        list(self.parse_mdl())  # consume to rid ourselves of lines
        # make data chunks by breaking on empty lines
        for _, chunk in groupby((line.strip() for line in self.lines), bool):
            record_line: str = next(chunk)
            if not record_line.startswith("> "):
                continue
            # something of the form `> ___<___>___` where `_` is anything
            record_name = (
                record_line.split("> ", 1)[1]
                .strip()
                .rsplit(">")[0]
                .split("<", 1)[1]
            )
            yield record_name, str.join("\n", chunk)

    def write(self, outh: TextIOWrapper) -> None:
        """Write an SDBlock to a file-like handle.

        Args:
            outh: file-like handle, such as what is returned by open(..., "w")

        """
        print(list(self.lines))
        for line in self.lines:
            print(line, file=outh)


def parse_sdf(sdfpath: Pathy) -> Iterator[SDBlock]:
    """Parse in SDBlocks from a path. Wrapper around SDBlock.from_lines.

    Args:
        sdfpath: file-like handle, such as what is returned by open(..., "r")

    Returns:
        An Iterator of SDBlocks parsed in from the lines in the file
    """
    return SDBlock.from_lines(open(sdfpath).readlines())
