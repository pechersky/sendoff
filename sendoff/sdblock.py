"""Separating an SDF into molecule "SDBlock"s of chemical data and metadata."""
from __future__ import annotations

import os
from collections import deque
from dataclasses import dataclass
from io import TextIOWrapper
from itertools import chain, groupby
from typing import Iterable, Iterator, Tuple, Union

Pathy = Union[str, bytes, os.PathLike]


@dataclass
class SDBlock:
    """Handle a single molecule block from an SD file."""

    title: str
    mdl: deque[str]
    metadata: deque[str]

    @classmethod
    def parse_mdl(cls, lines: Iterable[str]) -> Iterable[str]:
        """Parse in an MDL block from the ingested lines.

        Args:
            lines: an Iterable of str that comprises the MDL block

        Yields:
            str lines comprising the MDL block
        """
        for line in lines:
            yield line
            if line.startswith("M  END"):
                return

    @classmethod
    def parse_metadata(cls, lines: Iterable[str]) -> Iterable[str]:
        """Parse in the SDF metadata from the ingested lines.

        Args:
            lines: an Iterable of str that metadata, optionally
                including the $$$$ delimiter

        Yields:
            str lines comprising the metadata, excluding the $$$$ delimiter
        """
        for line in lines:
            if line.startswith("$$$$"):
                return
            yield line

    @classmethod
    def from_block_lines(cls, lines: Iterable[str]) -> SDBlock:
        """Parse in an SDBlock from a sequence of lines.

        Args:
            lines: an Iterable of str that
                optionally have whitespace that is stripped

        Returns:
            An SDBlock with a parsed in title and lines
        """
        iterlines = iter(lines)
        title = next(iterlines).strip()
        mdl = deque(cls.parse_mdl(iterlines))
        metadata = deque(cls.parse_metadata(iterlines))
        return SDBlock(title, mdl, metadata)

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

    def records(self) -> Iterable[Tuple[str, str]]:
        """Generate SD metadata records one by one.

        Yields:
            Tuples of str, str of record, value.
                The value includes any newlines if it is multiline
        """
        # make data chunks by breaking on empty lines
        for _, chunk in groupby((line.strip() for line in self.metadata), bool):
            record_line: str = next(chunk)
            if not record_line.startswith("> "):
                continue
            # something of the form `> ___<___>___` where `_` is anything
            record_name = (
                record_line.split("> ", 1)[1].strip().rsplit(">")[0].split("<", 1)[1]
            )
            yield record_name, str.join("\n", chunk)

    def write(self, outh: TextIOWrapper) -> None:
        """Write an SDBlock to a file-like handle.

        Args:
            outh: file-like handle, such as what is returned by open(..., "w")

        """
        print(self.title, file=outh)
        for line in chain(self.mdl, self.metadata):
            outh.write(line)
        outh.write("$$$$\n")

    def append_record(self, record_name: str, value: str) -> None:
        """Append a field and value to the SD data record.

        Args:
            record_name: str for field's name
            value: str for value, with no terminating newline

        """
        self.metadata.append(f"> <{record_name}>\n")
        self.metadata.append(value + "\n")
        self.metadata.append("\n")


def parse_sdf(sdfpath: Pathy) -> Iterator[SDBlock]:
    """Parse in SDBlocks from a path. Wrapper around SDBlock.from_lines.

    Args:
        sdfpath: file-like handle, such as what is returned by open(..., "r")

    Returns:
        An Iterator of SDBlocks parsed in from the lines in the file
    """
    return SDBlock.from_lines(open(sdfpath).readlines())
