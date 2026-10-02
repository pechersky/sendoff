"""Separating an SDF into molecule "SDBlock"s of chemical data and metadata."""

from __future__ import annotations

import os
from typing import Iterator, Union

import sendoff.native as native

Pathy = Union[str, bytes, "os.PathLike[str]", "os.PathLike[bytes]"]

SDBlock = native.SDBlock


def parse_sdf(sdfpath: Pathy) -> Iterator[SDBlock]:
    """Parse SDBlocks from a path, preserving text or byte input.

    Args:
        sdfpath: SDF path; str paths use text mode and bytes paths use binary mode

    Returns:
        An Iterator of SDBlocks parsed in from the lines in the file
    """
    path = os.fspath(sdfpath)
    mode = "rb" if isinstance(path, bytes) else "r"
    with open(path, mode) as sdf:
        yield from SDBlock.from_lines(sdf)
