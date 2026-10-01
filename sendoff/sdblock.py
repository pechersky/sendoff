"""Separating an SDF into molecule "SDBlock"s of chemical data and metadata."""

from __future__ import annotations

import os
from typing import Iterator, Union, cast

import sendoff.native as native

Pathy = Union[str, bytes, "os.PathLike[str]"]

SDBlock = native.SDBlock


def parse_sdf(sdfpath: Pathy) -> Iterator[SDBlock]:
    """Parse in SDBlocks from a path. Wrapper around SDBlock.from_lines.

    Args:
        sdfpath: file-like handle, such as what is returned by open(..., "r")

    Returns:
        An Iterator of SDBlocks parsed in from the lines in the file
    """
    return cast(Iterator[SDBlock], SDBlock.from_lines(open(sdfpath).readlines()))
