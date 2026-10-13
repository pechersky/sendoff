"""Separating an SDF into molecule "SDBlock"s of chemical data and metadata."""

from __future__ import annotations

import os
from typing import Union

import sendoff.native as native

Pathy = Union[str, bytes, "os.PathLike[str]", "os.PathLike[bytes]"]

MalformedRecordError = native.MalformedRecordError
SDBlock = native.SDBlock
parse_sdf = native.parse_sdf
write_sdf = native.write_sdf
