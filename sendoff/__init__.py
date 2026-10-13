"""Get metadata of SD Files in a lightweight way.

The `sendoff` library provides lightweight tools for inspecting the metadata
of chemical SD Files.  It strives to do as little modification of the input
files as possible.
"""

import sendoff.native as native
from sendoff.ctable import (
    CTable,
    CTableError,
    IndicesDuplicateError,
    IndicesMismatchError,
    IndicesOutOfOrderError,
    MalformedCTableError,
)
from sendoff.native import MalformedRecordError, SDBlock, parse_sdf, write_sdf

__all__ = [
    "CTable",
    "CTableError",
    "IndicesDuplicateError",
    "IndicesMismatchError",
    "IndicesOutOfOrderError",
    "MalformedCTableError",
    "MalformedRecordError",
    "SDBlock",
    "parse_sdf",
    "write_sdf",
]

__version__ = native.__version__
