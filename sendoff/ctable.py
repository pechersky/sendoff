"""Reading and validating connection tables from SDFs."""

from sendoff.native import CTable

__all__ = [
    "CTable",
    "CTableError",
    "IndicesDuplicateError",
    "IndicesMismatchError",
    "IndicesOutOfOrderError",
    "MalformedCTableError",
]


class CTableError(ValueError):
    """A connection table's contents break its format or its counts."""


class MalformedCTableError(CTableError):
    """A connection table line breaks the format; the message names the line.

    Lines are numbered from 1, the title being line 1.
    """


class IndicesMismatchError(CTableError):
    """The number of atom or bond lines does not match the counts line."""


class IndicesOutOfOrderError(CTableError):
    """The atom or bond lines are not in increasing one-based index order."""


class IndicesDuplicateError(CTableError):
    """Atom or bond lines share an index."""
