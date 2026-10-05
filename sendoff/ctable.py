"""Reading and validating connection tables from SDFs."""

from sendoff.native import CTable

__all__ = [
    "CTable",
    "IndicesMismatchError",
    "IndicesOutOfOrderError",
    "IndicesDuplicateError",
]


class IndicesMismatchError(Exception):
    """When then number of atom or bond lines does not match the counts line."""


class IndicesOutOfOrderError(Exception):
    """When then atom or bond lines are not in increasing index order."""


class IndicesDuplicateError(Exception):
    """When then atom or bond lines have duplicate indices."""
