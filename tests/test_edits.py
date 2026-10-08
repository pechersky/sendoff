# mypy: allow-untyped-decorators
"""Pin the exact SDF text each record edit writes.

A record the block writes itself has a ``> <name>`` header, with one space
and no other fields, then its value lines and one blank line, every line
ending as the block's title line does. Edits copy every line they do not
rewrite byte for byte, including other writers' header styles and lines
between records.
"""
from typing import Any, Callable

import pytest

from sendoff.sdblock import SDBlock

# Headers as RDKit writes them, with a record number and a registry number,
# a duplicate name, and a line that is not a record.
SOURCE = (
    "t\nM  END\n"
    ">  <A>  (1) \n1\n\n"
    "> 7 <B> (MFCD01)\nb1\nb2\n\n"
    "note\n\n"
    ">  <A>\ndup\n\n"
    "$$$$\n"
)
KEPT = (
    "t\nM  END\n>  <A>  (1) \n1\n\n> 7 <B> (MFCD01)\nb1\nb2\n\nnote\n\n>  <A>\ndup\n\n"
)
UNTERMINATED = "t\nM  END\n> <A>\n1\n$$$$\n"
SPACED = "t\nM  END\n> <A>\n1\n\n\n\n$$$$\n"
EMPTY = "t\nM  END\n$$$$\n"

FORMS: dict[str, Callable[[str], Any]] = {
    "str": lambda text: text,
    "bytes": str.encode,
    "crlf-bytes": lambda text: str.encode(text.replace("\n", "\r\n")),
}

CASES = [
    pytest.param(SOURCE, "to_sdf", (), SOURCE, id="untouched"),
    pytest.param(
        SOURCE,
        "set_record",
        ("A", "x"),
        "t\nM  END\n> <A>\nx\n\n> 7 <B> (MFCD01)\nb1\nb2\n\nnote\n\n$$$$\n",
        id="set-rewrites-the-first-and-drops-duplicates",
    ),
    pytest.param(
        SOURCE,
        "set_record",
        ("B", "x"),
        "t\nM  END\n>  <A>  (1) \n1\n\n> <B>\nx\n\nnote\n\n>  <A>\ndup\n\n$$$$\n",
        id="set-rewrites-a-numbered-multiline-record",
    ),
    pytest.param(
        SOURCE,
        "set_record",
        ("C", "x"),
        KEPT + "> <C>\nx\n\n$$$$\n",
        id="set-appends-a-new-name",
    ),
    pytest.param(
        SOURCE,
        "append_record",
        ("A", "x"),
        KEPT + "> <A>\nx\n\n$$$$\n",
        id="append-adds-a-duplicate",
    ),
    pytest.param(
        SOURCE,
        "append_record",
        ("C", "x\ny"),
        KEPT + "> <C>\nx\ny\n\n$$$$\n",
        id="append-multiline",
    ),
    pytest.param(
        UNTERMINATED,
        "append_record",
        ("C", "x"),
        "t\nM  END\n> <A>\n1\n\n> <C>\nx\n\n$$$$\n",
        id="append-ends-an-unterminated-record",
    ),
    pytest.param(
        UNTERMINATED,
        "set_record",
        ("A", "x"),
        "t\nM  END\n> <A>\nx\n\n$$$$\n",
        id="set-ends-an-unterminated-record",
    ),
    pytest.param(
        SPACED,
        "append_record",
        ("C", "x"),
        "t\nM  END\n> <A>\n1\n\n\n\n> <C>\nx\n\n$$$$\n",
        id="append-keeps-extra-blank-lines",
    ),
    pytest.param(
        EMPTY,
        "append_record",
        ("C", "x"),
        "t\nM  END\n> <C>\nx\n\n$$$$\n",
        id="append-to-no-records",
    ),
    pytest.param(
        EMPTY,
        "append_record",
        ("C", ""),
        "t\nM  END\n> <C>\n\n\n$$$$\n",
        id="append-an-empty-value",
    ),
    pytest.param(
        SOURCE,
        "remove_record",
        ("A",),
        "t\nM  END\n> 7 <B> (MFCD01)\nb1\nb2\n\nnote\n\n$$$$\n",
        id="remove-drops-duplicates",
    ),
    pytest.param(
        SOURCE,
        "remove_record",
        ("B",),
        "t\nM  END\n>  <A>  (1) \n1\n\nnote\n\n>  <A>\ndup\n\n$$$$\n",
        id="remove-keeps-neighbors",
    ),
    pytest.param(SOURCE, "remove_record", ("C",), SOURCE, id="remove-missing"),
    pytest.param(
        SOURCE,
        "replace_records",
        ([("A", "1"), ("B", "b1\nb2"), ("A", "dup")],),
        "t\nM  END\n> <A>\n1\n\n> <B>\nb1\nb2\n\n> <A>\ndup\n\n$$$$\n",
        id="replace-rewrites-every-header-and-drops-other-lines",
    ),
    pytest.param(
        SOURCE,
        "replace_records",
        ({"C": "x", "A": "y"},),
        "t\nM  END\n> <C>\nx\n\n> <A>\ny\n\n$$$$\n",
        id="replace-from-a-mapping",
    ),
    pytest.param(
        SOURCE, "replace_records", ([],), "t\nM  END\n$$$$\n", id="replace-with-none"
    ),
]


def convert(value: Any, form: Callable[[str], Any]) -> Any:
    """Convert ``value``, a str or a container of str, to ``form``.

    Args:
        value: a str, a mapping of str, or a list of str tuples
        form: how to convert each str

    Returns:
        ``value`` with each str converted
    """
    if isinstance(value, str):
        return form(value)
    if isinstance(value, dict):
        return {form(name): form(text) for name, text in value.items()}
    return [tuple(map(form, pair)) for pair in value]


@pytest.mark.parametrize("form", FORMS.values(), ids=FORMS.keys())
@pytest.mark.parametrize("source, method, args, expected", CASES)
def test_edit_writes_exactly(
    form: Callable[[str], Any],
    source: str,
    method: str,
    args: tuple[Any, ...],
    expected: str,
) -> None:
    """Each edit writes exactly the expected SDF, which reads back the same.

    Args:
        form: str, bytes, or bytes with CRLF line endings
        source: the block's SDF before the edit
        method: the edit to make
        args: the edit's arguments, as str
        expected: the block's SDF after the edit
    """
    block = next(SDBlock.from_sdf(form(source)))
    getattr(block, method)(*(convert(arg, form) for arg in args))
    written = block.to_sdf()
    assert written == form(expected)
    reread = next(SDBlock.from_sdf(written))
    assert reread.records() == block.records()
    assert reread.to_sdf() == written
