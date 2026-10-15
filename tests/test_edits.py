# mypy: allow-untyped-decorators
"""Pin the exact SDF text each record edit writes.

A header the block writes is ``> <C>``, or ``>  <C>`` with
``header="padded"``. With ``header="copy"`` it copies the text around the name
on the first header with a name, so after RDKit's ``>  <A>  (1)`` it writes
``>  <C>  (1)`` with the same trailing space; without such a header, or when
its name has no closing ``>``, it writes ``> <C>``. Value lines and one blank
line follow, every line ending as the block's title line does.
``set_record`` keeps the header line of the record it sets, and
``replace_records`` gives the k-th record with a name the k-th old header line
with that name. Edits copy every line they do not rewrite byte for byte,
including lines between records, which ``replace_records`` drops.
"""
from functools import partial
from textwrap import dedent
from typing import Any, Callable

import pytest

from sendoff.sdblock import SDBlock
from tests.fixtures import sdf_of

# A header as RDKit writes it, with a record number and a trailing space
# (written \x20 so editors keep it), a header with a registry number, a
# duplicate name, and a line that is not a record.
SOURCE = """\
t
M  END
>  <A>  (1)\x20
1

> 7 <B> (MFCD01)
b1
b2

note

>  <A>
dup

$$$$
"""

# The last record has no blank line after its value.
UNTERMINATED = """\
t
M  END
> <A>
1
$$$$
"""

SPACED = """\
t
M  END
> <A>
1



$$$$
"""

EMPTY = """\
t
M  END
$$$$
"""

# The only header's name has no closing ``>``.
UNCLOSED = """\
t
M  END
>  <A
1

$$$$
"""

FORMS: dict[str, Callable[[str], Any]] = {
    "str": lambda text: text,
    "bytes": str.encode,
    "crlf-bytes": lambda text: str.encode(text.replace("\n", "\r\n")),
}

# Each case is a source block, an edit and its arguments, and the SDF the
# block writes after the edit, dedented.
CASES = [
    pytest.param(SOURCE, "__str__", (), SOURCE, id="untouched"),
    pytest.param(
        SOURCE,
        "set_record",
        ("A", "x"),
        """\
        t
        M  END
        >  <A>  (1)\x20
        x

        > 7 <B> (MFCD01)
        b1
        b2

        note

        $$$$
        """,
        id="set-keeps-the-first-header-and-drops-duplicates",
    ),
    pytest.param(
        SOURCE,
        "set_record",
        ("B", "x"),
        """\
        t
        M  END
        >  <A>  (1)\x20
        1

        > 7 <B> (MFCD01)
        x

        note

        >  <A>
        dup

        $$$$
        """,
        id="set-keeps-a-registry-header",
    ),
    pytest.param(
        SOURCE,
        "set_record",
        ("C", "x"),
        """\
        t
        M  END
        >  <A>  (1)\x20
        1

        > 7 <B> (MFCD01)
        b1
        b2

        note

        >  <A>
        dup

        > <C>
        x

        $$$$
        """,
        id="set-appends-a-new-name",
    ),
    pytest.param(
        SOURCE,
        "append_record",
        ("A", "x"),
        """\
        t
        M  END
        >  <A>  (1)\x20
        1

        > 7 <B> (MFCD01)
        b1
        b2

        note

        >  <A>
        dup

        > <A>
        x

        $$$$
        """,
        id="append-adds-a-duplicate",
    ),
    pytest.param(
        SOURCE,
        "append_record",
        ("C", "x\ny"),
        """\
        t
        M  END
        >  <A>  (1)\x20
        1

        > 7 <B> (MFCD01)
        b1
        b2

        note

        >  <A>
        dup

        > <C>
        x
        y

        $$$$
        """,
        id="append-multiline",
    ),
    pytest.param(
        UNTERMINATED,
        "append_record",
        ("C", "x"),
        """\
        t
        M  END
        > <A>
        1

        > <C>
        x

        $$$$
        """,
        id="append-ends-an-unterminated-record",
    ),
    pytest.param(
        UNTERMINATED,
        "set_record",
        ("A", "x"),
        """\
        t
        M  END
        > <A>
        x

        $$$$
        """,
        id="set-ends-an-unterminated-record",
    ),
    pytest.param(
        SPACED,
        "append_record",
        ("C", "x"),
        """\
        t
        M  END
        > <A>
        1



        > <C>
        x

        $$$$
        """,
        id="append-keeps-extra-blank-lines",
    ),
    pytest.param(
        EMPTY,
        "append_record",
        ("C", "x"),
        """\
        t
        M  END
        > <C>
        x

        $$$$
        """,
        id="append-to-no-records",
    ),
    pytest.param(
        EMPTY,
        "append_record",
        ("C", ""),
        """\
        t
        M  END
        > <C>


        $$$$
        """,
        id="append-an-empty-value",
    ),
    pytest.param(
        SOURCE,
        "remove_record",
        ("A",),
        """\
        t
        M  END
        > 7 <B> (MFCD01)
        b1
        b2

        note

        $$$$
        """,
        id="remove-drops-duplicates",
    ),
    pytest.param(
        SOURCE,
        "remove_record",
        ("B",),
        """\
        t
        M  END
        >  <A>  (1)\x20
        1

        note

        >  <A>
        dup

        $$$$
        """,
        id="remove-keeps-neighbors",
    ),
    pytest.param(SOURCE, "remove_record", ("C",), SOURCE, id="remove-missing"),
    pytest.param(
        SOURCE,
        "replace_records",
        ([("A", "1"), ("B", "b1\nb2"), ("A", "dup")],),
        """\
        t
        M  END
        >  <A>  (1)\x20
        1

        > 7 <B> (MFCD01)
        b1
        b2

        >  <A>
        dup

        $$$$
        """,
        id="replace-keeps-headers-and-drops-other-lines",
    ),
    pytest.param(
        SOURCE,
        "replace_records",
        ([("B", "x"), ("A", "y"), ("A", "z"), ("A", "w")],),
        """\
        t
        M  END
        > 7 <B> (MFCD01)
        x

        >  <A>  (1)\x20
        y

        >  <A>
        z

        > <A>
        w

        $$$$
        """,
        id="replace-keeps-each-names-headers-in-order",
    ),
    pytest.param(
        SOURCE,
        "replace_records",
        ({"C": "x", "A": "y"},),
        """\
        t
        M  END
        > <C>
        x

        >  <A>  (1)\x20
        y

        $$$$
        """,
        id="replace-from-a-mapping",
    ),
    pytest.param(
        SOURCE,
        "replace_records",
        ([],),
        """\
        t
        M  END
        $$$$
        """,
        id="replace-with-none",
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


def assert_writes(block: SDBlock[Any], expected: Any) -> None:
    """Assert that ``block`` writes ``expected``, which reads back the same.

    Args:
        block: the edited block
        expected: its SDF, in the block's type
    """
    written = sdf_of(block)
    assert written == expected
    reread = next(SDBlock.from_sdf(written))
    assert reread.records() == block.records()
    assert sdf_of(reread) == written


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
        expected: the block's SDF after the edit, dedented
    """
    block = next(SDBlock.from_sdf(form(source)))
    getattr(block, method)(*(convert(arg, form) for arg in args))
    assert_writes(block, form(dedent(expected)))


# SOURCE with a record C appended under a {header} line.
APPENDED = """\
t
M  END
>  <A>  (1)\x20
1

> 7 <B> (MFCD01)
b1
b2

note

>  <A>
dup

{header}
x

$$$$
"""

# SOURCE with its SD data replaced by a record C.
REPLACED = """\
t
M  END
{header}
x

$$$$
"""

# The header line each style writes for C after SOURCE's headers; None passes
# no style.
STYLES = {
    None: "> <C>",
    "plain": "> <C>",
    "padded": ">  <C>",
    "copy": ">  <C>  (1)\x20",
}

# Each way to write a new header: a method, or the constructor, and its
# arguments after SOURCE, and the SDF it writes.
STYLED = [
    pytest.param("append_record", ("C", "x"), APPENDED, id="append"),
    pytest.param("set_record", ("C", "x"), APPENDED, id="set"),
    pytest.param("replace_records", ({"C": "x"},), REPLACED, id="replace"),
    pytest.param("SDBlock", ({"C": "x"},), APPENDED, id="construct"),
]


@pytest.mark.parametrize("form", FORMS.values(), ids=FORMS.keys())
@pytest.mark.parametrize(
    "header, line", STYLES.items(), ids=[str(style) for style in STYLES]
)
@pytest.mark.parametrize("method, args, expected", STYLED)
def test_new_headers_follow_the_style(
    form: Callable[[str], Any],
    header: str | None,
    line: str,
    method: str,
    args: tuple[Any, ...],
    expected: str,
) -> None:
    """Every way to write a new header line writes it in the asked style.

    Args:
        form: str, bytes, or bytes with CRLF line endings
        header: the style, or None for the default
        line: the header line the style writes
        method: the edit, or "SDBlock" for the constructor
        args: the edit's arguments after SOURCE, as str
        expected: the block's SDF after the edit, with a {header} line
    """
    styled: dict[str, Any] = {} if header is None else {"header": header}
    converted = [convert(arg, form) for arg in args]
    if method == "SDBlock":
        block = SDBlock(form(SOURCE), *converted, **styled)
    else:
        block = next(SDBlock.from_sdf(form(SOURCE)))
        getattr(block, method)(*converted, **styled)
    assert_writes(block, form(expected.format(header=line)))


@pytest.mark.parametrize("form", FORMS.values(), ids=FORMS.keys())
@pytest.mark.parametrize(
    "source, expected",
    [
        pytest.param(
            EMPTY,
            """\
            t
            M  END
            > <C>
            x

            $$$$
            """,
            id="no-records",
        ),
        pytest.param(
            UNCLOSED,
            """\
            t
            M  END
            >  <A
            1

            > <C>
            x

            $$$$
            """,
            id="unclosed-name",
        ),
    ],
)
def test_copy_without_a_named_header_writes_plain(
    form: Callable[[str], Any], source: str, expected: str
) -> None:
    """With no header to copy, ``header="copy"`` writes ``> <C>``.

    Args:
        form: str, bytes, or bytes with CRLF line endings
        source: the block's SDF before the edit
        expected: the block's SDF after the edit, dedented
    """
    block = next(SDBlock.from_sdf(form(source)))
    block.append_record(form("C"), form("x"), header="copy")
    assert_writes(block, form(dedent(expected)))


@pytest.mark.parametrize("method, args, expected", STYLED)
def test_other_styles_raise_before_editing(
    method: str, args: tuple[Any, ...], expected: str
) -> None:
    """Another ``header`` raises ValueError and leaves the block as it was.

    Args:
        method: the edit, or "SDBlock" for the constructor
        args: the edit's arguments after SOURCE
        expected: unused
    """
    block = next(SDBlock.from_sdf(SOURCE))
    edit = partial(SDBlock, SOURCE) if method == "SDBlock" else getattr(block, method)
    with pytest.raises(ValueError, match="'plain', 'padded', or 'copy', not 'Plain'"):
        edit(*args, header="Plain")  # type: ignore[arg-type]
    assert sdf_of(block) == SOURCE
