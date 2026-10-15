# mypy: allow-untyped-decorators
"""Test the SDBlock record, setter, comparison, copy, and pickle API."""

import copy
import io
import pickle
import typing

import pytest

from sendoff.sdblock import MalformedRecordError, SDBlock
from tests.fixtures import sdf_of


def test_bytes_block_requires_binary_append_values() -> None:
    """Appending to a bytes block cannot introduce implicitly encoded text."""
    block = next(SDBlock.from_lines([b"title\n", b"\n", b"M  END\n"]))

    with pytest.raises(TypeError, match="record names must be bytes"):
        block.append_record("field", b"value")  # type: ignore[arg-type]
    with pytest.raises(TypeError, match="values must be bytes"):
        block.append_record(b"field", "value")  # type: ignore[arg-type]

    block.append_record(b"field", b"value")
    assert block.records() == [(b"field", b"value")]


def test_iter_records_preserves_order_duplicates_and_multiline_values() -> None:
    """Iterate records on demand without changing list-returning records()."""
    block = next(
        SDBlock.from_lines(
            [
                "Title\n",
                "\n",
                "M  END\n",
                "> <Record>\n",
                "first\n",
                "second\n",
                "\n",
                "> <Record>\n",
                "third\n",
                "\n",
                "$$$$\n",
            ]
        )
    )
    iterator = block.iter_records()
    assert iter(iterator) is iterator
    first = next(iterator)
    second = next(iterator)

    assert first == ("Record", "first\nsecond")
    assert second == ("Record", "third")
    assert first == ("Record", "first\nsecond")
    records = block.records()
    assert isinstance(records, list)
    assert records == [first, second]


def test_iter_records_defers_malformed_field_error() -> None:
    """A later malformed record does not prevent reading earlier fields."""
    block = next(
        SDBlock.from_lines(
            [
                "Title\n",
                "\n",
                "M  END\n",
                "> <Good>\n",
                "value\n",
                "\n",
                "> missing name\n",
                "other\n",
                "\n",
                "$$$$\n",
            ]
        )
    )
    iterator = block.iter_records()

    assert next(iterator) == ("Good", "value")
    with pytest.raises(MalformedRecordError, match="missing a name"):
        next(iterator)
    with pytest.raises(StopIteration):
        next(iterator)
    with pytest.raises(MalformedRecordError, match="missing a name"):
        block.records()


def test_iter_records_preserves_bytes_and_owns_its_block() -> None:
    """Byte record iteration remains lazy and retains source data safely."""
    block = next(
        SDBlock.from_lines(
            [
                b"Title\n",
                b"\n",
                b"M  END\n",
                b"> <Record>\n",
                b"first\xfe\n",
                b"second\n",
                b"\n",
                b"> <Record>\n",
                b"third\n",
                b"\n",
                b"$$$$\n",
            ]
        )
    )
    iterator = block.iter_records()
    del block

    record = next(iterator)

    assert record == (b"Record", b"first\xfe\nsecond")
    assert list(iterator) == [(b"Record", b"third")]


def test_iter_records_raises_when_the_records_change() -> None:
    """Like a dict, iteration stops loudly rather than skipping records."""
    block = SDBlock(
        "Title\n\n  0  0  0  0  0  0  0  0  0  0999 V2000\nM  END\n"
        "> <A>\n1\n\n> <B>\n2\n\n$$$$\n"
    )
    iterator = block.iter_records()
    assert next(iterator) == ("A", "1")
    block.title = "Changed"
    block.remove_record("missing")
    assert next(iterator) == ("B", "2")

    for edit in (
        lambda: block.set_record("A", "10"),
        lambda: block.append_record("C", "3"),
        lambda: block.remove_record("C"),
        lambda: block.replace_records(block.records()),
    ):
        iterator = block.iter_records()
        next(iterator)
        edit()
        with pytest.raises(RuntimeError, match="changed during iteration"):
            next(iterator)
        assert list(iterator) == []

    for name, value in block.records():
        block.set_record(name, value * 2)
    assert block.records() == [("A", "1010"), ("B", "22")]


def test_get_and_to_dict_match_dict_of_records() -> None:
    """get() and to_dict() let later duplicate records win, like dict()."""
    text = "t\nM  END\n> <A>\n1\n\n> <B>\nx\ny\n\n> <A>\n2\n\n$$$$\n"
    block = next(SDBlock.from_lines(io.StringIO(text)))
    byte_block = next(SDBlock.from_lines(io.BytesIO(text.encode())))

    assert block.to_dict() == dict(block.records()) == {"A": "2", "B": "x\ny"}
    assert byte_block.to_dict() == {b"A": b"2", b"B": b"x\ny"}
    assert block.get("A") == "2"
    assert block.get("missing") is None
    assert block.get("missing", "default") == "default"
    assert byte_block.get(b"B") == b"x\ny"
    with pytest.raises(TypeError, match="must be bytes"):
        byte_block.get("A")  # type: ignore[call-overload]
    with pytest.raises(TypeError, match="must be str"):
        block.get(b"A")  # type: ignore[call-overload]


def test_dict_style_access_matches_the_record_methods() -> None:
    """block[name], assignment, del, and iteration follow get/set/remove/names."""
    text = "t\nM  END\n> <A>  (1)\n1\n\n> <B>\nx\ny\n\n> <A>\n2\n\n$$$$\n"
    block = next(SDBlock.from_lines(io.StringIO(text)))
    byte_block = next(SDBlock.from_lines(io.BytesIO(text.encode())))

    assert block["A"] == block.get("A") == "2"
    assert byte_block[b"B"] == b"x\ny"
    assert list(block) == block.names() == ["A", "B", "A"]
    assert list(byte_block) == [b"A", b"B", b"A"]
    with pytest.raises(KeyError) as missing:
        block["missing"]
    assert missing.value.args == ("missing",)
    with pytest.raises(TypeError, match="must be bytes"):
        byte_block["A"]  # type: ignore[index]
    with pytest.raises(TypeError):
        block[0]  # type: ignore[index]
    assert not isinstance(block, typing.Sequence)

    block["A"] = "new"
    block["C"] = "c"
    del block["B"]
    with pytest.raises(KeyError):
        del block["B"]
    assert sdf_of(block) == "t\nM  END\n> <A>  (1)\nnew\n\n> <C>\nc\n\n$$$$\n"

    names = iter(block)
    assert next(names) == "A"
    block["D"] = "d"
    with pytest.raises(RuntimeError, match="changed during iteration"):
        next(names)


def test_record_mutators_keep_order_and_other_lines() -> None:
    """set_record edits in place; remove_record and replace_records rebuild."""
    text = "t\nM  END\n> <A>  (1)\n1\n\nnote\n\n> <B>\nb\n\n> <A>\n2\n\n$$$$\n"
    block = next(SDBlock.from_lines(io.StringIO(text)))

    block.set_record("A", "new")
    block.set_record("C", "c")
    assert block.records() == [("A", "new"), ("B", "b"), ("C", "c")]
    assert block.remove_record("B") == 1
    assert block.remove_record("B") == 0
    assert str(block) == "t\nM  END\n> <A>  (1)\nnew\n\nnote\n\n> <C>\nc\n\n$$$$\n"

    block.replace_records({"X": "1", "Y": "2"})
    assert block.records() == [("X", "1"), ("Y", "2")]
    block.replace_records((name, value) for name, value in [("X", "a"), ("X", "b")])
    assert block.records() == [("X", "a"), ("X", "b")]
    block.replace_records([])
    assert block.records() == []
    assert block.molblock == "t\nM  END\n"


def test_record_mutators_check_types_and_headers() -> None:
    """Names and values must match the block's type; bad headers raise."""
    byte_block = next(SDBlock.from_lines([b"t\n", b"M  END\n", b"> <A>\n", b"1\n"]))
    byte_block.set_record(b"A", b"2")
    assert byte_block.to_dict() == {b"A": b"2"}
    with pytest.raises(TypeError, match="record names must be bytes"):
        byte_block.set_record("A", b"2")  # type: ignore[arg-type]
    with pytest.raises(TypeError, match="values must be bytes"):
        byte_block.replace_records({b"A": "2"})  # type: ignore[arg-type]
    with pytest.raises(TypeError, match="record names must be bytes"):
        byte_block.remove_record("A")  # type: ignore[arg-type]
    with pytest.raises(TypeError):
        byte_block.replace_records([b"AB"])  # type: ignore[list-item]

    bad = next(SDBlock.from_lines(["t\n", "M  END\n", "> no name\n", "1\n"]))
    with pytest.raises(MalformedRecordError, match="missing a name"):
        bad.set_record("A", "1")
    with pytest.raises(MalformedRecordError, match="missing a name"):
        bad.remove_record("A")
    bad.replace_records({"A": "1"})
    assert bad.records() == [("A", "1")]


def test_record_names_are_shared_between_blocks() -> None:
    """Repeated names reuse one object; evicted and long names stay exact."""

    def check(sdf: typing.AnyStr, names: list[typing.AnyStr]) -> None:
        first, second = SDBlock(sdf), SDBlock(sdf)
        assert [n for n, _ in first.records()] == names
        assert list(first.to_dict()) == [n for n, _ in second.records()]
        assert all(n + v[-1:] == v for n, v in second.iter_records())
        assert first.records()[-1][0] is second.records()[-1][0]

    names = [f"N{index}" for index in range(5000)] + ["L" * 65, "A"]
    sdf = "t\nM  END\n" + str.join("", [f"> <{name}>\n{name}v\n\n" for name in names])
    check(sdf, names)
    check(sdf.encode(), [name.encode() for name in names])


def test_record_names_follow_changing_orders() -> None:
    """Reordered, missing, repeated, and new names between blocks read exactly."""
    orders = [
        ["A", "B", "C", "D"],
        ["A", "C", "D"],
        ["D", "C", "B", "A"],
        ["A", "A", "E", "B"],
        [],
        ["B", "C", "D", "A"],
        ["A", "B", "C", "D"],
    ]
    sdf = str.join(
        "",
        [
            "t\nM  END\n"
            + str.join("", [f"> <{n}>\n{n}{i}\n\n" for n in order])
            + "$$$$\n"
            for i, order in enumerate(orders)
        ],
    )
    expected = [[(n, f"{n}{i}") for n in order] for i, order in enumerate(orders)]
    for _ in range(2):
        blocks = list(SDBlock.from_sdf(sdf))
        assert [block.records() for block in blocks] == expected
        assert [block.to_dict() for block in blocks] == [dict(e) for e in expected]
        assert [list(block.iter_records()) for block in blocks] == expected
    lazy = blocks[2].iter_records()
    assert next(lazy) == expected[2][0]
    assert blocks[3].to_dict() == dict(expected[3])
    assert list(lazy) == expected[2][1:]
    byte_blocks = list(SDBlock.from_sdf(sdf.encode()))
    assert [block.records() for block in byte_blocks] == [
        [(n.encode(), v.encode()) for n, v in e] for e in expected
    ]


def test_sdblock_is_generic_at_runtime() -> None:
    """SDBlock[str] and SDBlock[bytes] work in evaluated annotations."""
    alias = SDBlock[str]
    assert typing.get_origin(alias) is SDBlock
    assert typing.get_args(alias) == (str,)
    assert typing.get_args(SDBlock[bytes]) == (bytes,)


def test_mdl_is_the_ctab_text() -> None:
    """The mdl property is the CTAB through M  END as one value, without the title."""
    text = "title\r\n  prog\r\n\r\nM  END\r\n> <A>\r\n1\r\n"
    block = next(SDBlock.from_lines(io.StringIO(text, newline=None)))
    byte_block = next(SDBlock.from_lines(io.BytesIO(text.encode())))

    assert block.mdl == "  prog\n\nM  END\n"
    assert block.title + "\n" + block.mdl == block.molblock
    assert block.mdl.splitlines() == ["  prog", "", "M  END"]
    assert byte_block.mdl == b"  prog\r\n\r\nM  END\r\n"
    assert byte_block.mdl.splitlines() == [b"  prog", b"", b"M  END"]
    assert next(SDBlock.from_lines(["t\n", "M  END\n"])).mdl == "M  END\n"
    assert next(SDBlock.from_lines(["t\n"])).mdl == ""
    no_end = next(SDBlock.from_lines(["t\n", "  prog\n", "> <A>\n", "$$$$\n"]))
    assert no_end.mdl == "  prog\n> <A>\n"


def test_names_len_and_contains_skip_values() -> None:
    """names(), len() and in read record names without building values."""
    sdf = "t\nM  END\n> <A>\n1\n2\n\n> <B>\n\n> <A>\n3\n\n$$$$\n"
    (block,) = SDBlock.from_sdf(sdf)
    (byte_block,) = SDBlock.from_sdf(sdf.encode())

    assert block.names() == [name for name, _ in block.records()] == ["A", "B", "A"]
    assert byte_block.names() == [b"A", b"B", b"A"]
    assert len(block) == len(byte_block) == 3
    assert "A" in block and "B" in block and "C" not in block
    assert b"A" in byte_block and "" not in block
    with pytest.raises(TypeError, match="record names must be bytes"):
        assert "A" in byte_block  # type: ignore[operator]

    empty = SDBlock("t\n\nM  END\n")
    assert (empty.names(), len(empty), "A" in empty) == ([], 0, False)
    assert empty, "a block without records is still true"

    block.append_record("C", "4")
    block.remove_record("A")
    assert (block.names(), len(block), "A" in block) == (["B", "C"], 2, False)
    malformed = SDBlock("t\n\nM  END\n> <A>\n1\n\n> no name\n1\n\n$$$$\n")
    calls = (lambda: len(malformed), lambda: "A" in malformed, malformed.names)
    for call in (*calls, lambda: malformed.get("A")) * 2:
        with pytest.raises(MalformedRecordError):
            call()


@pytest.mark.parametrize(
    "read",
    [
        lambda block: block.records(),
        lambda block: list(block.iter_records()),
        lambda block: block.get("B"),
        lambda block: block["A"],
        lambda block: list(block),
        lambda block: block.to_dict(),
        lambda block: block.set_record("B", "2"),
        lambda block: block.remove_record("B"),
        lambda block: block.__setitem__("B", "2"),
        lambda block: block.__delitem__("A"),
    ],
    ids=[
        "records",
        "iter_records",
        "get",
        "getitem",
        "iter",
        "to_dict",
        "set_record",
        "remove_record",
        "setitem",
        "delitem",
    ],
)
def test_nameless_headers_raise_malformed_record_error(
    read: typing.Callable[[SDBlock[str]], object],
) -> None:
    """Reading or editing records past a header without `<name>` raises.

    Args:
        read: what reads the records
    """
    block = SDBlock("t\n\nM  END\n> <A>\n1\n\n> no name\n1\n\n$$$$\n")
    with pytest.raises(MalformedRecordError, match="^SD data field header is"):
        read(block)
    assert issubclass(MalformedRecordError, ValueError)
    assert sdf_of(block) == "t\n\nM  END\n> <A>\n1\n\n> no name\n1\n\n$$$$\n"


def test_appended_records_follow_a_blank_line() -> None:
    """Appending after a last record without its blank line adds a record."""
    for method in ("append_record", "set_record"):
        block = SDBlock("t\nM  END\n> <A>\n1")
        getattr(block, method)("B", "2")
        assert block.records() == [("A", "1"), ("B", "2")]
        assert sdf_of(block) == "t\nM  END\n> <A>\n1\n\n> <B>\n2\n\n$$$$\n"
    constructed = SDBlock("t\nM  END\n> <A>\n1\n$$$$\n", {"B": "2"})
    assert constructed.records() == [("A", "1"), ("B", "2")]


def test_record_edits_take_keywords() -> None:
    """Record edits name their arguments `name` and `value`."""
    block = SDBlock("t\nM  END\n$$$$\n")
    block.append_record(name="A", value="1")
    block.set_record(name="B", value="2")
    block.remove_record(name="B")
    assert str(block) == "t\nM  END\n> <A>\n1\n\n$$$$\n"


def test_replace_records_can_read_the_block_itself() -> None:
    """Records to replace with may come from, or change, the same block."""
    block = SDBlock("t\nM  END\n> <A>\n1\n\n> <B>\n2\n\n$$$$\n")
    block.replace_records(block.iter_records())
    assert block.records() == [("A", "1"), ("B", "2")]
    block.replace_records(
        (name, value) for name, value in block.iter_records() if name != "A"
    )
    assert block.records() == [("B", "2")]

    def retitling() -> typing.Iterator[tuple[str, str]]:
        block.title = "changed"
        yield "C", "3"

    block.replace_records(retitling())
    assert (block.title, block.records()) == ("changed", [("C", "3")])


def test_constructor_parses_one_block() -> None:
    """SDBlock(sdf) parses a molblock or one record and appends records."""
    molblock = "title\n  prog\n\nM  END\n"
    block = SDBlock(molblock)
    assert (block.title, block.molblock, block.records()) == ("title", molblock, [])
    assert sdf_of(block) == molblock + "$$$$\n"

    with_data = SDBlock(molblock + "> <A>\n1\n\n$$$$\n\n", {"B": "2"})
    assert with_data.molblock == molblock
    assert with_data.records() == [("A", "1"), ("B", "2")]
    pairs = SDBlock(molblock.encode(), [(b"A", b"1"), (b"A", b"2")])
    assert pairs.records() == [(b"A", b"1"), (b"A", b"2")]

    assert SDBlock(molblock.replace("\n", "\r\n")).molblock == molblock
    assert SDBlock(molblock.encode().replace(b"\n", b"\r\n")).mdl.endswith(
        b"M  END\r\n"
    )
    empty = SDBlock("")
    assert (empty.title, empty.mdl, empty.records()) == ("", "", [])

    with pytest.raises(ValueError, match="one block"):
        SDBlock(molblock + "$$$$\n" + molblock)
    with pytest.raises(TypeError, match="str or bytes"):
        SDBlock(1)  # type: ignore[type-var]
    with pytest.raises(TypeError, match="record names must be str"):
        SDBlock(molblock, [(b"A", b"1")])  # type: ignore[type-var]


def test_blocks_compare_by_content_and_are_unhashable() -> None:
    """== compares kind and content, and blocks are mutable so unhashable."""
    sdf = "t\n  prog\n\nM  END\n> <A>\n1\n\n$$$$\n"
    block = SDBlock(sdf)
    assert block == SDBlock(sdf) and not block != SDBlock(sdf)
    assert block != SDBlock(sdf.encode())
    assert SDBlock(sdf.encode()) == SDBlock(sdf.encode())
    assert block != sdf and block != object()

    edited = SDBlock(sdf)
    edited.set_record("A", "2")
    assert edited != block
    edited.set_record("A", "1")
    assert edited == block

    # The same text parses differently when SD data follows an MDL
    # without M  END, so equal text alone does not make equal blocks.
    appended = SDBlock("t\n  prog\n")
    appended.append_record("A", "1")
    reparsed = SDBlock(sdf_of(appended))
    assert sdf_of(reparsed) == sdf_of(appended)
    assert (reparsed.records(), appended.records()) == ([], [("A", "1")])
    assert reparsed != appended

    with pytest.raises(TypeError, match="unhashable"):
        hash(block)


def _state_edge_blocks() -> list[SDBlock[str]]:
    """Build blocks that SDBlock(str(block)) does not rebuild exactly.

    Returns:
        Edited and unusual blocks
    """
    no_end = SDBlock("t\n  prog\n")
    no_end.append_record("A", "1")
    return [no_end, SDBlock("")]


@pytest.mark.parametrize("protocol", range(pickle.HIGHEST_PROTOCOL + 1))
def test_blocks_pickle_their_exact_state(protocol: int) -> None:
    """Pickling keeps the kind and the state, even where text would not.

    Args:
        protocol: pickle protocol
    """
    sdf = "t\r\n  prog\r\n\r\nM  END\r\n> <A>\r\n1\r\n\r\n$$$$\r\n"
    blocks: list[SDBlock[str] | SDBlock[bytes]] = [
        SDBlock(sdf),
        SDBlock(sdf.encode()),
        SDBlock("caf\u00e9\nM  END\n"),
        SDBlock(b"\xff\nM  END\n> <\xfe>\n\xfd\n\n"),
        *_state_edge_blocks(),
    ]
    for block in blocks:
        restored = pickle.loads(pickle.dumps(block, protocol))
        assert restored == block and restored is not block
        assert type(restored.title) is type(block.title)
        assert sdf_of(restored) == sdf_of(block)
        assert restored.records() == block.records()
    assert pickle.loads(pickle.dumps(blocks, protocol)) == blocks


def test_block_state_is_validated() -> None:
    """Unpickling rejects state that no block has instead of panicking."""
    from_state = SDBlock._from_state
    assert from_state(b"t\n$$$$\n", 2, 2, True) == SDBlock("t\n")
    for state in [
        (b"", 0, 0, True),
        (b"t\n", 2, 2, True),
        (b"t\n$$$$", 2, 2, False),
        (b"t\n$$$$\n", 1, 1, False),
        (b"t\nx\n$$$$\n", 2, 3, False),
        (b"t\nx\n$$$$\n", 4, 4, False),
        (b"t\nM  END\n$$$$\n", 9, 2, False),
        (b"t\nM  END\n$$$$\n", 2, 99, False),
    ]:
        with pytest.raises(ValueError, match="invalid SDBlock state"):
            from_state(*state)
    assert from_state(b"\xff\nM  END\n$$$$\n", 2, 9, False).title == b"\xff"
    # Text blocks from str paths decode on access, so may hold any bytes.
    with pytest.raises(UnicodeDecodeError):
        from_state(b"\xff\nM  END\n$$$$\n", 2, 9, True).title


@pytest.mark.parametrize("copier", [copy.copy, copy.deepcopy])
def test_block_copies_are_independent(
    copier: typing.Callable[[SDBlock[str]], SDBlock[str]],
) -> None:
    """Copies, shallow or deep, are equal blocks that are edited separately.

    Args:
        copier: copy.copy or copy.deepcopy
    """
    for block in [SDBlock("t\nM  END\n> <A>\n1\n\n"), *_state_edge_blocks()]:
        duplicate = copier(block)
        assert duplicate == block and duplicate is not block
        duplicate.append_record("B", "2")
        duplicate.title = "other"
        assert duplicate != block and "B" not in block
    byte_block = SDBlock(b"t\nM  END\n")
    assert copier(byte_block) == byte_block  # type: ignore[arg-type]
    shared = SDBlock("t\n")
    first, second = copy.deepcopy([shared, shared])
    assert first is second and first is not shared


V3000 = b"""\
t
  x

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 9 9 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0 0 0 0
M  V30 END ATOM
M  V30 END CTAB
M  END
> <A>
1

$$$$
"""


EDITS: list[typing.Callable[[SDBlock[bytes]], object]] = [
    lambda block: block.__setitem__(b"A", b"2"),
    lambda block: block.__delitem__(b"A"),
    lambda block: setattr(block, "title", b"u"),
    lambda block: setattr(block, "molblock", b"u\n\n\nM  END\n"),
    lambda block: block.append_record(b"B", b"2"),
    lambda block: block.set_record(b"B", b"2"),
    lambda block: block.remove_record(b"A"),
    lambda block: block.replace_records({b"B": b"2"}),
    SDBlock.fix_counts,
    SDBlock.renumber_indices,
]


@pytest.mark.parametrize("edit", EDITS)
def test_exported_blocks_refuse_edits(
    edit: typing.Callable[[SDBlock[bytes]], object],
) -> None:
    """A bytes block exports its SDF until every view of it is released.

    Args:
        edit: an edit to the block
    """
    block = SDBlock(V3000)
    view = memoryview(block)
    assert view.readonly and view == sdf_of(block) == bytes(block) == V3000
    duplicate = copy.copy(block)
    edit(duplicate)
    part = view[1:]
    view.release()
    with pytest.raises(BufferError, match="exports"):
        edit(block)
    assert bytes(part) == V3000[1:]
    part.release()
    edit(block)
    assert block == duplicate
    out = io.BytesIO()
    assert out.write(block) == len(sdf_of(block))


@pytest.mark.parametrize("edit", EDITS)
def test_repeated_lookups_follow_edits(
    edit: typing.Callable[[SDBlock[bytes]], object],
) -> None:
    """Lookups after the first use an index of the records, which edits drop.

    Args:
        edit: an edit to the block
    """
    records = b"> <A>\n1\n\n> <B>\nx\r\ny\r\n\n> <A>\n3\n"
    block = SDBlock(V3000.replace(b"> <A>\n1\n", records))
    names = [b"A", b"B", b"C"]

    def lookups() -> list[object]:
        found = [block.get(name) for name in names]
        return [len(block), *found, *(name in block for name in names)]

    def walked() -> list[object]:
        values = block.to_dict()
        found = [values.get(name) for name in names]
        return [len(block.names()), *found, *(name in values for name in names)]

    expected = [3, b"3", b"x\ny", None, True, True, False]
    assert lookups() == lookups() == walked() == expected
    edit(block)
    assert lookups() == lookups() == walked()


def test_str_blocks_are_not_bytes_like() -> None:
    """A str block has no buffer, so binary files and bytes() refuse it."""
    block = SDBlock("t\nM  END\n")
    with pytest.raises(TypeError, match="not bytes-like"):
        memoryview(block)
    with pytest.raises(TypeError, match="not bytes-like"):
        io.BytesIO().write(block)


def test_block_repr_shows_the_title() -> None:
    """The repr shows the title in the block's type."""
    assert repr(SDBlock(" aspirin \nM  END\n")) == "<SDBlock title=' aspirin'>"
    assert repr(SDBlock(b"aspirin\nM  END\n")) == "<SDBlock title=b'aspirin'>"


def test_str_and_bytes_are_the_sdf_of_their_kind() -> None:
    """str() is a str block's SDF; a bytes block's is bytes(), as for bytes."""
    text = SDBlock("t\nM  END\n")
    raw = SDBlock(b"t\nM  END\n")
    assert str(text) == f"{text}" == "t\nM  END\n$$$$\n"
    assert bytes(raw) == b"t\nM  END\n$$$$\n"
    assert str(raw) == repr(raw)
    with pytest.raises(TypeError, match="str"):
        bytes(text)


def test_molblock_setter_replaces_the_structure_only() -> None:
    """Setting molblock swaps title and CTAB and keeps every record line."""
    text = "old\nsrc\nM  END\n> <A>  (1)\n1\n\nnote\n\n$$$$\n"
    block = next(SDBlock.from_lines(io.StringIO(text)))

    block.molblock = "new\r\n  prog\r\n\r\nM  END\r\n\r\n"
    assert block.title == "new"
    assert block.molblock == "new\n  prog\n\nM  END\n"
    assert str(block) == "new\n  prog\n\nM  END\n> <A>  (1)\n1\n\nnote\n\n$$$$\n"
    for bad, message in [
        ("t\nx\n", "no M  END"),
        ("t\nM  END\n> <A>\n1\n", "after its M  END"),
        ("t\n$$$$\nM  END\n", r"\$\$\$\$ line"),
        ("$$$$\nM  END\n", r"\$\$\$\$ line"),
    ]:
        with pytest.raises(ValueError, match=message):
            block.molblock = bad
    with pytest.raises(TypeError, match="molblocks must be str"):
        block.molblock = b"t\nM  END\n"  # type: ignore[assignment]
    assert block.molblock == "new\n  prog\n\nM  END\n"
    block.molblock = "last\r\nM  END"
    assert block.molblock == "last\nM  END\n"

    byte_block = next(SDBlock.from_lines([b"t\n", b"M  END\n", b"> <A>\n", b"1\n"]))
    byte_block.molblock = b"b\r\nM  END"
    assert byte_block.molblock == b"b\r\nM  END\r\n"
    byte_block.molblock = b"b\nM  END\r"
    assert byte_block.molblock == b"b\nM  END\r\n"
    byte_block.molblock = b"b\r\nM  END\r\n"
    assert byte_block.molblock == b"b\r\nM  END\r\n"
    assert byte_block.title == b"b"
    assert byte_block.to_dict() == {b"A": b"1"}


def test_title_setter_replaces_only_the_title() -> None:
    """Setting the title keeps the CTAB and records and rejects delimiters.

    A title is written as given, so one starting with $$$$ would read back
    as the end of the block; one starting with a space does not.
    """
    block = next(SDBlock.from_lines(["old\n", "src\n", "M  END\n", "> <A>\n", "1\n"]))

    block.title = "a much longer title"
    assert block.title == "a much longer title"
    assert block.molblock == "a much longer title\nsrc\nM  END\n"
    assert block.records() == [("A", "1")]
    block.title = " $$$$ "
    assert block.title == " $$$$"
    assert sdf_of(block) == " $$$$ \nsrc\nM  END\n> <A>\n1\n$$$$\n"
    assert [b.title for b in SDBlock.from_sdf(sdf_of(block))] == [" $$$$"]
    for bad, message in [
        ("two\nlines", "single line"),
        ("cr\r", "single line"),
        ("$$$$", "cannot start with"),
        ("$$$$ x", "cannot start with"),
    ]:
        with pytest.raises(ValueError, match=message):
            block.title = bad
    with pytest.raises(TypeError, match="titles must be str"):
        block.title = b"bytes"  # type: ignore[assignment]
    assert block.title == " $$$$"


def test_only_ascii_whitespace_strips() -> None:
    """Str and bytes blocks strip only ASCII whitespace, as RDKit does.

    Titles are right-stripped and values kept as read. A value line holding
    only Unicode whitespace is not blank, unlike str.strip() in pure-Python
    sendoff 0.1.x, so it does not end the record.
    """
    text_block = next(
        SDBlock.from_lines(
            [
                " \u2003Title\u3000 \n",
                "\n",
                "M  END\n",
                "> <Record>\u2003\n",
                "\u2003Text\x1f\n",
                "\n",
                "> <Blank>\n",
                "\u3000\n",
                "skipped\n",
                "\n",
                "$$$$\n",
            ]
        )
    )
    byte_block = next(
        SDBlock.from_lines(
            [
                b"\xc2\xa0Title\xc2\xa0\n",
                b"\n",
                b"M  END\n",
                b"> <Record>\n",
                b"\xc2\xa0Bytes\x1f\n",
                b"\n",
                b"$$$$\n",
            ]
        )
    )

    assert text_block.title == " \u2003Title\u3000"
    assert text_block.records() == [
        ("Record", "\u2003Text\x1f"),
        ("Blank", "\u3000\nskipped"),
    ]
    assert byte_block.title == b"\xc2\xa0Title\xc2\xa0"
    assert byte_block.records() == [(b"Record", b"\xc2\xa0Bytes\x1f")]
