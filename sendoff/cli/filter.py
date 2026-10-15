"""Keep the blocks that have a record, a record with a value, or a new value."""

from __future__ import annotations

import os
from collections.abc import Iterable, Iterator, Sequence

from sendoff import SDBlock, write_sdf
from sendoff.cli._common import (
    Blocks,
    add_file,
    failures,
    input_path,
    output,
    parser,
    run,
)


def unique(
    blocks: Iterable[SDBlock[bytes]], name: bytes, keep: bool
) -> Iterator[SDBlock[bytes]]:
    """Yield the first block with each value of a record, and those without it.

    Args:
        blocks: the blocks
        name: the record name
        keep: False to yield the other blocks instead

    Yields:
        The kept blocks
    """
    seen: set[bytes] = set()
    for block in blocks:
        value = block.get(name)
        if value is None:
            first = True
        else:
            first = value not in seen
            seen.add(value)
        if first is keep:
            yield block


def main(argv: Sequence[str] | None = None) -> int:
    """Write the matching blocks unchanged.

    Args:
        argv: arguments, or None for ``sys.argv``

    Returns:
        Exit status
    """
    p = parser(
        "filter",
        "Write the blocks that have a record, a record with a value, or a "
        "value not seen before, unchanged. The last record wins for "
        "duplicates.",
    )
    add_file(p)
    tests = p.add_mutually_exclusive_group(required=True)
    tests.add_argument("--has", metavar="NAME", help="keep blocks with the record")
    tests.add_argument(
        "--eq",
        nargs=2,
        metavar=("NAME", "VALUE"),
        help="keep blocks whose record is VALUE",
    )
    tests.add_argument(
        "--unique",
        metavar="NAME",
        help="keep the first block with each value, and blocks without NAME",
    )
    p.add_argument("-v", "--invert", action="store_true", help="keep the others")
    args = p.parse_intermixed_args(argv)
    path = input_path(p, args.file, args.input)
    blocks = Blocks(path)
    out = output(p, args.output, path)
    keep = not args.invert
    kept: Iterator[SDBlock[bytes]]
    if args.has is not None:
        name = os.fsencode(args.has)
        kept = (block for block in blocks if (name in block) is keep)
    elif args.unique is not None:
        kept = unique(blocks, os.fsencode(args.unique), keep)
    else:
        name, value = map(os.fsencode, args.eq)
        kept = (block for block in blocks if (block.get(name) == value) is keep)
    with failures(p, blocks):
        write_sdf(out, kept)
    return 0


if __name__ == "__main__":
    run(main)
