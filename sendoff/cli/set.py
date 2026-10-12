"""Set a record on each block."""

from __future__ import annotations

import os
import sys
from collections.abc import Iterator, Sequence

from sendoff import SDBlock, write_sdf
from sendoff.cli._common import (
    Blocks,
    add_file,
    add_sources,
    failures,
    pairs,
    parser,
    run,
)


def main(argv: Sequence[str] | None = None) -> int:
    """Write the blocks with a record set.

    Args:
        argv: arguments, or None for ``sys.argv``

    Returns:
        Exit status
    """
    p = parser(
        "set",
        "Set a record on each block: replace the first record with that name, "
        "keeping its header line and place, and remove the others; or append it.",
    )
    p.add_argument("name", metavar="NAME", help="record name")
    add_file(p)
    add_sources(p, "value")
    p.add_argument(
        "--append",
        action="store_true",
        help="always append, keeping records with that name",
    )
    p.add_argument(
        "--header",
        choices=("plain", "padded", "copy"),
        default="plain",
        help="new header lines: plain '> <NAME>', padded '>  <NAME>', or copy "
        "the block's first header line around the name (default: plain)",
    )
    args = p.parse_intermixed_args(argv)
    name = os.fsencode(args.name)
    blocks = Blocks(args.file)
    edit = SDBlock.append_record if args.append else SDBlock.set_record

    def edited() -> Iterator[SDBlock[bytes]]:
        for block, value in pairs(p, args, blocks):
            if value or not args.skip_empty:
                edit(block, name, value, header=args.header)
            yield block

    with failures(p, blocks):
        write_sdf(sys.stdout.buffer, edited())
    return 0


if __name__ == "__main__":
    run(main)
