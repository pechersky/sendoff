"""Keep the blocks that have a record, or a record with a value."""

from __future__ import annotations

import os
import sys
from collections.abc import Sequence

from sendoff import write_sdf
from sendoff.cli._common import Blocks, add_file, failures, parser, run


def main(argv: Sequence[str] | None = None) -> int:
    """Write the matching blocks unchanged.

    Args:
        argv: arguments, or None for ``sys.argv``

    Returns:
        Exit status
    """
    p = parser(
        "filter",
        "Write the blocks that have a record, or a record with a value, "
        "unchanged. The last record wins for duplicates.",
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
    p.add_argument("-v", "--invert", action="store_true", help="keep the others")
    args = p.parse_intermixed_args(argv)
    blocks = Blocks(args.file)
    keep = not args.invert
    if args.has is not None:
        name = os.fsencode(args.has)
        kept = (block for block in blocks if (name in block) is keep)
    else:
        name, value = map(os.fsencode, args.eq)
        kept = (block for block in blocks if (block.get(name) == value) is keep)
    with failures(p, blocks):
        write_sdf(sys.stdout.buffer, kept)
    return 0


if __name__ == "__main__":
    run(main)
