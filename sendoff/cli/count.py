"""Count the blocks of SD files."""

from __future__ import annotations

import os
import sys
from collections.abc import Sequence

from sendoff.cli._common import STDIN, Blocks, failures, parser, run


def main(argv: Sequence[str] | None = None) -> int:
    """Print the block count of each file, then a total for several.

    Args:
        argv: arguments, or None for ``sys.argv``

    Returns:
        Exit status
    """
    p = parser("count", "Count blocks. Several files print COUNT<TAB>FILE lines.")
    p.add_argument(
        "files", nargs="*", default=[STDIN], metavar="FILE", help="default: stdin"
    )
    args = p.parse_intermixed_args(argv)
    out = sys.stdout.buffer
    total = 0
    for path in args.files:
        blocks = Blocks(path)
        with failures(p, blocks):
            count = sum(1 for _ in blocks)
        total += count
        if len(args.files) == 1:
            out.write(b"%d\n" % count)
        else:
            out.write(b"%d\t%s\n" % (count, os.fsencode(path)))
    if len(args.files) > 1:
        out.write(b"%d\ttotal\n" % total)
    return 0


if __name__ == "__main__":
    run(main)
