"""Count the blocks of SD files."""

from __future__ import annotations

import os
from collections.abc import Sequence

from sendoff.cli._common import (
    STDIO,
    Blocks,
    add_io,
    failures,
    opened,
    output,
    parser,
    run,
)


def main(argv: Sequence[str] | None = None) -> int:
    """Print the block count of each file, then a total for several.

    Args:
        argv: arguments, or None for ``sys.argv``

    Returns:
        Exit status
    """
    p = parser("count", "Count blocks. Several files print COUNT<TAB>FILE lines.")
    p.add_argument("files", nargs="*", metavar="FILE", help="default: stdin")
    add_io(p, many=True)
    args = p.parse_intermixed_args(argv)
    paths = [*args.files, *(args.input or [])] or [STDIO]
    total = 0
    with opened(p, output(p, args.output, *paths)) as out:
        for path in paths:
            blocks = Blocks(path)
            with failures(p, blocks):
                count = sum(1 for _ in blocks)
            total += count
            if len(paths) == 1:
                out.write(b"%d\n" % count)
            else:
                out.write(b"%d\t%s\n" % (count, os.fsencode(path)))
        if len(paths) > 1:
            out.write(b"%d\ttotal\n" % total)
    return 0


if __name__ == "__main__":
    run(main)
