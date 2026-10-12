"""Print one record value, or the title, of each block."""

from __future__ import annotations

import os
import sys
from collections.abc import Sequence

from sendoff.cli._common import STDIN, Blocks, failures, parser, run


def main(argv: Sequence[str] | None = None) -> int:
    """Print a line per block, so several runs line up with ``paste``.

    Args:
        argv: arguments, or None for ``sys.argv``

    Returns:
        Exit status
    """
    p = parser(
        "get",
        "Print a record's value, or with -t the title, of each block, one line "
        "per block. Values are raw bytes; the last record wins for duplicates.",
    )
    p.usage = "%(prog)s [-h] [-z] [--missing STR] (NAME | -t) [FILE]"
    p.add_argument("name", nargs="?", metavar="NAME", help="record name")
    p.add_argument("file", nargs="?", metavar="FILE", help="SD file (default: stdin)")
    p.add_argument("-t", "--title", action="store_true", help="print titles")
    p.add_argument("-z", "--null", action="store_true", help="end lines with NUL")
    p.add_argument(
        "--missing",
        default="",
        metavar="STR",
        help="print for blocks without the record (default: empty)",
    )
    args = p.parse_intermixed_args(argv)
    if args.title:
        if args.file is not None:
            p.error("-t takes no NAME")
        path = args.name
    elif args.name is None:
        p.error("NAME or -t is required")
    else:
        path = args.file
    blocks = Blocks(STDIN if path is None else path)
    end = b"\0" if args.null else b"\n"
    out = sys.stdout.buffer
    with failures(p, blocks):
        if args.title:
            out.writelines(block.title + end for block in blocks)
        else:
            name = os.fsencode(args.name)
            missing = os.fsencode(args.missing)
            out.writelines(block.get(name, missing) + end for block in blocks)
    return 0


if __name__ == "__main__":
    run(main)
