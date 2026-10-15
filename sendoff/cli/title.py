"""Replace the title of each block."""

from __future__ import annotations

import os
from collections.abc import Iterator, Sequence

from sendoff import SDBlock, write_sdf
from sendoff.cli._common import (
    Blocks,
    add_file,
    add_sources,
    failures,
    input_path,
    output,
    pairs,
    parser,
    run,
)


def main(argv: Sequence[str] | None = None) -> int:
    """Write the blocks with new titles.

    Args:
        argv: arguments, or None for ``sys.argv``

    Returns:
        Exit status
    """
    p = parser("title", "Replace the title of each block.")
    add_file(p)
    add_sources(p, "title").add_argument(
        "--field",
        metavar="NAME",
        help="copy a record's value, as get NAME | title --from - would",
    )
    args = p.parse_intermixed_args(argv)
    path = input_path(p, args.file, args.input)
    blocks = Blocks(path)
    out = output(p, args.output, path, args.source)

    def titled() -> Iterator[SDBlock[bytes]]:
        if args.field is None:
            titles = pairs(p, args, blocks)
        else:
            name = os.fsencode(args.field)
            titles = ((block, block.get(name, b"")) for block in blocks)
        for block, title in titles:
            if title or not args.skip_empty:
                block.title = title
            yield block

    with failures(p, blocks):
        write_sdf(out, titled())
    return 0


if __name__ == "__main__":
    run(main)
