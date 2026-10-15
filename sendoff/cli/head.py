"""Write the first or all-but-last blocks of an SD file."""

from __future__ import annotations

from collections import deque
from collections.abc import Iterator, Sequence
from itertools import islice

from sendoff import SDBlock, write_sdf
from sendoff.cli._common import (
    Blocks,
    add_file,
    failures,
    input_path,
    output,
    parser,
    run,
    slice_spec,
)


def take(blocks: Blocks, count: int, from_end: bool) -> Iterator[SDBlock[bytes]]:
    """Yield the first `count` blocks, or all but the last `count`.

    Args:
        blocks: the input
        count: blocks to keep, or with `from_end` to drop from the end
        from_end: whether `-n` had a leading ``-``

    Yields:
        The kept blocks
    """
    if not from_end or count == 0:
        yield from islice(blocks, None if from_end else count)
        return
    window: deque[SDBlock[bytes]] = deque()
    for block in blocks:
        if len(window) == count:
            yield window.popleft()
        window.append(block)


def main(argv: Sequence[str] | None = None) -> int:
    """Write the first blocks of an SD file.

    Args:
        argv: arguments, or None for ``sys.argv``

    Returns:
        Exit status
    """
    p = parser(
        "head",
        "Write the first N blocks (default: 10), or with -n -K all but the "
        "last K blocks, unchanged.",
    )
    add_file(p)
    p.add_argument(
        "-n",
        "--lines",
        default="10",
        metavar="[-]N",
        help="keep the first N blocks, or with -N all but the last N (default: 10)",
    )
    args, count, from_end = slice_spec(p, argv)
    path = input_path(p, args.file, args.input)
    blocks = Blocks(path)
    out = output(p, args.output, path)
    with failures(p, blocks):
        write_sdf(out, take(blocks, count, from_end))
    return 0


if __name__ == "__main__":
    run(main)
