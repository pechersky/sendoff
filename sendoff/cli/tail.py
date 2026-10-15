"""Write the last blocks of an SD file, or from block K on."""

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


def take(blocks: Blocks, count: int, from_start: bool) -> Iterator[SDBlock[bytes]]:
    """Yield the last `count` blocks, or from 1-based block `count` on.

    Args:
        blocks: the input
        count: blocks to keep from the end, or 1-based start block
        from_start: whether `-n` had a leading ``+``

    Returns:
        The kept blocks
    """
    if from_start:
        return islice(blocks, max(count - 1, 0), None)
    return iter(deque(blocks, maxlen=count))


def main(argv: Sequence[str] | None = None) -> int:
    """Write the last blocks of an SD file, or from block K on.

    Args:
        argv: arguments, or None for ``sys.argv``

    Returns:
        Exit status
    """
    p = parser(
        "tail",
        "Write the last N blocks (default: 10), or with -n +K starting with "
        "1-based block K, unchanged.",
    )
    add_file(p)
    p.add_argument(
        "-n",
        "--lines",
        default="10",
        metavar="[+]N",
        help="keep the last N blocks, or with +N from block N on (default: 10)",
    )
    args, count, from_start = slice_spec(p, argv, sign="+")
    path = input_path(p, args.file, args.input)
    blocks = Blocks(path)
    out = output(p, args.output, path)
    with failures(p, blocks):
        write_sdf(out, take(blocks, count, from_start))
    return 0


if __name__ == "__main__":
    run(main)
