Shell helpers run as ``python -m sendoff.cli.<command>``, reading a file
(``FILE`` or ``-i``) or stdin in bytes mode and writing ``-o FILE`` or stdout:
``count`` prints the number of blocks, ``get`` prints a record's value or the
title of each block, one line per block, ``set`` sets a record and ``title``
the title from ``--value`` or from ``--from`` lines, one per block, ``title
--field`` copies a record into the title, and ``filter`` keeps the blocks that
have a record, a record with a value, or, with ``--unique``, a value not seen
before. ``filter --unique NAME`` replaces ``scripts/sdf-unique-record``.
