Shell helpers run as ``python -m sendoff.cli.<command>``, reading a file or
stdin in bytes mode: ``count`` prints the number of blocks, ``get`` prints a
record's value or the title of each block, one line per block, ``set`` sets a
record and ``title`` the title from ``--value`` or from ``--from`` lines, one
per block, and ``title --field`` copies a record into the title.
