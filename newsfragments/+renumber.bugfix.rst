``CTable.renumber_indices`` now updates ``num_atoms`` and ``num_bonds`` to the
renumbered counts, which used to keep the counts the table was read with.
The rewritten ``COUNTS`` line keeps the line ending it had, so a table built
from lines that end in ``\n`` no longer gets one line without it.
