``CTable.renumber_indices`` now updates ``num_atoms`` and ``num_bonds`` to the
renumbered counts, which used to keep the counts the table was read with.
The rewritten ``COUNTS`` line keeps the line ending it had, so a table built
from lines that end in ``\n`` no longer gets one line without it.
Renumbering keeps what follows a bond's atoms, such as ``CFG`` stereo and the
line ending, which it used to drop. It now raises ``NotImplementedError``
rather than silently corrupt a table: when it would change indices that
other data refers to (atom ``ATTCHORD``, bond ``ENDPTS``, ``SGROUP``,
``COLLECTION``, ``OBJ3D``, and ``LINKNODE`` lines), which it does not remap,
and for ``-`` continued atom or bond lines and CTABs nested in ``RGROUP`` or
``TEMPLATE`` blocks.
