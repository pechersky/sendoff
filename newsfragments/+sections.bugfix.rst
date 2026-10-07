Validating indices, renumbering, and ``fix_counts`` find atom and bond lines by
the root CTAB's ``BEGIN`` and ``END`` lines, so they no longer read a CTAB
nested in an ``RGROUP`` or ``TEMPLATE`` block as the root's. They raise
``MalformedCTableError`` when an ``ATOM`` or ``BOND`` section has no ``END``
line or appears twice.
Validating indices and ``fix_counts`` read a line ending in ``-`` and the next
line as one atom or bond line, as V3000 continues lines, rather than as two,
and raise ``MalformedCTableError`` when a section's last line ends in ``-``.
