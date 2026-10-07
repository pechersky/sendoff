Validating indices, renumbering, and ``fix_counts`` find atom and bond lines by
the root CTAB's ``BEGIN`` and ``END`` lines, so they no longer read a CTAB
nested in an ``RGROUP`` or ``TEMPLATE`` block as the root's. They raise
``ValueError`` when an ``ATOM`` or ``BOND`` section has no ``END`` line or
appears twice.
