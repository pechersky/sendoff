Records the block writes keep the file's header lines. ``set_record`` keeps
the header line of the record it sets, such as ``>  <A>  (1)``, rather than
rewriting it as ``> <A>``, and ``replace_records`` gives the k-th record with a
name the k-th old header line with that name. New header lines are still
``> <C>``; ``append_record``, ``set_record``, ``replace_records``, and
``SDBlock(sdf, records)`` take ``header="padded"`` for ``>  <C>``, or
``header="copy"`` to copy the first named header's line around the name, so
after RDKit's ``>  <A>  (1)`` they write ``>  <C>  (1)``.
