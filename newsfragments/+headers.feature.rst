Records the block writes keep the file's header style. ``set_record`` keeps
the header line of the record it sets, such as ``>  <A>  (1)``, rather than
rewriting it as ``> <A>``, and ``replace_records`` gives the k-th record with a
name the k-th old header line with that name. Other new headers copy the
first named header's line around the name, so after RDKit's ``>  <A>  (1)``
``append_record`` writes ``>  <C>  (1)``; a block without one still gets
``> <C>``.
