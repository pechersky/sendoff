Records the block writes keep the file's header style. ``set_record`` keeps
the header line of the record it sets, such as ``>  <A>  (1)``, rather than
rewriting it as ``> <A>``. ``append_record``, and ``set_record`` of a new name,
copy the first named header's line around the name, so after RDKit's
``>  <A>  (1)`` they write ``>  <C>  (1)``; a block without one still gets
``> <C>``.
