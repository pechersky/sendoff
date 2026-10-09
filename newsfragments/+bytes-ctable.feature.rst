``ctable``, ``renumber_indices``, and ``fix_counts`` work on bytes blocks.
Their CTable lines are decoded as Latin-1, so ``line.encode("latin-1")``
gives the bytes back and edits keep every other byte.
