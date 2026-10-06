An unterminated final block (no ``$$$$``) is no longer dropped. A block
without an ``M  END`` line no longer keeps its ``$$$$`` line in the CTAB, which
made ``write`` double the delimiter and round trips grow extra blocks.
