An unterminated final block (no ``$$$$``) is no longer dropped. A block
without an ``M  END`` line no longer keeps its ``$$$$`` line in the CTAB, which
made ``write`` double the delimiter and round trips grow extra blocks.
Blocks read with ``\r\n`` line endings, as from a ``bytes`` path, now write
them back unchanged. ``write`` used to end the title and ``$$$$`` lines with
``\n``, mixing line endings. Lines a block adds or ends, such as a file's
unterminated last line, now end like its title line.
