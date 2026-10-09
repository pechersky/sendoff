``block.title`` strips only trailing whitespace, and ``molblock``,
``to_sdf``, and ``write`` keep the title line as read. A title such as
`` $$$$`` used to be written stripped, as a ``$$$$`` line that read back as
the end of the block. Setting a ``title`` that starts with ``$$$$``, or a
``molblock`` whose title line does, raises ``ValueError``.
