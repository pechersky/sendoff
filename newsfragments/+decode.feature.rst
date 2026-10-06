A ``UnicodeDecodeError`` from a block read through a ``str`` path now carries
a note naming the record or part being read, the block's title, its index in
``write_sdf``, and how to read files that are not UTF-8.
