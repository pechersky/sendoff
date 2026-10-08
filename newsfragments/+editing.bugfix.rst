``append_record`` and ``set_record`` add a separate record after a last record
without a blank line after it, as at the end of an unterminated file. The new
record used to run on as part of the last value.
``replace_records`` accepts records read from the same block, such as
``block.iter_records()``, which used to raise a ``PanicException``.
