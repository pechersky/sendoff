Reading or editing records past a header line with no ``<name>`` raises
``sendoff.MalformedRecordError``, a ``ValueError``, instead of a bare
``ValueError``.
