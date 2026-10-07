CTable errors form one hierarchy under ``sendoff.ctable.CTableError``, a
``ValueError``. Malformed tables raise ``MalformedCTableError``, whose message
names the line, instead of ``IndexError``, ``KeyError``, or ``ValueError``.
``IndicesMismatchError``, ``IndicesOutOfOrderError``, and
``IndicesDuplicateError`` were plain ``Exception`` subclasses and are now
``CTableError`` subclasses, and every message is rewritten to name the line or
counts involved. Catch ``CTableError`` instead of ``IndexError``, ``KeyError``,
or ``ValueError``. Renumbering now rejects a bond line whose index is not an
integer.
