The ``SDBlock`` API changed. Migrate 0.1.x code as follows:

- ``str.join("", block.mdl)`` becomes ``block.mdl``, one ``str`` (or
  ``bytes``) of the CTAB after the title through ``M  END``, each line ending
  in ``\n``. Iterate lines with ``block.mdl.splitlines(keepends=True)`` and
  search with ``needle in block.mdl``. Code that iterated lines still runs but
  walks characters, so ``any(s in line for line in block.mdl)``,
  ``islice(block.mdl, ...)``, and ``str.join("\n", block.mdl)`` give wrong
  results without raising: check every ``.mdl`` use.
- ``block.title + "\n" + str.join("", block.mdl)`` becomes ``block.molblock``.
  ``mdl`` is read-only; assign ``block.molblock`` to replace the structure.
- ``block.metadata`` is gone. Read records with ``records()``, ``get``, or
  ``to_dict``; edit them with ``append_record``, ``set_record``,
  ``remove_record``, or ``replace_records``.
- ``records()`` returns a list; use ``iter_records()`` for ``next(...)``.
  ``dict(block.records()).get(name)`` becomes ``block.get(name)``,
  ``dict(block.records())`` becomes ``block.to_dict()``,
  ``len(list(block.records()))`` becomes ``len(block)``,
  ``name in dict(block.records())`` becomes ``name in block``, and
  ``[name for name, _ in block.records()]`` becomes ``block.names()``.
- ``SDBlock`` is no longer a dataclass. ``SDBlock(title, mdl, metadata)``
  becomes ``SDBlock(sdf, records)``, and ``SDBlock.from_block_lines(lines)``
  becomes ``SDBlock(str.join("", lines))``. ``parse_mdl`` and
  ``parse_metadata`` are removed. Blocks still compare by content and copy
  and pickle, but pickles from 0.1.x do not load, and ``repr`` shows only the
  title.
- ``SDBlock.from_lines`` no longer takes a ``str`` or ``bytes``, which it
  iterated character by character: ``from_lines(text.splitlines(True))``
  becomes ``SDBlock.from_sdf(text)``. ``from_lines`` also takes file objects.
- ``write(outh, with_newlines=...)`` becomes ``write(file)``, and
  ``append_record(record_name=..., value=...)`` becomes
  ``append_record(name=..., value=...)``, matching ``set_record``. Replace
  ``block.write(handle)`` loops with ``write_sdf(dest, blocks)``;
  ``block.to_sdf()`` returns the text.
- ``parse_sdf`` decodes ``str`` paths as UTF-8, not with the locale encoding.
  Pass a ``bytes`` path for other encodings to get ``bytes`` titles and values.
  Its ``sdfpath`` keyword is now ``source``.
- ``CTable`` loses ``title``, ``source``, ``comment``, ``counts``, ``format``,
  ``CTableFormat``, the ``parse_*`` static methods, and ``atomlines`` and
  ``bondlines``; ``lines`` is a list.
- ``block.num_atoms()`` and ``block.num_bonds()`` become properties,
  ``block.num_atoms`` and ``block.num_bonds``, as on ``CTable``. They read the
  counts line without building a ``CTable``.
- ``renumber_indices`` no longer rewrites the V3000 ``COUNTS`` line, so it no
  longer repairs counts after atom or bond lines are added or removed. Call
  ``fix_counts()`` for that.
- ``sendoff.sdblock.Pathy`` is gone; annotate paths as
  ``str | bytes | os.PathLike[str] | os.PathLike[bytes]``.
