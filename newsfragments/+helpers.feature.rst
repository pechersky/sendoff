New ``SDBlock`` helpers for common tasks:

- ``len(block)``, ``name in block``, and ``block.names()`` read record names
  without building values.
- ``get(name, default=None)`` and ``to_dict()`` read records without
  ``dict(block.records())``, and ``iter_records()`` iterates them lazily.
- ``molblock`` returns the title and CTAB for ``Chem.MolFromMolBlock``.
  Assigning it swaps the structure and keeps the SD data.
- ``title`` is assignable, and ``set_record``, ``remove_record``, and
  ``replace_records`` edit the SD data.
- ``SDBlock(sdf, records=None)`` parses one block, such as an RDKit molblock,
  and ``SDBlock.from_sdf(text)`` parses SDF held in memory. ``to_sdf()``
  returns a block as SDF text.
- ``sendoff.sdblock.write_sdf(dest, blocks)`` writes blocks to a path or a
  file object and returns how many it wrote.
- ``fix_counts()`` on ``SDBlock`` and ``CTable`` sets the V3000 ``COUNTS``
  line to the number of atom and bond lines.
- ``parse_sdf`` of a ``bytes`` path reads files that are not UTF-8 and gives
  ``bytes`` titles and values.
- The package is typed: ``SDBlock`` is ``Generic[AnyStr]``, so ``parse_sdf`` of
  a ``str`` path gives ``SDBlock[str]`` and of a ``bytes`` path gives
  ``SDBlock[bytes]``.
