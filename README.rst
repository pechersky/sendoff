sendoff
=======

.. image:: https://results.pre-commit.ci/badge/github/pechersky/sendoff/main.svg
   :target: https://results.pre-commit.ci/latest/github/pechersky/sendoff/main
   :alt: pre-commit.ci status

.. image:: https://github.com/pechersky/sendoff/actions/workflows/tox.yml/badge.svg
   :target: https://github.com/pechersky/sendoff/actions/workflows/tox.yml
   :alt: Tox status

The minimal SDF metadata parser. It reads and edits titles and SD records
without building molecules, so it is fast, and leaves the structure and other
records as written, so it handles files that chemistry toolkits reject or
rewrite.

Usage
-----

.. code-block:: python

   from collections import Counter
   from rdkit import Chem
   from sendoff import SDBlock, parse_sdf, write_sdf

   for block in parse_sdf("compounds.sdf"):
       block.title
       block["ID"]  # the last value named "ID"; KeyError if none
       block.get("ID")  # or None: blocks need not share record names
       block.to_dict()  # every record, later duplicates winning
       block.records()  # [(name, value), ...] in file order
       block.iter_records()  # the same, lazily
       "ID" in block, len(block), list(block)  # names only
       "V3000" in block.mdl  # the CTAB after the title, through M  END
       mol = Chem.MolFromMolBlock(block.molblock)

   Counter(name for block in parse_sdf("compounds.sdf") for name in block)

   def tag(block):
       block["Source"] = "vendor"  # replace in place, or append
       block.remove_record("Internal note")
       return block

   # Faster than a block.write loop. Write to a new path: parse_sdf
   # memory-maps local files, and truncating one mid-parse crashes.
   write_sdf("out.sdf", (tag(block) for block in parse_sdf("in.sdf")))

   block = SDBlock(Chem.MolToMolBlock(mol), {"ID": "X-1"})
   block.to_sdf()  # SDF text ending in $$$$
   SDBlock.from_sdf(sdf_text)  # blocks from SDF in memory

New headers are ``> <name>``; ``header="padded"`` writes ``>  <name>``, as
RDKit does, and ``header="copy"`` follows the block's first header. Names are
not checked: the CTfile format wants a letter, then letters, digits, or
underscores.

Non-UTF-8 files
---------------

str paths decode UTF-8. For other encodings, decode yourself:

.. code-block:: python

   with open("vendor.sdf", encoding="gbk") as f:  # read through Python
       write_sdf("vendor-utf8.sdf", parse_sdf(f))  # str blocks write UTF-8

   for block in parse_sdf(os.fsencode("vendor.sdf")):  # bytes, as written
       block[b"Product name"].decode("gbk", "replace")

Shell pipelines
---------------

.. code-block:: bash

   cat my.sdf | python -c '
   from sendoff import parse_sdf
   for block in parse_sdf("/dev/stdin"):
       print(block.get("Supplier"))
   ' | sort | uniq -c

Write SDF to a pipe with ``write_sdf(sys.stdout, blocks)``; for bytes, read
``b"/dev/stdin"`` and write ``sys.stdout.buffer``.
