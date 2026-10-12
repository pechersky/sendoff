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
       block["ID"]
       block.get("ID")
       block.to_dict()
       block.records()
       block.iter_records()
       "ID" in block, len(block), list(block)
       "V3000" in block.mdl
       mol = Chem.MolFromMolBlock(block.molblock)

   Counter(name for block in parse_sdf("compounds.sdf") for name in block)

   def tag(block):
       block["Source"] = "vendor"
       block.remove_record("Internal note")
       return block

   write_sdf("out.sdf", (tag(block) for block in parse_sdf("in.sdf")))

   block = SDBlock(Chem.MolToMolBlock(mol), {"ID": "X-1"})
   block.to_sdf()
   SDBlock.from_sdf(sdf_text)

``parse_sdf`` memory-maps local files, so write to a new path: truncating a
file while it is parsed crashes the process.

New headers are ``> <name>``; ``header="padded"`` writes ``>  <name>``, as
RDKit does, and ``header="copy"`` follows the block's first header. Names are
not validated.

Non-UTF-8 files
---------------

str paths decode UTF-8. Open other encodings yourself, or stay in bytes:

.. code-block:: python

   with open("vendor.sdf", encoding="gbk") as f:
       names = {block["CAS"]: block["Product name"] for block in parse_sdf(f)}

   in_stock = (
       block
       for block in parse_sdf(os.fsencode("vendor.sdf"))
       if block.get(b"Availability") == b"In-stock"
   )
   write_sdf("in-stock.sdf", in_stock)

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
