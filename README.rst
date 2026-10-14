sendoff
=======

.. image:: https://results.pre-commit.ci/badge/github/pechersky/sendoff/main.svg
   :target: https://results.pre-commit.ci/latest/github/pechersky/sendoff/main
   :alt: pre-commit.ci status

.. image:: https://github.com/pechersky/sendoff/actions/workflows/tox.yml/badge.svg
   :target: https://github.com/pechersky/sendoff/actions/workflows/tox.yml
   :alt: Tests

The minimal SDF metadata parser.

sendoff reads and edits SDF titles and data records without building
molecules, and writes back everything it doesn't edit byte for byte, even
malformed chemistry. Tests check it against RDKit on edge cases such as
multiline values and ``$$$$`` titles.

In time per operation, it is 27–59× faster than RDKit, 6.4–23× faster than
OEChem, and 8.2–48× faster than the pure-Python sendoff 0.1.6 (geomeans by
dataset; see the `benchmarks
<https://github.com/pechersky/sendoff/blob/main/benchmark-results.md>`_).

Install
-------

.. code-block:: bash

   pip install sendoff

Wheels cover CPython 3.11 and later on Linux (glibc and musl), macOS, and
Windows, on x86_64 and arm64.

Usage
-----

.. code-block:: python

   from collections import Counter
   from rdkit import Chem
   from sendoff import SDBlock, parse_sdf, sdf_writer, write_sdf

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

   with sdf_writer("active.sdf") as active, sdf_writer("rest.sdf") as rest:
       for block in parse_sdf("in.sdf"):
           (active if block.get("Active") == "1" else rest).write(block)

   block = SDBlock(Chem.MolToMolBlock(mol), {"ID": "X-1"})
   str(block)
   SDBlock.from_sdf(sdf_text)

``parse_sdf`` memory-maps local files, so write to a new path: truncating a
file while it is parsed crashes the process. ``parse_sdf(path, borrow=True)``
also skips copying each block: blocks read the file in place until edited,
so the same holds for as long as they live.

New headers are ``> <name>``; ``header="padded"`` writes ``>  <name>``, as
RDKit does, and ``header="copy"`` follows the block's first header.

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

   python -m sendoff.cli.count my.sdf
   zcat chembl.sdf.gz | python -m sendoff.cli.count
   python -m sendoff.cli.get Supplier my.sdf | sort | uniq -c
   paste <(python -m sendoff.cli.get -t my.sdf) <(python -m sendoff.cli.get CAS my.sdf)
   python -m sendoff.cli.get ID my.sdf | sed 's/^/TL-/' |
       python -m sendoff.cli.title --from - my.sdf > titled.sdf
   python -m sendoff.cli.filter --eq Availability In-stock vendor.sdf |
       python -m sendoff.cli.set Batch --value B7 -o batch.sdf
   python -m sendoff.cli.filter --unique InChIKey my.sdf > unique.sdf

The ``count``, ``get``, ``set``, ``title``, and ``filter`` commands read
``FILE`` or ``-i FILE``, else stdin, as bytes and write ``-o FILE``, else
stdout; ``-o`` cannot name an input. ``get`` prints one line per block, empty
where a record is missing, and ``--from`` reads one line per block; ``-z`` uses
NUL instead. Run any with ``-h``.

In Python, read compressed files with ``parse_sdf(gzip.open(path))``, and
write a pipe with ``write_sdf(sys.stdout, blocks)``, or ``sys.stdout.buffer``
for bytes blocks.
