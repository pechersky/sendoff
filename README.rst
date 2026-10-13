sendoff
=======

.. image:: https://results.pre-commit.ci/badge/github/pechersky/sendoff/main.svg
   :target: https://results.pre-commit.ci/latest/github/pechersky/sendoff/main
   :alt: pre-commit.ci status

.. image:: https://github.com/pechersky/sendoff/actions/workflows/tox.yml/badge.svg
   :target: https://github.com/pechersky/sendoff/actions/workflows/tox.yml
   :alt: Tox status

The minimal SDF metadata parser.
 
Often, SDFs have lots of useful metadata on them in the title and record fields/values.
However, reading a molecule (via rdkit, OpenEye toolkits, etc) can be slow because those
libraries also construct the molecules. Modifying the metadata, or filtering/sorting based
on the metadata also can induce non-idempotent differences in the file based on
opinionated approaches in chemical libraries.
 
This library strives to be able to handle SDF files even with malformed chemistry or
metadata. Since much debugging of our files and data deals with such files, having access
to simple fast tools to interrogate the files while not modifying the file is crucial.
 
This package also tries to document the "canonical" ways metadata is handled by the larger
packages. To wit, there are tests to monitor how, for example, rdkit deals with molecules
that have multiline record values, or a "$$$$" molecule title.


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

   python -m sendoff.cli.count my.sdf
   python -m sendoff.cli.get Supplier my.sdf | sort | uniq -c
   paste <(python -m sendoff.cli.get -t my.sdf) <(python -m sendoff.cli.get CAS my.sdf)
   python -m sendoff.cli.get ID my.sdf | sed 's/^/TL-/' |
       python -m sendoff.cli.title --from - my.sdf > titled.sdf
   python -m sendoff.cli.filter --eq Availability In-stock vendor.sdf |
       python -m sendoff.cli.set Batch --value B7 > batch.sdf

The ``count``, ``get``, ``set``, ``title``, and ``filter`` commands read a file
or stdin as bytes and write stdout. ``get`` prints one line per block, empty
where a record is missing, and ``--from`` reads one line per block; ``-z`` uses
NUL instead. Run any with ``-h``.

In Python, write a pipe with ``write_sdf(sys.stdout, blocks)``, or
``sys.stdout.buffer`` for bytes blocks.
