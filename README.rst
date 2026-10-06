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
to simple tools to interrogate the files while not modifying the file is crucial.

This package also tried to document the "canonical" ways metadata is handled by the larger
packages. To wit, there are tests to monitor how, for example, rdkit deals with molecules
that have multiline record values, or a "$$$$" molecule title.

Usage
-----

.. code-block:: python

   from rdkit import Chem
   from sendoff.sdblock import SDBlock, parse_sdf, write_sdf

   # str paths give SDBlock[str]; bytes paths give SDBlock[bytes], for feeds
   # that are not valid UTF-8.
   for block in parse_sdf("compounds.sdf"):
       block.title
       block.get("ID")  # the last value named "ID", or None
       block.to_dict()  # every record, later duplicates winning
       block.records()  # [(name, value), ...] in file order
       "ID" in block, len(block), block.names()  # names only, no values built
       "V3000" in block.mdl  # the CTAB after the title, through M  END
       block.mdl_contains("V3000")  # the same search, without building it
       mol = Chem.MolFromMolBlock(block.molblock)

Many SDFs do not carry the same records on every block: ADMET Predictor
exports, for one, have hundreds of record names, some on only part of the
blocks. Access is per block, so ``get`` and ``to_dict`` assume no fixed columns;
``Counter(n for b in parse_sdf(path) for n in b.names())`` surveys a file.

Edit blocks and write them back:

.. code-block:: python

   def tag(block):
       block.set_record("Source", "vendor")  # replace in place, or append
       block.remove_record("Internal note")
       return block

   write_sdf("out.sdf", (tag(block) for block in parse_sdf("in.sdf")))

   # Keep blocks whose CTAB has a charge line; the rest are skipped unbuilt.
   write_sdf("charged.sdf", parse_sdf("in.sdf", mdl_contains="M  CHG"))

   block = SDBlock(Chem.MolToMolBlock(mol), {"ID": "X-1"})
   block.molblock = Chem.MolToMolBlock(other)  # new structure, same SD data
   block.to_sdf()  # the block as SDF text ending in $$$$
   blocks = list(SDBlock.from_sdf(sdf_text))  # parse SDF held in memory

``write_sdf`` to a path writes in large chunks with the GIL released, so prefer
it to a ``block.write`` loop.

Native typing stubs
-------------------

``sendoff/native.pyi`` is generated: ``scripts/native-stub`` builds the
extension with the ``stubgen`` feature, takes the signatures, the type hints in
``#[pyo3(signature = ...)]``, and the Rust doc comments from PyO3 introspection,
then applies ``scripts/native-stub.patch`` for what introspection cannot
express. That covers ``SDBlock`` being generic over ``str`` and ``bytes``,
overloads, typed getters and iterators, and the file-object protocols.

* Change docstrings and simple types in the Rust source; the ``native-stub``
  pre-commit hook regenerates the stub.
* For typing that needs the patch, edit ``sendoff/native.pyi`` and run
  ``scripts/native-stub --update-patch``.
* ``scripts/native-stub --check`` fails if the committed stub is stale.

``tests/test_stubs.py`` runs ``mypy.stubtest`` against the built extension, so
a signature that drifts from the Rust API fails the test suite.
