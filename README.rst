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
   from sendoff import SDBlock, parse_sdf, write_sdf

   # str paths give SDBlock[str]; bytes paths give SDBlock[bytes], for feeds
   # that are not valid UTF-8.
   for block in parse_sdf("compounds.sdf"):
       block.title
       block["ID"]  # the last value named "ID"; KeyError if none
       block.get("ID")  # or None
       block.to_dict()  # every record, later duplicates winning
       block.records()  # [(name, value), ...] in file order
       block.iter_records()  # the same, lazily
       "ID" in block, len(block), list(block)  # names only, no values built
       "V3000" in block.mdl  # the CTAB after the title, through M  END
       mol = Chem.MolFromMolBlock(block.molblock)

Many SDFs do not carry the same records on every block: ADMET Predictor
exports, for one, have hundreds of record names, some on only part of the
blocks. Access is per block, so ``get`` and ``to_dict`` assume no fixed columns;
``Counter(n for b in parse_sdf(path) for n in b.names())`` surveys a file.

Edit blocks and write them back:

.. code-block:: python

   def tag(block):
       block["Source"] = "vendor"  # set_record: replace in place, or append
       block.remove_record("Internal note")
       return block

   write_sdf("out.sdf", (tag(block) for block in parse_sdf("in.sdf")))

   block = SDBlock(Chem.MolToMolBlock(mol), {"ID": "X-1"})
   block.molblock = Chem.MolToMolBlock(other)  # new structure, same SD data
   block.to_sdf()  # the block as SDF text ending in $$$$
   blocks = list(SDBlock.from_sdf(sdf_text))  # parse SDF held in memory

New header lines are ``> <name>``; pass ``header="padded"`` for ``>  <name>``,
as RDKit and Open Babel write, or ``header="copy"`` to follow the block's first
header, which is a bit slower. Existing header lines are kept. Names are written
unchecked; the CTfile format wants a letter, then letters, digits, or
underscores. Reading records past a header line with no ``<name>`` raises
``sendoff.MalformedRecordError``, a ``ValueError``.

``write_sdf`` to a path writes in large chunks with the GIL released, so prefer
it to a ``block.write`` loop. ``parse_sdf`` memory-maps regular files on local
Linux filesystems, so write results to a new path: truncating a file while it
is parsed crashes the process.

Non-UTF-8 files
---------------

A str path decodes text as UTF-8 when it is read, so a stray byte raises
``UnicodeDecodeError`` from the title, name, or value holding it. Only the
caller knows a feed's encoding:

.. code-block:: python

   # One known encoding: a text file, read through Python, not memory-mapped.
   # write_sdf writes str blocks as UTF-8, so this also transcodes the feed.
   with open("vendor.sdf", encoding="gbk") as f:
       write_sdf("vendor-utf8.sdf", parse_sdf(f))

   # Mixed or unknown: bytes, decoded per value as the feed needs.
   for block in parse_sdf(os.fsencode("vendor.sdf")):
       name = block[b"Product name"].decode("gbk", "replace")

Bytes blocks write back byte for byte. ``open(path, errors="replace")`` reads
on past bad bytes, at the cost of the text they held.

Shell pipelines
---------------

.. code-block:: bash

   sdf-record() {  # one value per block, newlines escaped
       python -c '
   import signal, sys
   from sendoff import parse_sdf
   signal.signal(signal.SIGPIPE, signal.SIG_DFL)
   for block in parse_sdf("/dev/stdin"):
       value = block.get(sys.argv[1])
       if value is not None:
           print(value.replace("\n", r"\n"))
   ' "$1"
   }

   cat my.sdf | sdf-record Supplier | sort | uniq -c

Write SDF to a pipe with ``write_sdf(sys.stdout, blocks)``; the path
``"/dev/stdout"`` is opened with truncation. For non-UTF-8 feeds, parse
``b"/dev/stdin"`` and write to ``sys.stdout.buffer``. ``parse_sdf`` also reads
file objects such as ``sys.stdin.buffer``.

CTable errors
-------------

``CTable`` and the ``SDBlock`` CTable methods raise the same exception for a
problem whichever method meets it. All but ``NotImplementedError`` are
``sendoff.CTableError``, a ``ValueError``:

- ``MalformedCTableError``: a line breaks the format. The message names the
  line, numbered from the title as line 1, as the block numbers it.
- ``IndicesMismatchError``: fewer or more atom or bond lines than the counts
  line declares.
- ``IndicesOutOfOrderError``: with ``strict``, an index is not its one-based
  position.
- ``IndicesDuplicateError``: an index repeats, or renumbering meets a bond to
  an atom index that several atoms share.
- ``NotImplementedError``: the operation does not support the table: V2000
  validation, renumbering, or ``fix_counts``, or constructs renumbering or
  ``fix_counts`` would corrupt.

Migrating from 0.1.x
--------------------

``SDBlock`` is no longer a dataclass of line deques. The common rewrites:

.. code-block:: python

   str.join("", block.mdl)  # block.mdl
   any(s in line for line in block.mdl)  # s in block.mdl
   block.title + "\n" + str.join("", block.mdl)  # block.molblock
   dict(block.records())  # block.to_dict()
   dict(block.records())[name]  # block[name]
   dict(block.records()).get(name)  # block.get(name)
   len(list(block.records()))  # len(block)
   name in dict(block.records())  # name in block
   [name for name, _ in block.records()]  # block.names()
   next(block.records())  # next(block.iter_records())
   SDBlock.from_block_lines(lines)  # SDBlock(str.join("", lines))
   SDBlock.from_lines(text.splitlines(True))  # SDBlock.from_sdf(text)
   block.write(handle, with_newlines=True)  # block.write(handle)

``block.mdl`` is now one string, so code that iterated it as lines still runs
but walks characters: ``any(s in line for line in block.mdl)`` is always
``False`` for a longer ``s``. Check every ``.mdl`` use.

``block.metadata`` is gone: edit records with ``set_record``,
``remove_record``, ``replace_records``, or ``append_record``. ``mdl`` is
read-only, so assign ``block.molblock`` to change the structure. The release
notes list every change.

Rust development checks
-----------------------

Ownership-focused Clippy lints are configured in ``Cargo.toml``. Pre-commit
and GitHub Actions check all targets and features, treating warnings as errors:

.. code-block:: console

   cargo clippy --locked --all-targets --all-features -- -D warnings

The shared VS Code settings make rust-analyzer run Clippy instead of
``cargo check``, using the same lint configuration.

Release wheels
--------------

Build release wheels with profile-guided optimization on top of the release
profile's fat LTO:

.. code-block:: console

   rustup component add llvm-tools
   scripts/pgo-build --out dist

It builds an instrumented wheel, runs it over a synthetic SDF it generates
(heterogeneous records, V2000 and V3000, str and bytes, CRLF), and rebuilds
with that profile; the result is reproducible. Real catalog and ADMET files
stay held out: on them record access (``get``, ``in``, ``len``, ``names``,
``to_dict``) got 5-18% faster and scan-bound work stayed within noise.
``maturin build --release`` still makes a correct, LTO-only wheel.

The tests normally load a coverage-instrumented debug build of the checkout.
To test a built wheel, put it ahead of the checkout; the report header names
the ``sendoff.native`` under test:

.. code-block:: console

   python -m zipfile -e dist/sendoff-*.whl /tmp/sendoff-wheel
   PYTHONPATH=/tmp/sendoff-wheel .venv/bin/pytest --import-mode=append --no-rust-coverage

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
