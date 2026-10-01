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

Native migration status
-----------------------

This feature-branch scaffold is **unreleased**, still version 0.1.6. Importing
``sendoff`` always loads the mandatory private ``sendoff._native`` extension.
There is no optional import or supported pure-Python installation. Existing
Python algorithms remain transitionally unchanged until their assigned native
leaves are integrated; this scaffold does **not** claim a completed Rust port.
Private algorithm exports currently raise explicit ``NotImplementedError``
messages and must not be used as implemented operations.

Runtime policy is GIL-enabled CPython 3.11 through the latest stable release on
Linux and macOS, not PyPy, Windows, musllinux or free-threaded Python. Intended
wheel baselines are manylinux2014/2.17 x86_64 and aarch64, macOS Intel 10.15 and
arm64 11.0. The A2 ABI feasibility probe succeeded on four native platforms at
https://github.com/pechersky/sendoff/actions/runs/36867567700; it is not production
package CI or evidence of execution on the oldest macOS versions. Local native
Linux wheels are development artifacts, not portable manylinux release wheels.
Production wheel/testing/release jobs remain separate work; the existing
Poetry-based ``tox.yml`` workflow is not yet migrated.

Development and artifacts
-------------------------

Rust/Cargo and Python >=3.11 are required for source builds. Use a virtual
environment and project-local scratch space::

    python3.11 -m venv .venv
    mkdir -p .cache/build-tmp
    export TMPDIR="$PWD/.cache/build-tmp"
    .venv/bin/python -m pip install -e '.[dev]'
    VIRTUAL_ENV="$PWD/.venv" .venv/bin/maturin develop --locked
    .venv/bin/python -m pytest --basetemp=.pytest_cache/dev
    .venv/bin/mypy sendoff tests --ignore-missing-imports --strict
    PRE_COMMIT_HOME="$PWD/.cache/pre-commit" .venv/bin/pre-commit run --all-files
    PYO3_PYTHON="$PWD/.venv/bin/python" cargo build --locked
    cargo fmt --all -- --check
    PYO3_PYTHON="$PWD/.venv/bin/python" cargo clippy --locked --all-targets -- -D warnings
    .venv/bin/maturin build --release --locked --compatibility linux --out dist
    .venv/bin/maturin sdist --out dist

Editable installation is mandatory for source-tree tests; rebuilding with
``maturin develop`` picks up leaf changes. ``tox -e py311`` (tox >=4.24) installs
via PEP 660 and runs pytest and mypy directly, without Poetry or
tox-poetry-installer. Uninstall the obsolete tox-poetry-installer when upgrading
an existing environment to tox 4. Other declared environments are py312, py313
and py314; their presence is not a claim they were run locally. Check wheel and
sdist installations in fresh environments away from repository imports, e.g.
``python -I -c 'import sendoff, sendoff._native'`` with the artifact installed.
Do not infer artifact validity from a source-directory import.

The ``test`` extra includes all existing RDKit probes and pytest plugins.
CPython 3.11 pins the A1 oracle environment's published RDKit 2022.9.5 and NumPy
1.23.5 (2022.9.2 disappeared from PyPI). Newer interpreters use RDKit 2026.3.6
and compatible NumPy >=1.26.4,<3 selected by pip. New toolkit versions are a
declared test policy, not a waiver: the same expected outcomes and strict
xfails must pass in production validation. Toolkit availability/OS floors are
separate from sendoff's runtime policy; RDKit 2026.3.6 has no CPython 3.14 macOS
Intel wheel. That environment needs a source-built toolkit or an explicit
upstream-supported solution, not skipped RDKit probes. No toolkit is a sendoff
runtime dependency.

The CPython 3.11 development/typing extra retains A1's mypy 1.7.0; newer
interpreters use mypy >=1.19. The unchanged baseline source also reports five
``StringIO``/``TextIOWrapper`` argument mismatches under mypy 1.19/1.20.
Modern whole-suite typing remains a documented integration risk, not a passing
gate or permission to alter frozen public annotations. The existing mypy 1.19
pre-commit hook remains unchanged and checks scoped edits.

Coverage.py/pytest-cov report **Python adapter statements and branches**, not
Rust source coverage. Native behavioral/regression checks remain Python tests;
Cargo build/fmt/Clippy are compilation/style/static checks, not behavior tests.
The pre-commit Git exclusion is anchored to ``.git/`` so it no longer silently
excludes ``.github/`` YAML. Rust fmt/Clippy hooks run for Cargo/Rust changes.

``Cargo.toml`` is the version source: maturin reads it into PEP 621 metadata,
Rust exports ``CARGO_PKG_VERSION``, and Python reads that export. tbump updates
Cargo and its own configuration, then runs ``cargo check`` to synchronize and
stage ``Cargo.lock`` before a version commit. Inspect
``tbump --dry-run --no-push --no-tag <version>`` before release work.
After a version bump, rebuild the editable extension with ``maturin develop``
before reading Python's version or running ``.venv/bin/towncrier --draft``;
``cargo check`` synchronizes the lock but does not reinstall the binary.
No 1.0.0 bump, tag or publishing is authorized by this scaffold.

Private leaf interface and ownership
------------------------------------

Main owns Cargo metadata/lock, ``rust/src/lib.rs``, ``sendoff/_native.pyi``,
``sendoff/__init__.py`` and both public facades. Each leaf implements only its
Rust file and its own Python tests. All four
``pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()>`` functions are
already called by ``lib.rs``. Workers compile and test directly via
``import sendoff._native as native`` without editing manifests, shared export
wiring or public Python files. Rebuild editable, then run the leaf's Python
tests. Main alone wires public facades after review.

All approved A2 arguments are required positional-or-keyword and all native
input operands are generic ``&Bound<'_, PyAny>``. Return contracts are typed
in ``sendoff/_native.pyi``. No shared behavioral helper, shadow native facade
object, algorithm or Rust behavioral test has been created.

* A4, ``rust/src/framing.rs``:
  ``_mdl_iter(lines)``, ``_metadata_iter(lines)``,
  ``_from_block_lines(cls, block_type, lines)``, ``_blocks_iter(cls, lines)``,
  ``_read_sdf_lines(sdfpath)``.
* A5, ``rust/src/sddata.rs``:
  ``_records_iter(block)``, ``_write(block, outh, with_newlines)``,
  ``_append_record(block, record_name, value)``.
* A6, ``rust/src/ctable.rs``:
  ``_ctable_init(table, lines, v3000)``, ``_parse_format(line, formats)``,
  ``_parse_v2000_counts(line)``, ``_parse_v3000_counts(line)``,
  ``_atomlines(table)``, ``_bondlines(table)``.
* A7, ``rust/src/indices.rs``:
  ``_valid_atom_indices(table, strict, v3000, errors)``,
  ``_valid_bond_indices(table, strict, v3000, errors)``,
  ``_renumber_ctable(table, v3000, duplicate_error)``.

``errors`` is the original Python (Mismatch, OutOfOrder, Duplicate) exception
class tuple, ``duplicate_error`` its existing duplicate exception class, and
``v3000`` the original ``CTableFormat.V3000`` member. Never recreate identities
in Rust or coerce operands to Rust strings/integers/bools. A1 classes, dataclass,
deques, annotations, subclass dispatch and pickles remain the oracle.

Later main-owned lazy adapters must use ``for item in native.factory(...):
yield item`` inside genuine Python generators, **not** ``yield from``.
Factories execute on first resume. Native iterators must preserve live deque
mutation errors, distinguish legitimate exhaustion from callback
StopIteration/PEP 479 failures, and implement GC traversal/clear for owned
``Py<PyAny>`` references that can form cycles. A6 raw access must return actual
``itertools.takewhile`` objects with deque iterators captured at call time.
These are handoff requirements, not algorithms implemented by A3.
Do not begin leaf implementation until parent review.
