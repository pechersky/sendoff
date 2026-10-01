1.0.0 (unreleased)
==================

The parser, SDData processing, count parsing, validation and renumbering now
use a Rust/PyO3 backend. Ordinary text and index processing use Rust's standard
library. Python classes, enums, exceptions, mutable deques and generator
interfaces remain compatible. There is no supported pure-Python fallback.
Existing documented parser defects remain unchanged; no measured speedup is
claimed by this compatibility stage.

Supported binaries
------------------

CPython 3.11 through the latest stable GIL-enabled release uses cp311-abi3
wheels for Linux manylinux2014 x86_64/aarch64 and macOS x86_64/arm64.
Windows, musllinux, PyPy and free-threaded Python are not supported.
Byte-oriented/non-UTF-8 APIs, additional lazy views and extended CTAB parsing
remain separate, release-gated stages.

Development and source builds
-----------------------------

Source builds require a current stable Rust toolchain and CPython 3.11 or
newer. Maturin is the build backend; Cargo.toml is the package version source.
The historical README is intentionally unchanged.

.. code-block:: sh

    python3.11 -m venv .venv
    .venv/bin/python -m pip install -e ".[dev]"
    .venv/bin/maturin develop --locked
    .venv/bin/python -m pytest
    .venv/bin/pre-commit run --all-files

Rebuild with maturin after changing Rust sources. To exercise another
supported interpreter, run tox under that interpreter with ``tox -e py``.
Behavioral checks remain Python tests; Rust formatting and Clippy run through
pre-commit.

Release boundary
----------------

Do not tag or publish without explicit approval. Prepare the version and
Cargo.lock together, rebuild the extension, and check the changelog before
creating an approved ``v<version>`` tag. A tag must match Cargo.toml exactly.
Branch pushes never publish.

The tag workflow runs the full Python suite and typing checks, then calls the
same four-platform wheel workflow used by development CI. Each wheel is ABI
audited and runs the standalone Python parity corpus on supported
interpreters. Publication consumes those four verified wheels and the sdist;
it does not rebuild distributions. The existing ``PYPI_API_TOKEN`` secret is
required for publication.
