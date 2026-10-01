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
