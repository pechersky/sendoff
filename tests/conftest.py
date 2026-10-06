"""Configure pytests."""

import sys
from importlib import import_module

import pytest
from _pytest.terminal import TerminalReporter

from tests.rust_coverage import report, start

coverage = (
    None
    if any(
        argument in sys.argv
        for argument in (
            "--help",
            "-h",
            "--version",
            "-V",
            "--collect-only",
            "--co",
            "--no-rust-coverage",
        )
    )
    else start()
)

pytest_plugins = ["tests.fixtures"]


def pytest_addoption(parser: pytest.Parser) -> None:
    """Add the option to test the importable extension.

    Args:
        parser: Pytest's command-line parser.
    """
    parser.addoption(
        "--no-rust-coverage",
        action="store_true",
        help="test the sendoff.native that imports normally, such as an "
        "installed release wheel, instead of a coverage build",
    )


def pytest_report_header() -> str:
    """Name the native extension under test.

    Returns:
        The extension's path.
    """
    return f"sendoff.native: {import_module('sendoff.native').__file__}"


def pytest_terminal_summary(terminalreporter: TerminalReporter) -> None:
    """Include Rust coverage collected by the public Python test workload.

    Args:
        terminalreporter: The active terminal reporter.
    """
    if coverage is not None:
        report(*coverage, terminalreporter)
