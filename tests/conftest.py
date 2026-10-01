"""Configure pytests."""

import sys

from _pytest.terminal import TerminalReporter

from tests.rust_coverage import report, start

coverage = (
    None
    if any(
        argument in sys.argv
        for argument in ("--help", "-h", "--version", "-V", "--collect-only", "--co")
    )
    else start()
)

pytest_plugins = ["tests.fixtures"]


def pytest_terminal_summary(terminalreporter: TerminalReporter) -> None:
    """Include Rust coverage collected by the public Python test workload.

    Args:
        terminalreporter: Pytest's terminal summary.
    """
    if coverage is not None:
        report(*coverage, terminalreporter)
