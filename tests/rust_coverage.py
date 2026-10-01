"""Collect Rust line/function coverage from the default public Python tests."""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from importlib import import_module
from importlib.machinery import ExtensionFileLoader
from importlib.util import module_from_spec, spec_from_loader
from pathlib import Path
from tempfile import TemporaryDirectory
from types import ModuleType

import pytest
from _pytest.terminal import TerminalReporter


def start() -> tuple[ModuleType, TemporaryDirectory[str], Path]:
    """Load a cached test-only instrumented extension before fixture imports.

    Returns:
        Native module, isolated session files, and compiler-matched LLVM tools.

    Raises:
        UsageError: If coverage cannot be initialized safely.
    """
    if "sendoff.native" in sys.modules:
        raise pytest.UsageError("Rust coverage must start before importing sendoff.")
    library_directory = Path(
        subprocess.check_output(
            ["rustc", "--print", "target-libdir"], text=True
        ).strip()
    )
    tools = library_directory.parent / "bin"
    if not all((tools / name).is_file() for name in ("llvm-profdata", "llvm-cov")):
        raise pytest.UsageError(
            "Rust coverage requires compiler-matched LLVM tools. "
            "Run: rustup component add llvm-tools-preview"
        )
    root = Path(__file__).resolve().parents[1]
    cache = root / ".cache" / "rust-coverage"
    cache.mkdir(parents=True, exist_ok=True)
    session = TemporaryDirectory(prefix="run-", dir=cache)
    directory = Path(session.name)
    environment = os.environ.copy()
    environment["RUSTFLAGS"] = (
        f"{environment.get('RUSTFLAGS', '')} -C instrument-coverage"
    )
    environment["PYO3_PYTHON"] = sys.executable
    environment["LLVM_PROFILE_FILE"] = str(directory / "%m-%p.profraw")
    build = subprocess.run(
        [
            "cargo",
            "build",
            "--locked",
            "--features",
            "coverage",
            "--target-dir",
            str(cache / "target"),
            "--message-format=json-render-diagnostics",
        ],
        cwd=root,
        env=environment,
        stdout=subprocess.PIPE,
        text=True,
        check=True,
    )
    artifacts = [
        Path(message["filenames"][0])
        for line in build.stdout.splitlines()
        if (message := json.loads(line))["reason"] == "compiler-artifact"
        and message["target"]["name"] == "native"
    ]
    if len(artifacts) != 1:
        raise pytest.UsageError("Cargo did not produce one native coverage library.")
    for profile in directory.glob("*.profraw"):
        profile.unlink()
    library = directory / "native.abi3.so"
    shutil.copyfile(artifacts[0], library)
    os.environ["LLVM_PROFILE_FILE"] = environment["LLVM_PROFILE_FILE"]
    loader = ExtensionFileLoader("sendoff.native", str(library))
    spec = spec_from_loader("sendoff.native", loader)
    if spec is None:
        raise pytest.UsageError("Cannot load the native coverage library.")
    module = module_from_spec(spec)
    sys.modules[spec.name] = module
    loader.exec_module(module)
    setattr(import_module("sendoff"), "native", module)
    return module, session, tools


def report(
    module: ModuleType,
    session: TemporaryDirectory[str],
    tools: Path,
    terminalreporter: TerminalReporter,
) -> None:
    """Dump and report coverage before removing isolated session artifacts.

    Args:
        module: The instrumented extension exercised by the Python tests.
        session: Profiles and library for this run alone.
        tools: LLVM tools belonging to the active Rust compiler.
        terminalreporter: Pytest's terminal summary.

    Raises:
        RuntimeError: If the instrumented workload produced no profile.
    """
    module.flush_coverage()
    directory = Path(session.name)
    profiles = sorted(directory.glob("*.profraw"))
    if not profiles:
        raise RuntimeError("The Python tests produced no Rust coverage profile.")
    merged = directory / "coverage.profdata"
    subprocess.run(
        [
            str(tools / "llvm-profdata"),
            "merge",
            "-sparse",
            *map(str, profiles),
            "-o",
            str(merged),
        ],
        check=True,
    )
    root = Path(__file__).resolve().parents[1]
    sources = sorted((root / "rust" / "src").rglob("*.rs"))
    arguments = [str(module.__file__), f"--instr-profile={merged}", *map(str, sources)]
    summary = subprocess.check_output(
        [
            str(tools / "llvm-cov"),
            "report",
            *arguments,
            "--show-region-summary=false",
            "--show-branch-summary=false",
        ],
        text=True,
    )
    html = directory.parent / "html"
    subprocess.run(
        [
            str(tools / "llvm-cov"),
            "show",
            *arguments,
            "--format=html",
            f"--output-dir={html}",
        ],
        check=True,
    )
    terminalreporter.write_sep("=", "Rust coverage from the Python tests")
    terminalreporter.write_line(summary.rstrip())
    terminalreporter.write_line(f"Rust HTML coverage: {html / 'index.html'}")
    session.cleanup()
