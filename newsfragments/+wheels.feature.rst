Wheels are built with profile-guided optimization for CPython 3.11 to 3.15,
one per version, on Linux (x86_64, aarch64; glibc 2.28 or later, and musl
1.2), macOS (x86_64, arm64), and Windows (x86_64, ARM64). An abi3 wheel covers
later CPython versions; per-version wheels are 1-3% faster, up to 10% when
building dicts.
