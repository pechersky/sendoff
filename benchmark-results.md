# Benchmark results

sendoff's native rewrite against 0.1.6 (the last pure-Python release), RDKit
2026.3.6 and OEChem 2026.1.0, measured 2026-10-14 on EC2 c6i.2xlarge (x86_64) and
c7g.2xlarge (arm64): Amazon Linux 2023, Python 3.11.16, sendoff from a
profile-guided release wheel (`scripts/pgo-build`) and a plain release wheel
(`sendoff plain`). sendoff is faster than 0.1.6, RDKit and OEChem in every
cell.

Every number is other ÷ sendoff, so above 1 means sendoff is faster. Cells read
`wall / op`:

- *wall*: Hyperfine's median of 5 randomized runs after a warmup, each a fresh
  interpreter, so startup and imports count.
- *op*: the median time of the operation itself, timed in process after imports.

† marks output that differs from sendoff's (see [Differences](#differences)),
*error* a run that raised, and - a run not made: 0.1.6 lacks the newer methods
and cannot read B, which is not UTF-8. Geomeans take each toolkit's fastest
mode per case, whether or not its output matched. Write cases are not compared,
as the toolkits reserialize records.

## At a glance

How many times faster sendoff is than each baseline, in operation time on
x86_64 / arm64, for a sample of the full matrix: a case or two per
[dataset](#datasets), including where the margin is smallest. Write cases are
not like for like, as the toolkits reserialize records.

| Dataset | Case | sendoff plain | 0.1.6 | RDKit default | RDKit fast | OEChem default | OEChem tokenizer | OEChem threaded |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| B full | `molecule_count` | 1.08× / 1.07× | - | 197× / 165× † | 92.7× / 81.5× | 90× / 73.8× † | 26.4× / 13.8× | 20.2× / 10.5× |
| A full | `title_extraction` | 1.09× / 1.02× | 30× / 27.4× | 107× / 91.7× † | 56.9× / 49.1× | 43.4× / 34.7× | 10.8× / 5.15× | 9.84× / 5.24× |
| C full | `field_get` | 1.21× / 1.07× | - | 90.3× / 72.4× | 44× / 36.2× | 39.2× / 29.3× | 32× / 25× | 23.9× / 20.6× |
| B full | `field_get` | 1.22× / 1.11× | - | 86.2× / 67.9× † | 41.3× / 34.3× | 40.4× / 30.5× † | 49.8× / 39.9× | 39.1× / 34.3× |
| A full | `record_dicts` | 1.11× / 1.08× | - | 53.6× / 47.3× † | 32.4× / 29.5× | 32× / 26.8× | 19.1× / 16× | 14.4× / 13.4× |
| B full | `record_name_set` | 1.12× / 1.1× | - | 66.5× / 54.2× | 35.9× / 30.7× | 49.1× / 39.8× | 32.7× / 27.5× | 25.8× / 23.6× |
| D full | `field_extraction` | 1.06× / 1.1× | 17.7× / 17.4× | 23.1× / 19.3× † | 19.6× / 16.4× | 21.6× / 20.1× | 16.7× / 15.5× | 16.5× / 15.3× |
| A full | `write_round_trip` | 1.06× / 1.01× | 24.5× / 25.4× | 115× / 99.4× | 73.1× / 62.1× | 52.9× / 55.4× | 7.66× / 4.6× | 5.4× / 3.11× |
| A full | `write_sdf_round_trip` | 1.07× / 1.02× | - | 188× / 150× | 119× / 93.5× | 86.4× / 83.4× | 12.5× / 6.92× | 8.89× / 4.69× |
| A full | `mdl_scan` | 1.05× / 1.02× | 14.7× / 15.2× | 2.82× / 2.76× | 2.79× / 2.76× | 6.44× / 3.45× | 7.78× / 4.41× | 4.36× / 2.57× |
| D full | `mdl_scan` | 0.977× / 1× | 63.1× / 68.6× | 2.73× / 2.68× | 2.75× / 2.67× | 19.9× / 8.39× | 21.4× / 9.67× | 21.1× / 9.06× |

On 100-record files, where start-up dominates, wall time is 3.0–3.5× faster than
0.1.6, 8.7–9.9× than RDKit and 10–12× than OEChem ([geomeans](#summary)).
Profile-guided optimization cuts operation time by 5–15% on representative and
full datasets (1.09× geomean on x86_64, 1.06× on arm64), and
`parse_sdf(path, borrow=True)` by another 7–16% on x86_64 and 3–12% on arm64
across B–D.

## Implementations

| Name | Reads with |
|---|---|
| sendoff | `parse_sdf` and `SDBlock` methods (PGO release wheel, `scripts/pgo-build`) |
| sendoff plain | `parse_sdf` and `SDBlock` methods (plain release wheel) |
| 0.1.6 | `parse_sdf` and `SDBlock` methods (pure Python) |
| RDKit default | `Chem.SDMolSupplier(path)`; skips records it cannot parse |
| RDKit fast | `Chem.SDMolSupplier(path, sanitize=False, removeHs=False, strictParsing=False)` |
| OEChem default | `oemolistream` and `OEReadMolecule` |
| OEChem tokenizer | `OEGetTokenizer` raw records, SD data parsed in Python |
| OEChem threaded | `oemolithread` raw records, SD data parsed in Python |

Toolkit runs compute what sendoff's do (counts, lists, sets), not more.

## Datasets

A, B and C are vendor catalogs; D is in-house property predictions with many
fields per record. sendoff reads B's representative and full tiers as bytes,
as their values are GBK.

| Dataset | Tier | Records | MB |
|---|---|---:|---:|
| A | smoke | 100 | 0.19 |
| A | representative | 10,000 | 18.3 |
| A | full | 151,527 | 279 |
| B | smoke | 100 | 0.22 |
| B | representative | 10,000 | 22.1 |
| B | full | 328,135 | 722 |
| C | smoke | 100 | 0.16 |
| C | representative | 10,000 | 15.5 |
| C | full | 11,251 | 17.7 |
| D | representative | 1,995 | 35.1 |
| D | full | 1,998 | 44.8 |

## Cases

`blocks` is `parse_sdf(path)`; output goes to `os.devnull`.

| Case | sendoff |
|---|---|
| `molecule_count` | `sum(1 for _ in blocks)` |
| `title_extraction` | `[b.title for b in blocks]` |
| `record_count_distribution` | `[sum(1 for _ in b.records()) for b in blocks]` |
| `field_extraction` | `[v for b in blocks for k, v in b.records() if k == name]` |
| `field_get` | `[v for b in blocks if (v := b.get(name)) is not None]` |
| `field_presence` | `sum(name in b for b in blocks)` |
| `record_len_distribution` | `[len(b) for b in blocks]` |
| `record_name_set` | `{n for b in blocks for n in b.names()}` |
| `record_dicts` | `sum(len(b.to_dict()) for b in blocks)` |
| `mdl_scan` | `sum("M  CHG" in b.mdl for b in blocks)` |
| `write_round_trip` | `b.write(out)` for each block |
| `append_record_write` | `b.append_record(name, str(i))`, then `b.write(out)` |
| `write_sdf_round_trip` | `write_sdf(os.devnull, blocks)` |

## Summary

### Geomean speedup, wall time

| Dataset | sendoff plain x86_64 | sendoff plain arm64 | 0.1.6 x86_64 | 0.1.6 arm64 | RDKit x86_64 | RDKit arm64 | OEChem x86_64 | OEChem arm64 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| A smoke | 0.988 | 1.03 | 3.41 | 3.07 | 9.92 | 8.99 | 11.9 | 10.2 |
| A representative | 1.05 | 1.03 | 13.5 | 12.8 | 27.5 | 23.1 | 12.7 | 9.89 |
| A full | 1.09 | 1.04 | 20.2 | 19.8 | 35.7 | 31.4 | 11.8 | 9.09 |
| B smoke | 0.989 | 0.995 | 3.49 | 3.1 | 9.9 | 8.72 | 12.1 | 10 |
| B representative | 1.09 | 1.03 | - | - | 39.6 | 32.8 | 16.8 | 12.9 |
| B full | 1.1 | 1.08 | - | - | 68.6 | 58 | 18.2 | 13.5 |
| C smoke | 0.979 | 1.01 | 3 | 3.08 | 8.99 | 9.06 | 10.6 | 10.4 |
| C representative | 1.06 | 1.03 | 11.4 | 11.4 | 26.4 | 23.4 | 12.4 | 10.7 |
| C full | 1.03 | 1.03 | 11.7 | 11.9 | 27 | 24.3 | 12.2 | 10.6 |
| D representative | 1.05 | 1.05 | 32.3 | 31.8 | 32.2 | 26.4 | 22.1 | 15.9 |
| D full | 1.01 | 1.06 | 35.5 | 37.2 | 33.7 | 28.8 | 22.2 | 16.5 |

### Geomean speedup, operation time

| Dataset | sendoff plain x86_64 | sendoff plain arm64 | 0.1.6 x86_64 | 0.1.6 arm64 | RDKit x86_64 | RDKit arm64 | OEChem x86_64 | OEChem arm64 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| A smoke | 1.02 | 1.05 | 10.8 | 11.5 | 34.5 | 31.4 | 8.66 | 7.09 |
| A representative | 1.1 | 1.05 | 21.3 | 19.7 | 40.1 | 32.5 | 12.7 | 8.83 |
| A full | 1.09 | 1.05 | 21.2 | 20.6 | 37.1 | 32.4 | 11.7 | 8.8 |
| B smoke | 1.03 | 1.04 | 15 | 15.7 | 36.3 | 32.3 | 11.9 | 9.19 |
| B representative | 1.15 | 1.07 | - | - | 78.5 | 59.1 | 20.7 | 13.6 |
| B full | 1.11 | 1.08 | - | - | 71.1 | 59.8 | 18.4 | 13.5 |
| C smoke | 1.01 | 1.04 | 8.81 | 10.2 | 35.5 | 36.1 | 7.77 | 7.41 |
| C representative | 1.13 | 1.05 | 20.9 | 20.1 | 46.3 | 38.3 | 13.8 | 10.3 |
| C full | 1.08 | 1.05 | 20.3 | 20.2 | 45.2 | 38.4 | 13.3 | 10.2 |
| D representative | 1.07 | 1.06 | 58.7 | 53.6 | 46.7 | 35.4 | 29.9 | 17.4 |
| D full | 1.01 | 1.06 | 57.1 | 58.5 | 44.7 | 36.7 | 28.2 | 18 |

### Peak RSS, MB, x86_64 (largest over cases)

| Dataset | sendoff | sendoff plain | 0.1.6 | RDKit default | RDKit fast | OEChem default | OEChem tokenizer | OEChem threaded |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| A smoke | 10 | 10 | 12 | 57 | 57 | 55 | 53 | 53 |
| A representative | 12 | 12 | 68 | 58 | 57 | 56 | 71 | 60 |
| A full | 26 | 26 | 864 | 68 | 67 | 64 | 332 | 71 |
| B smoke | 10 | 10 | 12 | 57 | 56 | 56 | 53 | 53 |
| B representative | 12 | 12 | - | 57 | 57 | 56 | 75 | 62 |
| B full | 44 | 44 | - | 84 | 82 | 76 | 792 | 83 |
| C smoke | 10 | 10 | 12 | 57 | 57 | 54 | 53 | 53 |
| C representative | 12 | 12 | 56 | 57 | 57 | 55 | 69 | 60 |
| C full | 12 | 12 | 62 | 57 | 57 | 55 | 71 | 60 |
| D representative | 12 | 11 | 158 | 58 | 57 | 55 | 87 | 86 |
| D full | 12 | 11 | 211 | 57 | 57 | 55 | 97 | 96 |

## Full results

### x86_64 (c6i.2xlarge)

| Dataset | Case | sendoff ms | sendoff plain | 0.1.6 | RDKit default | RDKit fast | OEChem default | OEChem tokenizer | OEChem threaded |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| A smoke | `append_record_write` | 8.5 / 0.277 | 0.983 / 1.03 | 3.39 / 10.8 | 10.9 / 67.9 | 10.1 / 47.4 | 12.4 / 26.4 | 11.8 / 8.42 | 11.5 / 7.76 |
| A smoke | `field_extraction` | 8.42 / 0.3 | 0.986 / 1.04 | 3.45 / 11.1 | 10.4 / 44.9 | 10 / 32.4 | 12.3 / 19.8 | 12 / 11.3 | 11.7 / 10.5 |
| A smoke | `field_get` | 8.28 / 0.191 | 0.987 / 1.04 | - | 10.5 / 60.3 | 9.89 / 39.8 | 12.6 / 19.2 | 12.8 / 17.9 | 12.5 / 16.8 |
| A smoke | `field_presence` | 8.39 / 0.183 | 0.973 / 1 | - | 10.3 / 63 | 9.9 / 41.9 | 12.1 / 20 | 12.1 / 18.6 | 12 / 17.4 |
| A smoke | `mdl_scan` | 8.47 / 0.311 | 0.998 / 1.01 | 3.37 / 7.08 | 8.9 / 2.24 | 8.94 / 2.21 | 11.9 / 3.59 | 11.8 / 4.49 | 11.7 / 3.75 |
| A smoke | `molecule_count` | 8.24 / 0.159 | 0.989 / 0.938 | 3.34 / 12.5 | 10.4 / 71.5 | 9.87 / 50.6 | 12.5 / 21.9 | 11.6 / 6.83 | 11.8 / 5.41 |
| A smoke | `record_count_distribution` | 8.39 / 0.317 | 0.995 / 1.05 | 3.46 / 10.5 | 10.5 / 42.2 | 9.79 / 29 | 12.2 / 18.5 | 12 / 10.5 | 12.1 / 10.1 |
| A smoke | `record_dicts` | 8.3 / 0.267 | 0.999 / 1.1 | - | 10.7 / 50.1 | 9.82 / 34.5 | 12.6 / 22.8 | 12.6 / 12.8 | 12.4 / 12.3 |
| A smoke | `record_len_distribution` | 8.29 / 0.177 | 0.985 / 1.05 | - | 10.4 / 70.9 | 9.99 / 49.7 | 12.6 / 31.9 | 11.7 / 18.9 | 12 / 18.1 |
| A smoke | `record_name_set` | 8.28 / 0.214 | 0.987 / 1.04 | - | 10.5 / 58.8 | 10.1 / 39.8 | 12.9 / 26.9 | 12.3 / 15.8 | 12.1 / 15.1 |
| A smoke | `title_extraction` | 8.28 / 0.155 | 0.984 / 0.993 | 3.4 / 12.8 | 10.4 / 74.3 | 9.9 / 52.2 | 12.2 / 22.7 | 12 / 5.43 | 12.2 / 6.13 |
| A smoke | `write_round_trip` | 8.38 / 0.239 | 0.996 / 1.03 | 3.45 / 12.3 | 11 / 78 | 10.3 / 54.4 | 12.7 / 29.7 | 11.9 / 5.01 | 11.5 / 3.97 |
| A smoke | `write_sdf_round_trip` | 8.3 / 0.215 | 0.982 / 1.01 | - | 11.2 / 88.2 | 10.5 / 60.2 | 12.5 / 33.2 | 12 / 5.59 | 12.1 / 4.55 |
| A representative | `append_record_write` | 23.1 / 14.8 | 1.04 / 1.06 | 14 / 20.1 | 65.2 / 96.7 | 41.4 / 59.4 | 32.4 / 43.9 | 12.5 / 12.8 | 9.72 / 8.64 |
| A representative | `field_extraction` | 26.9 / 18.5 | 1.05 / 1.09 | 12.8 / 17.1 | 36.4 / 48.9 † | 22.2 / 28.2 | 22.7 / 27.6 | 14.8 / 16.1 | 12.4 / 12.7 |
| A representative | `field_get` | 17.6 / 9.42 | 1.11 / 1.2 | - | 46.6 / 79.1 † | 25.6 / 40.1 | 23.1 / 33 | 22.9 / 32 | 19.1 / 25.2 |
| A representative | `field_presence` | 16.5 / 8.42 | 1.06 / 1.11 | - | 49.6 / 88.5 † | 27.3 / 44.8 | 24.9 / 37.1 | 24 / 35.4 | 20.4 / 28 |
| A representative | `mdl_scan` | 22.3 / 14.1 | 1.04 / 1.05 | 10.3 / 14.4 | 5.09 / 2.8 | 5.1 / 2.84 | 8.34 / 6.42 | 9.29 / 7.75 | 7.24 / 4.52 |
| A representative | `molecule_count` | 13.6 / 5.57 | 1.04 / 1.09 | 15.5 / 33.1 | 58.8 / 130 † | 32.3 / 65.4 | 29.1 / 52.8 | 13.3 / 14.7 | 11.9 / 11.2 |
| A representative | `record_count_distribution` | 29 / 20.8 | 1.06 / 1.08 | 11.8 / 15.2 | 33.6 / 43.2 † | 20.6 / 25.1 | 21.1 / 24.7 | 13.5 / 14 | 11.3 / 10.9 |
| A representative | `record_dicts` | 23.8 / 15.7 | 1.08 / 1.11 | - | 41.6 / 58.1 † | 25.6 / 34 | 26.4 / 33.6 | 16.8 / 19 | 14 / 15 |
| A representative | `record_len_distribution` | 16.2 / 8.21 | 1.06 / 1.12 | - | 56 / 102 † | 32.7 / 55.9 | 35.8 / 59.2 | 24.5 / 36.1 | 20.3 / 28.4 |
| A representative | `record_name_set` | 19.9 / 11.8 | 1.06 / 1.1 | - | 46.1 / 71.3 | 27.2 / 39.5 | 29.9 / 42.1 | 20.1 / 25.5 | 16.8 / 20.1 |
| A representative | `title_extraction` | 14 / 5.92 | 1.04 / 1.09 | 15.1 / 31.4 | 58.1 / 125 † | 32.2 / 63.9 | 28.2 / 50.3 | 11.9 / 11.4 | 11.5 / 10.7 |
| A representative | `write_round_trip` | 19 / 10.8 | 1.05 / 1.07 | 15.9 / 25.3 | 78.6 / 131 | 49.8 / 80.3 | 38.4 / 58.2 | 9.78 / 8.02 | 8.54 / 5.82 |
| A representative | `write_sdf_round_trip` | 14.9 / 6.68 | 1.03 / 1.07 | - | 100 / 212 | 63.2 / 130 | 49.1 / 94.9 | 12.4 / 13 | 10.8 / 9.43 |
| A full | `append_record_write` | 236 / 227 | 1.07 / 1.07 | 20.2 / 20.9 | 84.6 / 87.7 | 54 / 55.9 | 39.8 / 41 | 12.5 / 12.5 | 7.65 / 7.52 |
| A full | `field_extraction` | 290 / 281 | 1.08 / 1.09 | 16.8 / 17.3 | 43.5 / 44.7 † | 26.2 / 26.7 | 25.6 / 26 | 15.8 / 16 | 11.9 / 11.9 |
| A full | `field_get` | 157 / 148 | 1.16 / 1.17 | - | 64.8 / 68.3 † | 34.7 / 36.3 | 27.9 / 28.9 | 29.5 / 30.6 | 22.4 / 23 |
| A full | `field_presence` | 141 / 132 | 1.1 / 1.11 | - | 72.2 / 76.6 † | 38.3 / 40.4 | 30.9 / 32.3 | 32.6 / 34 | 24.5 / 25.4 |
| A full | `mdl_scan` | 222 / 212 | 1.05 / 1.05 | 14.2 / 14.7 | 3.05 / 2.82 | 3.02 / 2.79 | 6.66 / 6.44 | 7.99 / 7.78 | 4.66 / 4.36 |
| A full | `molecule_count` | 97.4 / 88.3 | 1.08 / 1.09 | 29.3 / 32 | 101 / 111 † | 53.6 / 58.3 | 42.2 / 45.4 | 13.8 / 14 | 10.4 / 10.3 |
| A full | `record_count_distribution` | 324 / 315 | 1.08 / 1.08 | 15 / 15.4 | 38.8 / 39.7 † | 23.4 / 23.8 | 23 / 23.3 | 13.9 / 14 | 10.4 / 10.4 |
| A full | `record_dicts` | 245 / 236 | 1.1 / 1.11 | - | 51.9 / 53.6 † | 31.6 / 32.4 | 31.2 / 32 | 18.9 / 19.1 | 14.3 / 14.4 |
| A full | `record_len_distribution` | 139 / 130 | 1.1 / 1.12 | - | 82.7 / 88.2 † | 47.6 / 50.5 | 50.7 / 53.7 | 32.9 / 34.4 | 24.6 / 25.6 |
| A full | `record_name_set` | 192 / 183 | 1.1 / 1.11 | - | 60.3 / 62.9 | 35 / 36.4 | 37.3 / 38.7 | 24.2 / 24.8 | 18.1 / 18.5 |
| A full | `title_extraction` | 103 / 93.8 | 1.09 / 1.09 | 27.7 / 30 | 98.5 / 107 † | 52.7 / 56.9 | 40.6 / 43.4 | 10.9 / 10.8 | 10 / 9.84 |
| A full | `write_round_trip` | 180 / 170 | 1.06 / 1.06 | 23.4 / 24.5 | 109 / 115 | 69.7 / 73.1 | 50.6 / 52.9 | 7.91 / 7.66 | 5.71 / 5.4 |
| A full | `write_sdf_round_trip` | 114 / 105 | 1.07 / 1.07 | - | 173 / 188 | 110 / 119 | 80.3 / 86.4 | 12.4 / 12.5 | 9.11 / 8.89 |
| B smoke | `append_record_write` | 8.47 / 0.254 | 0.983 / 1.02 | 3.45 / 14.6 | 11 / 76.3 | 10.3 / 51.2 | 14.6 / 96.9 | 11.9 / 11 | 11.9 / 10.3 |
| B smoke | `field_extraction` | 8.37 / 0.331 | 0.997 / 1.07 | 3.64 / 14 | 10.6 / 42 | 10.1 / 28.8 | 13 / 22.9 | 12.6 / 14.3 | 12.5 / 13.8 |
| B smoke | `field_get` | 8.37 / 0.18 | 0.977 / 1.02 | - | 10.3 / 63.4 | 9.79 / 41 | 12.1 / 21.9 | 12.5 / 26.7 | 12.2 / 25.7 |
| B smoke | `field_presence` | 8.18 / 0.159 | 0.999 / 1.08 | - | 10.6 / 71.6 | 9.94 / 47.7 | 12.3 / 24.7 | 12.8 / 30.2 | 12.5 / 29.8 |
| B smoke | `mdl_scan` | 8.51 / 0.297 | 0.992 / 0.967 | 3.37 / 8.65 | 8.77 / 2.2 | 8.66 / 2.21 | 11.5 / 4.28 | 11.4 / 5.36 | 11.7 / 4.72 |
| B smoke | `molecule_count` | 8.2 / 0.115 | 0.999 / 1.02 | 3.46 / 20.5 | 10.5 / 97 | 10.1 / 62.1 | 12.9 / 32.3 | 12.4 / 11.1 | 12 / 9.22 |
| B smoke | `record_count_distribution` | 8.42 / 0.351 | 0.994 / 1.06 | 3.61 / 13.1 | 10.4 / 39.5 | 9.87 / 27.5 | 12.5 / 21.3 | 12.6 / 13.4 | 12.2 / 12.7 |
| B smoke | `record_dicts` | 8.3 / 0.291 | 1 / 1.08 | - | 10.5 / 49.2 | 9.9 / 34.6 | 12.4 / 26.4 | 12.3 / 16.5 | 12.7 / 15.8 |
| B smoke | `record_len_distribution` | 8.3 / 0.16 | 0.982 / 1.06 | - | 10.7 / 79.8 | 10.1 / 54 | 13.2 / 44.6 | 12.8 / 29.9 | 12.7 / 28.8 |
| B smoke | `record_name_set` | 8.51 / 0.235 | 0.966 / 1.01 | - | 10.3 / 54.3 | 9.74 / 36.7 | 12.4 / 30.9 | 12.2 / 20.5 | 12.2 / 19.8 |
| B smoke | `title_extraction` | 8.19 / 0.119 | 0.991 / 1.03 | 3.41 / 19.8 | 10.7 / 95.7 | 9.75 / 61.2 | 12.5 / 31.6 | 11.9 / 8.65 | 11.8 / 9.2 |
| B smoke | `write_round_trip` | 8.46 / 0.204 | 0.987 / 0.991 | 3.51 / 17.7 | 11 / 93.9 | 10.2 / 62.5 | 14.4 / 119 | 11.9 / 6.61 | 11.5 / 5.51 |
| B smoke | `write_sdf_round_trip` | 8.23 / 0.192 | 0.99 / 1.05 | - | 11.4 / 99.5 | 10.4 / 66.3 | 15 / 126 | 12 / 7 | 12.1 / 5.9 |
| B representative | `append_record_write` | 22.2 / 13.9 | 1.09 / 1.11 | - | 74.5 / 113 | 44.8 / 66.2 | 106 / 162 | 15.1 / 17.1 | 11.7 / 11.6 |
| B representative | `field_extraction` | 29.7 / 21.3 | 1.12 / 1.15 | - | error | error | 27 / 32.9 † | 18.3 / 20.7 | 15.5 / 17 |
| B representative | `field_get` | 16.4 / 8.28 | 1.14 / 1.23 | - | 51.9 / 94.1 † | 27.3 / 45.3 | 28 / 43.8 † | 33.3 / 54 | 28.3 / 44.3 |
| B representative | `field_presence` | 15.8 / 7.63 | 1.1 / 1.17 | - | 54.2 / 102 † | 28.4 / 49.2 | 29 / 47.6 † | 34.3 / 57.9 | 29.1 / 47.6 |
| B representative | `mdl_scan` | 22.2 / 13.9 | 1.05 / 1.07 | - | error | error | 9.44 / 8.1 | 10.5 / 9.58 | 8.28 / 6.1 |
| B representative | `molecule_count` | 11.3 / 3.24 | 1.06 / 1.18 | - | 73.3 / 234 † | 38.3 / 112 | 39.5 / 108 † | 17.8 / 31.2 | 15.9 / 24.8 |
| B representative | `record_count_distribution` | 31.3 / 22.8 | 1.08 / 1.12 | - | error | error | 25.5 / 30.6 † | 16.9 / 18.8 | 14.4 / 15.5 |
| B representative | `record_dicts` | 25.4 / 17.3 | 1.12 / 1.16 | - | error | error | 32 / 41.2 † | 21.3 / 25.6 | 18.3 / 21.1 |
| B representative | `record_len_distribution` | 15.5 / 7.29 | 1.11 / 1.23 | - | 63.2 / 124 † | 36.3 / 67.1 | 48.8 / 90.6 † | 34.6 / 59.8 | 29.6 / 49.4 |
| B representative | `record_name_set` | 21.3 / 13.2 | 1.09 / 1.14 | - | 46.5 / 69.7 | 26.8 / 37.7 | 36.1 / 50.9 | 25.6 / 33.8 | 21.9 / 27.9 |
| B representative | `title_extraction` | 11.6 / 3.54 | 1.07 / 1.18 | - | 73.2 / 219 † | 38.5 / 106 | 38.5 / 98.4 † | 16.2 / 24.4 | 15.3 / 22.4 |
| B representative | `write_round_trip` | 17 / 8.56 | 1.03 / 1.06 | - | 95.6 / 181 | 57.7 / 106 | 137 / 260 | 12.2 / 12.5 | 10.6 / 9.3 |
| B representative | `write_sdf_round_trip` | 11.8 / 3.72 | 1.08 / 1.2 | - | 137 / 416 | 82.9 / 243 | 197 / 598 | 17.3 / 28.8 | 15.1 / 21.5 |
| B full | `append_record_write` | 484 / 474 | 1.06 / 1.06 | - | 106 / 108 | 61.9 / 63 | 152 / 155 | 16.6 / 16.7 | 9.99 / 9.96 |
| B full | `field_extraction` | 710 / 701 | 1.13 / 1.13 | - | error | error | 32.2 / 32.5 † | 20.6 / 20.6 | 16.1 / 16.2 |
| B full | `field_get` | 302 / 292 | 1.22 / 1.22 | - | 83.8 / 86.2 † | 40.2 / 41.3 | 39.4 / 40.4 † | 48.6 / 49.8 | 38.2 / 39.1 |
| B full | `field_presence` | 274 / 265 | 1.13 / 1.14 | - | 92.2 / 95.1 † | 44.3 / 45.6 | 43.4 / 44.5 † | 53 / 54.4 | 41.6 / 42.7 |
| B full | `mdl_scan` | 477 / 468 | 1.04 / 1.04 | - | error | error | 7.99 / 7.93 | 9.39 / 9.3 | 5.57 / 5.46 |
| B full | `molecule_count` | 133 / 124 | 1.08 / 1.08 | - | 184 / 197 † | 87 / 92.7 | 84.6 / 90 † | 25.5 / 26.4 | 19.6 / 20.2 |
| B full | `record_count_distribution` | 770 / 760 | 1.1 / 1.1 | - | error | error | 29.9 / 30.1 † | 18.6 / 18.7 | 14.4 / 14.4 |
| B full | `record_dicts` | 595 / 586 | 1.12 / 1.12 | - | error | error | 39.3 / 39.7 † | 24.5 / 24.6 | 19.3 / 19.4 |
| B full | `record_len_distribution` | 265 / 255 | 1.17 / 1.17 | - | 111 / 115 † | 60 / 61.9 | 81.9 / 84.5 † | 54.3 / 55.8 | 42.5 / 43.7 |
| B full | `record_name_set` | 455 / 446 | 1.12 / 1.12 | - | 65.3 / 66.5 | 35.4 / 35.9 | 48.3 / 49.1 | 32.3 / 32.7 | 25.4 / 25.8 |
| B full | `title_extraction` | 144 / 135 | 1.11 / 1.11 | - | 175 / 186 † | 83.9 / 89 | 79.6 / 84.2 † | 20.2 / 20.7 | 18.2 / 18.7 |
| B full | `write_round_trip` | 300 / 291 | 1.04 / 1.04 | - | 168 / 173 | 97.8 / 101 | 242 / 249 | 12 / 12 | 8.68 / 8.63 |
| B full | `write_sdf_round_trip` | 146 / 137 | 1.07 / 1.07 | - | 345 / 367 | 201 / 213 | 495 / 526 | 24.5 / 25.3 | 17.8 / 18.2 |
| C smoke | `append_record_write` | 8.47 / 0.24 | 1.01 / 1.06 | 3.1 / 8.65 | 9.98 / 67.2 | 9.47 / 46.8 | 11.2 / 24.6 | 11 / 7.84 | 10.6 / 6.85 |
| C smoke | `field_extraction` | 8.84 / 0.236 | 0.952 / 0.978 | 2.9 / 8.57 | 9.02 / 46.4 | 8.67 / 31.7 | 10.7 / 18.4 | 10.3 / 10.1 | 10.5 / 8.95 |
| C smoke | `field_get` | 8.44 / 0.148 | 0.976 / 1.07 | - | 9.45 / 66.7 | 9.38 / 44.3 | 10.9 / 20.2 | 10.8 / 16.3 | 10.8 / 14.6 |
| C smoke | `field_presence` | 8.81 / 0.15 | 0.94 / 0.95 | - | 9.09 / 65.8 | 8.64 / 43.5 | 10.5 / 19.9 | 10.6 / 16.2 | 10.3 / 14.4 |
| C smoke | `mdl_scan` | 8.69 / 0.274 | 0.982 / 0.999 | 2.93 / 5.43 | 8.03 / 2.1 | 8.01 / 2.12 | 10.4 / 3.29 | 10.6 / 4.3 | 10.4 / 3.41 |
| C smoke | `molecule_count` | 8.33 / 0.115 | 1.02 / 1.06 | 3.03 / 11.2 | 9.51 / 83.8 | 9.01 / 55.2 | 11.1 / 24.2 | 10.9 / 7.82 | 11 / 5.92 |
| C smoke | `record_count_distribution` | 8.46 / 0.246 | 0.99 / 1.02 | 3.04 / 8.42 | 9.39 / 44.5 | 8.89 / 30.4 | 11.2 / 17.4 | 10.8 / 9.65 | 10.9 / 8.6 |
| C smoke | `record_dicts` | 8.37 / 0.2 | 1 / 1.05 | - | 9.59 / 55.6 | 9.18 / 37.8 | 11.2 / 21.9 | 11.1 / 12.2 | 10.8 / 10.8 |
| C smoke | `record_len_distribution` | 8.67 / 0.145 | 0.949 / 0.999 | - | 9.27 / 72.7 | 8.94 / 48.9 | 10.8 / 28.3 | 10.7 / 16.6 | 10.7 / 14.6 |
| C smoke | `record_name_set` | 8.9 / 0.182 | 0.926 / 0.89 | - | 9.02 / 58.1 | 8.69 / 39.1 | 10.6 / 22.8 | 10.5 / 13.5 | 10.4 / 12 |
| C smoke | `title_extraction` | 8.3 / 0.116 | 0.983 / 0.989 | 3.03 / 11.2 | 9.48 / 84.3 | 9.11 / 55.5 | 11.1 / 24.6 | 10.9 / 5.85 | 11 / 6.17 |
| C smoke | `write_round_trip` | 8.72 / 0.204 | 0.965 / 0.938 | 2.99 / 9.74 | 9.8 / 78.1 | 9.22 / 54.6 | 10.9 / 28.1 | 10.5 / 4.79 | 10.3 / 3.61 |
| C smoke | `write_sdf_round_trip` | 8.26 / 0.172 | 1.02 / 1.11 | - | 10.2 / 92.6 | 9.81 / 65 | 11.8 / 34.9 | 11.4 / 5.74 | 10.9 / 4.31 |
| C representative | `append_record_write` | 21.4 / 13 | 1.06 / 1.08 | 11.1 / 16.2 | 59.5 / 92.3 | 36.7 / 54.7 | 29.4 / 41.1 | 11.7 / 12.1 | 9.2 / 8.03 |
| C representative | `field_extraction` | 20.5 / 12 | 1.06 / 1.1 | 11.2 / 17 | 37.4 / 58.1 | 21.9 / 31.6 | 23.1 / 31.8 | 14.7 / 17.4 | 12.1 / 12.9 |
| C representative | `field_get` | 14.8 / 6.42 | 1.11 / 1.24 | - | 44.8 / 91.9 | 24.2 / 44.8 | 23.8 / 40.2 | 20.6 / 32.8 | 16.9 / 24.6 |
| C representative | `field_presence` | 13.7 / 5.34 | 1.07 / 1.16 | - | 48.4 / 110 | 26.3 / 53.9 | 25.8 / 48.5 | 22.2 / 39.1 | 18.2 / 29.2 |
| C representative | `mdl_scan` | 19.8 / 11.4 | 1.04 / 1.05 | 8.63 / 12.7 | 5.44 / 3.3 | 5.4 / 3.31 | 8.46 / 6.7 | 9.45 / 8.27 | 7.35 / 4.75 |
| C representative | `molecule_count` | 11.5 / 3.23 | 1.05 / 1.14 | 13.4 / 39.4 | 56.1 / 177 | 30 / 84.6 | 28.8 / 74.1 | 14.1 / 21.2 | 12.3 / 15.3 |
| C representative | `record_count_distribution` | 22.3 / 13.8 | 1.05 / 1.1 | 10.2 / 14.7 | 34.5 / 50.7 | 20.2 / 27.6 | 21 / 27.2 | 13.2 / 14.9 | 11 / 11 |
| C representative | `record_dicts` | 18.1 / 9.9 | 1.08 / 1.13 | - | 42.8 / 71.1 | 25.2 / 38.9 | 26.2 / 38.6 | 16.7 / 21.2 | 13.8 / 15.9 |
| C representative | `record_len_distribution` | 13.6 / 5.22 | 1.04 / 1.12 | - | 53.2 / 125 | 30.4 / 65.4 | 33.1 / 68.2 | 22.1 / 39.8 | 18.2 / 29.6 |
| C representative | `record_name_set` | 15.6 / 7.21 | 1.07 / 1.14 | - | 46.5 / 90.9 | 26.6 / 47.7 | 29.1 / 50.2 | 19.4 / 29.1 | 16 / 21.9 |
| C representative | `title_extraction` | 11.9 / 3.62 | 1.05 / 1.15 | 12.9 / 35.5 | 55.4 / 163 | 29.9 / 78.5 | 28.4 / 67.7 | 12.4 / 15.1 | 11.8 / 13.8 |
| C representative | `write_round_trip` | 16.7 / 8.17 | 1.05 / 1.08 | 13.3 / 24 | 75.2 / 145 | 46.4 / 85.8 | 36.5 / 62.9 | 9.93 / 8.89 | 8.57 / 6.18 |
| C representative | `write_sdf_round_trip` | 12.2 / 4.07 | 1.06 / 1.15 | - | 103 / 291 | 63.3 / 172 | 50.2 / 127 | 13.5 / 17.9 | 11.7 / 12.4 |
| C full | `append_record_write` | 23.5 / 15 | 1.04 / 1.06 | 11.3 / 16.1 | 60.5 / 90.2 | 37.4 / 53.9 | 29.6 / 40 | 11.4 / 11.7 | 8.91 / 7.77 |
| C full | `field_extraction` | 22.2 / 13.6 | 1.06 / 1.09 | 11.4 / 16.7 | 38.9 / 58.2 | 22.5 / 31.5 | 23.1 / 31 | 14.9 / 17.4 | 12 / 12.8 |
| C full | `field_get` | 15.9 / 7.4 | 1.09 / 1.21 | - | 46.5 / 90.3 | 25 / 44 | 24 / 39.2 | 20.8 / 32 | 16.9 / 23.9 |
| C full | `field_presence` | 14.6 / 6.21 | 1.04 / 1.08 | - | 50.7 / 108 | 27.1 / 52.2 | 26.3 / 46.8 | 22.6 / 38 | 18.4 / 28.2 |
| C full | `mdl_scan` | 21.9 / 13.2 | 1.02 / 1.04 | 8.71 / 12.4 | 5.15 / 3.22 | 5.23 / 3.25 | 8.16 / 6.59 | 9.22 / 8.1 | 6.77 / 4.4 |
| C full | `molecule_count` | 12.3 / 3.86 | 0.996 / 1.02 | 13.9 / 37.6 | 58.5 / 168 | 30.9 / 80.2 | 29.7 / 70 | 13.8 / 20.1 | 12 / 14.6 |
| C full | `record_count_distribution` | 24.4 / 15.7 | 1.03 / 1.07 | 10.5 / 14.6 | 35.4 / 50.4 | 20.7 / 27.4 | 21 / 26.8 | 13.3 / 14.6 | 10.8 / 10.9 |
| C full | `record_dicts` | 20.1 / 11.5 | 1.04 / 1.08 | - | 43.3 / 69.5 | 25.2 / 37.9 | 25.8 / 37.2 | 16.4 / 20.6 | 13.4 / 15.3 |
| C full | `record_len_distribution` | 14.3 / 5.93 | 1.05 / 1.11 | - | 56.9 / 125 | 32 / 65 | 34.5 / 67.4 | 23 / 39.6 | 18.7 / 29.1 |
| C full | `record_name_set` | 16.7 / 8.32 | 1.05 / 1.08 | - | 48.7 / 89.3 | 27.7 / 46.9 | 29.8 / 48.8 | 19.8 / 28.6 | 16.2 / 21.3 |
| C full | `title_extraction` | 12.6 / 4.28 | 1.02 / 1.06 | 13.5 / 34 | 58.6 / 156 | 31.3 / 75.8 | 29.2 / 64.4 | 12.4 / 14.5 | 11.8 / 13.2 |
| C full | `write_round_trip` | 18 / 9.5 | 1.01 / 1.02 | 13.6 / 23.1 | 77.9 / 140 | 47.8 / 83.3 | 37.3 / 61.2 | 9.69 / 8.67 | 8.2 / 5.96 |
| C full | `write_sdf_round_trip` | 13.4 / 4.92 | 1.01 / 1.06 | - | 105 / 271 | 64.3 / 161 | 50.3 / 118 | 13 / 16.8 | 11 / 11.5 |
| D representative | `append_record_write` | 19 / 10.6 | 1.03 / 1.04 | 45 / 78.3 | 93.6 / 161 | 67.1 / 114 | 56.6 / 92.9 | 27.3 / 40 | 26.6 / 38.9 |
| D representative | `field_extraction` | 77.8 / 69.2 | 1.07 / 1.08 | 16.7 / 18.4 | 22.7 / 24.5 † | 18.4 / 19.7 | 21.5 / 22.9 | 16.6 / 17.3 | 16.5 / 17.2 |
| D representative | `field_get` | 22.5 / 14.2 | 1.08 / 1.11 | - | 42.1 / 61.8 † | 28.5 / 40.1 | 25.6 / 34.1 | 58.3 / 85.9 | 57.8 / 85.2 |
| D representative | `field_presence` | 22.5 / 14.3 | 1.06 / 1.09 | - | 42.1 / 61.4 † | 28.4 / 39.8 | 25.9 / 34.4 | 58.5 / 85.6 | 57.7 / 84.5 |
| D representative | `mdl_scan` | 18 / 9.53 | 1.02 / 1.03 | 31.1 / 55.9 | 5.36 / 2.83 | 5.42 / 2.84 | 14.8 / 18.1 | 15.6 / 19.6 | 15.5 / 19.3 |
| D representative | `molecule_count` | 12.3 / 3.86 | 0.989 / 1.02 | 44.4 / 135 | 76.6 / 226 † | 51.5 / 146 | 46.4 / 124 | 21.2 / 43.7 | 21.1 / 42.8 |
| D representative | `record_count_distribution` | 80.2 / 71.5 | 1.06 / 1.07 | 16 / 17.6 | 21.9 / 23.5 † | 17.7 / 18.9 † | 20.8 / 22 | 15.9 / 16.5 | 15.8 / 16.4 |
| D representative | `record_dicts` | 69.4 / 61 | 1.09 / 1.1 | - | 26.1 / 28.5 † | 21.1 / 22.8 † | 24.6 / 26.4 | 19 / 20.1 | 18.8 / 19.8 |
| D representative | `record_len_distribution` | 22.2 / 13.9 | 1.09 / 1.14 | - | 58.5 / 88.7 † | 43.8 / 65.1 † | 68.3 / 103 | 57.9 / 86.1 | 57.4 / 85.2 |
| D representative | `record_name_set` | 40.8 / 32.4 | 1.1 / 1.12 | - | 32.4 / 38.6 † | 24.6 / 28.7 † | 37.9 / 44.8 | 32.4 / 37.9 | 32.2 / 37.5 |
| D representative | `title_extraction` | 12.2 / 3.95 | 1.01 / 1.05 | 44.7 / 132 | 77.6 / 222 † | 52.4 / 144 | 47.7 / 124 | 21.2 / 41.6 | 21.2 / 42.1 |
| D representative | `write_round_trip` | 16.8 / 8.4 | 1.01 / 1.01 | 49.4 / 95.5 | 105 / 202 | 75.4 / 142 | 63.7 / 116 | 15.7 / 20.2 | 15.5 / 19.8 |
| D representative | `write_sdf_round_trip` | 14 / 5.74 | 1.02 / 1.02 | - | 126 / 294 | 90.7 / 209 | 76.2 / 170 | 18.9 / 29.6 | 18.6 / 29 |
| D full | `append_record_write` | 22.2 / 13.6 | 0.983 / 0.991 | 52.8 / 84.3 | 99.4 / 157 | 76.5 / 120 | 60.6 / 92.1 | 30.6 / 43.1 | 30 / 42 |
| D full | `field_extraction` | 111 / 102 | 1.05 / 1.06 | 16.5 / 17.7 | 21.8 / 23.1 † | 18.7 / 19.6 | 20.7 / 21.6 | 16.2 / 16.7 | 16 / 16.5 |
| D full | `field_get` | 29.6 / 20.8 | 1.04 / 1.07 | - | 39.9 / 53.3 † | 29.1 / 37.9 | 24.3 / 30 | 61.8 / 83.3 | 61.2 / 82.4 |
| D full | `field_presence` | 29.6 / 21 | 1.04 / 1.05 | - | 40 / 52.9 † | 29.1 / 37.6 | 24 / 29.5 | 61.8 / 82.6 | 61.3 / 81.8 |
| D full | `mdl_scan` | 20 / 11.5 | 0.981 / 0.977 | 37.6 / 63.1 | 5.07 / 2.73 | 5.14 / 2.75 | 16 / 19.9 | 17 / 21.4 | 16.9 / 21.1 |
| D full | `molecule_count` | 15.2 / 6.34 | 0.912 / 0.877 | 49.4 / 115 | 77.2 / 174 † | 56.6 / 124 | 46.6 / 97.2 | 21.2 / 36.3 | 21 / 35.4 |
| D full | `record_count_distribution` | 114 / 105 | 1.04 / 1.05 | 15.9 / 17 | 20.8 / 21.9 † | 17.9 / 18.7 | 19.8 / 20.6 | 15.5 / 15.9 | 15.3 / 15.7 |
| D full | `record_dicts` | 101 / 92.3 | 1.06 / 1.07 | - | 24.5 / 26.1 † | 20.9 / 22.2 | 23.1 / 24.3 | 18.1 / 18.8 | 18 / 18.7 |
| D full | `record_len_distribution` | 28.4 / 19.9 | 1.09 / 1.13 | - | 58.9 / 80.6 † | 47.7 / 64.4 | 72.6 / 99 | 63.1 / 85.3 | 62.3 / 84.2 |
| D full | `record_name_set` | 56.5 / 47.6 | 1.07 / 1.09 | - | 30.5 / 34.7 | 24.7 / 27.8 | 37.4 / 42.5 | 32.5 / 36.6 | 32.2 / 36.2 |
| D full | `title_extraction` | 15.5 / 6.54 | 0.897 / 0.852 | 48.3 / 110 | 76.3 / 170 † | 55.7 / 121 | 45.8 / 93.9 | 20.5 / 33.9 | 20.7 / 34.6 |
| D full | `write_round_trip` | 19.9 / 11.3 | 0.99 / 0.979 | 57.4 / 98.5 | 110 / 188 | 84.9 / 143 | 67.4 / 110 | 16.2 / 20.3 | 16.1 / 20 |
| D full | `write_sdf_round_trip` | 16.9 / 8.18 | 0.962 / 0.957 | - | 130 / 261 | 100 / 198 | 79.5 / 153 | 19.1 / 28.2 | 18.9 / 27.6 |

### arm64 (c7g.2xlarge)

| Dataset | Case | sendoff ms | sendoff plain | 0.1.6 | RDKit default | RDKit fast | OEChem default | OEChem tokenizer | OEChem threaded |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| A smoke | `append_record_write` | 10.2 / 0.332 | 1.01 / 1.02 | 3.05 / 11.6 | 9.82 / 60 | 9.14 / 42.1 | 10.7 / 27.5 | 9.97 / 6.89 | 10.1 / 6.7 |
| A smoke | `field_extraction` | 10.1 / 0.36 | 1.05 / 1.07 | 3.13 / 11.4 | 9.32 / 41 | 8.85 / 28.5 | 10.5 / 17.4 | 10.2 / 10.1 | 10.2 / 9.96 |
| A smoke | `field_get` | 9.76 / 0.237 | 1.05 / 1.06 | - | 9.47 / 53.5 | 8.97 / 36.3 | 10.5 / 15.9 | 10.6 / 15.6 | 10.4 / 15.1 |
| A smoke | `field_presence` | 9.79 / 0.222 | 1.06 / 1.09 | - | 9.35 / 57.7 | 8.94 / 39.3 | 10.5 / 17 | 10.5 / 16.7 | 10.4 / 16.1 |
| A smoke | `mdl_scan` | 10.3 / 0.352 | 1.03 / 1.02 | 2.91 / 7.67 | 7.74 / 2.19 | 7.75 / 2.2 | 9.65 / 2.12 | 9.76 / 2.87 | 9.9 / 2.68 |
| A smoke | `molecule_count` | 9.58 / 0.172 | 1.04 / 1.04 | 3.08 / 14.1 | 9.84 / 72 | 9.41 / 48.4 | 11.2 / 20.9 | 10.9 / 5.04 | 10.5 / 4.63 |
| A smoke | `record_count_distribution` | 9.79 / 0.379 | 1.05 / 1.08 | 3.22 / 10.8 | 9.6 / 38.8 | 9.08 / 26.9 | 10.8 / 16.5 | 10.4 / 9.53 | 10.5 / 9.26 |
| A smoke | `record_dicts` | 10.3 / 0.333 | 0.997 / 1.04 | - | 9.08 / 44.3 | 8.69 / 31.3 | 10.3 / 19.1 | 9.93 / 11.2 | 9.97 / 10.8 |
| A smoke | `record_len_distribution` | 9.84 / 0.222 | 1.03 / 1.08 | - | 9.39 / 62.5 | 9.18 / 43.1 | 10.6 / 26.8 | 10.4 / 16.6 | 10.5 / 16.1 |
| A smoke | `record_name_set` | 9.82 / 0.265 | 1.06 / 1.09 | - | 9.51 / 52.3 | 9.14 / 36.8 | 10.7 / 22.7 | 10.5 / 14.1 | 10.5 / 13.6 |
| A smoke | `title_extraction` | 9.84 / 0.184 | 1.01 / 1.02 | 3.01 / 13.2 | 9.27 / 68.5 | 8.99 / 47.3 | 10.4 / 19.6 | 10.2 / 3.25 | 10.2 / 4.5 |
| A smoke | `write_round_trip` | 9.95 / 0.28 | 1.03 / 1.04 | 3.1 / 13.4 | 9.86 / 69.9 | 9.28 / 48.9 | 10.8 / 31.1 | 9.97 / 3.37 | 9.97 / 3.08 |
| A smoke | `write_sdf_round_trip` | 9.72 / 0.265 | 1.02 / 1.02 | - | 10.2 / 74.4 | 9.56 / 51.6 | 11.2 / 33 | 10.3 / 3.63 | 10.3 / 3.23 |
| A representative | `append_record_write` | 30.1 / 19.5 | 1.01 / 1.03 | 13.4 / 19.2 | 52.8 / 77.4 | 33.7 / 48 | 31 / 42.7 | 9.57 / 9.62 | 8.1 / 7.42 |
| A representative | `field_extraction` | 34.8 / 24.6 | 1.05 / 1.07 | 12.1 / 16 | 31.8 / 41.8 † | 20.1 / 25.4 | 18.9 / 22.6 | 12.3 / 13.4 | 11 / 11.6 |
| A representative | `field_get` | 24.2 / 14.1 | 1.04 / 1.07 | - | 37.2 / 58.5 † | 21.4 / 31.1 | 17.4 / 22.9 | 17.9 / 23.7 | 16 / 20.4 |
| A representative | `field_presence` | 23 / 12.8 | 1.03 / 1.08 | - | 39.3 / 64.7 † | 22.4 / 34.3 | 18.4 / 25.2 | 18.6 / 25.9 | 16.7 / 22.3 |
| A representative | `mdl_scan` | 27.8 / 17.3 | 0.996 / 1.01 | 10.1 / 14.6 | 4.46 / 2.69 | 4.47 / 2.7 | 5.59 / 3.37 | 6.2 / 4.3 | 5.16 / 2.62 |
| A representative | `molecule_count` | 18.4 / 8.34 | 1.02 / 1.03 | 14.1 / 27.8 | 47.2 / 94.4 † | 27.3 / 50.7 | 21.8 / 36.3 | 8.59 / 7.12 | 7.76 / 5.38 |
| A representative | `record_count_distribution` | 37.6 / 27.2 | 1.05 / 1.07 | 11.2 / 14.5 | 29.2 / 37.4 † | 18.4 / 22.5 | 17.3 / 20.3 | 11.2 / 11.9 | 10 / 10.2 |
| A representative | `record_dicts` | 31.7 / 21.4 | 1.04 / 1.07 | - | 35.2 / 48.3 † | 21.7 / 28.5 | 21 / 26.5 | 13.6 / 15.5 | 12.1 / 13.3 |
| A representative | `record_len_distribution` | 22.8 / 12.6 | 1.03 / 1.07 | - | 44.3 / 74 † | 27.1 / 42.9 | 27.1 / 41.3 | 18.7 / 26 | 16.7 / 22.4 |
| A representative | `record_name_set` | 27.4 / 17.1 | 1.04 / 1.07 | - | 37.3 / 55.1 | 22.7 / 31.9 | 22.9 / 30.9 | 15.8 / 19.5 | 14.1 / 16.8 |
| A representative | `title_extraction` | 19.1 / 8.78 | 1 / 1.03 | 13.6 / 26.3 | 47.1 / 93.4 † | 26.9 / 49.4 | 21.2 / 34.9 | 7.47 / 5.05 | 7.44 / 5.07 |
| A representative | `write_round_trip` | 24.7 / 14.6 | 1.02 / 1.02 | 15.6 / 24.5 | 63.1 / 102 | 40.1 / 62.6 | 36.3 / 54.6 | 6.66 / 4.53 | 5.8 / 3.03 |
| A representative | `write_sdf_round_trip` | 19.8 / 9.79 | 1.03 / 1.03 | - | 78.7 / 151 | 50.1 / 93.3 | 45.3 / 81.4 | 8.41 / 6.73 | 7.37 / 4.64 |
| A full | `append_record_write` | 296 / 285 | 1.03 / 1.03 | 20 / 20.6 | 74.4 / 76.8 | 46.6 / 48 | 42.3 / 43.5 | 9.95 / 9.89 | 7.16 / 7.04 |
| A full | `field_extraction` | 373 / 363 | 1.07 / 1.08 | 16.3 / 16.7 | 39.4 / 40.2 † | 24.6 / 25.1 | 22.4 / 22.8 | 13.7 / 13.7 | 11.4 / 11.4 |
| A full | `field_get` | 218 / 208 | 1.06 / 1.07 | - | 54.1 / 56.3 † | 29.5 / 30.5 | 22.1 / 22.7 | 23.6 / 24.1 | 19.7 / 20.1 |
| A full | `field_presence` | 199 / 189 | 1.07 / 1.07 | - | 60.1 / 63 † | 32.7 / 34.1 | 24.2 / 24.9 | 25.7 / 26.4 | 21.4 / 22 |
| A full | `mdl_scan` | 265 / 255 | 1.02 / 1.02 | 14.7 / 15.2 | 2.98 / 2.76 | 2.98 / 2.76 | 3.75 / 3.45 | 4.73 / 4.41 | 2.91 / 2.57 |
| A full | `molecule_count` | 132 / 122 | 1.02 / 1.02 | 27.1 / 29.2 | 87.2 / 93.8 † | 45.9 / 49 | 34.3 / 36.3 | 7.66 / 7.28 | 6.02 / 5.6 |
| A full | `record_count_distribution` | 413 / 403 | 1.07 / 1.07 | 14.7 / 15 | 35.5 / 36.2 † | 22.1 / 22.4 | 20.3 / 20.5 | 12.1 / 12.1 | 10.1 / 10.1 |
| A full | `record_dicts` | 323 / 313 | 1.08 / 1.08 | - | 46 / 47.3 † | 28.8 / 29.5 | 26.3 / 26.8 | 15.9 / 16 | 13.3 / 13.4 |
| A full | `record_len_distribution` | 197 / 186 | 1.06 / 1.07 | - | 68.4 / 71.7 † | 39.5 / 41.2 | 40.1 / 41.6 | 25.8 / 26.6 | 21.6 / 22.1 |
| A full | `record_name_set` | 264 / 254 | 1.07 / 1.07 | - | 51.9 / 53.7 | 31.5 / 32.5 | 30.3 / 31.1 | 19.6 / 19.9 | 16.4 / 16.6 |
| A full | `title_extraction` | 140 / 130 | 1.02 / 1.02 | 25.6 / 27.4 | 85.6 / 91.7 † | 46.1 / 49.1 | 32.9 / 34.7 | 5.66 / 5.15 | 5.68 / 5.24 |
| A full | `write_round_trip` | 226 / 216 | 1.01 / 1.01 | 24.3 / 25.4 | 95.1 / 99.4 | 59.6 / 62.1 | 53.3 / 55.4 | 4.93 / 4.6 | 3.47 / 3.11 |
| A full | `write_sdf_round_trip` | 153 / 143 | 1.02 / 1.02 | - | 141 / 150 | 87.9 / 93.5 | 78.6 / 83.4 | 7.29 / 6.92 | 5.12 / 4.69 |
| B smoke | `append_record_write` | 10.2 / 0.309 | 1.03 / 1.04 | 3.16 / 15.6 | 9.8 / 65.8 | 9.14 / 44.7 | 12.4 / 86.9 | 10.1 / 8.85 | 10 / 8.54 |
| B smoke | `field_extraction` | 10.2 / 0.403 | 1.02 / 1.08 | 3.26 / 14.3 | 9.28 / 38.1 | 8.83 / 26.9 | 10.6 / 20.8 | 10.2 / 13.1 | 10.3 / 12.9 |
| B smoke | `field_get` | 10.2 / 0.231 | 1.01 / 1.06 | - | 8.92 / 54.5 | 8.53 / 36 | 10.2 / 18.3 | 10.2 / 23.2 | 10.2 / 22.7 |
| B smoke | `field_presence` | 9.95 / 0.21 | 1.02 / 1.09 | - | 9.26 / 60.2 | 8.77 / 39.1 | 10.4 / 20.1 | 10.5 / 25.5 | 10.4 / 25 |
| B smoke | `mdl_scan` | 10.6 / 0.332 | 0.973 / 1.01 | 2.89 / 9.57 | 7.62 / 2.24 | 7.54 / 2.18 | 9.45 / 2.5 | 9.48 / 3.34 | 9.47 / 3.11 |
| B smoke | `molecule_count` | 9.78 / 0.141 | 0.996 / 1.01 | 3.07 / 20.9 | 9.36 / 87.1 | 8.98 / 57.9 | 10.7 / 28.5 | 10.4 / 6.47 | 10.2 / 6.13 |
| B smoke | `record_count_distribution` | 10.4 / 0.433 | 0.983 / 1.07 | 3.19 / 13.3 | 9.12 / 35.6 | 8.72 / 25 | 10.5 / 19.3 | 10.1 / 12 | 10.1 / 11.9 |
| B smoke | `record_dicts` | 10.3 / 0.359 | 0.986 / 1.08 | - | 9.18 / 43.1 | 8.7 / 30.6 | 10.4 / 23.5 | 10.2 / 14.9 | 10.1 / 14.6 |
| B smoke | `record_len_distribution` | 10.3 / 0.214 | 0.971 / 1.04 | - | 9.07 / 65.6 | 8.66 / 44.3 | 10.7 / 36.5 | 10.1 / 24.8 | 10.2 / 24.4 |
| B smoke | `record_name_set` | 10.3 / 0.288 | 0.993 / 1.05 | - | 9.05 / 49.2 | 8.6 / 33.9 | 10.4 / 27.2 | 10.2 / 18.7 | 10.1 / 18.3 |
| B smoke | `title_extraction` | 9.87 / 0.144 | 0.987 / 1.02 | 3.07 / 20.5 | 9.25 / 86.9 | 8.88 / 57.3 | 10.4 / 28.3 | 10.1 / 4.74 | 10.1 / 6.28 |
| B smoke | `write_round_trip` | 10.3 / 0.241 | 1 / 1.02 | 3.1 / 19.7 | 9.69 / 84.4 | 9.03 / 56.8 | 12.3 / 110 | 9.76 / 4.2 | 9.7 / 3.98 |
| B smoke | `write_sdf_round_trip` | 10.2 / 0.258 | 0.975 / 0.976 | - | 9.87 / 79.1 | 9.1 / 52.4 | 12.4 / 104 | 9.81 / 3.94 | 9.82 / 3.68 |
| B representative | `append_record_write` | 29 / 18.6 | 1.01 / 1.04 | - | 59.7 / 88.9 | 36.2 / 52.3 | 88.2 / 132 | 11.6 / 12.6 | 9.79 / 9.88 |
| B representative | `field_extraction` | 38.8 / 28.4 | 1.03 / 1.06 | - | error | error | 23 / 27.7 † | 15.5 / 17.5 | 14.1 / 15.6 |
| B representative | `field_get` | 23.2 / 12.9 | 1.05 / 1.11 | - | 40.6 / 66.6 † | 22.3 / 33.9 | 21.2 / 30.2 † | 26 / 38.8 | 23.7 / 34.7 |
| B representative | `field_presence` | 21.9 / 11.8 | 1.06 / 1.13 | - | 43.6 / 73.9 † | 24.1 / 37.8 | 22.5 / 33 † | 27.6 / 42.3 | 25 / 37.6 |
| B representative | `mdl_scan` | 25.8 / 15.5 | 1.03 / 1.05 | - | error | error | 6.62 / 4.43 | 7.36 / 5.56 | 6.06 / 3.53 |
| B representative | `molecule_count` | 15.3 / 5.18 | 1.01 / 1.06 | - | 59.6 / 160 † | 33 / 81.6 | 31 / 71.7 † | 11.3 / 13.4 | 10.1 / 10.1 |
| B representative | `record_count_distribution` | 40.7 / 30.5 | 1.05 / 1.07 | - | error | error | 21.8 / 25.7 † | 14.6 / 16 | 13.2 / 14.2 |
| B representative | `record_dicts` | 33.9 / 23.7 | 1.08 / 1.12 | - | error | error | 26.6 / 33.8 † | 17.8 / 21.2 | 16.2 / 18.9 |
| B representative | `record_len_distribution` | 22 / 11.7 | 1.05 / 1.12 | - | 49.7 / 86.6 † | 29.2 / 48.1 | 37.8 / 62.2 † | 27.2 / 42.4 | 24.7 / 37.7 |
| B representative | `record_name_set` | 29 / 18.8 | 1.05 / 1.09 | - | 38 / 54.3 | 23.1 / 31.4 | 29.1 / 39.4 | 21 / 26.9 | 19.1 / 23.9 |
| B representative | `title_extraction` | 15.5 / 5.59 | 1.03 / 1.07 | - | 60.8 / 155 † | 32.9 / 77.2 | 30.8 / 67.1 † | 10.1 / 9.67 | 10.1 / 9.58 |
| B representative | `write_round_trip` | 21.9 / 11.7 | 0.997 / 1.02 | - | 78.8 / 141 | 47 / 81.1 | 115 / 207 | 8.2 / 6.61 | 7.06 / 4.55 |
| B representative | `write_sdf_round_trip` | 16 / 5.66 | 0.98 / 1.04 | - | 106 / 286 | 64.7 / 169 | 158 / 428 | 11.3 / 13.6 | 9.63 / 9.33 |
| B full | `append_record_write` | 609 / 599 | 1.05 / 1.05 | - | 90.1 / 91.5 | 52.4 / 53.1 | 133 / 135 | 13 / 13 | 9.3 / 9.28 |
| B full | `field_extraction` | 928 / 918 | 1.07 / 1.07 | - | error | error | 27.6 / 27.8 † | 17.9 / 17.9 | 15.4 / 15.4 |
| B full | `field_get` | 428 / 417 | 1.11 / 1.11 | - | 66.5 / 67.9 † | 33.6 / 34.3 | 30 / 30.5 † | 39.3 / 39.9 | 33.7 / 34.3 |
| B full | `field_presence` | 393 / 383 | 1.13 / 1.13 | - | 72 / 73.7 † | 36.6 / 37.4 | 32.7 / 33.3 † | 42.4 / 43.2 | 36.3 / 37 |
| B full | `mdl_scan` | 507 / 496 | 1.04 / 1.04 | - | error | error | 4.68 / 4.55 | 5.85 / 5.68 | 3.5 / 3.35 |
| B full | `molecule_count` | 174 / 164 | 1.07 / 1.07 | - | 156 / 165 † | 77.3 / 81.5 | 70.2 / 73.8 † | 13.9 / 13.8 | 10.6 / 10.5 |
| B full | `record_count_distribution` | 1004 / 993 | 1.07 / 1.07 | - | error | error | 25.6 / 25.8 † | 16.3 / 16.3 | 13.9 / 13.9 |
| B full | `record_dicts` | 776 / 766 | 1.11 / 1.11 | - | error | error | 33.7 / 34 † | 21.6 / 21.6 | 18.5 / 18.6 |
| B full | `record_len_distribution` | 385 / 375 | 1.12 / 1.13 | - | 85.5 / 87.6 † | 47.6 / 48.7 | 62.4 / 63.8 † | 42.8 / 43.6 | 36.8 / 37.5 |
| B full | `record_name_set` | 617 / 607 | 1.1 / 1.1 | - | 53.4 / 54.2 | 30.3 / 30.7 | 39.3 / 39.8 | 27.2 / 27.5 | 23.4 / 23.6 |
| B full | `title_extraction` | 189 / 179 | 1.06 / 1.06 | - | 149 / 156 † | 76.9 / 80.8 | 65.3 / 68.3 † | 9.99 / 9.74 | 9.76 / 9.68 |
| B full | `write_round_trip` | 387 / 377 | 1.03 / 1.03 | - | 138 / 142 | 80.7 / 82.7 | 204 / 209 | 6.87 / 6.67 | 4.71 / 4.54 |
| B full | `write_sdf_round_trip` | 186 / 176 | 1.06 / 1.07 | - | 285 / 301 | 167 / 176 | 425 / 449 | 14.3 / 14.3 | 9.7 / 9.64 |
| C smoke | `append_record_write` | 10.1 / 0.292 | 0.988 / 1 | 3.04 / 9.32 | 9.62 / 58.9 | 9.09 / 40.8 | 10.9 / 25.8 | 10.2 / 6.7 | 10.2 / 6.43 |
| C smoke | `field_extraction` | 9.74 / 0.252 | 1.02 / 1.07 | 3.12 / 10.2 | 9.44 / 49.4 | 9.1 / 34.1 | 10.9 / 18.3 | 10.4 / 10.2 | 10.4 / 9.96 |
| C smoke | `field_get` | 9.8 / 0.176 | 1.01 / 1.05 | - | 9.21 / 62.5 | 8.9 / 42.2 | 10.6 / 17.8 | 10.5 / 14.7 | 10.5 / 14.3 |
| C smoke | `field_presence` | 9.87 / 0.163 | 1 / 1.04 | - | 9.34 / 67.3 | 8.99 / 46.7 | 10.7 / 19.3 | 10.4 / 15.9 | 10.4 / 15.5 |
| C smoke | `mdl_scan` | 9.94 / 0.296 | 1.02 / 1.02 | 3.06 / 6.32 | 8.14 / 2.34 | 8.16 / 2.37 | 10.3 / 2.27 | 10.3 / 3.14 | 10.2 / 2.92 |
| C smoke | `molecule_count` | 9.89 / 0.127 | 0.989 / 1.03 | 2.98 / 13.1 | 9.24 / 83.6 | 8.96 / 56.9 | 10.4 / 22.9 | 10.2 / 6.46 | 10.2 / 5.53 |
| C smoke | `record_count_distribution` | 9.78 / 0.276 | 1.02 / 1.06 | 3.11 / 9.47 | 9.36 / 45.2 | 9.13 / 31.9 | 10.7 / 16.7 | 10.5 / 9.21 | 10.4 / 8.91 |
| C smoke | `record_dicts` | 9.77 / 0.233 | 1.01 / 1.05 | - | 9.43 / 53.8 | 9.09 / 37.3 | 10.7 / 20.1 | 10.6 / 11.1 | 10.5 / 10.8 |
| C smoke | `record_len_distribution` | 9.46 / 0.157 | 1.04 / 1.09 | - | 9.83 / 76 | 9.26 / 51.9 | 11.2 / 28 | 10.8 / 16.3 | 10.8 / 16 |
| C smoke | `record_name_set` | 9.77 / 0.191 | 1.01 / 1.05 | - | 9.44 / 63 | 9.02 / 43.1 | 10.8 / 23 | 10.7 / 13.7 | 10.6 / 13.2 |
| C smoke | `title_extraction` | 9.61 / 0.132 | 1.04 / 1.06 | 3.08 / 12.6 | 9.55 / 83.5 | 9.12 / 56.6 | 11 / 22.6 | 10.5 / 4.05 | 10.7 / 5.79 |
| C smoke | `write_round_trip` | 9.75 / 0.212 | 1.02 / 1.02 | 3.18 / 12.3 | 10 / 80 | 9.49 / 55.8 | 11.1 / 34.1 | 10.5 / 3.94 | 10.4 / 3.71 |
| C smoke | `write_sdf_round_trip` | 9.64 / 0.208 | 1 / 1 | - | 10.1 / 81.7 | 9.59 / 56.9 | 11.4 / 34.9 | 10.6 / 4.03 | 10.6 / 3.74 |
| C representative | `append_record_write` | 26.8 / 16.6 | 1.03 / 1.03 | 11.2 / 16.4 | 49.6 / 75.1 | 30.6 / 44.5 | 29.2 / 41 | 9.73 / 9.49 | 8.31 / 7.24 |
| C representative | `field_extraction` | 25.5 / 15.5 | 1.04 / 1.06 | 11 / 16.4 | 34.7 / 52.2 | 20.5 / 28.5 | 19.7 / 26 | 12.9 / 14.6 | 11.3 / 12.1 |
| C representative | `field_get` | 19.1 / 9.09 | 1.04 / 1.09 | - | 38.2 / 71.7 | 22.2 / 37.8 | 19.4 / 29.5 | 17.2 / 25.1 | 15.4 / 20.9 |
| C representative | `field_presence` | 17.8 / 7.92 | 1.04 / 1.07 | - | 41.4 / 82.9 | 24 / 43.8 | 20.9 / 34 | 18.5 / 28.5 | 16.3 / 23.6 |
| C representative | `mdl_scan` | 23.9 / 13.5 | 1 / 1.01 | 8.79 / 13.4 | 5.26 / 3.35 | 5.27 / 3.36 | 6.46 / 3.86 | 7.15 / 5.01 | 6.19 / 3.27 |
| C representative | `molecule_count` | 14.7 / 4.79 | 1.02 / 1.04 | 13 / 33.8 | 48.3 / 131 | 26.4 / 64.3 | 24 / 52.4 | 10.6 / 11 | 9.92 / 8.35 |
| C representative | `record_count_distribution` | 27.9 / 17.8 | 1.04 / 1.07 | 10.1 / 14.2 | 30.5 / 43.3 | 18.6 / 24.6 | 18 / 22.5 | 11.7 / 12.5 | 10.2 / 10.3 |
| C representative | `record_dicts` | 23.3 / 13.3 | 1.05 / 1.08 | - | 37.8 / 60.4 | 23.2 / 34.7 | 21.8 / 30.6 | 14.1 / 17.1 | 12.5 / 14.2 |
| C representative | `record_len_distribution` | 17.9 / 7.77 | 1.02 / 1.07 | - | 45.4 / 94.4 | 26 / 49.6 | 26.8 / 48.8 | 18.3 / 29 | 16.2 / 23.9 |
| C representative | `record_name_set` | 20.5 / 10.4 | 1.05 / 1.07 | - | 40.3 / 71.7 | 23.4 / 38.4 | 23.6 / 36.9 | 16.1 / 21.9 | 14.2 / 18.2 |
| C representative | `title_extraction` | 15.4 / 5.3 | 0.977 / 1.02 | 12.4 / 30.6 | 47.7 / 123 | 26.4 / 61.6 | 23.2 / 48.2 | 9.1 / 7.22 | 9.07 / 7.28 |
| C representative | `write_round_trip` | 20.3 / 10.3 | 1.03 / 1.03 | 14 / 24.8 | 63.9 / 118 | 39.5 / 69.9 | 37.1 / 63.1 | 7.95 / 5.68 | 7.19 / 3.93 |
| C representative | `write_sdf_round_trip` | 16.1 / 6.03 | 1.01 / 1.04 | - | 80.7 / 202 | 49.8 / 119 | 46.9 / 108 | 10 / 9.73 | 9.1 / 6.85 |
| C full | `append_record_write` | 29.1 / 18.6 | 1.02 / 1.04 | 11.7 / 16.7 | 50.9 / 75.1 | 31.5 / 44.9 | 29.8 / 41.1 | 9.72 / 9.56 | 8.14 / 7.28 |
| C full | `field_extraction` | 27.7 / 17.4 | 1.04 / 1.06 | 11.3 / 16.3 | 35.1 / 51.2 | 20.7 / 28.3 | 20 / 25.9 | 12.9 / 14.6 | 11.3 / 12 |
| C full | `field_get` | 20.3 / 10.3 | 1.05 / 1.07 | - | 40.7 / 72.4 | 22.4 / 36.2 | 19.9 / 29.3 | 17.7 / 25 | 15.5 / 20.6 |
| C full | `field_presence` | 19 / 8.89 | 1.02 / 1.08 | - | 43.6 / 84.3 | 23.9 / 42.1 | 21.2 / 34 | 18.8 / 28.7 | 16.4 / 23.6 |
| C full | `mdl_scan` | 25.5 / 15.1 | 1 / 1.02 | 9.22 / 13.6 | 5.14 / 3.36 | 5.1 / 3.38 | 6.25 / 3.88 | 7.03 / 5.04 | 5.89 / 3.17 |
| C full | `molecule_count` | 15.5 / 5.42 | 1.01 / 1.03 | 13.7 / 33.8 | 51.3 / 131 | 27.7 / 64.1 | 24.8 / 52 | 10.4 / 11 | 9.56 / 8.15 |
| C full | `record_count_distribution` | 30.4 / 20 | 1.04 / 1.07 | 10.4 / 14.3 | 31.7 / 44 | 19.7 / 25.9 | 18.1 / 22.4 | 11.6 / 12.5 | 10.1 / 10.2 |
| C full | `record_dicts` | 24.9 / 14.9 | 1.04 / 1.07 | - | 39.7 / 61 | 23.3 / 33.6 | 22.5 / 30.7 | 14.4 / 17.2 | 12.5 / 14.1 |
| C full | `record_len_distribution` | 18.4 / 8.65 | 1.05 / 1.08 | - | 50.2 / 97.5 | 29.3 / 53.1 | 28.9 / 49.4 | 19.5 / 29.4 | 17 / 24.2 |
| C full | `record_name_set` | 21.9 / 11.8 | 1.03 / 1.07 | - | 42.2 / 71.8 | 23.9 / 37.7 | 24.4 / 36.7 | 16.4 / 21.9 | 14.4 / 18 |
| C full | `title_extraction` | 15.9 / 5.94 | 1.02 / 1.04 | 13.3 / 30.8 | 51.3 / 124 | 29.4 / 65.3 | 24.5 / 48.2 | 9.2 / 7.3 | 9.24 / 7.65 |
| C full | `write_round_trip` | 21.8 / 11.7 | 0.996 / 1.01 | 14.6 / 24.8 | 66.6 / 118 | 41.2 / 70 | 38.4 / 62.9 | 7.69 / 5.68 | 6.7 / 3.9 |
| C full | `write_sdf_round_trip` | 16.8 / 6.75 | 0.995 / 1.03 | - | 86.2 / 203 | 53.2 / 120 | 49.8 / 109 | 10.2 / 9.82 | 8.94 / 6.81 |
| D representative | `append_record_write` | 25 / 14.6 | 0.995 / 0.997 | 45.1 / 75.1 | 77.9 / 128 | 56.2 / 90.5 | 59.5 / 94.7 | 21.5 / 29.7 | 20.9 / 28.5 |
| D representative | `field_extraction` | 99.4 / 89.2 | 1.09 / 1.09 | 16.3 / 17.9 | 19.2 / 20.5 † | 15.6 / 16.5 | 20.2 / 21.3 | 15.4 / 16 | 15.1 / 15.7 |
| D representative | `field_get` | 32.8 / 22.7 | 1.09 / 1.13 | - | 32.3 / 43.1 † | 22.1 / 28.4 | 20.7 / 25.4 | 47.2 / 63.6 | 46.4 / 62.2 |
| D representative | `field_presence` | 32.5 / 22.4 | 1.12 / 1.16 | - | 32.5 / 43.5 † | 22.2 / 28.7 | 21 / 25.9 | 47.5 / 64.2 | 46.7 / 63 |
| D representative | `mdl_scan` | 22.4 / 12 | 1.01 / 1.01 | 31.9 / 56.8 | 4.99 / 2.61 | 4.97 / 2.58 | 8.51 / 7.32 | 9.06 / 8.44 | 8.84 / 7.88 |
| D representative | `molecule_count` | 16.2 / 6.25 | 1.01 / 0.997 | 42.9 / 106 | 64.9 / 155 † | 44.6 / 102 | 41.7 / 91.3 | 12.2 / 14.9 | 11.9 / 13.7 |
| D representative | `record_count_distribution` | 105 / 95.3 | 1.09 / 1.09 | 15.3 / 16.7 | 17.9 / 19 † | 14.6 / 15.2 † | 18.9 / 19.8 | 14.2 / 14.6 | 14 / 14.4 |
| D representative | `record_dicts` | 91.3 / 81.4 | 1.08 / 1.08 | - | 21.4 / 23 † | 17.4 / 18.5 † | 22.3 / 23.7 | 17 / 17.7 | 16.7 / 17.4 |
| D representative | `record_len_distribution` | 32.5 / 22.4 | 1.09 / 1.14 | - | 42.6 / 58.1 † | 32 / 42.8 † | 55 / 75.2 | 46.8 / 63.1 | 45.9 / 61.9 |
| D representative | `record_name_set` | 56.1 / 46.2 | 1.1 / 1.11 | - | 25.5 / 29.1 † | 19.2 / 21.6 † | 33.3 / 38.1 | 27.8 / 31.5 | 27.3 / 30.9 |
| D representative | `title_extraction` | 16.3 / 6.38 | 1.03 / 1.01 | 42.7 / 105 | 64.6 / 153 † | 44.4 / 101 | 41.3 / 89.9 | 11.6 / 13.3 | 12.1 / 14 |
| D representative | `write_round_trip` | 22.3 / 12.1 | 0.994 / 1 | 50 / 89.9 | 86.8 / 153 | 62.5 / 109 | 66.6 / 115 | 8.89 / 7.89 | 8.69 / 7.32 |
| D representative | `write_sdf_round_trip` | 18.7 / 8.84 | 1.02 / 1.01 | - | 103 / 210 | 74.3 / 148 | 79.5 / 157 | 10.6 / 10.8 | 10.3 / 10 |
| D full | `append_record_write` | 28 / 17.7 | 1 / 1 | 55.3 / 85.8 | 87.1 / 133 | 66.7 / 101 | 68.9 / 103 | 25.3 / 34.1 | 24.5 / 32.6 |
| D full | `field_extraction` | 140 / 130 | 1.09 / 1.1 | 16.4 / 17.4 | 18.5 / 19.3 † | 15.8 / 16.4 | 19.4 / 20.1 | 15.2 / 15.5 | 14.9 / 15.3 |
| D full | `field_get` | 41.4 / 31.2 | 1.1 / 1.14 | - | 31.5 / 39.1 † | 23.2 / 28.2 | 20.5 / 23.9 | 52.4 / 66.1 | 51.5 / 64.9 |
| D full | `field_presence` | 41.1 / 30.8 | 1.14 / 1.17 | - | 31.7 / 39.5 † | 23.3 / 28.4 | 20.7 / 24.2 | 52.7 / 66.7 | 51.7 / 65.4 |
| D full | `mdl_scan` | 23.7 / 13.7 | 1 / 1 | 40.9 / 68.6 | 5.03 / 2.68 | 5.06 / 2.67 | 9.3 / 8.39 | 10.1 / 9.67 | 9.82 / 9.06 |
| D full | `molecule_count` | 18.1 / 7.88 | 1.01 / 1 | 53 / 118 | 71.3 / 153 † | 52.8 / 111 | 46.6 / 93.6 | 12.7 / 15.6 | 12.4 / 14.7 |
| D full | `record_count_distribution` | 148 / 138 | 1.09 / 1.09 | 15.4 / 16.3 | 17.4 / 18.1 † | 14.8 / 15.3 | 18.3 / 18.9 | 14.1 / 14.4 | 13.9 / 14.2 |
| D full | `record_dicts` | 127 / 117 | 1.08 / 1.09 | - | 20.8 / 22 † | 17.8 / 18.7 | 22 / 23 | 17.1 / 17.6 | 16.7 / 17.3 |
| D full | `record_len_distribution` | 41 / 30.7 | 1.1 / 1.14 | - | 42.9 / 54.5 † | 34.4 / 43.2 | 59.5 / 75.9 | 51.7 / 65.4 | 50.9 / 64.3 |
| D full | `record_name_set` | 75.5 / 65.1 | 1.1 / 1.12 | - | 24.1 / 26.7 | 19.4 / 21.3 | 33.2 / 36.9 | 28.9 / 31.9 | 28.4 / 31.3 |
| D full | `title_extraction` | 18 / 7.94 | 1.02 / 1.02 | 53.1 / 117 | 72.2 / 153 † | 53.2 / 110 | 46.8 / 93 | 12.2 / 14.2 | 12.5 / 14.8 |
| D full | `write_round_trip` | 24.9 / 14.7 | 1.01 / 1 | 61.8 / 102 | 97.5 / 159 | 74.9 / 121 | 76.9 / 123 | 9.4 / 8.55 | 8.98 / 7.98 |
| D full | `write_sdf_round_trip` | 21.5 / 11.3 | 0.991 / 0.99 | - | 113 / 207 | 86.6 / 157 | 89.7 / 160 | 10.7 / 11.1 | 10.4 / 10.3 |

## Differences

Every `sendoff` run is checked against a standalone `awk`/`rg`/`sed` baseline
before Hyperfine times it. Toolkit outputs are hashed and compared against
`sendoff`'s on the read cases; `†` marks a difference:

- RDKit default skips records it cannot parse: 1 in A's representative tier,
  17 in B's and 5 in D's.
- RDKit hides record names that start with `_` (D representative).
- RDKit decodes values as UTF-8 and raises on B's GBK bytes.
- OEChem default drops B's 12 zero-atom records per 10,000.
