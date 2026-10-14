# Benchmark results

sendoff `8e357ba` against `main` (`4d997bd`, pure Python), RDKit 2026.3.6 and
OEChem 2026.1.0, measured 2026-10-13 on EC2 c6i.2xlarge (x86_64) and
c7g.2xlarge (arm64): Amazon Linux 2023, Python 3.11.16, sendoff from a plain
release wheel. sendoff is faster in every cell.

Every number is other ÷ sendoff, so above 1 means sendoff is faster. Cells read
`wall / op`:

- *wall*: Hyperfine's median of 5 randomized runs after a warmup, each a fresh
  interpreter, so startup and imports count.
- *op*: the median time of the operation itself, timed in process after imports.

† marks output that differs from sendoff's (see [Differences](#differences)),
*error* a run that raised, and - a run not made: `main` lacks the newer methods
and cannot read B, which is not UTF-8. Geomeans take each toolkit's fastest
mode per case, whether or not its output matched. Write cases are not compared,
as the toolkits reserialize records.

## At a glance

How many times faster sendoff is than each toolkit, in operation time on
x86_64 / arm64, for a sample of the full matrix: a case or two per
[dataset](#datasets), including where the margin is smallest. Write cases are
not like for like, as the toolkits reserialize records.

| Dataset | Case | main | RDKit default | RDKit fast | OEChem default | OEChem tokenizer | OEChem threaded |
|---|---|---:|---:|---:|---:|---:|---:|
| B full | `molecule_count` | - | 159× / 122× † | 75× / 61× | 72× / 55× † | 21× / 10× | 16× / 7.8× |
| A full | `title_extraction` | 25× / 22× | 89× / 73× † | 48× / 39× | 36× / 28× | 9.0× / 4.1× | 8.2× / 4.2× |
| C full | `field_get` | - | 75× / 62× | 37× / 31× | 33× / 25× | 27× / 21× | 20× / 18× |
| B full | `field_get` | - | 74× / 59× † | 36× / 30× | 35× / 26× † | 43× / 34× | 34× / 30× |
| A full | `record_dicts` | - | 47× / 41× † | 28× / 25× † | 28× / 23× † | 17× / 14× † | 13× / 12× † |
| B full | `record_name_set` | - | 57× / 47× | 31× / 26× | 42× / 34× | 28× / 24× | 22× / 20× |
| D full | `field_extraction` | 16× / 15× | 20× / 16× † | 17× / 14× | 19× / 17× | 15× / 13× | 14× / 13× |
| A full | `write_round_trip` | 21× / 21× | 97× / 83× | 61× / 52× | 44× / 46× | 6.4× / 3.8× | 4.5× / 2.6× |
| A full | `write_sdf_round_trip` | - | 161× / 124× | 102× / 77× | 74× / 69× | 11× / 5.7× | 7.7× / 3.9× |
| A full | `mdl_scan` | 14× / 14× | 2.7× / 2.6× | 2.7× / 2.6× | 6.2× / 3.2× | 7.5× / 4.1× | 4.2× / 2.4× |
| D full | `mdl_scan` | 55× / 57× | 2.4× / 2.2× | 2.4× / 2.2× | 17× / 7.0× | 19× / 8.0× | 18× / 7.5× |

On 100-record files, where start-up dominates, wall time is 1.5× faster than
`main`, 4.4–4.6× than RDKit and 5.1–5.4× than OEChem ([geomeans](#summary)).
The perf round since `8e357ba` cut operation time a further 5% (x86_64) and 6%
(arm64), and `parse_sdf(path, borrow=True)` another 9–21% on B–D (`00a0e25`,
c6i/c7g.16xlarge); these tables predate both.

## Implementations

| Name | Reads with |
|---|---|
| sendoff, main | `parse_sdf` and `SDBlock` methods |
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

| Dataset | main x86_64 | main arm64 | RDKit x86_64 | RDKit arm64 | OEChem x86_64 | OEChem arm64 |
|---|---:|---:|---:|---:|---:|---:|
| A smoke | 1.52 | 1.56 | 4.42 | 4.56 | 5.29 | 5.18 |
| A representative | 8.09 | 8.69 | 16.2 | 15.6 | 7.53 | 6.69 |
| A full | 16.8 | 16.3 | 29.7 | 25.9 | 9.81 | 7.5 |
| B smoke | 1.55 | 1.56 | 4.41 | 4.38 | 5.37 | 5.05 |
| B representative | - | - | 21.2 | 20.8 | 9.59 | 8.51 |
| B full | - | - | 55.4 | 46.4 | 15.1 | 11.1 |
| C smoke | 1.46 | 1.5 | 4.38 | 4.41 | 5.17 | 5.08 |
| C representative | 6.87 | 7.15 | 15.5 | 14.5 | 7.29 | 6.66 |
| C full | 7.29 | 7.6 | 16.5 | 15.3 | 7.41 | 6.7 |
| D representative | 19.7 | 20.6 | 20 | 17.6 | 13.7 | 10.6 |
| D full | 23.8 | 24.8 | 23 | 19.7 | 15.2 | 11.3 |

### Geomean speedup, operation time

| Dataset | main x86_64 | main arm64 | RDKit x86_64 | RDKit arm64 | OEChem x86_64 | OEChem arm64 |
|---|---:|---:|---:|---:|---:|---:|
| A smoke | 8.91 | 11.1 | 27.7 | 29.8 | 6.96 | 6.73 |
| A representative | 17.6 | 17.2 | 33.3 | 28.4 | 10.5 | 7.73 |
| A full | 18.4 | 17.5 | 32.3 | 27.7 | 10.2 | 7.53 |
| B smoke | 11.4 | 13.6 | 27.1 | 27.9 | 8.89 | 7.93 |
| B representative | - | - | 58.4 | 49 | 16.2 | 11.5 |
| B full | - | - | 59.1 | 48.9 | 15.6 | 11.3 |
| C smoke | 8.16 | 8.92 | 32.3 | 31.1 | 7.06 | 6.39 |
| C representative | 16.7 | 16.5 | 37 | 31.7 | 11 | 8.52 |
| C full | 16.7 | 16.6 | 37.2 | 31.8 | 11 | 8.45 |
| D representative | 43.5 | 41.6 | 35.2 | 27.9 | 22.5 | 13.8 |
| D full | 47.7 | 46.3 | 37.1 | 29.4 | 23.4 | 14.4 |

### Peak RSS, MB, x86_64 (largest over cases)

| Dataset | sendoff | main | RDKit default | RDKit fast | OEChem default | OEChem tokenizer | OEChem threaded |
|---|---:|---:|---:|---:|---:|---:|---:|
| A smoke | 11 | 12 | 57 | 57 | 55 | 53 | 53 |
| A representative | 14 | 68 | 58 | 57 | 56 | 71 | 60 |
| A full | 24 | 864 | 68 | 67 | 64 | 332 | 71 |
| B smoke | 11 | 12 | 57 | 56 | 56 | 53 | 53 |
| B representative | 14 | - | 57 | 57 | 56 | 75 | 62 |
| B full | 36 | - | 84 | 82 | 76 | 792 | 83 |
| C smoke | 11 | 12 | 57 | 57 | 54 | 53 | 53 |
| C representative | 14 | 56 | 57 | 57 | 55 | 69 | 60 |
| C full | 14 | 62 | 57 | 57 | 55 | 71 | 60 |
| D representative | 14 | 158 | 58 | 57 | 55 | 87 | 86 |
| D full | 14 | 211 | 57 | 57 | 55 | 97 | 96 |

## Full results

### x86_64 (c6i.2xlarge)

| Dataset | Case | sendoff ms | main | RDKit default | RDKit fast | OEChem default | OEChem tokenizer | OEChem threaded |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| A smoke | `append_record_write` | 18.8 / 0.349 | 1.53 / 8.61 | 4.93 / 53.9 | 4.57 / 37.6 | 5.58 / 20.9 | 5.32 / 6.68 | 5.21 / 6.16 |
| A smoke | `field_extraction` | 19.2 / 0.38 | 1.51 / 8.75 | 4.57 / 35.5 | 4.39 / 25.6 | 5.39 / 15.6 | 5.23 / 8.96 | 5.1 / 8.32 |
| A smoke | `field_get` | 18.3 / 0.244 | - | 4.72 / 47.1 | 4.47 / 31.1 | 5.7 / 15 | 5.76 / 14 | 5.63 / 13.2 |
| A smoke | `field_presence` | 18.6 / 0.224 | - | 4.67 / 51.6 | 4.47 / 34.3 | 5.47 / 16.3 | 5.47 / 15.2 | 5.41 / 14.2 |
| A smoke | `mdl_scan` | 18.9 / 0.293 | 1.51 / 7.52 | 3.99 / 2.38 | 4 / 2.35 | 5.35 / 3.81 | 5.27 / 4.77 | 5.23 / 3.98 |
| A smoke | `molecule_count` | 18.9 / 0.201 | 1.46 / 9.84 | 4.54 / 56.5 | 4.3 / 39.9 | 5.47 / 17.3 | 5.07 / 5.4 | 5.13 / 4.27 |
| A smoke | `record_count_distribution` | 18.7 / 0.392 | 1.55 / 8.48 | 4.73 / 34.1 | 4.39 / 23.4 | 5.48 / 14.9 | 5.36 / 8.5 | 5.42 / 8.14 |
| A smoke | `record_dicts` | 18.8 / 0.339 | - | 4.69 / 39.4 | 4.32 / 27.2 | 5.53 / 18 | 5.53 / 10.1 | 5.46 / 9.7 |
| A smoke | `record_len_distribution` | 18.6 / 0.228 | - | 4.62 / 55.1 | 4.45 / 38.7 | 5.6 / 24.8 | 5.21 / 14.7 | 5.32 / 14 |
| A smoke | `record_name_set` | 18.7 / 0.277 | - | 4.65 / 45.4 | 4.45 / 30.7 | 5.69 / 20.8 | 5.44 / 12.2 | 5.38 / 11.7 |
| A smoke | `title_extraction` | 18.9 / 0.206 | 1.49 / 9.66 | 4.56 / 56.1 | 4.34 / 39.4 | 5.37 / 17.1 | 5.26 / 4.1 | 5.34 / 4.63 |
| A smoke | `write_round_trip` | 18.5 / 0.302 | 1.57 / 9.74 | 4.98 / 61.8 | 4.69 / 43.1 | 5.76 / 23.5 | 5.41 / 3.97 | 5.2 / 3.15 |
| A smoke | `write_sdf_round_trip` | 18.6 / 0.285 | - | 4.99 / 66.5 | 4.68 / 45.4 | 5.59 / 25 | 5.35 / 4.21 | 5.41 / 3.43 |
| A representative | `append_record_write` | 36.7 / 17.6 | 8.82 / 16.9 | 41 / 81.3 | 26 / 49.9 | 20.4 / 36.9 | 7.86 / 10.7 | 6.11 / 7.26 |
| A representative | `field_extraction` | 41.4 / 21.9 | 8.29 / 14.5 | 23.6 / 41.3 † | 14.4 / 23.8 | 14.8 / 23.3 | 9.61 / 13.6 | 8.07 / 10.7 |
| A representative | `field_get` | 30 / 11.3 | - | 27.3 / 65.9 † | 15 / 33.4 | 13.6 / 27.4 | 13.4 / 26.6 | 11.2 / 21 |
| A representative | `field_presence` | 28.7 / 9.94 | - | 28.5 / 74.9 † | 15.7 / 37.9 | 14.3 / 31.4 | 13.8 / 30 | 11.7 / 23.7 |
| A representative | `mdl_scan` | 34.5 / 15.2 | 6.67 / 13.3 | 3.3 / 2.59 | 3.3 / 2.62 | 5.4 / 5.92 | 6.01 / 7.15 | 4.69 / 4.17 |
| A representative | `molecule_count` | 26 / 7.21 | 8.1 / 25.6 | 30.8 / 101 † | 16.9 / 50.5 | 15.3 / 40.7 | 6.98 / 11.3 | 6.22 / 8.64 |
| A representative | `record_count_distribution` | 43.4 / 24.1 | 7.9 / 13.1 | 22.5 / 37.4 † | 13.8 / 21.7 | 14.1 / 21.4 | 9.06 / 12.1 | 7.57 / 9.44 |
| A representative | `record_dicts` | 37.6 / 18.5 | - | 26.4 / 49.4 † | 16.2 / 28.9 † | 16.7 / 28.5 † | 10.6 / 16.1 † | 8.9 / 12.8 † |
| A representative | `record_len_distribution` | 28.9 / 9.85 | - | 31.5 / 84.8 † | 18.4 / 46.6 | 20.2 / 49.3 | 13.8 / 30.1 | 11.4 / 23.6 |
| A representative | `record_name_set` | 33.2 / 14.1 | - | 27.6 / 59.7 | 16.2 / 33 | 17.9 / 35.2 | 12 / 21.4 | 10 / 16.8 |
| A representative | `title_extraction` | 26.4 / 7.62 | 8.01 / 24.4 | 30.9 / 97.4 † | 17.1 / 49.6 | 15 / 39.1 | 6.31 / 8.86 | 6.09 / 8.29 |
| A representative | `write_round_trip` | 33.1 / 13.7 | 9.07 / 20 | 44.9 / 103 | 28.5 / 63.3 | 22 / 45.9 | 5.6 / 6.33 | 4.88 / 4.59 |
| A representative | `write_sdf_round_trip` | 27.7 / 8.35 | - | 53.7 / 169 | 33.9 / 104 | 26.3 / 75.9 | 6.67 / 10.4 | 5.77 / 7.54 |
| A full | `append_record_write` | 287 / 267 | 16.6 / 17.7 | 69.6 / 74.4 | 44.4 / 47.4 | 32.7 / 34.7 | 10.2 / 10.6 | 6.29 / 6.38 |
| A full | `field_extraction` | 340 / 320 | 14.4 / 15.1 | 37.2 / 39.2 † | 22.3 / 23.5 | 21.8 / 22.8 | 13.5 / 14 | 10.2 / 10.5 |
| A full | `field_get` | 188 / 169 | - | 54.1 / 59.7 † | 29 / 31.7 | 23.3 / 25.3 | 24.7 / 26.7 | 18.7 / 20.1 |
| A full | `field_presence` | 170 / 150 | - | 60.1 / 67.5 † | 31.9 / 35.6 | 25.8 / 28.5 | 27.1 / 30 | 20.4 / 22.4 |
| A full | `mdl_scan` | 239 / 220 | 13.2 / 14.2 | 2.84 / 2.72 | 2.8 / 2.69 | 6.19 / 6.23 | 7.42 / 7.52 | 4.33 / 4.21 |
| A full | `molecule_count` | 128 / 108 | 22.4 / 26.2 | 77.1 / 90.5 † | 40.9 / 47.7 | 32.2 / 37.1 | 10.6 / 11.4 | 7.9 / 8.4 |
| A full | `record_count_distribution` | 372 / 353 | 13.1 / 13.7 | 33.8 / 35.4 † | 20.3 / 21.2 | 20 / 20.8 | 12.1 / 12.4 | 9.02 / 9.25 |
| A full | `record_dicts` | 288 / 269 | - | 44.3 / 47.1 † | 26.9 / 28.5 † | 26.7 / 28.1 † | 16.1 / 16.8 † | 12.2 / 12.6 † |
| A full | `record_len_distribution` | 164 / 145 | - | 70 / 78.7 † | 40.3 / 45.1 | 42.9 / 47.9 | 27.8 / 30.7 | 20.8 / 22.8 |
| A full | `record_name_set` | 229 / 209 | - | 50.7 / 55.1 | 29.5 / 31.8 | 31.4 / 33.8 | 20.3 / 21.7 | 15.3 / 16.2 |
| A full | `title_extraction` | 131 / 112 | 21.7 / 25.1 | 77 / 89.4 † | 41.2 / 47.5 | 31.8 / 36.2 | 8.55 / 9.03 | 7.84 / 8.22 |
| A full | `write_round_trip` | 222 / 203 | 19 / 20.6 | 88.8 / 96.7 | 56.5 / 61.4 | 41.1 / 44.4 | 6.42 / 6.44 | 4.63 / 4.54 |
| A full | `write_sdf_round_trip` | 141 / 121 | - | 140 / 161 | 88.9 / 102 | 65 / 74.4 | 10 / 10.7 | 7.36 / 7.65 |
| B smoke | `append_record_write` | 19.4 / 0.339 | 1.51 / 11 | 4.82 / 57.2 | 4.51 / 38.4 | 6.38 / 72.6 | 5.19 / 8.21 | 5.2 / 7.71 |
| B smoke | `field_extraction` | 19.2 / 0.446 | 1.59 / 10.4 | 4.61 / 31.2 | 4.4 / 21.4 | 5.65 / 17 | 5.51 / 10.7 | 5.46 / 10.2 |
| B smoke | `field_get` | 18.7 / 0.247 | - | 4.61 / 46.3 | 4.39 / 29.9 | 5.43 / 16 | 5.63 / 19.5 | 5.49 / 18.7 |
| B smoke | `field_presence` | 18.5 / 0.224 | - | 4.67 / 51 | 4.39 / 34 | 5.45 / 17.6 | 5.64 / 21.5 | 5.51 / 21.3 |
| B smoke | `mdl_scan` | 18.8 / 0.272 | 1.52 / 9.43 | 3.96 / 2.4 | 3.91 / 2.41 | 5.2 / 4.66 | 5.16 / 5.84 | 5.29 / 5.15 |
| B smoke | `molecule_count` | 18.6 / 0.179 | 1.53 / 13.3 | 4.63 / 62.7 | 4.45 / 40.2 | 5.7 / 20.9 | 5.47 / 7.18 | 5.27 / 5.97 |
| B smoke | `record_count_distribution` | 18.9 / 0.453 | 1.61 / 10.2 | 4.65 / 30.6 | 4.4 / 21.3 | 5.58 / 16.5 | 5.62 / 10.4 | 5.44 / 9.87 |
| B smoke | `record_dicts` | 18.8 / 0.384 | - | 4.66 / 37.2 | 4.37 / 26.2 | 5.5 / 20 | 5.44 / 12.5 | 5.6 / 12 |
| B smoke | `record_len_distribution` | 18.8 / 0.224 | - | 4.71 / 57.1 | 4.46 / 38.6 | 5.83 / 31.9 | 5.68 / 21.4 | 5.63 / 20.6 |
| B smoke | `record_name_set` | 18.6 / 0.304 | - | 4.7 / 42 | 4.47 / 28.3 | 5.69 / 23.8 | 5.6 / 15.8 | 5.6 / 15.3 |
| B smoke | `title_extraction` | 18.4 / 0.18 | 1.51 / 13.1 | 4.75 / 63.3 | 4.33 / 40.5 | 5.54 / 20.9 | 5.31 / 5.72 | 5.26 / 6.08 |
| B smoke | `write_round_trip` | 18.6 / 0.28 | 1.6 / 13 | 4.98 / 68.6 | 4.62 / 45.6 | 6.55 / 86.6 | 5.41 / 4.83 | 5.22 / 4.03 |
| B smoke | `write_sdf_round_trip` | 18.6 / 0.276 | - | 5.04 / 69.4 | 4.63 / 46.3 | 6.62 / 87.9 | 5.32 / 4.88 | 5.34 / 4.12 |
| B representative | `append_record_write` | 36.8 / 17.6 | - | 44.9 / 89.5 | 27 / 52.3 | 63.9 / 128 | 9.11 / 13.5 | 7.05 / 9.16 |
| B representative | `field_extraction` | 44 / 24.7 | - | error | error | 18.2 / 28.4 † | 12.3 / 17.9 | 10.5 / 14.6 |
| B representative | `field_get` | 29.4 / 10.5 | - | 29 / 74.3 † | 15.3 / 35.7 | 15.7 / 34.5 † | 18.6 / 42.7 | 15.8 / 35 |
| B representative | `field_presence` | 29 / 9.53 | - | 29.4 / 81.9 † | 15.4 / 39.4 | 15.8 / 38.1 † | 18.6 / 46.3 | 15.8 / 38.1 |
| B representative | `mdl_scan` | 34.2 / 14.8 | - | error | error | 6.12 / 7.58 | 6.8 / 8.97 | 5.37 / 5.71 |
| B representative | `molecule_count` | 23.7 / 4.89 | - | 35.1 / 155 † | 18.3 / 73.8 | 18.9 / 71.4 † | 8.5 / 20.7 | 7.62 / 16.4 |
| B representative | `record_count_distribution` | 46.1 / 26.7 | - | error | error | 17.3 / 26.1 † | 11.5 / 16 | 9.81 / 13.2 |
| B representative | `record_dicts` | 39.1 / 19.9 | - | error | error | 20.8 / 35.9 † | 13.8 / 22.3 | 11.9 / 18.4 |
| B representative | `record_len_distribution` | 28 / 9.18 | - | 35 / 98.8 † | 20.1 / 53.3 | 27.1 / 71.9 † | 19.2 / 47.5 | 16.4 / 39.2 |
| B representative | `record_name_set` | 35.3 / 16 | - | 28.1 / 57.5 | 16.2 / 31.1 | 21.9 / 42 | 15.5 / 27.9 | 13.2 / 23 |
| B representative | `title_extraction` | 24.6 / 5.35 | - | 34.6 / 145 † | 18.2 / 70.1 | 18.2 / 65.1 † | 7.69 / 16.2 | 7.26 / 14.9 |
| B representative | `write_round_trip` | 31.7 / 12.4 | - | 51.1 / 125 | 30.9 / 73.2 | 73.2 / 180 | 6.54 / 8.68 | 5.69 / 6.43 |
| B representative | `write_sdf_round_trip` | 24.1 / 5.32 | - | 67.4 / 291 | 40.7 / 170 | 96.5 / 418 | 8.5 / 20.1 | 7.43 / 15 |
| B full | `append_record_write` | 587 / 567 | - | 87.5 / 90.4 | 51 / 52.7 | 125 / 130 | 13.7 / 13.9 | 8.23 / 8.34 |
| B full | `field_extraction` | 826 / 806 | - | error | error | 27.7 / 28.2 † | 17.7 / 17.9 | 13.8 / 14 |
| B full | `field_get` | 358 / 339 | - | 70.6 / 74.4 † | 33.9 / 35.6 | 33.3 / 34.8 † | 41 / 42.9 | 32.2 / 33.7 |
| B full | `field_presence` | 323 / 304 | - | 78.4 / 83 † | 37.7 / 39.8 | 36.9 / 38.9 † | 45.1 / 47.5 | 35.4 / 37.3 |
| B full | `mdl_scan` | 503 / 483 | - | error | error | 7.57 / 7.68 | 8.9 / 9.01 | 5.28 / 5.29 |
| B full | `molecule_count` | 173 / 154 | - | 142 / 159 † | 66.9 / 74.7 | 65.1 / 72.5 † | 19.6 / 21.3 | 15.1 / 16.3 |
| B full | `record_count_distribution` | 902 / 883 | - | error | error | 25.5 / 26 † | 15.9 / 16.1 | 12.3 / 12.4 |
| B full | `record_dicts` | 672 / 653 | - | error | error | 34.8 / 35.7 † | 21.7 / 22.1 | 17.1 / 17.4 |
| B full | `record_len_distribution` | 313 / 294 | - | 93.9 / 99.7 † | 50.8 / 53.9 | 69.4 / 73.6 † | 45.9 / 48.5 | 36 / 38 |
| B full | `record_name_set` | 536 / 517 | - | 55.4 / 57.3 | 30 / 31 | 41 / 42.3 | 27.4 / 28.2 | 21.6 / 22.2 |
| B full | `title_extraction` | 185 / 166 | - | 136 / 151 † | 65.1 / 72.1 | 61.8 / 68.3 † | 15.6 / 16.8 | 14.1 / 15.1 |
| B full | `write_round_trip` | 409 / 389 | - | 123 / 129 | 71.9 / 75.2 | 177 / 186 | 8.82 / 8.96 | 6.38 / 6.45 |
| B full | `write_sdf_round_trip` | 188 / 168 | - | 269 / 299 | 156 / 174 | 386 / 429 | 19.1 / 20.6 | 13.8 / 14.9 |
| C smoke | `append_record_write` | 17.6 / 0.274 | 1.49 / 7.59 | 4.8 / 59 | 4.56 / 41.1 | 5.39 / 21.6 | 5.28 / 6.88 | 5.12 / 6.01 |
| C smoke | `field_extraction` | 17.6 / 0.248 | 1.46 / 8.16 | 4.53 / 44.2 | 4.36 / 30.2 | 5.39 / 17.5 | 5.19 / 9.62 | 5.26 / 8.53 |
| C smoke | `field_get` | 17.7 / 0.185 | - | 4.51 / 53.1 | 4.48 / 35.3 | 5.19 / 16.1 | 5.16 / 13 | 5.18 / 11.7 |
| C smoke | `field_presence` | 17.3 / 0.157 | - | 4.64 / 62.9 | 4.41 / 41.6 | 5.37 / 19.1 | 5.42 / 15.5 | 5.24 / 13.8 |
| C smoke | `mdl_scan` | 17.7 / 0.223 | 1.44 / 6.68 | 3.95 / 2.58 | 3.94 / 2.61 | 5.12 / 4.05 | 5.2 / 5.28 | 5.11 / 4.19 |
| C smoke | `molecule_count` | 17.4 / 0.137 | 1.45 / 9.33 | 4.55 / 70.1 | 4.31 / 46.1 | 5.32 / 20.2 | 5.22 / 6.54 | 5.27 / 4.95 |
| C smoke | `record_count_distribution` | 17.9 / 0.27 | 1.44 / 7.67 | 4.44 / 40.5 | 4.2 / 27.7 | 5.3 / 15.9 | 5.11 / 8.79 | 5.17 / 7.84 |
| C smoke | `record_dicts` | 17.4 / 0.22 | - | 4.6 / 50.4 | 4.4 / 34.3 | 5.36 / 19.8 | 5.33 / 11 | 5.2 / 9.79 |
| C smoke | `record_len_distribution` | 17.5 / 0.158 | - | 4.59 / 66.8 | 4.43 / 45 | 5.34 / 26 | 5.28 / 15.3 | 5.29 / 13.4 |
| C smoke | `record_name_set` | 18.1 / 0.196 | - | 4.44 / 53.9 | 4.27 / 36.3 | 5.2 / 21.1 | 5.14 / 12.5 | 5.1 / 11.2 |
| C smoke | `title_extraction` | 17.3 / 0.139 | 1.45 / 9.38 | 4.55 / 70.6 | 4.37 / 46.5 | 5.32 / 20.6 | 5.24 / 4.9 | 5.27 / 5.16 |
| C smoke | `write_round_trip` | 17.4 / 0.23 | 1.5 / 8.64 | 4.91 / 69.3 | 4.62 / 48.5 | 5.46 / 25 | 5.27 / 4.25 | 5.17 / 3.2 |
| C smoke | `write_sdf_round_trip` | 17.5 / 0.206 | - | 4.81 / 77.3 | 4.64 / 54.3 | 5.58 / 29.2 | 5.39 / 4.8 | 5.14 / 3.6 |
| C representative | `append_record_write` | 33.8 / 16 | 7.02 / 13.2 | 37.7 / 75.3 | 23.2 / 44.6 | 18.6 / 33.5 | 7.4 / 9.84 | 5.82 / 6.55 |
| C representative | `field_extraction` | 31.7 / 14 | 7.27 / 14.7 | 24.3 / 49.9 | 14.2 / 27.2 | 15 / 27.3 | 9.49 / 15 | 7.82 / 11.1 |
| C representative | `field_get` | 25.5 / 7.91 | - | 26 / 74.6 | 14.1 / 36.4 | 13.8 / 32.6 | 11.9 / 26.6 | 9.82 / 20 |
| C representative | `field_presence` | 24.1 / 6.57 | - | 27.4 / 89.5 | 14.9 / 43.8 | 14.6 / 39.4 | 12.5 / 31.8 | 10.3 / 23.7 |
| C representative | `mdl_scan` | 30.2 / 12.4 | 5.67 / 11.7 | 3.57 / 3.05 | 3.54 / 3.06 | 5.56 / 6.2 | 6.2 / 7.64 | 4.83 / 4.39 |
| C representative | `molecule_count` | 22.3 / 4.79 | 6.88 / 26.6 | 28.8 / 119 | 15.4 / 57.1 | 14.8 / 50 | 7.25 / 14.3 | 6.35 / 10.3 |
| C representative | `record_count_distribution` | 33.6 / 15.9 | 6.78 / 12.7 | 22.9 / 44 | 13.4 / 24 | 13.9 / 23.6 | 8.78 / 12.9 | 7.27 / 9.54 |
| C representative | `record_dicts` | 29.4 / 11.9 | - | 26.4 / 59.4 | 15.5 / 32.5 | 16.1 / 32.2 | 10.3 / 17.7 | 8.46 / 13.3 |
| C representative | `record_len_distribution` | 23.9 / 6.36 | - | 30.3 / 103 | 17.3 / 53.7 | 18.8 / 56 | 12.6 / 32.7 | 10.4 / 24.3 |
| C representative | `record_name_set` | 26.4 / 8.86 | - | 27.5 / 73.9 | 15.7 / 38.8 | 17.2 / 40.8 | 11.5 / 23.7 | 9.47 / 17.8 |
| C representative | `title_extraction` | 22.6 / 5.13 | 6.83 / 25 | 29.2 / 115 | 15.8 / 55.4 | 15 / 47.8 | 6.52 / 10.6 | 6.23 / 9.75 |
| C representative | `write_round_trip` | 28.3 / 10.5 | 7.86 / 18.6 | 44.3 / 112 | 27.4 / 66.7 | 21.5 / 48.9 | 5.85 / 6.91 | 5.05 / 4.8 |
| C representative | `write_sdf_round_trip` | 23.2 / 5.72 | - | 54.1 / 207 | 33.2 / 122 | 26.3 / 90.4 | 7.06 / 12.7 | 6.15 / 8.82 |
| C full | `append_record_write` | 36 / 18.1 | 7.41 / 13.3 | 39.6 / 74.7 | 24.4 / 44.6 | 19.3 / 33.1 | 7.43 / 9.69 | 5.83 / 6.43 |
| C full | `field_extraction` | 33.4 / 15.8 | 7.61 / 14.5 | 25.9 / 50.3 | 15 / 27.2 | 15.4 / 26.8 | 9.91 / 15 | 8 / 11.1 |
| C full | `field_get` | 26.5 / 8.87 | - | 27.9 / 75.3 | 15 / 36.7 | 14.4 / 32.7 | 12.5 / 26.7 | 10.1 / 19.9 |
| C full | `field_presence` | 24.7 / 7.36 | - | 29.9 / 90.9 | 16 / 44.1 | 15.5 / 39.5 | 13.3 / 32.1 | 10.9 / 23.8 |
| C full | `mdl_scan` | 31.9 / 14.1 | 5.99 / 11.7 | 3.54 / 3.03 | 3.59 / 3.06 | 5.61 / 6.21 | 6.33 / 7.63 | 4.65 / 4.15 |
| C full | `molecule_count` | 22.9 / 5.33 | 7.5 / 27.3 | 31.4 / 122 | 16.6 / 58.2 | 15.9 / 50.8 | 7.4 / 14.6 | 6.45 / 10.6 |
| C full | `record_count_distribution` | 35.8 / 18 | 7.15 / 12.8 | 24.1 / 44 | 14.1 / 24 | 14.3 / 23.4 | 9.03 / 12.8 | 7.37 / 9.5 |
| C full | `record_dicts` | 31 / 13.2 | - | 28.1 / 60.5 | 16.4 / 33 | 16.8 / 32.4 | 10.6 / 18 | 8.7 / 13.3 |
| C full | `record_len_distribution` | 25 / 7.13 | - | 32.4 / 104 | 18.2 / 54.1 | 19.7 / 56.1 | 13.1 / 33 | 10.7 / 24.3 |
| C full | `record_name_set` | 27.3 / 9.9 | - | 29.8 / 75 | 16.9 / 39.5 | 18.2 / 41 | 12.1 / 24 | 9.9 / 17.9 |
| C full | `title_extraction` | 23.2 / 5.79 | 7.38 / 25.1 | 31.9 / 116 | 17.1 / 56 | 15.9 / 47.6 | 6.75 / 10.7 | 6.46 / 9.73 |
| C full | `write_round_trip` | 30 / 12 | 8.2 / 18.3 | 46.9 / 111 | 28.8 / 65.8 | 22.5 / 48.4 | 5.83 / 6.85 | 4.94 / 4.71 |
| C full | `write_sdf_round_trip` | 24 / 6.49 | - | 58.4 / 205 | 35.9 / 122 | 28.1 / 89.7 | 7.24 / 12.7 | 6.12 / 8.72 |
| D representative | `append_record_write` | 31.3 / 13.2 | 27.3 / 62.8 | 56.8 / 129 | 40.7 / 91.3 | 34.3 / 74.6 | 16.5 / 32.1 | 16.1 / 31.2 |
| D representative | `field_extraction` | 99.5 / 80.9 | 13.1 / 15.7 | 17.8 / 21 † | 14.4 / 16.8 | 16.8 / 19.6 | 13 / 14.8 | 12.9 / 14.7 |
| D representative | `field_get` | 36.5 / 18.6 | - | 26 / 47.2 † | 17.5 / 30.7 | 15.8 / 26.1 | 35.9 / 65.7 | 35.6 / 65.1 |
| D representative | `field_presence` | 36.6 / 18.4 | - | 25.8 / 47.6 † | 17.4 / 30.9 | 15.9 / 26.6 | 35.9 / 66.3 | 35.5 / 65.4 |
| D representative | `mdl_scan` | 29.6 / 11.6 | 18.9 / 45.8 | 3.26 / 2.32 | 3.29 / 2.33 | 8.98 / 14.8 | 9.48 / 16 | 9.43 / 15.8 |
| D representative | `molecule_count` | 24.7 / 6.58 | 22.2 / 79.2 | 38.2 / 133 † | 25.7 / 85.8 | 23.2 / 72.8 | 10.6 / 25.7 | 10.5 / 25.1 |
| D representative | `record_count_distribution` | 101 / 83.2 | 12.7 / 15.1 | 17.3 / 20.2 † | 14.1 / 16.2 † | 16.5 / 18.9 | 12.6 / 14.2 | 12.5 / 14.1 |
| D representative | `record_dicts` | 90.3 / 72.3 | - | 20 / 24.1 † | 16.2 / 19.2 † | 18.9 / 22.3 | 14.6 / 17 | 14.5 / 16.8 |
| D representative | `record_len_distribution` | 36.4 / 18.5 | - | 35.8 / 66.6 † | 26.8 / 48.8 † | 41.8 / 77.2 | 35.4 / 64.6 | 35.1 / 64 |
| D representative | `record_name_set` | 58.8 / 40.4 | - | 22.5 / 31 † | 17.1 / 23 † | 26.3 / 36 | 22.5 / 30.4 | 22.3 / 30.1 |
| D representative | `title_extraction` | 25.1 / 6.72 | 21.7 / 77.3 | 37.7 / 130 † | 25.4 / 84.5 | 23.2 / 72.7 | 10.3 / 24.4 | 10.3 / 24.7 |
| D representative | `write_round_trip` | 29.4 / 11.4 | 28.1 / 70.6 | 59.9 / 149 | 42.9 / 105 | 36.2 / 85.6 | 8.96 / 15 | 8.83 / 14.7 |
| D representative | `write_sdf_round_trip` | 26.8 / 8.51 | - | 65.6 / 198 | 47.4 / 141 | 39.9 / 114 | 9.86 / 20 | 9.74 / 19.6 |
| D full | `append_record_write` | 33.8 / 15.7 | 34.7 / 72.8 | 65.4 / 136 | 50.3 / 103 | 39.9 / 79.5 | 20.1 / 37.2 | 19.7 / 36.3 |
| D full | `field_extraction` | 135 / 117 | 13.5 / 15.5 | 17.9 / 20.2 † | 15.3 / 17.2 | 17 / 19 | 13.3 / 14.7 | 13.2 / 14.5 |
| D full | `field_get` | 43.4 / 25.5 | - | 27.3 / 43.6 † | 19.9 / 31 | 16.6 / 24.6 | 42.2 / 68.2 | 41.8 / 67.4 |
| D full | `field_presence` | 42.8 / 25 | - | 27.6 / 44.4 † | 20.1 / 31.6 | 16.6 / 24.8 | 42.7 / 69.4 | 42.3 / 68.7 |
| D full | `mdl_scan` | 31.3 / 13.3 | 24.1 / 54.7 | 3.24 / 2.37 | 3.29 / 2.38 | 10.3 / 17.2 | 10.9 / 18.6 | 10.8 / 18.3 |
| D full | `molecule_count` | 26.7 / 8.36 | 28.1 / 86.9 | 43.9 / 132 † | 32.2 / 94.3 | 26.5 / 73.7 | 12.1 / 27.5 | 12 / 26.9 |
| D full | `record_count_distribution` | 138 / 120 | 13.2 / 15 | 17.3 / 19.3 † | 14.8 / 16.5 | 16.4 / 18.1 | 12.8 / 14 | 12.7 / 13.9 |
| D full | `record_dicts` | 123 / 104 | - | 20.2 / 23.2 † | 17.3 / 19.7 | 19.1 / 21.6 | 15 / 16.7 | 14.9 / 16.6 |
| D full | `record_len_distribution` | 43.7 / 25.5 | - | 38.3 / 63 † | 31 / 50.3 | 47.2 / 77.3 | 41 / 66.6 | 40.5 / 65.7 |
| D full | `record_name_set` | 74.9 / 56.7 | - | 23 / 29.1 | 18.6 / 23.3 | 28.2 / 35.6 | 24.5 / 30.7 | 24.3 / 30.4 |
| D full | `title_extraction` | 26.3 / 8.29 | 28.4 / 87 | 44.9 / 134 † | 32.8 / 95.3 | 26.9 / 74.1 | 12 / 26.8 | 12.2 / 27.3 |
| D full | `write_round_trip` | 31.7 / 13.8 | 36.1 / 80.8 | 69.5 / 154 | 53.4 / 117 | 42.4 / 90.4 | 10.2 / 16.7 | 10.1 / 16.4 |
| D full | `write_sdf_round_trip` | 28.2 / 10.5 | - | 78.1 / 203 | 60 / 155 | 47.6 / 119 | 11.5 / 22 | 11.3 / 21.5 |

### arm64 (c7g.2xlarge)

| Dataset | Case | sendoff ms | main | RDKit default | RDKit fast | OEChem default | OEChem tokenizer | OEChem threaded |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| A smoke | `append_record_write` | 19.7 / 0.338 | 1.57 / 11.4 | 5.06 / 59 | 4.71 / 41.3 | 5.51 / 27 | 5.14 / 6.77 | 5.19 / 6.58 |
| A smoke | `field_extraction` | 19.5 / 0.398 | 1.61 / 10.3 | 4.81 / 37.1 | 4.56 / 25.8 | 5.4 / 15.8 | 5.24 / 9.15 | 5.27 / 9 |
| A smoke | `field_get` | 19.5 / 0.254 | - | 4.74 / 50 | 4.49 / 33.9 | 5.27 / 14.9 | 5.33 / 14.6 | 5.23 / 14.2 |
| A smoke | `field_presence` | 19.5 / 0.236 | - | 4.7 / 54.4 | 4.49 / 37 | 5.26 / 16 | 5.26 / 15.8 | 5.24 / 15.2 |
| A smoke | `mdl_scan` | 19.7 / 0.287 | 1.52 / 9.43 | 4.03 / 2.69 | 4.04 / 2.7 | 5.03 / 2.61 | 5.09 / 3.53 | 5.16 / 3.3 |
| A smoke | `molecule_count` | 19.3 / 0.188 | 1.53 / 12.9 | 4.89 / 66.2 | 4.68 / 44.4 | 5.58 / 19.2 | 5.43 / 4.63 | 5.22 / 4.25 |
| A smoke | `record_count_distribution` | 19.6 / 0.426 | 1.6 / 9.65 | 4.79 / 34.5 | 4.53 / 24 | 5.36 / 14.7 | 5.19 / 8.49 | 5.22 / 8.25 |
| A smoke | `record_dicts` | 19.5 / 0.355 | - | 4.77 / 41.6 | 4.56 / 29.3 | 5.41 / 17.9 | 5.22 / 10.5 | 5.24 / 10.1 |
| A smoke | `record_len_distribution` | 19.6 / 0.234 | - | 4.72 / 59.1 | 4.62 / 40.8 | 5.35 / 25.3 | 5.21 / 15.7 | 5.27 / 15.3 |
| A smoke | `record_name_set` | 19.5 / 0.292 | - | 4.79 / 47.5 | 4.6 / 33.4 | 5.37 / 20.6 | 5.3 / 12.8 | 5.3 / 12.4 |
| A smoke | `title_extraction` | 19.4 / 0.196 | 1.53 / 12.3 | 4.7 / 64.2 | 4.56 / 44.3 | 5.28 / 18.3 | 5.16 / 3.05 | 5.15 / 4.22 |
| A smoke | `write_round_trip` | 19.7 / 0.309 | 1.57 / 12.2 | 4.98 / 63.4 | 4.69 / 44.3 | 5.46 / 28.2 | 5.04 / 3.05 | 5.03 / 2.79 |
| A smoke | `write_sdf_round_trip` | 19.5 / 0.279 | - | 5.08 / 70.8 | 4.76 / 49.1 | 5.56 / 31.4 | 5.12 / 3.45 | 5.12 / 3.08 |
| A representative | `append_record_write` | 41.8 / 20.8 | 9.67 / 18.1 | 38.1 / 72.7 | 24.3 / 45.1 | 22.4 / 40.1 | 6.9 / 9.03 | 5.84 / 6.96 |
| A representative | `field_extraction` | 48.5 / 28.6 | 8.69 / 13.8 | 22.8 / 35.9 † | 14.5 / 21.8 | 13.6 / 19.5 | 8.86 / 11.5 | 7.93 / 9.95 |
| A representative | `field_get` | 36.6 / 15.8 | - | 24.6 / 52.1 † | 14.1 / 27.7 | 11.5 / 20.4 | 11.8 / 21.1 | 10.6 / 18.2 |
| A representative | `field_presence` | 34 / 14.4 | - | 26.6 / 57.3 † | 15.2 / 30.4 | 12.4 / 22.3 | 12.6 / 22.9 | 11.3 / 19.7 |
| A representative | `mdl_scan` | 38.4 / 18.2 | 7.3 / 13.9 | 3.23 / 2.56 | 3.24 / 2.57 | 4.05 / 3.21 | 4.49 / 4.09 | 3.74 / 2.5 |
| A representative | `molecule_count` | 29.9 / 10.1 | 8.69 / 23 | 29 / 78.1 † | 16.8 / 41.9 | 13.4 / 30 | 5.28 / 5.89 | 4.77 / 4.44 |
| A representative | `record_count_distribution` | 51.2 / 31.4 | 8.24 / 12.5 | 21.4 / 32.4 † | 13.5 / 19.5 | 12.7 / 17.6 | 8.24 / 10.3 | 7.36 / 8.82 |
| A representative | `record_dicts` | 44.6 / 24.5 | - | 24.9 / 42.2 † | 15.4 / 24.9 † | 14.9 / 23.2 † | 9.63 / 13.5 † | 8.61 / 11.7 † |
| A representative | `record_len_distribution` | 33.8 / 14.2 | - | 29.9 / 65.9 † | 18.3 / 38.2 | 18.3 / 36.8 | 12.6 / 23.1 | 11.2 / 19.9 |
| A representative | `record_name_set` | 38.9 / 19.2 | - | 26.3 / 49.3 | 16 / 28.5 | 16.1 / 27.6 | 11.1 / 17.4 | 9.92 / 15 |
| A representative | `title_extraction` | 30.8 / 10.7 | 8.42 / 21.7 | 29.2 / 77 † | 16.7 / 40.7 | 13.1 / 28.8 | 4.63 / 4.16 | 4.6 / 4.18 |
| A representative | `write_round_trip` | 38.1 / 17.2 | 10.1 / 20.7 | 40.9 / 85.8 | 26 / 52.9 | 23.5 / 46.2 | 4.31 / 3.83 | 3.75 / 2.56 |
| A representative | `write_sdf_round_trip` | 32.5 / 11.7 | - | 47.9 / 127 | 30.5 / 78.2 | 27.6 / 68.2 | 5.12 / 5.64 | 4.49 / 3.89 |
| A full | `append_record_write` | 334 / 312 | 17.7 / 18.8 | 65.9 / 70.2 | 41.3 / 43.9 | 37.5 / 39.7 | 8.82 / 9.04 | 6.34 / 6.43 |
| A full | `field_extraction` | 452 / 431 | 13.4 / 14 | 32.5 / 33.9 † | 20.3 / 21.1 | 18.5 / 19.2 | 11.3 / 11.5 | 9.43 / 9.63 |
| A full | `field_get` | 261 / 240 | - | 45.2 / 48.9 † | 24.7 / 26.5 | 18.5 / 19.7 | 19.7 / 21 | 16.5 / 17.5 |
| A full | `field_presence` | 240 / 218 | - | 50 / 54.6 † | 27.2 / 29.5 | 20.1 / 21.6 | 21.4 / 22.9 | 17.8 / 19.1 |
| A full | `mdl_scan` | 294 / 272 | 13.2 / 14.2 | 2.68 / 2.58 | 2.69 / 2.59 | 3.38 / 3.23 | 4.26 / 4.13 | 2.62 / 2.41 |
| A full | `molecule_count` | 176 / 154 | 20.4 / 23.1 | 65.7 / 74.3 † | 34.5 / 38.8 | 25.8 / 28.7 | 5.77 / 5.76 | 4.53 / 4.43 |
| A full | `record_count_distribution` | 494 / 473 | 12.2 / 12.7 | 29.7 / 30.8 † | 18.5 / 19.1 | 16.9 / 17.5 | 10.1 / 10.3 | 8.42 / 8.56 |
| A full | `record_dicts` | 386 / 364 | - | 38.5 / 40.6 † | 24.1 / 25.3 † | 22 / 23 † | 13.3 / 13.7 † | 11.1 / 11.5 † |
| A full | `record_len_distribution` | 238 / 215 | - | 56.7 / 62.1 † | 32.7 / 35.6 | 33.2 / 36 | 21.4 / 23 | 17.9 / 19.2 |
| A full | `record_name_set` | 310 / 288 | - | 44.2 / 47.3 | 26.9 / 28.6 | 25.9 / 27.4 | 16.7 / 17.5 | 13.9 / 14.6 |
| A full | `title_extraction` | 184 / 162 | 19.5 / 21.9 | 65.2 / 73.3 † | 35.1 / 39.3 | 25.1 / 27.7 | 4.31 / 4.11 | 4.32 / 4.19 |
| A full | `write_round_trip` | 280 / 258 | 19.7 / 21.2 | 76.9 / 83 | 48.2 / 51.9 | 43.1 / 46.3 | 3.99 / 3.84 | 2.8 / 2.6 |
| A full | `write_sdf_round_trip` | 196 / 174 | - | 110 / 124 | 68.9 / 77 | 61.6 / 68.7 | 5.72 / 5.7 | 4.01 / 3.86 |
| B smoke | `append_record_write` | 20.8 / 0.357 | 1.54 / 13.5 | 4.79 / 57 | 4.47 / 38.7 | 6.07 / 75.3 | 4.94 / 7.66 | 4.9 / 7.4 |
| B smoke | `field_extraction` | 19.9 / 0.474 | 1.67 / 12.1 | 4.76 / 32.4 | 4.53 / 22.9 | 5.42 / 17.6 | 5.25 / 11.1 | 5.28 / 11 |
| B smoke | `field_get` | 20.6 / 0.276 | - | 4.44 / 45.7 | 4.25 / 30.2 | 5.08 / 15.3 | 5.1 / 19.5 | 5.08 / 19.1 |
| B smoke | `field_presence` | 20.5 / 0.252 | - | 4.49 / 50.2 | 4.26 / 32.6 | 5.02 / 16.7 | 5.11 / 21.2 | 5.07 / 20.8 |
| B smoke | `mdl_scan` | 20.4 / 0.28 | 1.49 / 11.3 | 3.94 / 2.64 | 3.9 / 2.58 | 4.89 / 2.95 | 4.9 / 3.95 | 4.9 / 3.68 |
| B smoke | `molecule_count` | 20 / 0.177 | 1.5 / 16.6 | 4.57 / 69.4 | 4.39 / 46.2 | 5.2 / 22.7 | 5.07 / 5.15 | 5 / 4.88 |
| B smoke | `record_count_distribution` | 20.1 / 0.506 | 1.65 / 11.4 | 4.7 / 30.5 | 4.5 / 21.4 | 5.41 / 16.5 | 5.22 / 10.3 | 5.2 / 10.2 |
| B smoke | `record_dicts` | 20.1 / 0.416 | - | 4.71 / 37.2 | 4.47 / 26.4 | 5.36 / 20.3 | 5.22 / 12.8 | 5.19 / 12.6 |
| B smoke | `record_len_distribution` | 19.8 / 0.244 | - | 4.71 / 57.6 | 4.5 / 38.8 | 5.55 / 32 | 5.26 / 21.7 | 5.31 / 21.4 |
| B smoke | `record_name_set` | 20.4 / 0.343 | - | 4.56 / 41.2 | 4.34 / 28.4 | 5.25 / 22.8 | 5.15 / 15.7 | 5.11 / 15.4 |
| B smoke | `title_extraction` | 20.3 / 0.188 | 1.49 / 15.7 | 4.5 / 66.6 | 4.32 / 43.9 | 5.08 / 21.7 | 4.9 / 3.63 | 4.92 / 4.82 |
| B smoke | `write_round_trip` | 20.3 / 0.303 | 1.57 / 15.7 | 4.91 / 67.2 | 4.58 / 45.2 | 6.21 / 87.8 | 4.94 / 3.35 | 4.91 / 3.17 |
| B smoke | `write_sdf_round_trip` | 20.4 / 0.285 | - | 4.94 / 71.7 | 4.56 / 47.5 | 6.21 / 94.2 | 4.91 / 3.57 | 4.92 / 3.33 |
| B representative | `append_record_write` | 41.2 / 21 | - | 42.1 / 78.7 | 25.6 / 46.3 | 62.2 / 117 | 8.16 / 11.1 | 6.91 / 8.74 |
| B representative | `field_extraction` | 53.3 / 33 | - | error | error | 16.7 / 23.8 † | 11.3 / 15.1 | 10.3 / 13.4 |
| B representative | `field_get` | 34.7 / 14.8 | - | 27.1 / 58 † | 14.9 / 29.6 | 14.2 / 26.4 † | 17.4 / 33.8 | 15.9 / 30.3 |
| B representative | `field_presence` | 34 / 13.7 | - | 28 / 63.8 † | 15.5 / 32.7 | 14.5 / 28.5 † | 17.8 / 36.6 | 16.1 / 32.5 |
| B representative | `mdl_scan` | 36 / 15.7 | - | error | error | 4.75 / 4.37 | 5.28 / 5.48 | 4.35 / 3.48 |
| B representative | `molecule_count` | 26.7 / 6.78 | - | 34.1 / 122 † | 18.8 / 62.3 | 17.7 / 54.8 † | 6.45 / 10.3 | 5.75 / 7.74 |
| B representative | `record_count_distribution` | 55.5 / 35.2 | - | error | error | 15.9 / 22.3 † | 10.7 / 13.9 | 9.68 / 12.3 |
| B representative | `record_dicts` | 48.2 / 27.9 | - | error | error | 18.7 / 28.7 † | 12.5 / 17.9 | 11.4 / 16 |
| B representative | `record_len_distribution` | 33.3 / 13.3 | - | 32.8 / 76.2 † | 19.3 / 42.3 | 25 / 54.7 † | 18 / 37.3 | 16.4 / 33.2 |
| B representative | `record_name_set` | 41.7 / 21.3 | - | 26.4 / 48 | 16.1 / 27.8 | 20.2 / 34.8 | 14.6 / 23.7 | 13.2 / 21.1 |
| B representative | `title_extraction` | 27.3 / 7.3 | - | 34.7 / 118 † | 18.8 / 59 | 17.5 / 51.3 † | 5.75 / 7.4 | 5.73 / 7.33 |
| B representative | `write_round_trip` | 35.3 / 14.9 | - | 49 / 111 | 29.2 / 63.8 | 71.7 / 163 | 5.1 / 5.2 | 4.39 / 3.58 |
| B representative | `write_sdf_round_trip` | 27.2 / 7.2 | - | 62.5 / 225 | 38.1 / 133 | 92.8 / 337 | 6.62 / 10.7 | 5.67 / 7.33 |
| B full | `append_record_write` | 710 / 688 | - | 77.3 / 79.6 | 45 / 46.2 | 114 / 117 | 11.2 / 11.3 | 7.98 / 8.07 |
| B full | `field_extraction` | 1097 / 1075 | - | error | error | 23.3 / 23.7 † | 15.1 / 15.3 | 13 / 13.2 |
| B full | `field_get` | 507 / 484 | - | 56.2 / 58.6 † | 28.4 / 29.5 | 25.3 / 26.3 † | 33.2 / 34.4 | 28.5 / 29.6 |
| B full | `field_presence` | 467 / 445 | - | 60.6 / 63.4 † | 30.9 / 32.2 | 27.5 / 28.6 † | 35.7 / 37.2 | 30.6 / 31.9 |
| B full | `mdl_scan` | 531 / 510 | - | error | error | 4.46 / 4.43 | 5.58 / 5.53 | 3.34 / 3.26 |
| B full | `molecule_count` | 243 / 221 | - | 112 / 122 † | 55.4 / 60.6 | 50.3 / 54.8 † | 9.94 / 10.3 | 7.57 / 7.8 |
| B full | `record_count_distribution` | 1172 / 1150 | - | error | error | 21.9 / 22.2 † | 13.9 / 14 | 11.9 / 12 |
| B full | `record_dicts` | 930 / 907 | - | error | error | 28.2 / 28.7 † | 18 / 18.3 | 15.4 / 15.7 |
| B full | `record_len_distribution` | 456 / 435 | - | 72.2 / 75.5 † | 40.2 / 41.9 | 52.7 / 55 † | 36.2 / 37.6 | 31 / 32.3 |
| B full | `record_name_set` | 726 / 704 | - | 45.4 / 46.7 | 25.8 / 26.5 | 33.4 / 34.3 | 23.1 / 23.7 | 19.9 / 20.4 |
| B full | `title_extraction` | 257 / 236 | - | 109 / 119 † | 56.6 / 61.3 | 48 / 51.8 † | 7.35 / 7.39 | 7.18 / 7.35 |
| B full | `write_round_trip` | 500 / 479 | - | 107 / 112 | 62.5 / 65.1 | 158 / 164 | 5.32 / 5.26 | 3.64 / 3.58 |
| B full | `write_sdf_round_trip` | 251 / 230 | - | 211 / 230 | 123 / 134 | 315 / 344 | 10.6 / 10.9 | 7.18 / 7.37 |
| C smoke | `append_record_write` | 20.3 / 0.326 | 1.51 / 8.35 | 4.78 / 52.8 | 4.52 / 36.6 | 5.43 / 23.1 | 5.05 / 6 | 5.07 / 5.76 |
| C smoke | `field_extraction` | 20.2 / 0.299 | 1.5 / 8.61 | 4.55 / 41.6 | 4.38 / 28.7 | 5.25 / 15.4 | 5.01 / 8.56 | 5.02 / 8.37 |
| C smoke | `field_get` | 19.9 / 0.211 | - | 4.53 / 52.1 | 4.38 / 35.2 | 5.21 / 14.9 | 5.19 / 12.2 | 5.19 / 11.9 |
| C smoke | `field_presence` | 20 / 0.189 | - | 4.62 / 58 | 4.45 / 40.2 | 5.28 / 16.6 | 5.13 / 13.7 | 5.13 / 13.3 |
| C smoke | `mdl_scan` | 20.2 / 0.249 | 1.51 / 7.51 | 4.01 / 2.79 | 4.02 / 2.82 | 5.06 / 2.7 | 5.08 / 3.73 | 5.03 / 3.47 |
| C smoke | `molecule_count` | 20.1 / 0.164 | 1.46 / 10.1 | 4.54 / 64.7 | 4.4 / 44.1 | 5.11 / 17.7 | 5.03 / 5 | 5.01 / 4.28 |
| C smoke | `record_count_distribution` | 20.3 / 0.321 | 1.5 / 8.14 | 4.51 / 38.9 | 4.4 / 27.4 | 5.14 / 14.4 | 5.04 / 7.91 | 5.02 / 7.66 |
| C smoke | `record_dicts` | 20 / 0.263 | - | 4.61 / 47.7 | 4.45 / 33.1 | 5.23 / 17.8 | 5.18 / 9.8 | 5.14 / 9.6 |
| C smoke | `record_len_distribution` | 20 / 0.187 | - | 4.65 / 63.7 | 4.38 / 43.5 | 5.3 / 23.5 | 5.1 / 13.7 | 5.13 / 13.4 |
| C smoke | `record_name_set` | 20 / 0.228 | - | 4.6 / 52.9 | 4.4 / 36.2 | 5.28 / 19.4 | 5.2 / 11.5 | 5.18 / 11.1 |
| C smoke | `title_extraction` | 20 / 0.164 | 1.48 / 10.1 | 4.58 / 66.9 | 4.37 / 45.4 | 5.27 / 18.2 | 5.05 / 3.25 | 5.14 / 4.64 |
| C smoke | `write_round_trip` | 20.1 / 0.262 | 1.55 / 9.98 | 4.87 / 64.9 | 4.61 / 45.3 | 5.38 / 27.7 | 5.09 / 3.2 | 5.07 / 3.01 |
| C smoke | `write_sdf_round_trip` | 20.2 / 0.247 | - | 4.85 / 68.7 | 4.58 / 47.9 | 5.44 / 29.4 | 5.06 / 3.38 | 5.07 / 3.14 |
| C representative | `append_record_write` | 39.4 / 18.6 | 7.65 / 14.7 | 33.7 / 67 | 20.8 / 39.7 | 19.9 / 36.6 | 6.62 / 8.47 | 5.66 / 6.46 |
| C representative | `field_extraction` | 38.7 / 18 | 7.28 / 14 | 22.9 / 44.7 | 13.5 / 24.4 | 13 / 22.3 | 8.49 / 12.5 | 7.48 / 10.4 |
| C representative | `field_get` | 31.2 / 10.8 | - | 23.5 / 60.2 | 13.6 / 31.7 | 11.9 / 24.8 | 10.6 / 21 | 9.46 / 17.5 |
| C representative | `field_presence` | 29.4 / 9.33 | - | 25.1 / 70.4 | 14.5 / 37.2 | 12.7 / 28.8 | 11.2 / 24.2 | 9.9 / 20 |
| C representative | `mdl_scan` | 34.8 / 14.3 | 6.02 / 12.7 | 3.6 / 3.16 | 3.61 / 3.17 | 4.42 / 3.63 | 4.89 / 4.73 | 4.23 / 3.08 |
| C representative | `molecule_count` | 27 / 6.72 | 7.06 / 24.1 | 26.3 / 93.5 | 14.4 / 45.8 | 13 / 37.3 | 5.76 / 7.82 | 5.4 / 5.95 |
| C representative | `record_count_distribution` | 41.2 / 20.7 | 6.81 / 12.2 | 20.6 / 37.2 | 12.6 / 21.1 | 12.2 / 19.3 | 7.89 / 10.7 | 6.93 / 8.82 |
| C representative | `record_dicts` | 35.7 / 15.6 | - | 24.7 / 51.6 | 15.1 / 29.6 | 14.2 / 26.1 | 9.22 / 14.6 | 8.15 / 12.1 |
| C representative | `record_len_distribution` | 29.4 / 9.18 | - | 27.7 / 79.9 | 15.8 / 42 | 16.3 / 41.3 | 11.1 / 24.5 | 9.84 / 20.2 |
| C representative | `record_name_set` | 32.9 / 12.2 | - | 25 / 61.2 | 14.5 / 32.8 | 14.7 / 31.5 | 10 / 18.7 | 8.85 / 15.5 |
| C representative | `title_extraction` | 27.4 / 7.21 | 6.97 / 22.5 | 26.8 / 90.4 | 14.8 / 45.2 | 13 / 35.3 | 5.12 / 5.3 | 5.1 / 5.34 |
| C representative | `write_round_trip` | 33.5 / 13 | 8.47 / 19.6 | 38.6 / 93.3 | 23.9 / 55.3 | 22.4 / 49.9 | 4.81 / 4.5 | 4.35 / 3.11 |
| C representative | `write_sdf_round_trip` | 28.3 / 7.92 | - | 46 / 154 | 28.3 / 91 | 26.7 / 82.3 | 5.72 / 7.42 | 5.18 / 5.22 |
| C full | `append_record_write` | 41.3 / 20.9 | 8.24 / 14.9 | 35.8 / 66.8 | 22.2 / 39.9 | 21 / 36.6 | 6.84 / 8.51 | 5.73 / 6.47 |
| C full | `field_extraction` | 40.8 / 20.4 | 7.68 / 14 | 23.8 / 43.8 | 14.1 / 24.2 | 13.6 / 22.1 | 8.75 / 12.5 | 7.68 / 10.3 |
| C full | `field_get` | 32.3 / 12.1 | - | 25.5 / 61.6 | 14 / 30.8 | 12.5 / 25 | 11.1 / 21.3 | 9.73 / 17.6 |
| C full | `field_presence` | 30.8 / 10.6 | - | 27 / 71 | 14.8 / 35.5 | 13.1 / 28.7 | 11.6 / 24.2 | 10.1 / 19.9 |
| C full | `mdl_scan` | 36.7 / 16.2 | 6.41 / 12.7 | 3.57 / 3.13 | 3.54 / 3.15 | 4.35 / 3.61 | 4.88 / 4.7 | 4.09 / 2.96 |
| C full | `molecule_count` | 28 / 7.63 | 7.59 / 24 | 28.4 / 93.4 | 15.3 / 45.6 | 13.7 / 37 | 5.78 / 7.83 | 5.29 / 5.79 |
| C full | `record_count_distribution` | 44.1 / 23.3 | 7.15 / 12.3 | 21.9 / 37.9 | 13.6 / 22.3 | 12.5 / 19.3 | 7.99 / 10.7 | 6.98 / 8.81 |
| C full | `record_dicts` | 37.9 / 17.4 | - | 26.1 / 52.3 | 15.3 / 28.8 | 14.7 / 26.3 | 9.48 / 14.8 | 8.22 / 12.1 |
| C full | `record_len_distribution` | 30.6 / 10.3 | - | 30.2 / 81.9 | 17.6 / 44.6 | 17.4 / 41.5 | 11.7 / 24.7 | 10.2 / 20.3 |
| C full | `record_name_set` | 34.5 / 13.7 | - | 26.8 / 61.5 | 15.2 / 32.3 | 15.5 / 31.5 | 10.4 / 18.8 | 9.1 / 15.4 |
| C full | `title_extraction` | 28.7 / 8.16 | 7.38 / 22.5 | 28.5 / 90.2 | 16.3 / 47.5 | 13.6 / 35.1 | 5.11 / 5.32 | 5.13 / 5.57 |
| C full | `write_round_trip` | 35.3 / 14.8 | 9.03 / 19.5 | 41.3 / 92.7 | 25.5 / 55.2 | 23.8 / 49.6 | 4.76 / 4.48 | 4.15 / 3.07 |
| C full | `write_sdf_round_trip` | 29.4 / 8.89 | - | 49.4 / 154 | 30.4 / 91.3 | 28.5 / 82.6 | 5.83 / 7.46 | 5.12 / 5.17 |
| D representative | `append_record_write` | 37 / 16.7 | 30.5 / 66 | 52.7 / 112 | 38 / 79.6 | 40.3 / 83.2 | 14.5 / 26.1 | 14.1 / 25.1 |
| D representative | `field_extraction` | 129 / 108 | 12.6 / 14.7 | 14.8 / 16.9 † | 12 / 13.6 | 15.5 / 17.6 | 11.8 / 13.1 | 11.6 / 12.9 |
| D representative | `field_get` | 49 / 28.6 | - | 21.6 / 34.3 † | 14.8 / 22.6 | 13.9 / 20.2 | 31.6 / 50.6 | 31 / 49.5 |
| D representative | `field_presence` | 49.1 / 28.3 | - | 21.5 / 34.4 † | 14.7 / 22.7 | 13.9 / 20.5 | 31.4 / 50.8 | 30.9 / 49.9 |
| D representative | `mdl_scan` | 35.1 / 14.5 | 20.4 / 47.3 | 3.19 / 2.17 | 3.17 / 2.15 | 5.44 / 6.09 | 5.78 / 7.02 | 5.64 / 6.55 |
| D representative | `molecule_count` | 29.9 / 9.77 | 23.2 / 68.2 | 35.1 / 99 † | 24.1 / 65.5 | 22.6 / 58.5 | 6.59 / 9.56 | 6.45 / 8.78 |
| D representative | `record_count_distribution` | 134 / 114 | 12 / 13.9 | 14 / 15.8 † | 11.4 / 12.7 † | 14.8 / 16.5 | 11.2 / 12.2 | 11 / 12 |
| D representative | `record_dicts` | 114 / 93.1 | - | 17.2 / 20.1 † | 13.9 / 16.2 † | 17.9 / 20.7 | 13.6 / 15.5 | 13.4 / 15.2 |
| D representative | `record_len_distribution` | 49 / 28.5 | - | 28.2 / 45.6 † | 21.2 / 33.6 † | 36.5 / 59 | 31 / 49.5 | 30.4 / 48.6 |
| D representative | `record_name_set` | 75.3 / 54.9 | - | 19 / 24.5 † | 14.3 / 18.2 † | 24.8 / 32.1 | 20.7 / 26.5 | 20.4 / 26 |
| D representative | `title_extraction` | 29.8 / 9.87 | 23.4 / 67.7 | 35.4 / 98.9 † | 24.3 / 65.5 | 22.7 / 58.1 | 6.33 / 8.58 | 6.61 / 9.04 |
| D representative | `write_round_trip` | 35.5 / 14.9 | 31.4 / 72.8 | 54.4 / 124 | 39.2 / 88.2 | 41.8 / 92.8 | 5.58 / 6.4 | 5.45 / 5.93 |
| D representative | `write_sdf_round_trip` | 31.9 / 11.8 | - | 60.6 / 157 | 43.6 / 111 | 46.7 / 118 | 6.24 / 8.09 | 6.06 / 7.52 |
| D full | `append_record_write` | 40.4 / 20 | 38.3 / 76.1 | 60.3 / 118 | 46.2 / 89.4 | 47.7 / 91.4 | 17.5 / 30.2 | 17 / 28.9 |
| D full | `field_extraction` | 175 / 154 | 13.1 / 14.7 | 14.8 / 16.3 † | 12.7 / 13.9 | 15.6 / 17 | 12.2 / 13.1 | 12 / 12.9 |
| D full | `field_get` | 59.6 / 38.8 | - | 21.8 / 31.4 † | 16.1 / 22.7 | 14.2 / 19.2 | 36.3 / 53.1 | 35.7 / 52.1 |
| D full | `field_presence` | 59.2 / 38.4 | - | 21.9 / 31.7 † | 16.2 / 22.8 | 14.4 / 19.4 | 36.5 / 53.5 | 35.9 / 52.5 |
| D full | `mdl_scan` | 37.1 / 16.5 | 26.1 / 56.9 | 3.21 / 2.22 | 3.23 / 2.22 | 5.94 / 6.97 | 6.47 / 8.03 | 6.27 / 7.53 |
| D full | `molecule_count` | 32.5 / 11.8 | 29.4 / 78.6 | 39.6 / 102 † | 29.3 / 73.9 | 25.9 / 62.5 | 7.04 / 10.4 | 6.9 / 9.83 |
| D full | `record_count_distribution` | 184 / 163 | 12.4 / 13.9 | 14 / 15.3 † | 11.9 / 13 | 14.8 / 16.1 | 11.4 / 12.2 | 11.2 / 12 |
| D full | `record_dicts` | 151 / 130 | - | 17.5 / 19.7 † | 15 / 16.8 | 18.5 / 20.6 | 14.4 / 15.8 | 14.1 / 15.5 |
| D full | `record_len_distribution` | 60 / 38.8 | - | 29.3 / 43.2 † | 23.5 / 34.2 | 40.6 / 60.1 | 35.3 / 51.7 | 34.8 / 50.9 |
| D full | `record_name_set` | 98.4 / 77.2 | - | 18.5 / 22.5 | 14.9 / 17.9 | 25.5 / 31.1 | 22.2 / 26.9 | 21.8 / 26.4 |
| D full | `title_extraction` | 32.3 / 11.9 | 29.6 / 77.7 | 40.2 / 102 † | 29.7 / 73.4 | 26.1 / 61.9 | 6.82 / 9.48 | 6.95 / 9.87 |
| D full | `write_round_trip` | 38.5 / 17.9 | 39.9 / 84.1 | 62.9 / 131 | 48.3 / 99.2 | 49.6 / 101 | 6.07 / 7.02 | 5.8 / 6.55 |
| D full | `write_sdf_round_trip` | 35.2 / 14.5 | - | 68.9 / 161 | 52.8 / 122 | 54.7 / 125 | 6.55 / 8.61 | 6.35 / 8.05 |

## Differences

- RDKit default skips records it cannot parse: 1 in A's representative tier,
  17 in B's and 5 in D's.
- RDKit hides record names that start with `_` (D representative).
- RDKit decodes values as UTF-8 and raises on B's GBK bytes.
- OEChem default drops B's 12 zero-atom records per 10,000.
- A `record_dicts`: in text mode, sendoff treats a value line holding only
  U+3000 (ideographic space) as blank, ending the value. The toolkits and
  sendoff's bytes mode treat only ASCII whitespace as blank.
