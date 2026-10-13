Parsing, record access, and writing are compiled and much faster; installing
from a source distribution needs a Rust toolchain. Against 0.1.x, on an ADMET
Predictor export of 1,998 blocks with about 560 records each, counting blocks
went from 502 ms to 9 ms, reading every record from 1,139 ms to 105 ms, and
rewriting the file (0.1.x: a ``block.write`` loop; now ``write_sdf``) from
781 ms to 43 ms. On a 151,527-block vendor catalog the same tasks went from
2.0 s to 95 ms, 3.3 s to 275 ms, and 3.2 s to 281 ms. ``parse_sdf`` now
memory-maps files on local Linux filesystems, which parses them about 6%
faster again; pipes and network filesystems are still read in chunks.
A block indexes its records on its second lookup, so further ``get``,
``block[name]``, ``in``, and ``len`` calls skip parsing until it is edited.
