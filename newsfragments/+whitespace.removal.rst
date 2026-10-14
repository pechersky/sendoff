Str blocks strip only ASCII whitespace, as bytes blocks, RDKit, and OEChem
do, so a file reads the same in either mode: a value line of Unicode spaces,
such as U+3000, no longer ends the value, and ``title`` keeps trailing
Unicode spaces and ``\x1c`` to ``\x1f``.
