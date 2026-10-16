# Out-of-spec behavior

Where sendoff knowingly departs from the CTfile Formats specification.

## `$$$$` lines always end the record

This includes header lines and data values, where the grammar allows
`$$$$` as content. A molecule or value with such a line is cut there.
Following the grammar made reading 1.5–4% slower on typical files and up to
17% slower on others (EC2 c6i, compared with 1.0.0).

## `S  SKPnnn` lines skip nothing

An `M  END` or `$$$$` line after one still ends the molfile or record. Few
tools write or honor these lines
([Molfile S  SKP](http://www.dalkescientific.com/writings/diary/archive/2020/09/21/molfile_s__skp.html)).

## `M  END` header lines end the molfile

A program or comment line starting with `M  END` ends the molfile, just as a
`$$$$` header line ends the record.
