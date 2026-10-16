# Out-of-spec behavior

Where sendoff knowingly departs from the CTfile Formats specification.

## `$$$$` lines always end the record

This includes header lines and data values, where the grammar allows
`$$$$` as content. Following the grammar made reading 1.5–4% slower on
typical files and up to 17% slower on others (EC2 c6i, compared with
1.0.0).

```text
methane
  program
$$$$ comment
  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 1 0 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0 0 0 0
M  V30 END ATOM
M  V30 END CTAB
M  END
> <NOTE>
first
$$$$ second
third

$$$$
```

sendoff reads three records: `methane`, cut after its program line; one
titled with the counts line, whose `NOTE` is `first`; and one titled
`third`.

## `S  SKPnnn` lines skip nothing

An `M  END` or `$$$$` line after one still ends the molfile or record. Only
V2000 molfiles have these lines, and few tools write or honor them
([Molfile S  SKP](http://www.dalkescientific.com/writings/diary/archive/2020/09/21/molfile_s__skp.html)).

```text
methane
  program

  1  0  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
S  SKP  1
M  END
M  END
$$$$
```

sendoff ends the molfile at the first `M  END`, which the specification
skips.

## `M  END` header lines end the molfile

A program or comment line starting with `M  END` ends the molfile, just as a
`$$$$` header line ends the record.

```text
methane
  program
M  END comment
  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 1 0 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0 0 0 0
M  V30 END ATOM
M  V30 END CTAB
M  END
$$$$
```

sendoff ends the molfile at the comment line, so its connection table does
not parse.
