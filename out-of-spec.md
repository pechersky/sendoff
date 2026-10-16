# Out-of-spec behavior

Where sendoff knowingly reads SDfiles differently from the CTfile Formats
specification, and why.

sendoff finds where a record ends by searching for lines that start with
`$$$$`. It finds where a molfile ends by searching for the first line after
the title that starts with `M  END`. It does not read the header or follow
the grammar to decide whether such a line is content. Following the grammar
would slow down every record, even though writers rarely produce these
lines.

The costs below were measured on an x86-64 EC2 c6i instance. A build that
follows the grammar was compared with sendoff 1.0.0 on 10,000-record files
from three suppliers, on a 279 MB file, and on two files derived from one of
the 10,000-record files: one without data items, and one without the blank
line after the last item.

## A `$$$$` line always ends the record

The specification says a line beginning with `$$$$` ends a record, but the
grammar has room for such lines elsewhere:

- The title line must not hold `$$$$`, but some writers write one anyway.
  The program and comment lines are free text, so they may start with
  `$$$$`.
- A data item's value runs to the blank line that ends the item, so a value
  line may start with `$$$$`.

sendoff ends the record at any line that starts with `$$$$`:

- A molecule with such a header line is cut at that line, and the rest of
  the molecule reads as a record of its own. When the title holds `$$$$`,
  the part before the cut is an empty record.
- A data value with such a line is cut before it, and the record ends
  there. The rest of the record reads as a record of its own, with the line
  after the `$$$$` line as its title.

Following the grammar for data values made reading 1.5–4.3% slower (4–21 ns
per record) on files whose items each end with a blank line. On the files
with no data items or without the blank line after the last item, it made
reading 5–11% slower (19–32 ns per record). There, it has to look ahead for
the next record's counts line. Following the grammar for header lines as
well made reading 2–9% and 7–17% slower. To tell a program or comment line
from the next record's title, it had to search back for the line breaks
before every `$$$$` line.

## `S  SKPnnn` lines skip nothing

The specification says an `S  SKPnnn` line in the properties block makes
readers ignore the next `nnn` lines. sendoff reads those lines like any
others, so an `M  END` line after `S  SKPnnn` ends the molfile and a `$$$$`
line after it ends the record.

Almost no tool writes these lines, and few honor them
([Molfile S  SKP](http://www.dalkescientific.com/writings/diary/archive/2020/09/21/molfile_s__skp.html)).
Skipped lines read the same either way unless they hold an `M  END` line or
a line starting with `$$$$`. Honoring them would mean searching every record
for them, a cost every file would pay for a feature that is barely used.

## An `M  END` header line ends the molfile

The program and comment lines are free text, so they may start with
`M  END`. sendoff ends the molfile at the first line after the title that
starts with `M  END`. When the header has such a line, the molfile ends
there and the rest of the molecule reads as the record's data. This matches
how sendoff reads a `$$$$` header line: a header line is read as the line it
looks like.
