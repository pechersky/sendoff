"""Keep legacy API probes independent of chemistry-toolkit fixture generation."""

V2000 = (
    "  literal title  \n"
    "  literal source  \n"
    " literal comment \n"
    "  2  1  0  0  0  0  0  0  0  0999 V2000\n"
    "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    1.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "  1  2  1  0  0  0  0\n"
    "M  END\n"
)

V3000 = (
    "  literal title  \n"
    "  literal source  \n"
    " literal comment \n"
    "  0  0  0  0  0  0  0  0  0  0999 V3000\n"
    "M  V30 BEGIN CTAB\n"
    "M  V30 COUNTS 2 1 0 0 0\n"
    "M  V30 BEGIN ATOM\n"
    "M  V30 1 C 0 0 0 0 CHG=1\n"
    "M  V30 2 O 1 0 0 0\n"
    "M  V30 END ATOM\n"
    "M  V30 BEGIN BOND\n"
    "M  V30 1 1 1 2 CFG=1\n"
    "M  V30 END BOND\n"
    "M  V30 END CTAB\n"
    "M  END\n"
)
