#!/usr/bin/env python3
#
# To be used to remove duplicates from raw data file before input processing.
# Removes duplicates based on (CHROM, POS, REF, ALT)

import sys

ifs = "\t"

CHROM_COL=0
POS_COL=1
REF_COL=2
ALT_COL=3

header = sys.stdin.readline()
sys.stdout.write(header)

printed = {}
for line in sys.stdin:
    toks = line.split(ifs)
    if (toks[CHROM_COL], toks[POS_COL], toks[REF_COL], toks[ALT_COL]) in printed:
        continue

    sys.stdout.write(line)
    printed[(toks[CHROM_COL], toks[POS_COL], toks[REF_COL], toks[ALT_COL])] = 1
