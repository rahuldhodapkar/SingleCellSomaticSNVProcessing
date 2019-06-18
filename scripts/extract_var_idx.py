#!/usr/bin/env python3
#
# Process output of Knight annotation pipeline to extract
# (CHROM, POS, REF, ALT) tuple for joining with VCF data
#
# Applied to the annotation data before loading into postgres
#

import re
import sys

ifs = "\t"
ofs = "\t"

CHROM_POS_COL = 2
CHROM_POS_REGEX = "(.+):(.+)"
re_chrom_pos = re.compile(CHROM_POS_REGEX)

CHANGE_COL= 3
CHANGE_REGEX = ".+:(.+)>(.+)"
re_change = re.compile(CHANGE_REGEX)

# process header

header_line = sys.stdin.readline()
nitems = len(header_line.split(ifs))

new_header_line = ofs.join(["CHROM", "POS", "REF", "ALT", header_line])
sys.stdout.write(new_header_line)

# process file
for line in sys.stdin:
    toks = line.split(ifs)

    if(len(toks) != nitems):
        print("Failed to correctly split line", file = sys.stderr)
        print(line, file = sys.stderr)
        sys.exit(1)

    m = re_chrom_pos.match(toks[CHROM_POS_COL])
    chrom = m.group(1)
    pos = m.group(2)

    if (m.start() != 0 or m.end() != len(toks[CHROM_POS_COL])):
        print("Invalid chrom:pos regex", file = sys.stderr)
        print(line, file = sys.stderr)
        sys.exit(1)

    m = re_change.match(toks[CHANGE_COL])
    ref = m.group(1)
    alt = m.group(2)

    if (m.start() != 0 or m.end() != len(toks[CHANGE_COL])):
        print("Invalid change detection regex", file = sys.stderr)
        print(line, file = sys.stderr)
        sys.exit(1)

    sys.stdout.write(ofs.join([chrom, pos, ref, alt] + toks))













