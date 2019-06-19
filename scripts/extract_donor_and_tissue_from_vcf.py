#!/usr/bin/env python3
#
# extract DonorNum and TissueType information from the merged_vcf data
#

import re
import sys

ifs = ","
ofs = ","

CELL_NAME_COL = 0
CELL_NAME_REGEX = '^"(\\d+)([BF])'
re_cell_name = re.compile(CELL_NAME_REGEX)

# process header
header_line = sys.stdin.readline()
nitems = len(header_line.split(ifs))
print(nitems, file=sys.stderr)

new_header_line = ofs.join(["DonorNum", "TissueType", header_line])
sys.stdout.write(new_header_line)

# process file
for line in sys.stdin:
    toks = line.split(ifs)

    ''' NOTE that this will fail
    if(len(toks) != nitems):
        print("Failed to correctly split line", file = sys.stderr)
        print("nitems {} != len(toks) {}".format(
            nitems, len(toks)), file=sys.stderr)
        print(line, file = sys.stderr)
        sys.exit(1)
    '''

    m = re_cell_name.match(toks[CELL_NAME_COL])

    donor_num = m.group(1)
    tissue_type = m.group(2)

    sys.stdout.write(ofs.join([donor_num, tissue_type, line]))

