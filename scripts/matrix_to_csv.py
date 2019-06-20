#!/usr/bin/env python3
#
# Convert from .matrix expression format to a database_friendly CSV format.
#
# Usage:
#           ./matrix_to_csv.py *.matrix
#
# NOTE: special handling for ERCC spike-in expression.
#       this script is not very portable, and suitable *only* for processing Alex Kitz'
#       Fat Treg data in its current form.
#

import sys
import os
import re

FS="\t"
OFS="\t"

CELL_NAME_REGEX = r'^"scratch60/ercc/\w+/([a-zA-Z0-9]+)_'
CELL_INFO_FROM_NAME_REGEX = r'^(\d+)([BF])(?:Treg|Index)(\d+)'
GENE_STRING_REGEX = r'^"(ENSG\d+)_(\w+)'
ERCC_REGEX = r'^"ERCC'

cell_name_re = re.compile(CELL_NAME_REGEX)
cell_info_from_name_re = re.compile(CELL_INFO_FROM_NAME_REGEX)
gene_string_re = re.compile(GENE_STRING_REGEX)
ercc_re = re.compile(ERCC_REGEX)

#################################################################
## HELPER FUNCTIONS
#################################################################

def extract_gene_and_ensembl(s):
    """Extract gene name and ensembl code from gene string"""

    is_ensembl = True
    m = gene_string_re.search(s)

    if not m:
        if not ercc_re.search(s):
            print("ERROR: unrecognized gene string \"{}\" encountered"
                                .format(s), file=sys.stderr)
            sys.exit(1)
        
        is_ensembl = False

    if not is_ensembl:
        # strip quotes for ERCC 
        ensembl_gene = '.'
        gene = re.sub(r'^"|"$', '', s)            
    else:
        ensembl_gene = m.group(1)
        gene = m.group(2)

    return (ensembl_gene, gene)

def output_line(cell=None, donor_num='.', tissue_type='.', cell_num='.',
                                     ensembl_gene='.', gene='.', expression="0.0"):
    """Process raw data for output to OFS-delim format"""

    if (float(expression) == 0):
        print("WARN: trying to write output with expression == 0", file=sys.stderr)

    if (expression.endswith("\n")):
        expression = expression[:-1] # trim trailing newline

    print(OFS.join([
        cell, donor_num, tissue_type, cell_num, ensembl_gene, gene, expression
    ]))

def extract_cell_name_data(cell_file_name):
    """Extract cell name data """

    m = cell_name_re.search(cell_file_name)

    if not m:
        print("ERROR: trying to extract bad cell name {}"
                    .format(cell_file_name), file=sys.stderr)
        sys.exit(1)

    cell_name = m.group(1)

    cell_name_m = cell_info_from_name_re.search(cell_name)
    ( donor_num, tissue_type, cell_num ) = cell_name_m.groups()

    return (cell_name, donor_num, tissue_type, str(int(cell_num)))

def process_file(f):
    """Process a and print formatted output lines"""

    header = f.readline()
    if header.startswith("\t"):
        header = header[1:]
    cell_file_names = header.split(FS)

    map(list, zip(*[(1, 2), (3, 4), (5, 6)]))

    [cell_names, donor_nums, tissue_types, cell_nums] = map(list, zip(
        *[ extract_cell_name_data(x) 
           for x in cell_file_names
        ]
    ))

    for line in f:

        toks = line.split(FS)
        gene_string = toks[0]

        (ensembl_gene, gene) = extract_gene_and_ensembl(gene_string)

        expr_vals = toks[1:len(toks)]

        for i in range(len(expr_vals)):
            if float(expr_vals[i]) == 0:
                continue

            # non-zero value
            output_line(cell=cell_names[i],
                        donor_num=donor_nums[i],
                        tissue_type=tissue_types[i],
                        cell_num=cell_nums[i],
                        ensembl_gene=ensembl_gene,
                        gene=gene,
                        expression=expr_vals[i])

#################################################################
## MAIN
#################################################################

# Print output file header.
print(OFS.join([
    "CellName", "DonorNum", "TissueType", "CellNum", "EnsemblGene", "Gene", "Expression"
]))

files = sys.argv[1:len(sys.argv)]

for filename in files:
    with open(filename, 'r') as f:
        process_file(f)

