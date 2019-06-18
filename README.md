# Annotated Variant Exploration

This repository contains code for preparing extracted data from variant calling
for easy downstream analysis.


## Required Tools

- [PostgreSQL](https://www.postgresql.org/)
- [csvsql](https://csvkit.readthedocs.io/en/1.0.2/scripts/csvsql.html)
- [pgloader](https://github.com/dimitri/pgloader)

## Raw Data Sources
### anno_coding.csv

A TSV-formatted file output from Jim Knight's variant annotation pipeline.
Contains all coding variants. Has been processed by two scripts:
`extract_var_idx.py` and `anno_remove_duplicates.py` to extract variant 
information and to remove duplicates.

### anno_noncoding.csv

A TSV-formatted file output from Jim Knight's variant annotation pipeline.
Contains all noncoding variants. The two output files from Jim Knight's
pipeline have the same headers and will be combined in the database-resident
description of the annotations.

### merged_vcf_data.csv

Created by `PrepareAnnoGenes.R`, this file contains a flattened representation
of all VCFs called by the variant calling pipeline from scRNAseq. Some donor
metadata has been extracted from the filenames and included in the data frame.
NAs are represented by "NA" and the data is CSV-formatted.

### pt_metadata.csv

Manually-modified CSV-formatted file containing patient metadata. NAs are
represented by the empty string.

## Orchestration and Database
Starts up a PostgreSQL database with a local data directory and logfile,
loads data into the database, and builds indexes.

In order to load data appropriately into Postgres, schemas were generated
for each of the raw data sets gathered. From these automatically generated
schemas, some modifications were made manually to correct for errors made
during the autogeneration process due to sampling error in the schema
generation program (in the `csvkit` toolset).

