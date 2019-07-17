#!/usr/bin/env bash
#
# After schema generation, copy raw data to server for analysis
#

set -e
set -x

CUR_DIR=`pwd`

DB_URI=postgresql://localhost:5432/somatic

echo "===== COPY DATA ====="

psql $DB_URI <<HERE
\\copy donor FROM '$CUR_DIR/raw_data/pt_metadata.csv' DELIMITER ',' CSV HEADER
HERE

psql $DB_URI <<HERE
\\copy annotation FROM '$CUR_DIR/raw_data/anno_all.csv' DELIMITER E'\\t' CSV HEADER NULL as '.'
HERE

psql $DB_URI <<HERE
\\copy variant FROM '$CUR_DIR/raw_data/merged_vcf_data.csv' DELIMITER ',' CSV HEADER NULL as 'NA'
HERE

psql $DB_URI <<HERE
\\copy expression FROM '$CUR_DIR/raw_data/expression.tsv' DELIMITER E'\\t' CSV HEADER NULL as '.'
HERE