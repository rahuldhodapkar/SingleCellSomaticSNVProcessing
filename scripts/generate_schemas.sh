#!/usr/bin/env bash
#
# Load Raw Data into PostgreSQL database
#

set -e
set -x

if [ ! -e schema ]
then
    mkdir schema
fi

csvsql -i postgresql --tables "donor"\
            ./raw_data/pt_metadata.csv \
            > ./schema/donor.sql

csvsql -i postgresql --tabs --tables "annotation" --no-constraints \
            <(head -n 10000 ./raw_data/all_merged_anno.coding.xls) \
            > ./schema/annotation.sql

csvsql -i postgresql --tables "variant" --no-constraints \
            <(head -n 1000 ./raw_data/merged_vcf_data.csv) \
            > ./schema/variant.sql
