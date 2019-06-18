#!/usr/bin/env bash
#
# Build indexes on the tables created with raw data in Postgres
#

set -e
set -x

DB_URI=postgresql://localhost:5432/ms

echo "===== 'donor' indexes ====="

psql $DB_URI <<HERE
CREATE UNIQUE INDEX 
    ON donor (
        "Sample"
    );
HERE

echo "===== 'annotation' indexes ====="

psql $DB_URI <<HERE
CREATE UNIQUE INDEX 
    ON annotation (
        "CHROM",
        "POS",
        "REF",
        "ALT"
    );
HERE

echo "===== 'variant' indexes ====="

psql $DB_URI <<HERE
CREATE UNIQUE INDEX 
    ON variant (
        "CHROM",
        "POS",
        "REF",
        "ALT",
        "CellName"
    );
HERE

