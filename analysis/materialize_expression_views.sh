#!/usr/bin/env bash
#
# Builds additional indexes and materialized views for downstream analysis.
#

set -e
set -x

DB_URI=postgresql://localhost:5432/ms

echo "===== 'expression_donor' ====="

psql $DB_URI <<HERE
CREATE MATERIALIZED VIEW expression_donor_unclean
AS
    SELECT * FROM
        expression
    INNER JOIN
        donor
    ON
        expression."DonorNum" = donor."Sample"
HERE

psql $DB_URI <<HERE
CREATE MATERIALIZED VIEW expression_donor
AS
    SELECT * FROM
        (SELECT * FROM
          expression
        INNER JOIN
          donor
        ON
          expression."DonorNum" = donor."Sample") as a
    WHERE
        "DonorNum" NOT IN (9,10)
HERE

psql $DB_URI <<HERE
CREATE INDEX ON expression_donor (
    "EnsemblGene",
    "Group",
    "TissueType"
)
HERE

echo "===== 'expression_by_gene' ====="

psql $DB_URI <<HERE
CREATE MATERIALIZED VIEW expression_by_gene
AS
    SELECT 
        "EnsemblGene",
        MAX("Gene") as "HRGeneName",
        "Group",
        "TissueType",
        COUNT(DISTINCT "CellName") as "NumCells",
        COUNT(DISTINCT "DonorNum") as "NumDonors"
    FROM
        expression_donor
    GROUP BY
        "EnsemblGene",
        "Group",
        "TissueType"
HERE