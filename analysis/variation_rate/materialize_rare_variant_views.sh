#!/usr/bin/env bash
#
# Builds additional indexes and materialized views for downstream analysis.
#

set -e
set -x

DB_URI=postgresql://localhost:5432/ms


echo "===== 'rare_annotation' ====="

psql $DB_URI <<HERE
CREATE MATERIALIZED VIEW rare_annotation
AS
    SELECT * FROM annotation
    WHERE
        "Top Consequence" IN (
            'frameshift_variant',
            'incomplete_terminal_codon_variant',
            'protein_altering_variant',
            'splice_acceptor_variant',
            'splice_donor_variant',
            'start_lost',
            'stop_gained',
            'stop_lost',
            'missense_variant'
        )
        AND "Max. Population Frequency" = 0;
HERE

psql $DB_URI <<HERE
CREATE INDEX
    ON rare_annotation (
        "CHROM",
        "POS",
        "REF",
        "ALT"
    );
HERE

echo "===== 'hq_rare_variant' ====="

psql $DB_URI <<HERE
CREATE MATERIALIZED VIEW hq_rare_variant
AS
    SELECT * FROM hq_variant
    INNER JOIN rare_annotation
    USING
        ("CHROM", "POS", "REF", "ALT");
HERE

echo "===== 'hq_rare_variant_donor' ====="

psql $DB_URI <<HERE
CREATE MATERIALIZED VIEW hq_rare_variant_donor
AS
    SELECT * FROM hq_rare_variant
    INNER JOIN donor
    ON
        hq_rare_variant."DonorNum" = donor."Sample";
HERE

psql $DB_URI <<HERE
CREATE INDEX ON hq_rare_variant_donor (
        "Gene",
        "Group",
        "TissueType"
    );
HERE

psql $DB_URI <<HERE
CREATE MATERIALIZED VIEW rare_variants_by_gene
AS
    SELECT 
        "Gene",
        "Group",
        "TissueType",
        COUNT(DISTINCT "CellName") as "NumCells",
        COUNT(DISTINCT "DonorNum") as "NumDonors"
    FROM 
        hq_rare_variant_donor
    GROUP BY
        "Gene",
        "Group",
        "TissueType"
HERE






