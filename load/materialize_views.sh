#!/usr/bin/env bash
#
# Builds additional indexes and materialized views for downstream analysis.
#

set -e
set -x

DB_URI=postgresql://localhost:5432/ms

echo "===== 'hq_variant' ====="

psql $DB_URI <<HERE
CREATE INDEX
    ON variant (
        "QD",
        "FS",
        "SOR",
        "MQ"
    );
HERE

# https://gatkforums.broadinstitute.org/gatk/discussion/6925/understanding-and-adapting-the-generic-hard-filtering-recommendations
QD_QUALITY_SCORE_CUTOFF=2
FS_QUALITY_SCORE_CUTOFF=60
SOR_QUALITY_SCORE_CUTOFF=3
MQ_QUALITY_SCORE_CUTOFF=40
psql $DB_URI <<HERE
CREATE MATERIALIZED VIEW hq_variant
AS
    SELECT * FROM variant
    WHERE
        "QD" > $QD_QUALITY_SCORE_CUTOFF
        AND "FS" < $FS_QUALITY_SCORE_CUTOFF
        AND "SOR" < $SOR_QUALITY_SCORE_CUTOFF
        AND "MQ" > $MQ_QUALITY_SCORE_CUTOFF;
HERE

psql $DB_URI <<HERE
CREATE INDEX
    ON hq_variant (
        "CHROM",
        "POS",
        "REF",
        "ALT"
    );
HERE

echo "===== 'deleterious_annotation' ====="

psql $DB_URI <<HERE
CREATE MATERIALIZED VIEW deleterious_annotation
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
            'stop_lost'
        )
        AND "Max. Population Frequency" < 0.05;
HERE

psql $DB_URI <<HERE
CREATE INDEX
    ON deleterious_annotation (
        "CHROM",
        "POS",
        "REF",
        "ALT"
    );
HERE

echo "===== 'hq_deleterious_variant' ====="

psql $DB_URI <<HERE
CREATE MATERIALIZED VIEW hq_deleterious_variant
AS
    SELECT * FROM hq_variant
    INNER JOIN deleterious_annotation
    USING
        ("CHROM", "POS", "REF", "ALT");
HERE

echo "===== 'hq_deleterious_variant_donor' ====="

psql $DB_URI <<HERE
CREATE MATERIALIZED VIEW hq_deleterious_variant_donor
AS
    SELECT * FROM hq_deleterious_variant
    INNER JOIN donor
    ON
        hq_deleterious_variant."DonorNum" = donor."Sample";
HERE

psql $DB_URI <<HERE
CREATE INDEX ON hq_deleterious_variant_donor (
        "Gene",
        "Group",
        "TissueType"
    );
HERE

psql $DB_URI <<HERE
CREATE MATERIALIZED VIEW variants_by_gene
AS
    SELECT 
        "Gene",
        "Group",
        "TissueType",
        COUNT(DISTINCT "CellName") as "NumCells",
        COUNT(DISTINCT "DonorNum") as "NumDonors"
    FROM 
        hq_deleterious_variant_donor
    GROUP BY
        "Gene",
        "Group",
        "TissueType"
HERE


psql $DB_URI <<HERE
CREATE MATERIALIZED VIEW all_merged
AS
    SELECT * FROM
      (SELECT * FROM annotation
        INNER JOIN hq_variant
        USING
          ("CHROM", "POS", "REF", "ALT")) as k
    INNER JOIN
      donor
    ON
      k."DonorNum" = donor."Sample";
HERE

psql $DB_URI <<HERE
CREATE INDEX ON all_merged (
    "Gene",
    "Top Consequence"
)
HERE




