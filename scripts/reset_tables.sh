#!/usr/bin/env bash
#
# Convenience script for resetting tables

DB_URI=postgresql://localhost:5432/ms

psql $DB_URI <<HERE
DROP TABLE donor;
DROP TABLE annotation;
DROP TABLE variant;
HERE
