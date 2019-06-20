#!/usr/bin/env bash
#
# create tables

set -e
set -x

DB_URI=postgresql://localhost:5432/ms

echo "===== CREATE TABLES ====="

psql $DB_URI -f ./schema/donor.sql
psql $DB_URI -f ./schema/annotation.sql
psql $DB_URI -f ./schema/variant.sql
psql $DB_URI -f ./schema/expression.sql