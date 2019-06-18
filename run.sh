#!/usr/bin/env bash
#
# Runs all scripts to start and build database server
#

set -e
set -x

echo "===== Cleaning Database Directories ====="
./load/clean_server.sh

echo "===== Starting Database Service ====="
./load/start_server.sh

echo "===== Listing Database Ports ====="
lsof -i :5432

echo "===== Loading Raw Data ====="
./load/copy_data.sh

echo "===== Define Database Schema ====="
./load/define_schemas.sh

echo "===== Building Database Indexes ====="
./load/build_indexes.sh


