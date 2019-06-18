#!/usr/bin/env bash
#
# Starts the PostgreSQL database server locally
#

set -e
set -x

LOCAL_DATA_DIR=db/data
LOCAL_LOG_DIR=db/log

mkdir -p $LOCAL_DATA_DIR
mkdir -p $LOCAL_LOG_DIR

initdb -D $LOCAL_DATA_DIR
postgres -D $LOCAL_DATA_DIR >$LOCAL_LOG_DIR/db.log 2>&1 &

# wait for the server to start
sleep 5

createdb ms