#!/usr/bin/env bash

set -e
set -x

postgres -D db/data >>db/log/db.log 2>&1 &
