#!/usr/bin/env bash

PG_PORT=:5432
PG_PROCS=`lsof -i $PG_PORT | grep LISTEN | awk '{ print $2 }' | uniq`

if [ ! -z $PG_PROCS ]
then
    kill $PG_PROCS
    sleep 5
fi

if [ -e db ]
then
    rm -rf db
fi
