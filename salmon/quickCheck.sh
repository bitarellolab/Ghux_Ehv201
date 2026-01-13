#!/usr/bin/bash

comp2=$1
notcomp=$(sacct|grep ${comp2}|grep RUNNING|grep -v batch|wc -l)

while [ $notcomp -gt 2 ]
    do
    comp=$(sacct|grep ${comp2}|grep COMPLETED|grep -v batch|wc -l)
    notcomp=$(sacct|grep ${comp2}|grep RUNNING|grep -v batch|wc -l)
    echo "${comp}" "have completed"
    echo "${notcomp}" "have NOT completed"
    sleep 10
    comp=$(sacct|grep 13668|grep RUNNING|grep -v batch|wc -l)
done


