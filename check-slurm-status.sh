#!/usr/bin/bash

##Usage: ./check-slurm-status <jobid>

job=$1

#job=13748
#sacct|grep ${job}


compl=$(sacct|grep ${job}|grep COMPLETED|grep -v batch|wc -l)
runn=$(sacct|grep ${job}|grep RUNNING|grep -v batch|wc -l)
pend=$(sacct|grep ${job}|grep PENDING|grep -v batch|wc -l)


echo "Completed: ${compl}"
echo "Running: ${runn}"
echo "Pending: ${pend}"

#sacct|grep ${job}|grep COMPLETED|grep -v batch|wc



