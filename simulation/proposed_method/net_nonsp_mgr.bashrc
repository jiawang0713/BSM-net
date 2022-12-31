#!/bin/bash

for loop in `seq 1 50`;
do
qsub net_nonsp_worker.PBS -v "loop=$loop"
done

