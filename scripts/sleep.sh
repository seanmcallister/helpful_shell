#!/bin/bash

#SBATCH --job-name=res
#SBATCH --ntasks=24
#SBATCH --mem=400000
#SBATCH --nodelist=biomix25

hostname
date

sleep 10000000

[ $? -eq 0 ] || echo "JOB FAILURE: code $?"

date
echo "done!"
