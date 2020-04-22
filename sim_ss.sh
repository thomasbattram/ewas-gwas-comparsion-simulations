#!/bin/bash
#
#
#PBS -t 1-10
#PBS -l nodes=1:ppn=2,walltime=200:00:00
# Define working directory
export WORK_DIR=$HOME
# Change into working directory
cd $WORK_DIR
# Execute code
#Set -e allows you to test the script
set -e

echo "Running on ${HOSTNAME}"
if [ -n "${1}" ]; then
  echo "${1}"
  PBS_ARRAYID=${1}
fi

i=${PBS_ARRAYID}

# 
module rm languages/R-3.3.3-ATLAS
module load languages/R-3.5.1-ATLAS-gcc-6.1

splits=("1,100" "101,200" "201,300" "301,400" "401,500" "501,600" "601,700" "701,800" "801,900" "901,1000")
echo ${splits[$i]}

time Rscript simulations.R ${splits[$i]}