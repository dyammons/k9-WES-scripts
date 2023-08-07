#!/bin/bash

#SBATCH --job-name=WESFastQC

#SBATCH --nodes=1
#SBATCH --ntasks=12       # modify this number to reflect how many cores you want to use (up to 24)
#SBATCH --time=24:00:00   # set time; default = 4 hours

#SBATCH --mail-type=END   # Keep these two lines of code if you want an e-mail sent to you when it is complete.
#SBATCH --mail-user=adh91@colostate.edu

#SBATCH --output=WESFastQC_%j.txt  # this will capture all output in a logfile with %j as the job #


# purge all existing modules
module purge

#Set path
DATADIR=/scratch/summit/adh91@colostate.edu/WES/01_input_raw_data

#Record date
today=$(date +"%Y-%m-%d")

#FASTQC version
fastqc --version

#Run Fastqc
fastqc *.fastq
