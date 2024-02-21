#!/bin/bash

#SBATCH --job-name=SAMTFlag

#SBATCH --nodes=1
#SBATCH --ntasks=4       # modify this number to reflect how many cores you want to use (up to 24)
#SBATCH --time=20:00:00   # set time; default = 4 hours
#SBATCH --mail-type=END   # Keep these two lines of code if you want an e-mail sent to you when it is complete.
#SBATCH --mail-user=adh91@colostate.edu

#SBATCH --output=SAMTFlag_%j.txt  # this will capture all output in a logfile with %j as the job #

#export PATH=/projects/evcon@colostate.edu/bin:$PATH

pthread=$1
echo $pthead


# purge all existing modules
module purge

#record data
 start=`date +%s`

source activate samtools-env
 
#Perform after sorting
for file in *.bam; do name=$(basename ${file} .bam); echo ${name};
samtools flagstat ${name}.bam > ${name}.fs.txt
done


