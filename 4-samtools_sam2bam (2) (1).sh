#!/bin/bash

#SBATCH --job-name=SAM2BAM

#SBATCH --nodes=1
#SBATCH --ntasks=4       # modify this number to reflect how many cores you want to use (up to 24)
#SBATCH --time=04:00:00   # set time; default = 4 hours

#SBATCH --mail-type=END   # Keep these two lines of code if you want an e-mail sent to you when it is complete.
#SBATCH --mail-user=adh91@colostate.edu

#SBATCH --output=SAM2BAM_%j.txt  # this will capture all output in a logfile with %j as the job #


#export PATH=/projects/evcon@colostate.edu/bin:$PATH

pthread=$1
echo $pthead


# purge all existing modules
module purge

conda activate samtools-env

#record data
 start=`date +%s`
 
#turn sam into bam and sort
for file in *.sam; do name=$(basename ${file} .sam); echo ${name};
samtools view -S -b ${name}.sam > ${name}.bam
done


#SORT AFTER BAM is CREATED
