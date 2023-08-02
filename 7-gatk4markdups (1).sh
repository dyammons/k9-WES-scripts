#!/bin/bash

#SBATCH --job-name=MarkDups

#SBATCH --nodes=1
#SBATCH --ntasks=4       # modify this number to reflect how many cores you want to use (up to 24)
#SBATCH --time=24:00:00   # set time; default = 4 hour
#SBATCH --mail-type=END   # Keep these two lines of code if you want an e-mail sent to you when it is complete.
#SBATCH --mail-user=adh91@colostate.edu

#SBATCH --output=MarkDups_%j.txt  # this will capture all output in a logfile with %j as the job #
export PATH=/projects/evcon@colostate.edu/bin:$PATH

pthread=$1
echo $pthead

# purge all existing modules
module purge

#record data
 start=`date +%s`
 
#turn sam into bam and sort
for file in *.bam; do name=$(basename ${file} .bam); echo ${name};
gatk MarkDuplicates -I ${name}.bam \
-O ${name}.markdup.bam \
--TMP_DIR tmp \
--METRICS_FILE ${name}_markdups.txt
done
