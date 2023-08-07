#!/bin/bash

#SBATCH --job-name=BQSR

#SBATCH --nodes=1
#SBATCH --ntasks=4       # modify this number to reflect how many cores you want to use (up to 24)
#SBATCH --time=2-00:00:00   # set time; default = 4 hours
#SBATCH --partition=amilan
#SBATCH --qos=long
#SBATCH --mem=64GB
#SBATCH --mail-type=END   # Keep these two lines of code if you want an e-mail sent to you when it is complete.
#SBATCH --mail-user=adh91@colostate.edu

#SBATCH --output=BQSR_%j.txt  # this will capture all output in a logfile with %j as the job #

export PATH=/projects/evcon@colostate.edu/bin:$PATH

pthread=$1
echo $pthead

# purge all existing modules
module purge

#record data
 start=`date +%s`
 
#BQSR
for file in *.bam; do name=$(basename ${file} .bam); echo ${name};
gatk ApplyBQSR \
-I ${name}.bam \
-R /projects/adh91@colostate.edu/references/Ensembl.CanFam3.1/Canis_lupus_familiaris.CanFam3.1.dna.toplevel.fa \
--bqsr-recal-file ${name}.recal_data.table \
--tmp-dir tmp \
-O ${name}.bqsr.bam
done
