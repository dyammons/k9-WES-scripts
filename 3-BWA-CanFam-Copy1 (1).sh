#!/bin/bash

#SBATCH --job-name=BAMalignWES

#SBATCH --nodes=1
#SBATCH --ntasks=4       # modify this number to reflect how many cores you want to use (up to 24)
#SBATCH --time=10:00:00   # set time; default = 4 hours
#SBATCH --mem=12GB
#SBATCH --mail-type=END   # Keep these two lines of code if you want an e-mail sent to you when it is complete.
#SBATCH --mail-user=adh91@colostate.edu

#SBATCH --output=BAMalignWES_%j.txt  # this will capture all output in a logfile with %j as the job #


# purge all existing modules
module purge

source activate mycustomenv

#record data
 start=`date +%s`
 
#align files 
for file in *_R1_001_val_1.fq.gz; do name=$(basename ${file} _R1_001_val_1.fq.gz); echo ${name};
bwa mem -t 6 \
/projects/adh91@colostate.edu/references/Ensembl.CanFam3.1/Canis_lupus_familiaris.CanFam3.1.dna.toplevel.fa \
${name}_R1_001_val_1.fq.gz ${name}_R2_001_val_2.fq.gz > ${name}_cfam3.1.sam
done
