#!/bin/bash

#SBATCH --job-name=Mutect2Norm
#SBATCH --nodes=1
#SBATCH --ntasks=4       # modify this number to reflect how many cores you want to use (up to 24)
#SBATCH --time=2-02:00:00   # set time; default = 4 hours
#SBATCH --mem=64GB
#SBATCH --partition=amilan  # modify this to reflect which queue you want to use. Either 'shas' or 'shas-testing'
#SBATCH --qos=long      # modify this to reflect which queue you want to use. Options are 'normal' and 'testing'

#SBATCH --mail-type=END   # Keep these two lines of code if you want an e-mail sent to you when it is complete.
#SBATCH --mail-user=adh91@colostate.edu

#SBATCH --output=Mutect2Norm_%j.txt  # this will capture all output in a logfile with %j as the job #

export PATH=/projects/evcon@colostate.edu/bin:$PATH

pthread=$1
echo $pthead

# purge all existing modules
module purge

#record data
 start=`date +%s`
 
#Step 1. Run Mutect2 in tumor-only mode for each normal sample.
for file in *_cfam3.1_sorted.markdup.RG.bqsr.bam; do name=$(basename ${file} _cfam3.1_sorted.markdup.RG.bqsr.bam); echo ${name};
gatk Mutect2 \
-R /projects/adh91@colostate.edu/references/Ensembl.CanFam3.1/Canis_lupus_familiaris.CanFam3.1.dna.toplevel.fa \
--germline-resource /scratch/alpine/adh91@colostate.edu/WES/03_GATK/rename.broad_umass_canid_variants.1.2.vcf.gz \
-max-mnp-distance 0 \
-I ${name}_cfam3.1_sorted.markdup.RG.bqsr.bam \
-O ${name}.vcf.gz
done
