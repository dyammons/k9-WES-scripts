#!/bin/bash

#SBATCH --job-name=TrimWES

#SBATCH --nodes=1
#SBATCH --ntasks=4       # modify this number to reflect how many cores you want to use (up to 24)
#SBATCH --time=04:00:00   # set time; default = 4 hours
#SBATCH --mem=2GB
#SBATCH --mail-type=END   # Keep these two lines of code if you want an e-mail sent to you when it is complete.
#SBATCH --mail-user=adh91@colostate.edu

#SBATCH --output=TrimgaloreWES_%j.txt  # this will capture all output in a logfile with %j as the job #


# purge all existing modules
module purge

#activate conda environment
source activate mycustomenv

#Print TrimGalore version
trim_galore --version

#run Trimgalore
trim_galore \
    --paired BL45144_38_S38_L001_R1_001.fastq BL45144_38_S38_L001_R2_001.fastq \
    -q 30 \
    --gzip \
   # --fastqc \
