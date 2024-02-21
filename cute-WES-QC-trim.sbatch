#!/usr/bin/env bash

#SBATCH --job-name=WES_QC_trim
#SBATCH --ntasks=8
#SBATCH --nodes=1
#SBATCH --time=04:00:00

#SBATCH --partition=amilan
#SBATCH --qos=normal

#SBATCH --mail-type=END
#SBATCH --mail-user=adh91@colostate.edu ### change to your email ###

#SBATCH --output=WES_QC_trim_%j.log  #modify as desired - will output a log file where the "%j" inserts the job ID number

######### Instructions ###########
#below, you don't need to change $SLURM_NTASKS. It will automatically populate whatever you put in --ntasks=# above.

#run WES-QC-trim script
### NOTE: make sure samples.txt is a file that lists all of the sample names present in the input directory that you wish to run!
#the following command should generate the list if run in the scripts folder: (cd ../01_input && ls *.fastq) | cut -d"_" -f1 | sort -u > samples.txt
bash WES-QC-trim.sh samples.txt $SLURM_NTASKS
