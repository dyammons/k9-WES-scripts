# k9-WES-scripts

howdy, here is a summary of how the scripts are designed to work.  
(Note: currently the scripts will only get you you through bwa mem and subsequent samtools steps.)

## Step 1: Organize, QC, and trim data

The first step uses the `WES-QC-trim.sh` script. This file should be run as a jobs with 8 threads with an estimated run time of 15 minutes per sample.  
So, with 100 samples you're looking at a 3.5 hr run time when you factor in 8 samples running at once.  
(Note: this is a rough estimate that should be updated as you test/use the script)

To execute this script, 
- update the *MODIFY THIS SECTION* in the `WES-QC-trim.sh` script
- generate a samples.txt metadat file that lists all the sample names (from the scripts folder `(cd ../01_input && ls *.fastq) | cut -d"_" -f1 | sort -u > samples.txt` should generate the desired list)
- check that the looks correct

Expected output:
```sh
$ head samples.txt 
sample_1
sample_2
sample_3
sample_4
sample_5
sample_6
sample_7
sample_8
sample_9
sample_10
```

Then run with `sbatch cute-WES-QC-trim.sbatch`.

## Step 2: Align the data

After files are trimmed and QC reports generated, the data are ready for alignment to the genome.  
This is completed using the `mkBatch-align.sh` which will divide samples into groups and provide code to copy-and-paste into your terminal to get multiple jobs queued at one.  

To run, set the preferences in the `mkBatch-align.sh file:
```sh
##### set user preferences #####
samplesPerJob=10            ### number of samples to groups together for a batch
numNode=1                   ### this will likely always be set to 1
numTasks=20                 ### number of threads - I have seen anywhere from 20 to 64 reccomended- I belive 20 will yield adequate job efficeincy
email="adh91@colostate.edu"
runTime="24:00:00"          ### runtime allocation for each batch of jobs -- this may need to be modififed
partition="amilan"          ### partition to use
```

Then, in your terminal, run:
```sh
bash mkBatch-align.sh samples.txt
```

This will generate multiple sbatch scripts and a file called `jobList.txt`.  
The contents in `jobList.txt` can be copied directly into a terminal to que all of your alignemnt jobs at once.

## Step 3: Mutect for downstream fun

Pending...

