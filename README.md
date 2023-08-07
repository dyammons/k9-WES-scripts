# k9-WES-scripts

howdy, here is a summary of how the scripts are designed to work.
(Note: currently the scripts will only get you you through bwa mem and subsequent samtools steps.)

## Step 1: Organize, QC, and trim data

The first step uses the `WES-QC-trim.sh` script. This file should be run as a jobs with 8 threads with an estimated run time of 15 minutes per sample.  
So, with 100 samples you're looking at a 3.5 hr run time when you factor in 8 samples running at once.
(this is a rough estimate that should be updated as you test/use the script)

To execute this script
