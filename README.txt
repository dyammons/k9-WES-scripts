#This document was generated to descrive the variant calling pipeline. 

#File structure:
All files are processed in a directory devoted to WES as such
/scratch/alpine/$USER/WES

Subfolders:
01_rawdata
02_alignment
03_Mutect


Scripts 1-2 are run in 01_rawdata folder

Post-trimming, trimmed fastq files are moved into the 02_alignment folder.

Scripts 3-5 are then performed. Script 5 ends with sorted BAM files

The sorted BAM files are then moved into a subfolder in 03_Mutect and scripts 6-9 are then performed sequentially.

Eventually, this creates BAM files that are ready to be processed in Mutect (the variant calling software throught the Broad). 




