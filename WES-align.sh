#!/usr/bin/env bash

####### MODIFY THIS SECTION #############

#this should be the same as the trim step for seemless useage
EXPERIMENT="creativeName"

### Path to input directory (modify if the files live in subdirectories) - can be absolute or relative
### Input files will be the trimmed files
inputdir="/pwd/to/input/"

### Path to tmp directory (dir will be made if does not already exist)- can be absolute or relative
tmpDir="/scratch/alpine/$USER/tmp/"

### This will be populated by the mkBatch-align.sh script
sampleList=$1

### Path to reference genome
genomeFA="/projects/adh91@colostate.edu/references/Ensembl.CanFam3.1/Canis_lupus_familiaris.CanFam3.1.dna.toplevel.fa"


### Set the output directory
outputdir="../03_output/"$EXPERIMENT"_output/"

########## DONE MODIFYING ###############



########## BEGIN CODE ###############

echo -e ">>> INITIATING analyzer with command:\n\t$0 $@"

# stash number of threads to used in sbatch job submission:
pthread=$2

# create main output directory
echo -e ">>> MAKING main output directory"
echo -e "\tmkdir $outputdir"
mkdir -p $outputdir
mkdir -p $tmpDir #also make tmp dir if not already in existance


####### META DATA #############

# this is the nickname to give the files
names=( $(cut -f1 --output-delimiter=' ' $sampleList) )


####### PIPELINE ##############

# update the user
echo -e ">>> INPUT: This script will process files from the metafile:\n\t$metadata"
echo -e ">>> PLAN: This script will process the sample files into the following names: "
echo -e "NAMES"

for (( counter=0; counter < ${#names[@]}; counter++ ))
do
    echo -e "${names[$counter]}"
done


# BWA to align to the genome
echo -e "\n>>> BWA: aligning each sample to the genome"
outBWA=$outputdir"03_bwa/"
mkdir -p $outBWA

for (( counter=0; counter < ${#names[@]}; counter++ ))
do
    samplename=${names[$counter]}

    ## execute BWA
    cmd3="bwa mem -t $pthread \
    $genomeFA \
    ${samplename}_val_1.fq ${samplename}_val_2.fq > ${outBWA}${samplename}_cfam3.1.sam"

    echo -e "\t$ ${cmd3}"
    time eval $cmd3

done


# SAMTOOLS and BAMCOVERAGE: to convert .sam output to uploadable .bam and .wg files
echo -e "\n>>> SAMTOOLS/BAMCOVERAGE: to convert files to uploadable _sort.bam and _sort.bam.bai files:"
samout=$outputdir"04_samtools/"
mkdir -p $samout

for seqname in ${names[@]}
do
    # echo
    echo -e "\tSamtools and BamCoverage convert: ${seqname}"
    
    # Samtools: compress .sam -> .bam
    cmd4="samtools view --threads $pthread -bS ${outBWA}${seqname}_cfam3.1.sam > ${samout}${seqname}.bam"
    
    echo -e "\t$ ${cmd4}"
    time eval $cmd4

    # Samtools: sort .bam -> _sort.bam
    cmd5="samtools sort --threads $pthread -o ${samout}${seqname}_sort.bam --reference $genomefa ${samout}${seqname}.bam"
    
    echo -e "\t$ ${cmd5}"
    time eval $cmd5

    #is this needed?
    #samtools flagstat --threads $pthread -o ${samout}${seqname}.fs.txt --reference ${samout}${seqname}.bam
    
done

######## VERSIONS #############
echo -e "\n>>> VERSIONS:"
echo -e "\n>>> BWA VERSION:"
$bwa --version
echo -e "\n>>> SAMTOOLS VERSION:"
$samtools --version

echo -e ">>> END: Analayzer complete."




