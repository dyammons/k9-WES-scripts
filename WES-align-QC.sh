#!/usr/bin/env bash

####### MODIFY THIS SECTION #############

EXPERIMENT="creativeName"

### Path to input directory (modify if the files live in subdirectories) - can be absolute or relative
inputdir="/pwd/to/input/"

### Path to tmp directory (dir will be made if does not already exist)- can be absolute or relative
tmpDir="/scratch/alpine/$USER/tmp/"


### This is the path to the metadata -- can specify or it will be pulled from the sbatch file when you run this script
metadata=$1

### Path to reference genome
genomeFA="/projects/adh91@colostate.edu/references/Ensembl.CanFam3.1/Canis_lupus_familiaris.CanFam3.1.dna.toplevel.fa"

#do you want to run QC before trimming files?
runPREQC=FALSE

#if samples are split across two lanes, this will concatinate the files if set to TRUE
multiLane=TRUE

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


# this is lane 1 for each sample
sample_list1=( $(cut -f1 -d',' --output-delimiter=' ' $metadata) )

# this is lane 2 for each sample
sample_list2=( $(cut -f2 -d',' --output-delimiter=' ' $metadata) )

# this is the nickname to give the files
names=( $(cut -f3 -d',' --output-delimiter=' ' $metadata) )


####### PIPELINE ##############

# update the user
echo -e ">>> INPUT: This script will process files from the metafile:\n\t$metadata"
echo -e ">>> PLAN: This script will process the sample files into the following names: "
echo -e "\tSAMPLE1\tSAMPLE2\tNAMES"

for (( counter=0; counter < ${#sample_list1[@]}; counter++ ))
do
    echo -e "\t${sample_list1[$counter]}\t${sample_list2[$counter]}\t${names[$counter]}"
done


#### STEP0: before starting, concatinate files based on metadata file
if [ $runPREQC == TRUE ]
then
  echo -e "\n>>> cat: concatinating files"
  for (( counter=0; counter < ${#sample_list1[@]}; counter++ )) ### <<<<< add an echo step of the actual command
  do
    samplename=${names[$counter]}
    sample1=${sample_list1[$counter]}
    sample2=${sample_list2[$counter]}
  
    cmd0="cat ${inputdir}/${sample1} ${inputdir}/${sample2} > ${inputdir}/${samplename}.fastq"
  
    echo -e "\t$ ${cmd0}"
    time eval $cmd0
  done
  
else
  echo -e "\n>>> skipping cat as experiment was not split across lanes"
fi


if [ $runPREQC == TRUE ]
then
    #### STEP1: fastqc to determine quality & trim reads -- be aware this will run all files in the input directory; it pretty quick though
    echo -e "\n>>> FASTQC: generating quality report"
    mkdir -p $outputdir"01_fastqc_pre"    

    # execute fastqc pre trim
    cmd1="fastqc -o ${outputdir}01_fastqc_pre -t $pthread ${inputdir}/*.fastq"

    echo -e "\t$ ${cmd1}"
    time eval $cmd1
else

    echo -e "\n>>> skipping pre FASTQC and moving to trim"

fi



#### STEP2: Trimgalore to trim reads
echo -e "\n>>> trim_galore: Trimgalore reads to size"
mkdir -p $outputdir"02_trim_galore"

for (( counter=0; counter < ${#sample_list1[@]}; counter++ ))
do

    samplename=${names[$counter]}
    sample1=${sample_list1[$counter]}
    sample2=${sample_list2[$counter]}
    
    # execute trim_galore
    cmd2="trim_galore --paired $inputdir/$sample1 $inputdir/$sample2 \
    -o $outputdir"02_trim_galore" \
    --basename $samplename \
    -q 30 &"
    
    echo -e "\t$ ${cmd2}"
    time eval $cmd2
done

wait
# note: the ampersand at the end of cmd2 and the 'wait' command enable parallelization


# BWA to align to the genome
echo -e "\n>>> BWA: aligning each sample to the genome"
outBWA=$outputdir"03_bwa/"
mkdir -p $outBWA

for (( counter=0; counter < ${#samples1[@]}; counter++ ))
do
    samplename=${names[$counter]}
    sample1=${samples1[$counter]}
    sample2=${samples2[$counter]}


    ## execute BWA
    cmd3="bwa mem -t $pthread \
    $genomeFA \
    ${samplename}_val_1.fq ${samplename}_val_2.fq > ${outBWA}${samplename}_cfam3.1.sam"

    echo -e "\t$ ${cmd3}"
    time eval $cmd3

done


# SAMTOOLS and BAMCOVERAGE: to convert .sam output to uploadable .bam and .wg files
echo -e "\n>>> SAMTOOLS/BAMCOVERAGE: to convert files to uploadable _sort.bam and _sort.bam.bai files:"
samout=$outputdir"03_samtools/"
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
echo -e "\n>>> FASTQC VERSION:"
$fastqc --version
echo -e "\n>>> trim_galore VERSION:"
$trim_galore --version
echo -e "\n>>> BWA VERSION:"
$bwa --version
echo -e "\n>>> SAMTOOLS VERSION:"
$samtools --version

echo -e ">>> END: Analayzer complete."




