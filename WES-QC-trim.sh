#!/usr/bin/env bash

####### MODIFY THIS SECTION #############
EXPERIMENT="creativeName"

### Path to input directory (modify if the files live in subdirectories) - can be absolute or relative
inputdir="/pwd/to/input/"

### Path to tmp directory (dir will be made if does not already exist)- can be absolute or relative
tmpDir="/scratch/alpine/$USER/tmp/"

### Ensure this point to a list of all sample names (Note: $1 is the first term after the script when calling `bash`.
sampleList=$1

#do you want to run QC before trimming files? (reccomend FALSE)
runPREQC=FALSE

#do you want to run QC after trimming files? (reccomend TRUE)
runPOSTQC=TRUE

### Set the output directory
outputdir="../03_output/"$EXPERIMENT"_output/"

########## DONE MODIFYING ###############


########## BEGIN CODE ###############

echo -e ">>> INITIATING QC and trim with command:\n\t$0 $@"

# stash number of threads to used in sbatch job submission:
pthread=$2

# create main output directory
echo -e ">>> MAKING main output directory"
echo -e "\tmkdir $outputdir"
mkdir -p $outputdir
mkdir -p $tmpDir #also make tmp dir if not already in existance


####### META DATA #############

# this is a list of the sample names that will be used for looping
names=( $(cut -f1 --output-delimiter=' ' $sampleList) )


####### PIPELINE ##############
# update the user
echo -e ">>> PLAN: This script will process the following samples: "
echo -e "NAMES"

for (( counter=0; counter < ${#names[@]}; counter++ ))
do
    echo -e "${names[$counter]}"
done


#### STEP0: before starting, concatinate (or rename) files based on metadata file
echo -e "\n>>> cat: concatinating files"
outCat=$outputdir"00_catFiles/"
mkdir -p $outCat

for (( counter=0; counter < ${#names[@]}; counter++ ))
do
	#stash sample name
	samplename=${names[$counter]}

 	#check if there are more than 1 R1 or R2 files for each sample
	testNum1=( $(ls ${inputdir}/${samplename}*R1* | wc -l ))
 	testNum2=( $(ls ${inputdir}/${samplename}*R2* | wc -l ))
	if [ $testNum1 > 1 ] || [ $testNum2 > 1 ]
	then
    		#if more than one R1 or R2 files for a given sample, then cat and rename them
      		cmd01="cat ${inputdir}/${samplename}_*R1* > ${outCat}/${samplename}_R1.fastq"
    		cmd02="cat ${inputdir}/${samplename}_*R2* > ${outCat}/${samplename}_R2.fastq"
	else
		#if one R1 or R2 file for a given sample, then mv and rename them
  		cmd01="mv ${inputdir}/${samplename}_*R1* ${outCat}/${samplename}_R1.fastq"
		cmd02="mv ${inputdir}/${samplename}_*R2* ${outCat}/${samplename}_R2.fastq"
	fi
	  
    	echo -e "\t$ ${cmd01}"
    	echo -e "\t$ ${cmd02}"

done


if [ $runPREQC == TRUE ]
then
    #### STEP1: fastqc to determine quality & trim reads
    echo -e "\n>>> FASTQC: generating quality report"
    mkdir -p $outputdir"01_fastqc_pre"    

    # execute fastqc pre trim
    cmd1="fastqc -o ${outputdir}01_fastqc_pre -t $pthread ${outCat}*.fastq"

    echo -e "\t$ ${cmd1}"
    time eval $cmd1
else

    echo -e "\n>>> skipping pre FASTQC and moving to trim"

fi


#### STEP2: Trimgalore to trim reads +/- QC
echo -e "\n>>> trim_galore: Trimgalore reads to size"
mkdir -p $outputdir"02_trim_galore"

for (( counter=0; counter < ${#names[@]}; counter++ ))
do

    sample1=${names[$counter]}"_R1.fastq"
    sample2=${names[$counter]}"_R2.fastq"
    
    # execute trim_galore +/- fastqc
    if [ $runPOSTQC == TRUE ]
    then
    	cmd2="trim_galore --paired $outCat/$sample1 $outCat/$sample2 \
     	-o $outputdir"02_trim_galore" \
	--basename $samplename \
	--fastqc \
	-q 30 &" 
 
    else
    	cmd2="trim_galore --paired $outCat/$sample1 $outCat/$sample2 \
	-o $outputdir"02_trim_galore" \
	--basename $samplename \
	-q 30 &"
    fi
    
    echo -e "\t$ ${cmd2}"
    time eval $cmd2

   #limit number of background jobs to the number of threads used
   if [[ $(jobs -r -p | wc -l) -ge $pthread ]]
   then
   	wait -n
   fi

done

wait
# note: the ampersand at the end of cmd2 and the 'wait' command enable parallelization

echo -e "\n>>> Data are now trimmed and ready for alignment"

######## VERSIONS #############
echo -e "\n>>> VERSIONS:"
echo -e "\n>>> FASTQC VERSION:"
$fastqc --version
echo -e "\n>>> trim_galore VERSION:"
$trim_galore --version


echo -e ">>> END: WES-QC-trim complete."



