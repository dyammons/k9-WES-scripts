#!/usr/bin/env bash

########################################################################
#  Function: automate job script creation for multiple WES samples!  #
#                                                                      #
#  Useage: bash mkBatch-align.sh samples.txt                                           #
#                                                                      #
#                                                                      #
#  Requirments: update the user preferences below                      #
#                                                                      #
#  Created: August 7, 2023 - by DA                                             #
#  Updated: NA                                     #
########################################################################


##### set user preferences #####
samplesPerJob=10
numNode=1
numTasks=20
email="adh91@colostate.edu"
runTime="24:00:00"
partition="amilan"

##### make the files #####
#retrieve samples based on what is present in the sample list
split -l $samplesPerJob $1 WES_metadata_

batch=$(ls -lh | grep "WES_metadata_" | awk '{print $9}')
declare -a StringArray=($batch)

> jobList.txt
for val in "${StringArray[@]}"; do

	#create sbatch file
	> cute_bwa_align_$val.sbatch
	echo "#!/usr/bin/env bash" >> cute_bwa_align_$val.sbatch
	echo "" >> cute_bwa_align_$val.sbatch
	echo "#SBATCH --job-name=cnt_$val" >> cute_bwa_align_$val.sbatch
	echo "" >> cute_bwa_align_$val.sbatch
	echo "#SBATCH --nodes=$numNode" >> cute_bwa_align_$val.sbatch                       
	echo "#SBATCH --ntasks=$numTasks" >> cute_bwa_align_$val.sbatch
	echo "#SBATCH --time=$runTime" >> cute_bwa_align_$val.sbatch
	echo "" >> cute_bwa_align_$val.sbatch
	echo "#SBATCH --partition=$partition" >> cute_bwa_align_$val.sbatch
	echo "#SBATCH --qos=normal" >> cute_bwa_align_$val.sbatch
	echo "" >> cute_bwa_align_$val.sbatch
	echo "#SBATCH --mail-type=END" >> cute_bwa_align_$val.sbatch
	echo "#SBATCH --mail-user=$email" >> cute_bwa_align_$val.sbatch
	echo "#SBATCH --output=cellRngr_$val-%j.log" >> cute_bwa_align_$val.sbatch
	echo "" >> cute_bwa_align_$val.sbatch
	echo "" >> cute_bwa_align_$val.sbatch
	echo "##### Call bash script #####" >> cute_bwa_align_$val.sbatch
	echo "" >> cute_bwa_align_$val.sbatch
	echo "bash WES-align.sh $val" >> cute_bwa_align_$val.sbatch
	echo "" >> cute_bwa_align_$val.sbatch
	
	#create list of cmds    
	echo "sbatch cute_bwa_align_$val.sbatch" >> jobList.txt
done 
