################# maybe split here ########################

#confirm input for each cmd should be sequential?

# GATK: begin mutant calling workflow
echo -e "\n>>> GATK: run MarkDuplicates:"
gatkOut=$outputdir"04_gatk/"
mkdir -p $gatkOut

for seqname in ${names[@]}
do
#run 20 threads for about 3.5 hrs per sample
#queryname grouped input unmodified from the aligne is faster
#https://hpc.nih.gov/training/gatk_tutorial/markdup.html#benchmarks-of-markduplicatesspark
  gatk MarkDuplicates -I ${samout}${seqname}_sort.bam \
  -O ${gatkOut}${seqname}_markdup.bam \
  --TMP_DIR ${tmpDir} \
  --METRICS_FILE ${gatkOut}${name}_markdups.txt
done

#find estimate run time
gatk AddOrReplaceReadGroups \
-I ${samout}${seqname}_markdup.bam  \
-O ${gatkOut}${seqname}_RG.bam \
-SORT_ORDER coordinate \
-RGID ${seqname} \
-RGLB ${seqname} \
-RGPU unit1 \
-RGPL ILLUMINA \
-RGSM ${name} \
--TMP_DIR ${tmpDir} \
-CREATE_INDEX True

#run with only 2 thread and no mem call for 24hrs
#https://hpc.nih.gov/training/gatk_tutorial/bqsr.html
#fix known-sites
gatk BaseRecalibrator \
-I ${samout}${seqname}_RG.bam \
-R ${genomeFA} \
--known-sites /scratch/alpine/adh91@colostate.edu/WES/03_GATK/rename.broad_umass_canid_variants.1.2.vcf.gz \
--known-sites /projects/adh91@colostate.edu/references/Ensembl.CanFam3.1/canis_lupus_familiaris.vcf.gz \
-O ${gatkOut}${seqname}.recal_data.table \
--tmp-dir ${tmpDir}

#run with only 2 thread and no mem call for 24hrs
#https://hpc.nih.gov/training/gatk_tutorial/bqsr.html
gatk ApplyBQSR \
-I ${samout}${seqname}_RG.bam \
-R ${genomeFA} \
--bqsr-recal-file ${gatkOut}${seqname}.recal_data.table \
--tmp-dir ${tmpDir} \
-O ${gatkOut}${seqname}_bqsr.bam

#should you be running this on the normal samples?#
#HaplotypeCaller might be more appropriate for the norms??

gatk Mutect2 \
-R ${genomeFA} \
--germline-resource /scratch/alpine/adh91@colostate.edu/WES/03_GATK/rename.broad_umass_canid_variants.1.2.vcf.gz \
-max-mnp-distance 0 \
-I ${gatkOut}${seqname}_bqsr.bam \
-O ${gatkOut}${seqname}.vcf.gz




