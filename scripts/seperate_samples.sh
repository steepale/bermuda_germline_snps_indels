#!/bin/bash -login
### Job name
### Resources
#PBS -l nodes=1:ppn=1,walltime=00:03:59:00,mem=10gb
### Send email if the job encounters an error
#PBS –m a
### Output files to where you submitted your batch file
#PBS -e ./jobs/${PBS_JOBNAME}_${PBS_JOBID}.err
#PBS -o ./jobs/${PBS_JOBNAME}_${PBS_JOBID}.log
#PBS -j oe

# Load modules
module load GATK/3.5.0

# Change to working directory
cd ${PBS_O_WORKDIR}

# Set variable
unzipped_file=`echo ${Var} | sed 's/.gz//'`
cohort=$(basename ${Var} "_hard_filtered_snps.g.vcf.gz")

# Determine samples (unzip file for python script)
gzip -d -c ${Var} > \
${unzipped_file}

python ./scripts/determine_samples.py \
${unzipped_file} \
./data/hard_filtered_variants/${cohort}_hard_filtered_samples.txt

# Iterate through sample name files to apply GATK to each sample
while read sample_name
do
  # Select individual samples and filter out filtered calls and all non-variants
  java -Xmx8g -cp $GATK -jar $GATK/GenomeAnalysisTK.jar \
  -T SelectVariants \
  -R ./data/Galgal5/galgal5.fa \
  -V ${Var} \
  -o ./data/hard_filtered_variants/${sample_name}_hard_filtered_snps.vcf.gz \
  -sn ${sample_name} \
  -env \
  -ef

  # Index the seperate sample vcf files
  tabix -f -p vcf ./data/hard_filtered_variants/${sample_name}_hard_filtered_snps.vcf.gz

done<./data/hard_filtered_variants/${cohort}_hard_filtered_samples.txt
# env: Dont include non-variant sites
# ef: Don't include filtered sites

# Job report
qstat -f ${PBS_JOBID}
