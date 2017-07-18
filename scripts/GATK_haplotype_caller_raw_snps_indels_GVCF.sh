#!/bin/bash -login
### Job name
### Resources
#PBS -l nodes=1:ppn=1,walltime=06:00:00:00,mem=120gb
### Send email if the job encounters an error
#PBS –m a
### Output files to where you submitted your batch file
#PBS -e ./jobs/$PBS_JOBNAME_$PBS_JOBID.err
#PBS -o ./jobs/$PBS_JOBNAME_$PBS_JOBID.log
#PBS -j oe

# Change to working directory
cd ${PBS_O_WORKDIR}

# Load modules
module load GATK/3.5.0

sample_name=$(basename ${Var} "_bwa_rg_dedupped_realigned_bqsr.bam")

java -Xmx120g -cp $GATK -jar $GATK/GenomeAnalysisTK.jar \
-T HaplotypeCaller \
-R ./data/Galgal5/galgal5.fa \
-I ${Var} \
-ERC GVCF \
-o ./data/raw_snps_indels/${sample_name}_raw_snps_indels.g.vcf

# Collect stats on run
qstat -f ${PBS_JOBID}

