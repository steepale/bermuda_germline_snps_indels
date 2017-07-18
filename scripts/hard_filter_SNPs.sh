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

# Variables
sample_name=$(basename ${Var} "_raw_snps_extracted.g.vcf.gz")

# Perform hard filtering
java -Xmx10g -cp $GATK -jar $GATK/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R ./data/Galgal5/galgal5.fa \
-V ${Var} \
--filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -2.5" \
--filterName "SNP_HARD_FILTER" \
-o ./data/hard_filtered_variants/${sample_name}_hard_filtered_snps.g.vcf.gz

# Job report
qstat -f ${PBS_JOBID}
