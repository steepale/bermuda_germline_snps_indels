#!/bin/bash -login
### Job name
### Resources
#PBS -l nodes=1:ppn=1,walltime=00:72:00:00,mem=60gb
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

java -Xmx60g -cp $GATK -jar $GATK/GenomeAnalysisTK.jar \
-T GenotypeGVCFs \
-R ./data/Galgal5/galgal5.fa \
-o ./data/raw_snps_indels/germline_raw_snps_indels_genotyped.g.vcf.gz

qstat -f ${PBS_JOBID}
