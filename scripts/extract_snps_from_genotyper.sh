#!/bin/bash -login
### Job name
### Resources
#PBS -l nodes=1:ppn=1,walltime=00:03:59:00,mem=10gb
### Send email if the job encounters an error
#PBS –m a
### Output files to where you submitted your batch file
#PBS -e ./jobs/${PBS_JOBNAME}_${Var}_${PBS_JOBID}.err
#PBS -o ./jobs/${PBS_JOBNAME}_${Var}_${PBS_JOBID}.log
#PBS -j oe

# Load modules
module load GATK/3.5.0
module load tabix/0.2.6

# Change to working directory
cd ${PBS_O_WORKDIR}

# Set variables
sample_name=$(basename ${Var} "_raw_snps_indels_genotyped.g.vcf.gz")

# Extract the SNPs from the call set
java -Xmx10g -cp $GATK -jar $GATK/GenomeAnalysisTK.jar \
-T SelectVariants \
-R ./data/Galgal5/genome.fa \
-V ${Var} \
-selectType SNP \
-o ./data/raw_snps_indels/${sample_name}_raw_snps_extracted.g.vcf.gz

# Job report
qstat -f ${PBS_JOBID}
