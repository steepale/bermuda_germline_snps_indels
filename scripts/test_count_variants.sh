#!/bin/bash -login
### Job name
### Resources
#PBS -l nodes=1:ppn=1,walltime=00:00:30:00,mem=10gb
### Send email if the job encounters an error
#PBS –m a
### Output files to where you submitted your batch file
#PBS -e ./jobs/${PBS_JOBNAME}_${Var}_${PBS_JOBID}.err
#PBS -o ./jobs/${PBS_JOBNAME}_${Var}_${PBS_JOBID}.log
#PBS -j oe

# Load modules
module load GATK/3.5.0
module load bcftools/1.2

# Change to working directory
cd ${PBS_O_WORKDIR}

# Select individual samples and filter out filtered calls and all non-variants
java -Xmx10g -cp $GATK -jar $GATK/GenomeAnalysisTK.jar \
-T SelectVariants \
-R ./data/Galgal5/galgal5.fa \
-V ./data/raw_snps_indels/germline_raw_snps_extracted.g.vcf.gz \
-o ./data/raw_snps_indels/${Var}_raw_snps_extracted.vcf.gz \
-sn ${Var} \
-env \
-ef

# Perform variant counts with bcftools
bcftools stats ./data/raw_snps_indels/${Var}_raw_snps_extracted.vcf.gz \
> ./analysis/variant_stats/snp_counts_${Var}.txt

# Job report
qstat -f ${PBS_JOBID}
