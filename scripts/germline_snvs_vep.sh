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

# Change to working directory
cd ${PBS_O_WORKDIR}

# Variables
sample=`echo ${Var} | xargs -i basename {} | sed 's/_hard_filtered_snps.vcf.gz//'`

# Annotate with VEP
perl ${HOME}/Apps/ensembl-tools-release-87/scripts/variant_effect_predictor/variant_effect_predictor.pl \
-i ${Var} \
-o ./data/hard_filtered_variants/${sample}_hard_filtered_snps_vep.vcf.gz \
--vcf \
--cache \
--species gallus_gallus \
--force_overwrite \
--protein \
--hgvs \
--domains \
--ccds \
--uniprot \
--tsl \
--appris \
--sift b

# Job report
qstat -f ${PBS_JOBID}
