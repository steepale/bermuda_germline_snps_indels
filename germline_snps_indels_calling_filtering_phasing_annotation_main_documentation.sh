#===============================================================================
#
#         FILE: /mnt/scratch/steepale/birdman/bermuda/germline_snps_indels/germline_snps_indels_calling_filtering_phasing_annotation_main_documentation.sh
#
#        USAGE: For development and documentation purposes, scripts inside
#
#  DESCRIPTION:  This script serves as a step by step documentation script and development script for germline variant
#                calls (SNPs and Indels), filtering, phasing, and annotation 
# REQUIREMENTS:  ---
#        NOTES:  ---
#       AUTHOR:  Alec Steep, steepale@msu.edu
#  AFFILIATION:  Michigan State University (MSU), East Lansing, MI, United States
#				         USDA ARS Avian Disease and Oncology Lab (ADOL), East Lansing, MI, United States
#				         Technical University of Munich (TUM), Weihenstephan, Germany
#      CREATED:  2017.06.11
#     REVISION:  
#===============================================================================

# PROJECT DIRECTORY (MSU Cluster)
PROJ_DIR='/mnt/scratch/steepale/birdman/bermuda/germline_snps_indels' 
cd ${PROJ_DIR}

# Calling germline variants with GATK
# https://www.broadinstitute.org/gatk/guide/article?id=2803

# make appopriate directories
mkdir -p ./{data,scripts,analysis,jobs}

mkdir ./analysis/variant_stats
mkdir ./data/filter_stats
mkdir ./data/Galgal5
mkdir ./data/bams
mkdir ./data/raw_snps_indels
mkdir ./analysis/raw_snps_indels
mkdir ./data/hard_filtered_variants
mkdir ./analysis/hard_filtered_variants
mkdir ./analysis/hard_filtered_variants/R_plots
mkdir ./data/truth_files
mkdir ./analysis/liftover
mkdir ./data/truth_files/divide_ma
mkdir ./scripts/compare_snps
mkdir ./data/vsqr
mkdir ./scripts/picard
mkdir ./data/truth_files/divide_ma/final_vcf
mkdir ./scripts/snpEff
mkdir ./data/germline_snps
mkdir ./analysis/snpEff
mkdir ./analysis/vsqr
mkdir ./data/germline_snps/indv_samples
mkdir ./scripts/indv_samples
mkdir ./data/arrays
mkdir ./data/germline_snps/indv_samples/tumors
mkdir ./data/germline_snps/indv_samples/normals
mkdir ./data/germline_snps/indv_samples/parents
mkdir ./data/germline_snps/indv_samples/private_snps_divided
mkdir ./data/haplosaurus
mkdir ./data/transmission_phased
mkdir ./data/GenVisR
mkdir ./data/proteinseqs

### Dependencies
# Move all the BAM files to the data directory, don't forget to move them back!
mv /mnt/scratch/steepale/birdman/bermuda/fastq2bam/data/final_bam/* ./data/bams/
# Copy the galgal5 reference genome and annotation files
cp /mnt/scratch/steepale/birdman/bermuda/fastq2bam/data/Galgal5/* ./data/Galgal5/


# Call Germline raw variants (SNPs and Indels) on DNA Sequencing samples
for sample in `find ./data/bams -name "*_bwa_rg_dedupped_realigned_bqsr.bam"`
do
echo $sample
qsub -v Var=$sample -N "GATK_haplotype_caller_raw_snps_indels_GVCF_"$sample ./scripts/GATK_haplotype_caller_raw_snps_indels_GVCF.sh
done

# ./scripts/GATK_haplotype_caller_raw_snps_indels_GVCF.sh
################################
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

##################################

# Perform Joint Genotyping
# https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_variantutils_GenotypeGVCFs.php

# Make sure th files are in uncompressed format if need be
# find ./data/raw_snps_indels -name "*.g.vcf.gz" |
# xargs -i sh -c "gunzip {}"

# Compress the files with bgzip (this process will take very long and you may need to edit the script)

for sample in `find ./data/raw_snps_indels -name "*_raw_snps_indels.g.vcf"`
do
echo $sample
qsub -v Var=$sample -N "compress_and_index_vcfs_"$sample ./scripts/compress_and_index_vcf.sh
done

# ./scripts/compress_and_index_vcf.sh
####################################
#!/bin/bash -login
### Job name
### Resources
#PBS -l nodes=1:ppn=1,walltime=00:03:00:00,mem=10gb
### Send email if the job encounters an error
#PBS –m a
### Output files to where you submitted your batch file
#PBS -e ./jobs/${PBS_JOBNAME}_${Var}_${PBS_JOBID}.err
#PBS -o ./jobs/${PBS_JOBNAME}_${Var}_${PBS_JOBID}.log
#PBS -j oe

# Change to working directory
cd ${PBS_O_WORKDIR}

# Load modules
module load tabix/0.2.6

# Compress
bgzip ${Var}

# Index
tabix -p vcf ${Var}'.gz'

# Remove previously made indexes
rm ${Var}'.idx'

# Collect stats on run
qstat -f ${PBS_JOBID}
#####################################


# Write the joint genotyping script before it can be submitted
cat <<Header_input > ./scripts/joint_genotype_germline_gvcfs.sh
#!/bin/bash -login
### Job name
### Resources
#PBS -l nodes=1:ppn=1,walltime=00:72:00:00,mem=60gb
### Send email if the job encounters an error
#PBS –m a
### Output files to where you submitted your batch file
#PBS -e ./jobs/\${PBS_JOBNAME}_\${PBS_JOBID}.err
#PBS -o ./jobs/\${PBS_JOBNAME}_\${PBS_JOBID}.log
#PBS -j oe

# Load modules
module load GATK/3.5.0

# Change to working directory
cd \${PBS_O_WORKDIR}

java -Xmx60g -cp \$GATK -jar \$GATK/GenomeAnalysisTK.jar \\
-T GenotypeGVCFs \\
-R ./data/Galgal5/galgal5.fa \\
Header_input

find ./data/raw_snps_indels -name "*_raw_snps_indels.g.vcf.gz" | \
xargs -i printf '%s\n' '-V {} \' >> ./scripts/joint_genotype_germline_gvcfs.sh

cat <<Footer_input >> ./scripts/joint_genotype_germline_gvcfs.sh
-o ./data/raw_snps_indels/germline_raw_snps_indels_genotyped.g.vcf.gz

qstat -f \${PBS_JOBID}
Footer_input

# # Perform Joint Genotyping on multiple samples
qsub -N "compress_and_index_vcfs" ./scripts/joint_genotype_germline_gvcfs.sh

# Perform Joint genotyping on a single sample:
# Line 6
#find `pwd` -name "002683_Line-6_raw_snps_indels.g.vcf" |\
#xargs -i echo 'qsub ./scripts/joint_genotype_single_sample.sh -v Var='{} |sh

# ./scripts/joint_genotype_single_sample.sh
################################
# #!/bin/bash -login
# #PBS -l nodes=1:ppn=1,walltime=24:00:00,mem=60gb
# #PBS -j oe
# set -e
# set -u
# set -o pipefail
# 
# module load GATK/3.5.0
# 
# cd /mnt/research/ADOL/OutsideCollaborations/20160201_Cheng_Steep_Xu_Zhang/germline_snps_indels
# 
# sample_name=$(basename ${Var} "_raw_snps_indels.g.vcf")
# 
# java -Xmx60g -cp $GATK -jar $GATK/GenomeAnalysisTK.jar \
# -T GenotypeGVCFs \
# -R ./data/Galgal5/galgal5.fa \
# -V ${Var} \
# -o ./data/raw_snps_indels/${sample_name}_raw_snps_indels_genotyped.g.vcf \
# 2> ./analysis/raw_snps_indels/${sample_name}_raw_snps_indels_genotyped.log
# 
# if test -f ./data/raw_snps_indels/${sample_name}_raw_snps_indels_genotyped.g.vcf
# then
#   echo "Joint Genotyping Complete: ./data/raw_snps_indels/${sample_name}_raw_snps_indels_genotyped.g.vcf $(date +%F)$" >> \
#   ./analysis/raw_snps_indels/${sample_name}_raw_snps_indels_genotyped.log
# fi
# 
# qstat -f ${PBS_JOBID}
##################################

###########################################################################################################

# Apply hard filters to the call set 
# https://www.broadinstitute.org/gatk/guide/article?id=2806

# We apply hard filters very stringently
# The strignent calls will account for the confident Germline variants already present 

# https://www.broadinstitute.org/gatk/guide/article?id=2806
# 1. Extract the SNPs from the call set
# 2. Determine parameters for filtering SNPs
# 3. Apply the filter to the SNP call set
# 4. Extract the Indels from the call set
# 5. Determine parameters for filtering indels
# 6. Apply the filter to the Indel call set

# 1. Extract the SNPs from the call set
find ./data/raw_snps_indels -name "germline_raw_snps_indels_genotyped.g.vcf.gz" |\
xargs -i echo 'qsub -N "extract_snps_from_genotyper" ./scripts/extract_snps_from_genotyper.sh -v Var='{} |sh

# ./scripts/extract_snps_from_genotyper.sh
################################
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
-R ./data/Galgal5/galgal5.fa \
-V ${Var} \
-selectType SNP \
-o ./data/raw_snps_indels/${sample_name}_raw_snps_extracted.g.vcf.gz

# Job report
qstat -f ${PBS_JOBID}
#################################

# Count all the SNPs per sample to make sure sample calls worked
# Strategy is to use the individually called samples from their pre-genotyped files to use as samples to run through script
find ./dev -name "*_raw_snps_indels.g.vcf.gz" | \
xargs -i basename {} | \
sed 's/_raw_snps_indels.g.vcf.gz//' | \
sort | uniq | \
xargs -i echo 'qsub -N "test_count_variants" ./scripts/test_count_variants.sh -v Var='{} |sh

# ./scripts/test_count_variants.sh
##############################
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
##############################

# 2. Determine parameters for filtering SNPs
# Extract stats to visualize in R
for sample_zip in `find ./data/raw_snps_indels -name "*_raw_snps_extracted.g.vcf.gz"`
do
sample_unzip=`echo ${sample_zip} | sed 's/.gz//'`
bgzip -d -c ${sample_zip} > ${sample_unzip}
python ./scripts/extract_stats_for_filtering.py ${sample_unzip}
done

# ./scripts/extract_stats_for_filtering.py
###################################
import sys
import os
import re

# infile
infile = sys.argv[1]

# Variables
sample = infile.split('./data/raw_snps_indels/')[1].split('_raw_snps_extracted.g.vcf')[0]

# outfile
outfile = open('./data/filter_stats/'+sample+'_raw_snps_filter_stats.txt', 'w')

# vcf format
# CHROM POS ID  REF ALT QUAL  FILTER  INFO  FORMAT  6x7-F1
# AC=2;AF=1.00;AN=2;DP=3;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=40.00;QD=32.76;SOR=2.833

# Print outfile header
outfile.write('QD' + '\t' + 'FS' + '\t' + 'MQ' + '\t' + 'MQRankSum' + '\t' + 'ReadPosRankSum' + '\n')

# Iterate through lines of vcf and collect ststs
for line in open(infile):
  if not line.startswith('#'):
    line = line.rstrip()
    col = line.split('\t')
    INFO = col[7]
    # Extract MQ Stats
    if re.search('MQ=', INFO):
      MQ = INFO.split('MQ=')[1].split(';')[0]
    else:
      MQ = '0.00'
    # Extract QD stats
    if re.search('QD=', INFO):
      QD = INFO.split('QD=')[1].split(';')[0]
    else:
      QD = '0.00'
    # Extract FS stats
    if re.search('FS=', INFO):
      FS = INFO.split('FS=')[1].split(';')[0]
    else:
      FS = '0.00'
    # Extract MQRankSum stats
    if re.search('MQRankSum=', INFO):
      MQRankSum = INFO.split('MQRankSum=')[1].split(';')[0]
    else:
      MQRankSum = '0.00'
    # Extract ReadPosRankSum stats
    if re.search('ReadPosRankSum=', INFO):
      ReadPosRankSum = INFO.split('ReadPosRankSum=')[1].split(';')[0]
    else:
      ReadPosRankSum = '0.00'
    # Print to outfile
    outfile.write(QD + '\t' + FS + '\t' + MQ + '\t' + MQRankSum + '\t' + ReadPosRankSum + '\n')
outfile.close()

print('Fin')
######################################

# Generate plots of stats associated with hard filtering with R. Used to make decisions with hard filtering.
module load R/3.2.0

for sample_file in `find ./data/filter_stats -name "*_raw_snps_filter_stats.txt"`
do
sample=`echo ${sample_file} | xargs -i basename {} | sed 's/_raw_snps_filter_stats.txt//'`
Rscript --vanilla ./scripts/hard_filter_plots.R ${sample_file} ${sample}
done

# ./scripts/hard_filter_plots.R
##########################################
#!/usr/bin/env Rscript

# Capture command line args
args = commandArgs(trailingOnly=TRUE)

library("ggplot2")

# Read in data and variables
df = read.table(args[1], header=TRUE)
sample = args[2]

# Variables for outfiles
out_dir="./analysis/hard_filtered_variants/R_plots/" 
end_string="_density_snps.png"
QD_out = paste0(out_dir, sample, "_QD", end_string)
FS_out = paste0(out_dir, sample, "_FS", end_string)
MQ_out = paste0(out_dir, sample, "_MQ", end_string)
MQ_40_out = paste0(out_dir, sample, "_MQ_40", end_string)
MQ_60_out = paste0(out_dir, sample, "_MQ_60", end_string)
MQRankSum_out = paste0(out_dir, sample, "_MQRankSum", end_string)
ReadPosRankSum_out = paste0(out_dir, sample, "_ReadPosRankSum", end_string)

## Generate the statistics plots
# QD Plot
qplot(QD, data = df[,1, drop = FALSE], geom = "density") + xlab("QD") + ylab("Density")
ggsave(QD_out)

# FS Plot (x-axis is logged)
qplot(FS, data = df[,2, drop = FALSE], geom = "density", log="x") + xlab("FS (x-axis logged)") + ylab("Density")
ggsave(FS_out)
# Cutoff right after right peak

# Generate the MQ Plot
qplot(MQ, data = df[,3, drop = FALSE], geom = "density") + xlab("MQ") + ylab("Density")
ggsave(MQ_out)

# Zoom in on the MQ peak at 40
qplot(MQ, data = df[,3, drop = FALSE], geom = "density") + xlab("MQ") + ylab("Density") + xlim(39.0, 41.0)
ggsave(MQ_40_out)
# 
# Zoom in on the MQ peak at 60
qplot(MQ, data = df[,3, drop = FALSE], geom = "density") + xlab("MQ") + ylab("Density") + xlim(59.0, 61.0)
ggsave(MQ_60_out)
# Discard anything that is not 60

# Generate the MQRankSum Plot
qplot(MQRankSum, data = df[,4, drop = FALSE], geom = "density") + xlab("MQRankSum") + ylab("Density")
ggsave(MQRankSum_out)
# Remove anything less than -2
# 
# Generate the ReadPosRankSum Plot
qplot(ReadPosRankSum, data = df[,5, drop = FALSE], geom = "density") + xlab("ReadPosRankSum") + ylab("Density")
ggsave(ReadPosRankSum_out)

##########################################

# 3. Apply the filter to the SNP call set

# Websites with recommendations on how to read plots and apply filters:
# https://www.broadinstitute.org/gatk/guide/article?id=6925
# https://www.broadinstitute.org/gatk/guide/article?id=2806

# Apply hard filters
for sample_file in `find ./data/raw_snps_indels -name "*_raw_snps_extracted.g.vcf.gz"`
do
sample=`echo ${sample_file} | xargs -i basename {} | sed 's/_raw_snps_extracted.g.vcf.gz//'`
qsub  -v Var=${sample_file} -N "hard_filter_SNPs_"$sample ./scripts/hard_filter_SNPs.sh
done

# ./scripts/hard_filter_SNPs.sh
################################
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
#################################

# Seperate samples from filtered SNPs

# Seperate samples from gvcf file into respective vcf files
# Collective germline samples
for sample_file in `find ./data/hard_filtered_variants/ -name "*_hard_filtered_snps.g.vcf.gz"`
do
sample=`echo ${sample_file} | xargs -i basename {} | sed 's/_hard_filtered_snps.g.vcf.gz//'`
qsub  -v Var=${sample_file} -N "seperate_samples_"$sample ./scripts/seperate_samples.sh
done

# ./scripts/seperate_samples.sh
##############################
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
##############################

# ./scripts/determine_samples.py
##############################
import os
import sys

# infile
infile = sys.argv[1]

# outfile
outfile = open(sys.argv[2], 'w')


# Grab the samples from the hard filtered gvcf file
for line in open(infile):
  if line.startswith("#CHROM"):
    line = line.rstrip()
    col = line.split('\t')
    for n in range(9,len(col)):
      outfile.write(col[n] + '\n')
outfile.close()
##############################

# VEP annotation for gallus gallus will be performed on the TUM cluster. If need be, we can always configure the MSU HPCC for later runs.

# Transfer files over (request from TUM cluster)
rsync -avp \
steepale@rsync.hpcc.msu.edu:/mnt/scratch/steepale/birdman/bermuda/germline_snps_indels/data/hard_filtered_variants/*_hard_filtered_snps.vcf* \
/scratch/steep/bermuda/germline_snps_indels/data/hard_filtered_variants/

# In TUM cluster
cd /scratch/steep/bermuda/germline_snps_indels/

# Run VEP annotation for germline snvs
for sample_file in `find ./data/hard_filtered_variants/ -name "*_hard_filtered_snps.vcf.gz"`
do
qsub -N "germline_snvs_vep" -v Var=${sample_file} ./scripts/germline_snvs_vep.sh
done

# ./scripts/germline_snvs_vep.sh
##############################
#!/bin/bash -login
### Job name
### Resources
#PBS -l nodes=1:ppn=1,walltime=04:00:00:00,mem=10gb
### Send email if the job encounters an error
#PBS –m a
### Output files to where you submitted your batch file
#PBS -e ./jobs/${PBS_JOBNAME}_${PBS_JOBID}.err
#PBS -o ./jobs/${PBS_JOBNAME}_${PBS_JOBID}.log
#PBS -j oe

# Change to working directory
#cd ${PBS_O_WORKDIR}
cd /scratch/steep/bermuda/germline_snps_indels

# Variables
#Var='./data/hard_filtered_variants/P4806_134_hard_filtered_snps.vcf.gz'
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
#######################################








