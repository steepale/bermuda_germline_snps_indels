# bermuda_germline_snps_indels
This pipeline aims to call hard-filtered germline snps (indels if necessary) with GATK from whole genome sequencing of Bermuda chickens.
Pipeline is an extension of the bermuda_fastq2bam repository.

See germline_snps_indels_calling_filtering_phasing_annotation_main_documentation.sh for step-by-step documnetation.

This pipeline is by no means a fool-proof method for calling germline SNPs, results will contain more false positives than a machine learning approach.

Outline of pipeline:
-Bam files were called for germline snps and indels with GATK's haplotype caller
-Haplotype calls were then joint genotyped by cohort (only one all inclusive cohort in this case)
-SNPs were extracted from the variant pool
-Useful stats for filtering were extarcted and visualized with R plots in order to determine hard filtering parameters
-SNPs were hard-filtered based on GATK's recommendations and R plots
-Samples were then seperated and annotated with Ensembl's Variant Effect Predictor (VEP)

All scripts have been written for reproducibility and portability.
