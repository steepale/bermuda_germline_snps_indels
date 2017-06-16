# bermuda_germline_snps_indels
This pipeline aims to call hard-filtered germline snps (indels if necessary) with GATK from whole genome sequencing of Bermuda chickens.
Pipeline is an extension of the bermuda_fastq2bam repository.

See germline_snps_indels_calling_filtering_phasing_annotation_main_documentation.sh for step-by-step documnetation.

This pipeline is by no means a fool-proof method for calling germline SNPs, results will contain more false positives than a machine learning approach.

Outline of pipeline:
1. Bam files were called for germline snps and indels with GATK's haplotype caller
2. Haplotype calls were then joint genotyped by cohort (only one all inclusive cohort in this case)
3. SNPs were extracted from the variant pool
4. Useful stats for filtering were extarcted and visualized with R plots in order to determine hard filtering parameters
5. SNPs were hard-filtered based on GATK's recommendations and R plots
6. Samples were then seperated and annotated with Ensembl's Variant Effect Predictor (VEP)

All scripts have been written for reproducibility and portability.
