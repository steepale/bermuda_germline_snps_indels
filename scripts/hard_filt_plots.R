#!/usr/bin/env Rscript

#setwd()

# Capture command line args
args = commandArgs(trailingOnly=TRUE)

library("ggplot2")

# In R, Generate the statistics plots
df = read.table(args[1], header=TRUE)
sample = args[2]

# In R, Generate the statistics plots

out_dir="./analysis/hard_filtered_variants/R_plots/" 
end_string="_density_snps.png"

cat("Current working dir: ")
