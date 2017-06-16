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

