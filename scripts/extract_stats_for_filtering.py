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
