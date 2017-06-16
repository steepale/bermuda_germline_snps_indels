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
