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
