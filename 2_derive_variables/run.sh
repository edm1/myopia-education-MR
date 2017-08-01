#!/bin/sh
#PBS -l nodes=1:ppn=1
#PBS -l walltime=00:06:00:00
#PBS -N n-AS-ukbbvis

# Change to work dir if set
if [ ! -z ${PBS_O_WORKDIR+x} ]; then
	cd $PBS_O_WORKDIR
fi

mkdir -p output

# Derive allele scores
#Rscript 1_derive_allele_scores.R

# Derive phenotypes
Rscript 2_derive_phenotypes.R

# Merge phenotypes and allele scores
Rscript 3_merge_phenotypes_allelescores.R

# Create a SNPTEST sample file
Rscript 4_make_snptest_sample.R
