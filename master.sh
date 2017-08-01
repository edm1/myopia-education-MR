#!/usr/bin/env bash
#

#
# Extract genotypes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

cd 1_extract_instruments
mkdir -p output
qsub run_education.sh
qsub run_myopia.sh
cd ..

#
# Derive variables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

cd 2_derive_variables
mkdir -p output
Rscript 1_derive_allele_scores.R
Rscript 2_derive_phenotypes.R
Rscript 3_merge_phenotypes_allelescores.R
Rscript 4_make_snptest_sample.R
cd ..

#
# Generate descriptives ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

cd 3_descriptives
mkdir -p output
Rscript 1_describe_data.R
cd ..

#
# Run MR ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

cd 4_MR
mkdir -p output
Rscript 1_MR.R
cd ..

#
# Produce variant plots ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

cd 5_variant_plots
mkdir -p output
Rscript eduAttainInstruments.R
Rscript myopiaInstruments.R
cd ..

#
# Run bias analysis ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

cd 6_bias-analysis
Stata -s do ds_1_bias_plots.do
cd ..

#
# Run MR-base ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

# Calculate split sample weights
cd 8_MRBase/using_split_sample/1_split_sample_weights
mkdir -p output
Rscript eduAttainInstruments_sampleA.R
Rscript eduAttainInstruments_sampleB.R
Rscript myopiaInstruments_sampleA.R
Rscript myopiaInstruments_sampleB.R
cd ..

# Run MR-base (check which samples used)
mkdir -p output
Rscript EducationA_RefractiveErrorB_MRBase.R
Rscript RefractiveErrorA_EducationB_MRBase.R
cd ../..
