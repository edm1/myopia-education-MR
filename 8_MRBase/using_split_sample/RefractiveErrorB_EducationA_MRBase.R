#!/usr/bin/env Rscript
#

library("TwoSampleMR")

setwd("/mnt/seconddrive/data/phd/ukbiobank/myopia_EA_project/8_MRBase/using_split_sample/")

# Args
in_exposure = "1_split_sample_weights/output/myopiaInstr_stats_splitB.tsv"
in_outcome =  "1_split_sample_weights/output/myopiaInstr_stats_splitA.tsv"
outname = "output/RefractiveErrorB_EducationA"

# Load education exposure data (Okbay SNPs)
exp_data = read_exposure_data(filename=in_exposure,
                              sep="\t",
                              snp_col="rsid",
                              beta_col="b_mse",
                              se_col="se_mse",
                              effect_allele_col="alleleB",
                              other_allele_col="alleleA",
                              eaf_col="eaf",
                              pval_col="p_mse",
                              clump=F)

exp_data$exposure = "aveMSE"


# Load refractive error outcome data (Pickrell SNPs with betas calculated in UKBB)
out_data = read_outcome_data(snps=exp_data$SNP,
                             filename=in_outcome,
                             sep="\t",
                             snp_col="rsid",
                             beta_col="b_ea",
                             se_col="se_ea",
                             effect_allele_col="alleleB",
                             other_allele_col="alleleA",
                             eaf_col="eaf",
                             pval_col="p_ea")
out_data$outcome = "eduyears"

# Harmonise data
dat = harmonise_data(exposure_dat=exp_data,
                     outcome_dat=out_data)

# Run MR
res = mr(dat)

# Make report
mr_report(dat, output_path=outname)
