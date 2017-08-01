#!/usr/bin/env Rscript
#

setwd("/mnt/seconddrive/data/phd/ukbiobank/myopia_EA_project/2_derive_variables/")

# Args
indata = "output/phenotypes_alleleScores_170608.Rdata"
insample = "../data/genetics/data.sample"
outsample  = "output/sample_170630.sample"

# Load data
load(indata)

# Load sample order
sampledf = read.table(insample, sep=" ", as.is=TRUE)
sample_names = sampledf[3:nrow(sampledf), 1]

# Keep only samples from sample file
data = data[sample_names, ]

# Make snptest sample
colnames(data)
header1 = c("ave_MSE_clean",
            "isMyopic",
            "isMyopicSevere",
            "isAmblyope_clean",
            "ave_logmar_clean_log",
            "best_logmar_clean_log",
            "worst_logmar_clean_log",
            "worst_logmar_clean",
            "rl_MSE_diff_clean",
            "eduyears_clean",
            "ave_logmar_clean_cc0",
            "best_logmar_clean_cc0",
            "worst_logmar_clean_cc0",
            "worst_covMSE",
            "tdi_log",
            "age",
            "sex_genetic",
            "genechip",
            paste0("PC", 1:10)
            )
header2 = c("P",
            "B",
            "B",
            "B",
            "P",
            "P",
            "P",
            "P",
            "P",
            "P",
            "B",
            "B",
            "B",
            "C",
            "C",
            "C",
            "D",
            "D",
            rep("C", 10)
            )
# Add sample name and missing columns
outdata = cbind(ID_1=sample_names, ID_2=sample_names, missing=0, data[, header1])
header2 = c(0, 0, 0, header2)
# Write headers and data
write(colnames(outdata), file=outsample, sep=" ", ncolumns=ncol(outdata))
write(header2, file=outsample, sep=" ", ncolumns=length(header2), append=T)
write.table(outdata, file=outsample, sep=" ", quote=F, col.names=F, row.names=F, append=T)
