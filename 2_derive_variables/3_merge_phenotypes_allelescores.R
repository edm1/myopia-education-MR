#!/usr/bin/env Rscript
#

library(tidyverse)

setwd("/mnt/seconddrive/data/phd/ukbiobank/myopia_EA_project/2_derive_variables/")

# Args
inpheno = "output/phenotypes_170728.Rdata"
inalscr = "output/alleleScores_170728.Rdata"
outpref  = "output/phenotypes_alleleScores_170728"

# Load
load(inpheno)
load(inalscr)

# Remove X prefix from allelescore names
rownames(alleleScr) = gsub("X", "", rownames(alleleScr))

# Merge
mergedf = cbind(data, alleleScr[match(rownames(data), rownames(alleleScr)), ])

# Write Rdata
data = mergedf
save(data, file=paste0(outpref, ".Rdata"))

# Write TSV for Neil
datasubset = data %>% filter(coreSampleGeno & !genoExcl) %>% select(eid,
                             ave_MSE_clean,
                             eduyears_clean,
                             ea_AS_dosage,
                             myopia_AS_dosage,
                             tdi_log,
                             age,
                             sex_genetic,
                             northing,
                             easting,
                             breastfed,
                             birthweight,
                             PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10)
write.table(datasubset, file=paste0(outpref, ".forNeil.tsv"), sep="\t", row.names=F, quote=F)

# Write full TSV
write.table(data, file=paste0(outpref, ".all.tsv"), sep="\t", row.names=F, quote=F)
