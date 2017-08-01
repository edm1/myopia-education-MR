#!/usr/bin/env Rscript
#

library("rio")

setwd("/mnt/seconddrive/data/phd/ukbiobank/myopia_EA_project/2_derive_variables/")

#
# Args
#

# Extracted myopia dosages
inMyopiaInstr = "../1_extract_instruments/output/Myopia_UKBB_Pickrell.dosage"
# Extracted EA dosages
inEaInstr = "../1_extract_instruments/output/EA_UKBB_Okbay.dosage"
# Extracted RE dosages
inReInstr = "../1_extract_instruments/output/RE_UKBB_Cream.dosage"
# Extracted VA dosages
inVaInstr = "../1_extract_instruments/output/VA_UKBB_ALSPAC.dosage"
# Out prefix
outprefix = "output/alleleScores_170728"
# Column at which sample dosages start
geno_start = 8
# Set threshold for hard calling best guess geno
hard_call_threshold = 0.1

# Test data
#inMyopiaInstr = "test_data/Myopia_test15.dosage"
#inMyopiaInstr = "test_data/Myopia_test1000.dosage"
#inEaInstr = "test_data/EA_test1000.dosage"
#setwd("/mnt/seconddrive/data/phd/ukbiobank/myopia_EA_project/2_derive_variables/")

#
# Do Myopia ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

# Load
df = rio::import(inMyopiaInstr, format="\t", header=T)
#df = read.csv(inMyopiaInstr, sep="\t", as.is=T)

#
# Create allele scores using imputed genotypes
#
colnames(df)
# Extract dosage columns
dosages = df[, geno_start:ncol(df)]
# Multiply by weights to give weighted dosages
weighted_dosages = sweep(dosages, MARGIN=1, df$weightB, "*")
# Sum rows to give per sample allele score
myopia_AS_dosage = colSums(weighted_dosages)
# Print to log
cat("Myopia dosage allele score made up of ", nrow(weighted_dosages), " variants\n",
    file=paste0(outprefix, ".log"), append=F)


#
# Create allele scores using bestguess
#

# # Define function to hard call
# hard_call = function(x, thresh) {
#   if(x <= (0 + thresh)) { 0 }
#   else if((1 - thresh) <= x & x <= (1 + thresh)) { 1 }
#   else if(x >= (2 - thresh)) { 2 }
#   else { NA }
# }
# # Hard call the dosages
# bestguess = apply(dosages, c(1,2), hard_call, thresh=hard_call_threshold)
# # Multiply by weights
# weighted_bestguess = sweep(bestguess, MARGIN=1, df$weightB, "*")
# # Remove variants that are missing in any sample
# weighted_bestguess = weighted_bestguess[complete.cases(weighted_bestguess), ]
# # Sum rows to give per sample allele score
# myopia_AS_bestguess = colSums(weighted_bestguess)
# # Print to log
# cat("Myopia bestguess allele score made up of ", nrow(weighted_bestguess), " variants\n",
#     file=paste0(outprefix, ".log"), append=T)

#
# Do EA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

# Load
df = rio::import(inEaInstr, format="\t", header=T)

#
# Create allele scores using imputed genotypes
#

# Extract dosage columns
dosages = df[, geno_start:ncol(df)]
# Multiply by weights to give weighted dosages
weighted_dosages = sweep(dosages, MARGIN=1, df$weightB, "*")
# Sum rows to give per sample allele score
ea_AS_dosage = colSums(weighted_dosages)
# Print to log
cat("EA bestguess allele score made up of ", nrow(weighted_dosages), " variants\n",
    file=paste0(outprefix, ".log"), append=T)

# #
# # Create allele scores using bestguess
# #
#
# # Hard call the dosages
# bestguess = apply(dosages, c(1,2), hard_call, thresh=hard_call_threshold)
# # Multiply by weights
# weighted_bestguess = sweep(bestguess, MARGIN=1, df$weightB, "*")
# # Remove variants that are missing in any sample
# weighted_bestguess = weighted_bestguess[complete.cases(weighted_bestguess), ]
# nrow(weighted_bestguess)
# # Sum rows to give per sample allele score
# ea_AS_bestguess = colSums(weighted_bestguess)
# # Print to log
# cat("EA bestguess allele score made up of ", nrow(weighted_bestguess), " variants\n",
#     file=paste0(outprefix, ".log"), append=T)

#
# Do RE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

# Load
df = rio::import(inReInstr, format="\t", header=T)

#
# Create allele scores using imputed genotypes
#

# Extract dosage columns
dosages = df[, geno_start:ncol(df)]
# Multiply by weights to give weighted dosages
weighted_dosages = sweep(dosages, MARGIN=1, df$weightB, "*")
# Sum rows to give per sample allele score
re_AS_dosage = colSums(weighted_dosages)
# Print to log
cat("RE bestguess allele score made up of ", nrow(weighted_dosages), " variants\n",
    file=paste0(outprefix, ".log"), append=T)

# #
# # Create allele scores using bestguess
# #
#
# # Hard call the dosages
# bestguess = apply(dosages, c(1,2), hard_call, thresh=hard_call_threshold)
# # Multiply by weights
# weighted_bestguess = sweep(bestguess, MARGIN=1, df$weightB, "*")
# # Remove variants that are missing in any sample
# weighted_bestguess = weighted_bestguess[complete.cases(weighted_bestguess), ]
# nrow(weighted_bestguess)
# # Sum rows to give per sample allele score
# re_AS_bestguess = colSums(weighted_bestguess)
# # Print to log
# cat("RE bestguess allele score made up of ", nrow(weighted_bestguess), " variants\n",
#     file=paste0(outprefix, ".log"), append=T)

#
# Do VA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

# Load
df = rio::import(inVaInstr, format="\t", header=T)

#
# Create allele scores using imputed genotypes
#

# Extract dosage columns
dosages = df[, geno_start:ncol(df)]
# Multiply by weights to give weighted dosages
weighted_dosages = sweep(dosages, MARGIN=1, df$weightB, "*")
# Sum rows to give per sample allele score
va_AS_dosage = colSums(weighted_dosages)
# Print to log
cat("VA bestguess allele score made up of ", nrow(weighted_dosages), " variants\n",
    file=paste0(outprefix, ".log"), append=T)

#
# Create allele scores using bestguess
#

# # Hard call the dosages
# bestguess = apply(dosages, c(1,2), hard_call, thresh=hard_call_threshold)
# # Multiply by weights
# weighted_bestguess = sweep(bestguess, MARGIN=1, df$weightB, "*")
# # Remove variants that are missing in any sample
# weighted_bestguess = weighted_bestguess[complete.cases(weighted_bestguess), ]
# nrow(weighted_bestguess)
# # Sum rows to give per sample allele score
# va_AS_bestguess = colSums(weighted_bestguess)
# # Print to log
# cat("VA bestguess allele score made up of ", nrow(weighted_bestguess), " variants\n",
#     file=paste0(outprefix, ".log"), append=T)

#
# Make output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

# Make df of results
alleleScr = data.frame(myopia_AS_dosage,
                       ea_AS_dosage,
                       re_AS_dosage,
                       va_AS_dosage)

# Save
save(alleleScr, file=paste0(outprefix, ".Rdata"))
