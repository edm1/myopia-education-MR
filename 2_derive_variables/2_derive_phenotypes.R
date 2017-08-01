#!/usr/bin/env Rscript
#

library("tidyverse")
library("psych")

setwd("/mnt/seconddrive/data/phd/ukbiobank/myopia_EA_project/2_derive_variables/")

#
# Args
#

in8434 = "../data/phenotypes/data.8434.Rdata"
in4263 = "../data/phenotypes/data.4263.Rdata"
in4688 = "../data/phenotypes/data.4688.Rdata"
#in8434 = "../data/phenotypes/data.8434.subset50k.Rdata"
#in4263 = "../data/phenotypes/data.4263.subset50k.Rdata"
#in4688 = "../data/phenotypes/data.4688.subset50k.Rdata"
outdata = "output/phenotypes_170728.Rdata"
outlog = "output/phenotypes_170728.log"

#
# Load data ####################################################################
#

# Load batch 8434
load(in8434)
df8434 = bd
rownames(df8434) = df8434$f.eid
rm(bd)

# Load batch 4263
load(in4263)
df4263 = bd
rownames(df4263) = df4263$f.eid
rm(bd)

# Load batch 4688
load(in4688)
df4688 = bd
rownames(df4688) = df4688$f.eid
rm(bd)

# Start log
cat("Derive phenotype log file\n\n", file=outlog, sep=" ",append=F)

#
# Define functions #############################################################
#

# Return columns using wildcards
grep_cols = function(df, query) {
  return(df[, grep(glob2rx(query), names(df), value=TRUE), drop=F])
}

# Convert 1 column df to a named vector. There must be an easier way to do this
to_named_vector = function(df) {
  # Asset 1 column
  stopifnot(ncol(df) == 1)
  # convert
  v = df[, 1]
  names(v) = rownames(df)
  return(v)
}

# Translate, log transform and centre a vector
log_transform = function(x) {
  # Translate, so that all are positive
  trans_x = x - min(x, na.rm=T) + 1
  # Log
  log_x = log(trans_x)
  # Centre
  centre_x = scale(log_x, scale=F)
  return(to_named_vector(centre_x))
}

# Merges multiple columns, taking the first if not NA
merge_instances = function(df) {
  while (ncol(df) > 1) {
    # Replace col 1 with col 2 if it doesn't exist
    to_replace = is.na(df[, 1]) & !is.na(df[, 2])
    df[to_replace, 1] = df[to_replace, 2]
    # Remove col 2
    df = df[, -c(2), drop=F]
  }
  df
}

# Processes UKBB multiple categorical variable. Df has a column of bools for
# each option
process_multicategory = function(df) {
  # Get all factor names
  all_facs = c()
  for (i in 1:ncol(df)) {
    all_facs = unique(c(all_facs, levels(df[, i])))
  }
  # For each factor, test for presence
  res = list()
  for (fac in all_facs) {
    res[[fac]] = apply(df == fac, 1, function(row) { any(row, na.rm=T) } )
  }
  # Make output dataframe
  resdf = data.frame(do.call("cbind", res))
  # Set missing as NA
  isMissing = apply(is.na(df), 1, function(row) { all(row) } )
  resdf[isMissing, ] = NA

  return(resdf)
}


#
# Derive educational attainment ################################################
#

# Note: People who say they have college or uni degree where not asked what age
#       they left school. We therefore set it to 21

#
# Derive using Neils method
#

# Get eduyears
eduyears = to_named_vector(merge_instances(grep_cols(df4263, "f.845.*.0")))

# Set any negative values to NA:
#-3 (Prefer not to answer), -2 (Never went to school), -1 (Do not know)
eduyears[eduyears == -3 | eduyears == -1] = NA
eduyears[eduyears == -2] = 15


# Get qualifications
qualif = grep_cols(df4263, "f.6138.*.*")
# Find if any have university or college
haveDegree = apply(qualif == "College or University degree", 1, function(row) { any(row, na.rm=T) } )

# Set anyone who has a degree and eduyears is NA to 21
eduyears[haveDegree & is.na(eduyears)] = 21

# Observe cumulative proportions
cumsum(table(eduyears)) / sum(!is.na(eduyears))

# Set anyone <15 to 15
eduyears[eduyears < 15 & !is.na(eduyears)] = 15
# Set anyone >21 to 21
eduyears[eduyears > 21 & !is.na(eduyears)] = 21

# Make same variable without college/university (for sensitivity)
eduyearsNo21 = eduyears
eduyearsNo21[haveDegree] = NA

# Make a binary variable of more/less that 15 years
eduyearsMore16 = NA
eduyearsMore16[eduyears > 16] = 1
eduyearsMore16[eduyears <= 16] = 0
eduyearsMore16 = factor(eduyearsMore16)
names(eduyearsMore16) = names(eduyears)

#
# Derive using Okbay method
#

eduyearsOkbay = NA
# None of the above (7)
tochange = apply(qualif == "None of the above", 1, function(row) { any(row, na.rm=T) } )
eduyearsOkbay[tochange] = 7
# CSEs or equivalent (10)
tochange = apply(qualif == "CSEs or equivalent", 1, function(row) { any(row, na.rm=T) } )
eduyearsOkbay[tochange] = 10
# O levels/GCSEs or equivalent (10)
tochange = apply(qualif == "O levels/GCSEs or equivalent", 1, function(row) { any(row, na.rm=T) } )
eduyearsOkbay[tochange] = 10
# A levels/AS levels or equivalent (13)
tochange = apply(qualif == "A levels/AS levels or equivalent", 1, function(row) { any(row, na.rm=T) } )
eduyearsOkbay[tochange] = 13
# Other professional qualifications eg: nursing, teaching (15)
tochange = apply(qualif == "Other professional qualifications eg: nursing, teaching", 1, function(row) { any(row, na.rm=T) } )
eduyearsOkbay[tochange] = 15
# NVQ or HND or HNC or equivalent (19)
tochange = apply(qualif == "NVQ or HND or HNC or equivalent", 1, function(row) { any(row, na.rm=T) } )
eduyearsOkbay[tochange] = 19
# College or University degree (20)
tochange = apply(qualif == "College or University degree", 1, function(row) { any(row, na.rm=T) } )
eduyearsOkbay[tochange] = 20
names(eduyearsOkbay) = rownames(qualif)

#
# Identify individuals who weren't born in england
#

# Country of birth for UK residents
birthCountryUK = grep_cols(df4263, "f.1647.*.*")
bornEngland = apply(birthCountryUK == "England", 1, function(row) { any(row, na.rm=T) } )
bornEngWales = apply(birthCountryUK == "England" | birthCountryUK == "Wales", 1, function(row) { any(row, na.rm=T) } )

# Country of birth for other - Will be NA if from UK
birthCountryOther = to_named_vector(grep_cols(df4263, "f.20115.*.*"))
isUK = is.na(birthCountryOther)

# Keep those true for both
isEngWales = bornEngWales & isUK

# Remove non england from eduyears
eduyears_clean = eduyears
eduyears_clean[!isEngWales] = NA
# And for Okbay version
eduyearsOkbay_clean = eduyearsOkbay
eduyearsOkbay_clean[!isEngWales] = NA
# And for version without college graduates
eduyearsNo21_clean = eduyearsNo21
eduyearsNo21_clean[!isEngWales] = NA

# Calculate cleaned cum proportions
cumsum(table(eduyears_clean)) / sum(!is.na(eduyears_clean))

#
# Descriptives
#

cat("\nEduyears descriptives:\n\n", file=outlog, sep=" ",append=TRUE)

# Cleaned
cat("eduyears:\n", file=outlog, sep=" ",append=TRUE)
write.table(describe(eduyears), file=outlog, append=TRUE, quote=F)
cat("eduyears_clean:\n", file=outlog, sep=" ",append=TRUE)
write.table(describe(eduyears_clean), file=outlog, append=TRUE, quote=F)

#
# Derive refractive error ######################################################
#

# Note: RE was measured at 2 instances, with 10 measurements per instance.
#       Unreliable autorefractor measures are marked with an E or e

#
# Do left RE for instance 0
#

# Get errors
l_refUnrel_0 = grep_cols(df8434, "f.5090.0.*")
# Create logical matrix
is_err = apply(l_refUnrel_0, 2, function(v) { grepl("e|E", v) })

# Get spherical power
l_sph_0 = grep_cols(df8434, "f.5085.0.*")
# Set errors to NA
l_sph_0 = as.matrix(l_sph_0)
l_sph_0[is_err] = NA
# Cal mean, sd, and n for spherical power
l_sph_mean_0 = apply(l_sph_0, 1, function(row) { mean(row, na.rm=T) })
l_sph_sd_0 = apply(l_sph_0, 1, function(row) { sd(row, na.rm=T) })
l_sph_n_0 = apply(l_sph_0, 1, function(row) { sum(!is.na(row)) })

# Get cylindrial power
l_cyl_0 = grep_cols(df8434, "f.5086.0.*")
# Set errors to NA
l_cyl_0 = as.matrix(l_cyl_0)
l_cyl_0[is_err] = NA
# Cal mean, sd, and n for spherical power
l_cyl_mean_0 = apply(l_cyl_0, 1, function(row) { mean(row, na.rm=T) })
l_cyl_sd_0 = apply(l_cyl_0, 1, function(row) { sd(row, na.rm=T) })
l_cyl_n_0 = apply(l_cyl_0, 1, function(row) { sum(!is.na(row)) })

# Calc mean spherical equivalence
l_MSE_0 = l_sph_mean_0 + (0.5 * l_cyl_mean_0)

#
# Do right RE for instance 0
#

# Get errors
r_refUnrel_0 = grep_cols(df8434, "f.5091.0.*")
# Create logical matrix
is_err = apply(r_refUnrel_0, 2, function(v) { grepl("e|E", v) })

# Get spherical power
r_sph_0 = grep_cols(df8434, "f.5084.0.*")
# Set errors to NA
r_sph_0 = as.matrix(r_sph_0)
r_sph_0[is_err] = NA
# Cal mean, sd, and n for spherical power
r_sph_mean_0 = apply(r_sph_0, 1, function(row) { mean(row, na.rm=T) })
r_sph_sd_0 = apply(r_sph_0, 1, function(row) { sd(row, na.rm=T) })
r_sph_n_0 = apply(r_sph_0, 1, function(row) { sum(!is.na(row)) })

# Get cylindrial power
r_cyl_0 = grep_cols(df8434, "f.5087.0.*")
# Set errors to NA
r_cyl_0 = as.matrix(r_cyl_0)
r_cyl_0[is_err] = NA
# Cal mean, sd, and n for spherical power
r_cyl_mean_0 = apply(r_cyl_0, 1, function(row) { mean(row, na.rm=T) })
r_cyl_sd_0 = apply(r_cyl_0, 1, function(row) { sd(row, na.rm=T) })
r_cyl_n_0 = apply(r_cyl_0, 1, function(row) { sum(!is.na(row)) })

# Calc mean spherical equivalence
r_MSE_0 = r_sph_mean_0 + (0.5 * r_cyl_mean_0)

#
# Calculate average MSE, CYL for instance 0
#

ave_MSE_0 = rowMeans(cbind(l_MSE_0, r_MSE_0), na.rm=FALSE)
ave_cyl_0 = rowMeans(cbind(l_cyl_mean_0, r_cyl_mean_0), na.rm=FALSE)

#
# Do left RE for instance 1
#

# Get errors
l_refUnrel_1 = grep_cols(df8434, "f.5090.1.*")
# Create logical matrix
is_err = apply(l_refUnrel_1, 2, function(v) { grepl("e|E", v) })

# Get spherical power
l_sph_1 = grep_cols(df8434, "f.5085.1.*")
# Set errors to NA
l_sph_1 = as.matrix(l_sph_1)
l_sph_1[is_err] = NA
# Cal mean, sd, and n for spherical power
l_sph_mean_1 = apply(l_sph_1, 1, function(row) { mean(row, na.rm=T) })
l_sph_sd_1 = apply(l_sph_1, 1, function(row) { sd(row, na.rm=T) })
l_sph_n_1 = apply(l_sph_1, 1, function(row) { sum(!is.na(row)) })

# Get cylindrial power
l_cyl_1 = grep_cols(df8434, "f.5086.1.*")
# Set errors to NA
l_cyl_1 = as.matrix(l_cyl_1)
l_cyl_1[is_err] = NA
# Cal mean, sd, and n for spherical power
l_cyl_mean_1 = apply(l_cyl_1, 1, function(row) { mean(row, na.rm=T) })
l_cyl_sd_1 = apply(l_cyl_1, 1, function(row) { sd(row, na.rm=T) })
l_cyl_n_1 = apply(l_cyl_1, 1, function(row) { sum(!is.na(row)) })

# Calc mean spherical equivalence
l_MSE_1 = l_sph_mean_1 + (0.5 * l_cyl_mean_1)

#
# Do right RE for instance 1
#

# Get errors
r_refUnrel_1 = grep_cols(df8434, "f.5091.1.*")
# Create logical matrix
is_err = apply(r_refUnrel_1, 2, function(v) { grepl("e|E", v) })

# Get spherical power
r_sph_1 = grep_cols(df8434, "f.5084.1.*")
# Set errors to NA
r_sph_1 = as.matrix(r_sph_1)
r_sph_1[is_err] = NA
# Cal mean, sd, and n for spherical power
r_sph_mean_1 = apply(r_sph_1, 1, function(row) { mean(row, na.rm=T) })
r_sph_sd_1 = apply(r_sph_1, 1, function(row) { sd(row, na.rm=T) })
r_sph_n_1 = apply(r_sph_1, 1, function(row) { sum(!is.na(row)) })

# Get cylindrial power
r_cyl_1 = grep_cols(df8434, "f.5087.1.*")
# Set errors to NA
r_cyl_1 = as.matrix(r_cyl_1)
r_cyl_1[is_err] = NA
# Cal mean, sd, and n for spherical power
r_cyl_mean_1 = apply(r_cyl_1, 1, function(row) { mean(row, na.rm=T) })
r_cyl_sd_1 = apply(r_cyl_1, 1, function(row) { sd(row, na.rm=T) })
r_cyl_n_1 = apply(r_cyl_1, 1, function(row) { sum(!is.na(row)) })

# Calc mean spherical equivalence
r_MSE_1 = r_sph_mean_1 + (0.5 * r_cyl_mean_1)

#
# Calculate average MSE, CYL for instance 1
#

ave_MSE_1 = rowMeans(cbind(l_MSE_1, r_MSE_1), na.rm=FALSE)
ave_cyl_1 = rowMeans(cbind(l_cyl_mean_1, r_cyl_mean_1), na.rm=FALSE)

#
# Take instance 0 if it exists, and instance 1 if it doesn't
#

# MSE
ave_MSE = ave_MSE_0
to_replace = is.na(ave_MSE_0) & !is.na(ave_MSE_1)
ave_MSE[to_replace] = ave_MSE_1[to_replace]
sum(!is.na(ave_MSE))

#
# Make excluions based on history
#

cat("Refractive error exclusions:\n", file=outlog, sep=" ",append=TRUE)
re_isOutlier = rep(FALSE, length(ave_MSE))

# 6148: Eye problems - combining across both instances
excl_eyeProb = grep_cols(df8434, "f.6148.*.*")
# Set as outlier if has history of Cataract
hasCat1 = apply(excl_eyeProb == "Cataract", 1, function(row) { any(row, na.rm=T) } )
re_isOutlier[hasCat1] = TRUE
cat("Num with cataracts (Eye problems multicategory):", sum(hasCat1), "\n", file=outlog, sep=" ",append=TRUE)

# 5324: Cataract surgery
excl_cataract = grep_cols(df8434, "f.5324.*.*")
# Set as outlier if not No. Exclude Don't knows (82)
hasCat2 = apply(excl_cataract != "No" & !is.na(excl_cataract), 1, function(row) { any(row, na.rm=T) } )
re_isOutlier[hasCat2] = TRUE
cat("Num with cataracts (surgery):", sum(hasCat2), "\n", file=outlog, sep=" ",append=TRUE)

# 5441: History of cataracts
excl_catHist = grep_cols(df8434, "f.5441.*.*")
hasCat3 = apply(!is.na(excl_catHist), 1, function(row) { any(row, na.rm=T) } )
re_isOutlier[hasCat3] = TRUE
cat("Num with cataract history:", sum(hasCat3), "\n", file=outlog, sep=" ",append=TRUE)

sum(hasCat1)
sum(hasCat2)
sum(hasCat3)
sum(hasCat1 | hasCat2 | hasCat3)

# 5325: Laser eye surgery
excl_laser = grep_cols(df8434, "f.5325.*.*")
# Set as outlier if not No. Exclude Don't knows (229)
to_excl = apply(excl_laser != "No" & !is.na(excl_laser), 1, function(row) { any(row, na.rm=T) } )
re_isOutlier[to_excl] = TRUE
cat("Num with laser eye surgery:", sum(to_excl), "\n", file=outlog, sep=" ",append=TRUE)

# 5419: History of trauma - any answer should be excluded
excl_trauma = grep_cols(df8434, "f.5419.*.*")
to_excl = apply(!is.na(excl_trauma), 1, function(row) { any(row, na.rm=T) } )
re_isOutlier[to_excl] = TRUE
cat("Num with eye trauma:", sum(to_excl), "\n", file=outlog, sep=" ",append=TRUE)

# 5328: Had corneal graft surgery
excl_cornGraft = grep_cols(df8434, "f.5328.*.*")
# Set as outlier if not No. Exclude Don't knows (229)
to_excl = apply(excl_cornGraft != "No" & !is.na(excl_cornGraft), 1, function(row) { any(row, na.rm=T) } )
re_isOutlier[to_excl] = TRUE
cat("Num with corneal graft surgery:", sum(to_excl), "\n", file=outlog, sep=" ",append=TRUE)

#
# Clean ave_MSE
#

ave_MSE_clean = ave_MSE
ave_MSE_clean[re_isOutlier] = NA

#
# Descriptives
#

cat("\nRefractive error descriptives:\n\n", file=outlog, sep=" ",append=TRUE)

# Cleaned
cat("ave_MSE_clean:\n", file=outlog, sep=" ",append=TRUE)
write.table(describe(ave_MSE_clean), file=outlog, append=TRUE, quote=F)

# Raw
cat("\nave_MSE:\n", file=outlog, sep=" ",append=TRUE)
write.table(describe(ave_MSE), file=outlog, append=TRUE, quote=F)

# Work out how many valid MSE and edu scores
sum(!is.na(ave_MSE_clean[names(eduyears_clean)] & eduyears_clean))

#
# Derive myopia from refractive error ##########################################
#

cat("\nMyopic samples:\n", file=outlog, sep=" ",append=TRUE)

# Myopia SE < -0.75
isMyopic = NA
isMyopic[!is.na(ave_MSE_clean)] = 0
isMyopic[ave_MSE_clean <= -0.75] = 1
names(isMyopic) = names(ave_MSE_clean)
isMyopic = factor(isMyopic)
cat("\nNum myopic (SE <= -0.75:\n", sum(isMyopic == 1), file=outlog, sep=" ",append=TRUE)

# Myopia SE < -6
isMyopicSevere = NA
isMyopicSevere[!is.na(ave_MSE_clean)] = 0
isMyopicSevere[ave_MSE_clean <= -6] = 1
names(isMyopicSevere) = names(ave_MSE_clean)
isMyopicSevere = factor(isMyopicSevere)
cat("\nNum severe myopic (SE <= -6:\n", sum(isMyopicSevere == 1), file=outlog, sep=" ",append=TRUE)

#
# Derive visual acuity logMAR ##################################################
#

#
# Get left and right final logMAR, there are 2 instances (0 and 1)
#

# Get left logMAR final
l_logmar = to_named_vector(merge_instances(grep_cols(df8434, "f.5208.*.0")))

# Get right logMAR final
r_logmar = to_named_vector(merge_instances(grep_cols(df8434, "f.5201.*.0")))

#
# Calc average, worst, best eyes
#

# Average
ave_logmar = rowMeans(cbind(l_logmar, r_logmar), na.rm=FALSE)
# Worst (max logmar)
worst_logmar = apply(cbind(l_logmar, r_logmar), 1, max, na.rm=FALSE)
# Best (min logmar)
best_logmar = apply(cbind(l_logmar, r_logmar), 1, min, na.rm=FALSE)

# Which is worst
which_worse = NA
which_worse[l_logmar > r_logmar] = "left"
which_worse[r_logmar > l_logmar] = "right"
which_worse[r_logmar == l_logmar] = "neither"

# Get MSE so to use as covariate
l_MSE = ifelse(!is.na(l_MSE_0), l_MSE_0, l_MSE_1)
r_MSE = ifelse(!is.na(r_MSE_0), r_MSE_0, r_MSE_1)
worst_covMSE = NA
worst_covMSE[which_worse == "left" & !is.na(which_worse)] = l_MSE[which_worse == "left" & !is.na(which_worse)]
worst_covMSE[which_worse == "right" & !is.na(which_worse)] = r_MSE[which_worse == "right" & !is.na(which_worse)]
worst_covMSE[which_worse == "neither" & !is.na(which_worse)] = ave_MSE[which_worse == "neither" & !is.na(which_worse)]
names(worst_covMSE) = names(ave_logmar)

# Code whether amblyopic or not
isAmblyope = factor(ifelse(abs(best_logmar - worst_logmar) > 0.2, 1, 0))
names(isAmblyope) = names(ave_logmar)

# Get right - left difference in RE
rl_MSE_diff = r_MSE-l_MSE
rl_MSE_diff_clean = rl_MSE_diff
rl_MSE_diff_clean[re_isOutlier] = NA

#
# Identify outliers
#

cat("\nVisual acuity (logMAR) exclusions:\n", file=outlog, sep=" ",append=TRUE)

# Same as refractive error + more
logmar_isOutlier = re_isOutlier

# 6148: Eye problems - anything that isn't "None of the above"
excl_eyeProb = grep_cols(df8434, "f.6148.*.*")
# Set as outlier if has history of Cataract
to_excl = apply(excl_eyeProb != "None of the above" & !is.na(excl_eyeProb), 1, function(row) { any(row, na.rm=T) } )
logmar_isOutlier[to_excl] = TRUE
cat("Num not 'None of the above' (Eye problems multicategory):", sum(to_excl), "\n", file=outlog, sep=" ",append=TRUE)

# 5326: Surgery for glaucoma or high eye pressure
excl_surgGluc = grep_cols(df8434, "f.5326.*.*")
to_excl = apply(excl_surgGluc != "No" & !is.na(excl_surgGluc), 1, function(row) { any(row, na.rm=T) } )
logmar_isOutlier[to_excl] = TRUE
cat("Num with glaucoma surgery:", sum(to_excl), "\n", file=outlog, sep=" ",append=TRUE)

# 5327: Had laser treatment for glaucoma or high eye pressure
excl_laserGluc = grep_cols(df8434, "f.5327.*.*")
to_excl = apply(excl_laserGluc != "No" & !is.na(excl_laserGluc), 1, function(row) { any(row, na.rm=T) } )
logmar_isOutlier[to_excl] = TRUE
cat("Num with glaucoma laser surgery:", sum(to_excl), "\n", file=outlog, sep=" ",append=TRUE)

# 6074: Has glasses but is not wearing them (right)
excl_noGlassR = grep_cols(df8434, "f.6074.*.*")
to_excl = apply(excl_noGlassR == "elsewhere", 1, function(row) { any(row, na.rm=T) } )
logmar_isOutlier[to_excl] = TRUE
cat("Num with glasses needed but not worn (right):", sum(to_excl), "\n", file=outlog, sep=" ",append=TRUE)

# 6074: Has glasses but is not wearing them (left)
excl_noGlassL = grep_cols(df8434, "f.6075.*.*")
to_excl = apply(excl_noGlassL == "elsewhere", 1, function(row) { any(row, na.rm=T) } )
logmar_isOutlier[to_excl] = TRUE
cat("Num with glasses needed but not worn (left):", sum(to_excl), "\n", file=outlog, sep=" ",append=TRUE)

#
# Clean variables
#

# Average
ave_logmar_clean = ave_logmar
ave_logmar_clean[logmar_isOutlier] = NA
# Worst
worst_logmar_clean = worst_logmar
worst_logmar_clean[logmar_isOutlier] = NA
# Best
best_logmar_clean = best_logmar
best_logmar_clean[logmar_isOutlier] = NA
# Amblyope
isAmblyope_clean = isAmblyope
isAmblyope_clean[logmar_isOutlier] = NA

#
# Transform continuous
#

# Log transform
ave_logmar_clean_log = log_transform(ave_logmar_clean)
best_logmar_clean_log  = log_transform(best_logmar_clean)
worst_logmar_clean_log = log_transform(worst_logmar_clean)

# Dichotomise to match Gen R
# Average
ave_logmar_clean_cc0 = NA
ave_logmar_clean_cc0[ave_logmar_clean <= 0] = 0
ave_logmar_clean_cc0[ave_logmar_clean > 0] = 1
ave_logmar_clean_cc0 = factor(ave_logmar_clean_cc0)
names(ave_logmar_clean_cc0) = names(ave_logmar_clean)
# Best
best_logmar_clean_cc0 = NA
best_logmar_clean_cc0[best_logmar_clean <= 0] = 0
best_logmar_clean_cc0[best_logmar_clean > 0] = 1
best_logmar_clean_cc0 = factor(best_logmar_clean_cc0)
names(best_logmar_clean_cc0) = names(best_logmar_clean)
# Worst
worst_logmar_clean_cc0 = NA
worst_logmar_clean_cc0[worst_logmar_clean <= 0] = 0
worst_logmar_clean_cc0[worst_logmar_clean > 0] = 1
worst_logmar_clean_cc0 = factor(worst_logmar_clean_cc0)
names(worst_logmar_clean_cc0) = names(worst_logmar_clean)

#
# Descriptives
#

cat("\nVisual acuity logMAR descriptives:\n\n", file=outlog, sep=" ",append=TRUE)

# Cleaned
cat("ave_logmar_clean:\n", file=outlog, sep=" ",append=TRUE)
write.table(describe(ave_logmar_clean), file=outlog, append=TRUE, quote=F)
cat("best_logmar_clean:\n", file=outlog, sep=" ",append=TRUE)
write.table(describe(best_logmar_clean), file=outlog, append=TRUE, quote=F)
cat("worst_logmar_clean:\n", file=outlog, sep=" ",append=TRUE)
write.table(describe(worst_logmar_clean), file=outlog, append=TRUE, quote=F)

# Raw
cat("ave_logmar:\n", file=outlog, sep=" ",append=TRUE)
write.table(describe(ave_logmar), file=outlog, append=TRUE, quote=F)
cat("best_logmar:\n", file=outlog, sep=" ",append=TRUE)
write.table(describe(best_logmar), file=outlog, append=TRUE, quote=F)
cat("worst_logmar:\n", file=outlog, sep=" ",append=TRUE)
write.table(describe(worst_logmar), file=outlog, append=TRUE, quote=F)

#
# Derive non-genetic covariates ################################################
#

# Get TDI and log transform
tdi = to_named_vector(grep_cols(df4263, "f.189.0.0"))
tdi_log = log_transform(tdi)

# Get age attended assessment centre
age = to_named_vector(merge_instances(grep_cols(df4263, "f.21003.*.0")))

# Birthweight
birthweight = to_named_vector(merge_instances(grep_cols(df4263, "f.20022.*.*")))

# Breastfed
bfdata = to_named_vector(merge_instances(grep_cols(df4263, "f.1677.*.*")))
breastfed = NA
breastfed[bfdata == "Yes"] = 1
breastfed[bfdata == "No"] = 0
breastfed = factor(breastfed)
names(breastfed) = names(bfdata)

# Northing
northing = to_named_vector(merge_instances(grep_cols(df8434, "f.129.*.*")))
# Easting
easting = to_named_vector(merge_instances(grep_cols(df8434, "f.130.*.*")))

# Load sex
sex = to_named_vector(grep_cols(df4263, "f.31.0.0"))

#
# Derive genetic covariates ####################################################
#

# Load data_covariates.txt
ukbb_cov_df = read_delim("../data/genetics/data.covariates.txt", " ", col_names=F)
ukbb_cov_df$X3 = ukbb_cov_df$X3 %>% recode_factor(M="Male", F="Female")
ukbb_cov_df$X4 = factor(ukbb_cov_df$X4)

length(ukbb_cov_df$X1)
length(unique(ukbb_cov_df$X1))

# Load genetic sex
sex_genetic = ukbb_cov_df$X3
names(sex_genetic) = ukbb_cov_df$X1

# Chip genotyped on
genechip = ukbb_cov_df$X4
names(genechip) = ukbb_cov_df$X1

# Load PC data
ukbb_pc_df = read_delim("../data/genetics/data.pca1-10.txt", " ", col_names=F)
pcs = data.frame(ukbb_pc_df[, -c(1,2)])
head(pcs)
rownames(pcs) = ukbb_pc_df$X1
colnames(pcs) = paste0("PC", 1:10)

length(ukbb_pc_df$X1)
length(unique(ukbb_pc_df$X1))

#
# Derive questionnaires Qs for predicting myopia ###############################
#

# 2207	Wears glasses or contact lenses
q_wearsGlasses = to_named_vector(merge_instances(grep_cols(df8434, "f.2207.*.0")))

# 2217	Age started wearing glasses or contact lenses
q_glassesAge = to_named_vector(merge_instances(grep_cols(df8434, "f.2217.*.0")))
# Set -3 (Prefer not to say) and -1 (Do no know) to NA
q_glassesAge[q_glassesAge == -3 & !is.na(q_glassesAge)] = NA
q_glassesAge[q_glassesAge == -1 & !is.na(q_glassesAge)] = NA

# 6147	Reason for glasses/contact lenses
q_glassesReason = grep_cols(df8434, "f.6147.*.*")
q_glassesReason_bools = process_multicategory(q_glassesReason)
# Keep only myopia, hypermetropia, presbyopia, astigmatism, strabismusm, amblyopia
q_glassesReason_bools = q_glassesReason_bools[, c(3, 4, 5, 6, 7, 8)]
colnames(q_glassesReason_bools) = paste0("q_glassesReason_", c("myopia", "hypermetropia", "presbyopia", "astigmatism", "strabismusm", "amblyopia"))

# 5843	Which eye(s) affected by myopia (short sight)
#q_whichEye_myopia = to_named_vector(merge_instances(grep_cols(df8434, "f.5843.*.0")))
# 5832	Which eye(s) affected by hypermetropia (long sight)
#q_whichEye_hypermetropia = to_named_vector(merge_instances(grep_cols(df8434, "f.5832.*.0")))
# 5610	Which eye(s) affected by presbyopia
#q_whichEye_presbyopia = to_named_vector(merge_instances(grep_cols(df8434, "f.5610.*.0")))


#
# Output variables #############################################################
#

# Get list of all sample names
sample_names = unique(c(names(ave_MSE_clean),
                        names(isMyopic),
                        names(isMyopicSevere),
                        names(isAmblyope_clean),
                        names(ave_logmar_clean_log),
                        names(best_logmar_clean_log),
                        names(worst_logmar_clean_log),
                        names(worst_logmar_clean),
                        names(rl_MSE_diff_clean),
                        names(eduyears_clean),
                        names(tdi_log),
                        names(age),
                        names(sex),
                        names(sex_genetic),
                        names(genechip),
                        names(eduyearsOkbay_clean),
                        names(eduyearsNo21_clean),
                        names(eduyearsMore16),
                        names(birthweight),
                        names(breastfed),
                        names(northing),
                        names(easting),
                        names(q_wearsGlasses),
                        names(q_glassesAge),
                        rownames(q_glassesReason_bools)
                      ))

# Make dataframe
data = data.frame(ave_MSE_clean = ave_MSE_clean[match(sample_names, names(ave_MSE_clean))],
                  isMyopic = isMyopic[match(sample_names, names(isMyopic))],
                  isMyopicSevere = isMyopicSevere[match(sample_names, names(isMyopicSevere))],
                  isAmblyope_clean = isAmblyope_clean[match(sample_names, names(isAmblyope_clean))],
                  ave_logmar_clean_log = ave_logmar_clean_log[match(sample_names, names(ave_logmar_clean_log))],
                  best_logmar_clean_log = best_logmar_clean_log[match(sample_names, names(best_logmar_clean_log))],
                  worst_logmar_clean_log = worst_logmar_clean_log[match(sample_names, names(worst_logmar_clean_log))],
                  worst_logmar_clean = worst_logmar_clean[match(sample_names, names(worst_logmar_clean))],
                  rl_MSE_diff_clean = rl_MSE_diff_clean[match(sample_names, names(rl_MSE_diff_clean))],
                  ave_logmar_clean_cc0 = ave_logmar_clean_cc0[match(sample_names, names(ave_logmar_clean_cc0))],
                  best_logmar_clean_cc0 = best_logmar_clean_cc0[match(sample_names, names(best_logmar_clean_cc0))],
                  worst_logmar_clean_cc0 = worst_logmar_clean_cc0[match(sample_names, names(worst_logmar_clean_cc0))],
                  worst_covMSE = worst_covMSE[match(sample_names, names(worst_covMSE))],
                  eduyears_clean = eduyears_clean[match(sample_names, names(eduyears_clean))],
                  tdi_log = tdi_log[match(sample_names, names(tdi_log))],
                  age = age[match(sample_names, names(age))],
                  sex = sex[match(sample_names, names(sex))],
                  sex_genetic = sex_genetic[match(sample_names, names(sex_genetic))],
                  genechip = genechip[match(sample_names, names(genechip))],
                  eduyearsOkbay_clean = eduyearsOkbay_clean[match(sample_names, names(eduyearsOkbay_clean))],
                  eduyearsNo21_clean = eduyearsNo21_clean[match(sample_names, names(eduyearsNo21_clean))],
                  eduyearsMore16 = eduyearsMore16[match(sample_names, names(eduyearsMore16))],
                  birthweight = birthweight[match(sample_names, names(birthweight))],
                  breastfed = breastfed[match(sample_names, names(breastfed))],
                  northing = northing[match(sample_names, names(northing))],
                  easting = easting[match(sample_names, names(easting))],
                  q_wearsGlasses = q_wearsGlasses[match(sample_names, names(q_wearsGlasses))],
                  q_glassesAge = q_glassesAge[match(sample_names, names(q_glassesAge))]
                  )
rownames(data) = sample_names

# Merge principle components dataframe
eid = sample_names
data2 = cbind(eid,
              data,
              pcs[match(sample_names, rownames(pcs)), ],
              q_glassesReason_bools[match(sample_names, rownames(q_glassesReason_bools)), ]
              )
data = data2

#
# Label core sample and split sample
#

# Core sample
data$coreSample = ifelse(apply(data[, c("ave_MSE_clean", "eduyears_clean", "age", "sex")], 1, function(row) {all(!is.na(row))} ), T, F)
sum(data$coreSample)

# Randomly assign to A or B
set.seed(123456)
data[data$coreSample, "splitSample"] = sample(c("A", "B"), size=sum(data$coreSample), replace=T)

# Core with genetics
data$coreSampleGeno = ifelse(data$coreSample == T & !is.na(data$genechip), T, F)
sum(data$coreSampleGeno)
table(data[data$coreSampleGeno == T, "splitSample"])

#
# Make inclusions / exclusions #################################################
#

# Load
excl_combined = read.table("../data/genetics/exclusions/combined_exclusions_170728.txt", sep=" ", as.is=T)[,1]
# Find samples to exclude
toExclude = sample_names %in% excl_combined
sum(toExclude)
# Exclude by setting to NA
data$genoExcl = toExclude

# Total sample size
data %>% filter(coreSampleGeno & !genoExcl) %>% dim

#
# Write output
#


save(data, file=outdata)

