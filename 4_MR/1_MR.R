#!/usr/bin/env Rscript
#
# Run MR in R using Stephen Burgess code from:
#   https://github.com/sb452/mr-code/blob/mr-code/ivrcode.pdf
#

setwd("/mnt/seconddrive/data/phd/ukbiobank/myopia_EA_project/4_MR/")

library("tidyverse")
library("sem")
library("ivpack")
library("sandwich")
library("lmtest")

# Args
indata = "../2_derive_variables/output/phenotypes_alleleScores_170728.Rdata"
outfile = "output/MR_results_170724.txt"

# Load
load(indata)

#
# Make inclusions / exclusions #################################################
#

data = data %>% filter(coreSampleGeno & !genoExcl)

#
# Define MR function ###########################################################
#

# Run an MR with code: https://github.com/sb452/mr-code/blob/mr-code/ivrcode.pdf
#  data = a data.frame
#  g, x, y = column names in data
#  c = vector of column names in data
run.MR = function(data, gcol, xcol, ycol, wcol=NULL, ccols=NULL) {

  ### Initiate results list
  res = list()
  
  ### Get complete cases
  data_cc = data[, c(gcol, xcol, ycol, wcol, ccols)]
  data_cc = data_cc[complete.cases(data_cc), ]
  # Define instruments (g) and num samples
  g = data_cc[, gcol]
  N = nrow(data_cc)
  # Save info to res
  res["g"] = gcol
  res["x"] = xcol
  res["y"] = ycol
  res["w"] = wcol
  res["covars"] = paste(ccols, collapse=" + ")
  res["N"] = N
  
  ### Regress out covariates
  if (is.null(ccols)) {
    x = data_cc[, xcol]
    y = data_cc[, ycol]
  } else {
    f = formula(paste0(xcol, " ~ ", paste(ccols, collapse=" + ")))
    x = as.numeric(lm(f, data=data_cc)$residuals)
    f = formula(paste0(ycol, " ~ ", paste(ccols, collapse=" + ")))
    y = as.numeric(lm(f, data=data_cc)$residuals)
  }
  
  # Do regressions
  if (is.null(wcol)) {
    xmod = lm(x ~ g)
    ymod = lm(y ~ g)
    ivmodel = ivreg(y~x|g, x=TRUE)
  } else {
    w = data_cc[, wcol]
    xmod = lm(x ~ g, weights=w)
    ymod = lm(y ~ g, weights=w)
    ivmodel = ivreg(y~x|g, weights=w, x=TRUE)
  }
  
  ### Estimate ratio manually (weighting both sides of the regression)
  bx = xmod$coef[2]
  bxse = summary(xmod)$coef[2,2]
  by = ymod$coef[2]
  byse = summary(ymod)$coef[2,2]
  beta_ratio = by / bx
  res["bx"] = paste0(bx, " ± ", bxse, " SE")
  res["by"] = paste0(by, " ± ", byse, " SE")
  res["ratio"] = beta_ratio
  
  ### Using ivreg
  summ = summary(ivmodel, diagnostics=T, vcov=sandwich)
  res["beta_tsls_ivreg"] = summ$coef[2,1]
  res["se_tsls_ivreg"] = summ$coef[2,2]
  res["p_tsls_ivreg"] = summ$coef[2,4]
  res["weak_instr_diagnostic"] = summ$diagnostics[1,4]
  res["wu-hausman_diagnostic"] = summ$diagnostics[2,4]
  res["sargan_diagnostic"] = summ$diagnostics[3,4]
  res["Anderson.Rubin.CI"] = anderson.rubin.ci(ivmodel)$confidence.interval
  
  ### Return results object
  res
}

#
# Calculate inverse probability weights
#

propAge = cumsum(table(data$eduyears_clean)) / sum(!is.na(data$eduyears_clean))
data = data %>% mutate(ipw = ifelse(eduyears_clean == 15, (0.33 / propAge["15"]), 1))


#
# Run MR:  #####################################################################
#

cat("MR results\n==========\n", file=outfile, append=F)
all_results = list()

#
# Main analysis
#

# Exposure/outcome = Refractive error <-> Educational Attainment
# Instruments      = Pickrell et al Myopia Instruments

cat("\n### Exposure/outcome = Refractive error <-> Educational Attainment\n", file=outfile, append=T)
cat("### Instruments      = Pickrell et al Myopia Instruments\n", file=outfile, append=T)

# RE -> EA - age + sex (not weighted)
cat("\n# RE -> EA\n", file=outfile, append=T)
res = run.MR(data, "myopia_AS_dosage", "ave_MSE_clean", "eduyears_clean", ccols=c("sex_genetic", "age"))
all_results[[length(all_results) + 1]] = res
write.table(t(data.frame(res)), col.names=F, quote=F, sep=":\t", file=outfile, append=T)
# RE -> EA - age + sex (weighted)
cat("\n# RE -> EA\n", file=outfile, append=T)
res = run.MR(data, "myopia_AS_dosage", "ave_MSE_clean", "eduyears_clean", wcol="ipw", ccols=c("sex_genetic", "age"))
all_results[[length(all_results) + 1]] = res
write.table(t(data.frame(res)), col.names=F, quote=F, sep=":\t", file=outfile, append=T)

# EA -> RE - age + sex (weighted)
cat("\n# EA -> RE\n", file=outfile, append=T)
res = run.MR(data, "ea_AS_dosage", "eduyears_clean", "ave_MSE_clean", ccols=c("sex_genetic", "age"))
all_results[[length(all_results) + 1]] = res
write.table(t(data.frame(res)), col.names=F, quote=F, sep=":\t", file=outfile, append=T)
# EA -> RE (not weighted)
cat("\n# EA -> RE\n", file=outfile, append=T)
res = run.MR(data, "ea_AS_dosage", "eduyears_clean", "ave_MSE_clean", wcol="ipw", ccols=c("sex_genetic", "age"))
all_results[[length(all_results) + 1]] = res
write.table(t(data.frame(res)), col.names=F, quote=F, sep=":\t", file=outfile, append=T)

#
# Sensitiviy analyses 
#

# Excluding university/college educated
#

cat("\n### Exposure/outcome = Refractive error <-> Educational Attainment (no college)\n", file=outfile, append=T)

# RE -> EA - age + sex (not weighted)
cat("\n# RE -> EA\n", file=outfile, append=T)
res = run.MR(data, "myopia_AS_dosage", "ave_MSE_clean", "eduyearsNo21_clean", ccols=c("sex_genetic", "age"))
all_results[[length(all_results) + 1]] = res
write.table(t(data.frame(res)), col.names=F, quote=F, sep=":\t", file=outfile, append=T)
# RE -> EA - age + sex (weighted)
cat("\n# RE -> EA\n", file=outfile, append=T)
res = run.MR(data, "myopia_AS_dosage", "ave_MSE_clean", "eduyearsNo21_clean", wcol="ipw", ccols=c("sex_genetic", "age"))
all_results[[length(all_results) + 1]] = res
write.table(t(data.frame(res)), col.names=F, quote=F, sep=":\t", file=outfile, append=T)

# EA -> RE - age + sex (weighted)
cat("\n# EA -> RE\n", file=outfile, append=T)
res = run.MR(data, "ea_AS_dosage", "eduyearsNo21_clean", "ave_MSE_clean", ccols=c("sex_genetic", "age"))
all_results[[length(all_results) + 1]] = res
write.table(t(data.frame(res)), col.names=F, quote=F, sep=":\t", file=outfile, append=T)
# EA -> RE (not weighted)
cat("\n# EA -> RE\n", file=outfile, append=T)
res = run.MR(data, "ea_AS_dosage", "eduyearsNo21_clean", "ave_MSE_clean", wcol="ipw", ccols=c("sex_genetic", "age"))
all_results[[length(all_results) + 1]] = res
write.table(t(data.frame(res)), col.names=F, quote=F, sep=":\t", file=outfile, append=T)

# Excluding using binary education variable: >=16 vs <16
#

cat("\n### Exposure/outcome = Refractive error <-> Educational Attainment (binary)\n", file=outfile, append=T)
colnames(data)

# EA -> RE - age + sex (not weighted)
xmod = glm(eduyearsMore16~ea_AS_dosage+sex_genetic+age, data=data, family="binomial")
bx   = xmod$coef[2]
bxse = coeftest(xmod, vcov = sandwich)[2,2]
ymod = lm(ave_MSE_clean~ea_AS_dosage+sex_genetic+age, data=data)
by   = ymod$coef[2]
byse = coeftest(ymod, vcov = sandwich)[2,2]
beta_ratio = by/bx
N = nrow(data)
res = list()
res["g"] = "ea_AS_dosage"
res["x"] = "eduyearsMore16"
res["y"] = "ave_MSE_clean"
res["covars"] = "sex_genetic + age"
res["N"] = N
res["bx"] = bx
res["by"] = by
res["ratio"] = beta_ratio
res["bxse"] = bxse
res["byse"] = byse
se_ratio_approx_1 = byse/bx
se_ratio_approx_2 = sqrt(byse^2/bx^2 + by^2*bxse^2/bx^4 - 2*1*by/bx^3)
res["se_ratio_approx_1"] = se_ratio_approx_1
res["se_ratio_approx_2"] = se_ratio_approx_2
f0 = by^2 - qt(0.975, N)^2 * byse^2
f1 = bx^2 - qt(0.975, N)^2 * bxse^2
f2 = by*bx
D = f2^2 - f0*f1
if(D>0) {
  r1 = (f2-sqrt(D))/f1
  r2 = (f2+sqrt(D))/f1
  if(f1>0) { res["Fiellers.Theorem.CI"] = paste0("Confidence interval is a closed interval [a,b]: a=", r1, ", b=", r2, sep="") }
  if(f1<0) { res["Fiellers.Theorem.CI"] = paste0("Confidence interval is the union of two open intervals (-Inf, a], [b, +Inf): a=", r2, ", b=", r1, sep="") }
}
if(D<0|D==0) { res["Fiellers.Theorem.CI"] = paste0("No finite confidence interval exists other than the entire real line.") }
cat("\n# EA -> RE - age + sex (not weighted)\n", file=outfile, append=T)
all_results[[length(all_results) + 1]] = res
write.table(t(data.frame(res)), col.names=F, quote=F, sep=":\t", file=outfile, append=T)

# EA -> RE - age + sex (weighted)
xmod = glm(eduyearsMore16~ea_AS_dosage+sex_genetic+age, data=data, weights=ipw, family="binomial")
bx   = xmod$coef[2]
bxse = coeftest(xmod, vcov = sandwich)[2,2]
ymod = lm(ave_MSE_clean~ea_AS_dosage+sex_genetic+age, weights=ipw, data=data)
by   = ymod$coef[2]
byse = coeftest(ymod, vcov = sandwich)[2,2]
beta_ratio = by/bx
N = nrow(data)
res = list()
res["g"] = "ea_AS_dosage"
res["x"] = "eduyearsMore16"
res["y"] = "ave_MSE_clean"
res["covars"] = "sex_genetic + age"
res["w"] = "ipw"
res["N"] = N
res["bx"] = bx
res["by"] = by
res["ratio"] = beta_ratio
res["bxse"] = bxse
res["byse"] = byse
se_ratio_approx_1 = byse/bx
se_ratio_approx_2 = sqrt(byse^2/bx^2 + by^2*bxse^2/bx^4 - 2*1*by/bx^3)
res["se_ratio_approx_1"] = se_ratio_approx_1
res["se_ratio_approx_2"] = se_ratio_approx_2
f0 = by^2 - qt(0.975, N)^2 * byse^2
f1 = bx^2 - qt(0.975, N)^2 * bxse^2
f2 = by*bx
D = f2^2 - f0*f1
if(D>0) {
  r1 = (f2-sqrt(D))/f1
  r2 = (f2+sqrt(D))/f1
  if(f1>0) { res["Fiellers.Theorem.CI"] = paste0("Confidence interval is a closed interval [a,b]: a=", r1, ", b=", r2, sep="") }
  if(f1<0) { res["Fiellers.Theorem.CI"] = paste0("Confidence interval is the union of two open intervals (-Inf, a], [b, +Inf): a=", r2, ", b=", r1, sep="") }
}
if(D<0|D==0) { res["Fiellers.Theorem.CI"] = paste0("No finite confidence interval exists other than the entire real line.") }
cat("\n# EA -> RE - age + sex (not weighted)\n", file=outfile, append=T)
all_results[[length(all_results) + 1]] = res
write.table(t(data.frame(res)), col.names=F, quote=F, sep=":\t", file=outfile, append=T)

# RE -> EA - age + sex (not weighted)
colnames(data)
xmod = lm(ave_MSE_clean ~ myopia_AS_dosage + sex_genetic + age, data=data)
bx   = xmod$coef[2]
bxse = coeftest(xmod, vcov = sandwich)[2,2]
ymod = glm(eduyearsMore16 ~ myopia_AS_dosage + sex_genetic + age, data=data, family="binomial")
by   = ymod$coef[2]
byse = coeftest(ymod, vcov = sandwich)[2,2]
beta_ratio = by/bx
N = nrow(data)
res = list()
res["g"] = "myopia_AS_dosage"
res["x"] = "ave_MSE_clean"
res["y"] = "eduyearsMore16"
res["covars"] = "sex_genetic + age"
res["N"] = N
res["bx"] = bx
res["by"] = by
res["ratio"] = beta_ratio
res["bxse"] = bxse
res["byse"] = byse
se_ratio_approx_1 = byse/bx
se_ratio_approx_2 = sqrt(byse^2/bx^2 + by^2*bxse^2/bx^4 - 2*1*by/bx^3)
res["se_ratio_approx_1"] = se_ratio_approx_1
res["se_ratio_approx_2"] = se_ratio_approx_2
f0 = by^2 - qt(0.975, N)^2 * byse^2
f1 = bx^2 - qt(0.975, N)^2 * bxse^2
f2 = by*bx
D = f2^2 - f0*f1
if(D>0) {
  r1 = (f2-sqrt(D))/f1
  r2 = (f2+sqrt(D))/f1
  if(f1>0) { res["Fiellers.Theorem.CI"] = paste0("Confidence interval is a closed interval [a,b]: a=", r1, ", b=", r2, sep="") }
  if(f1<0) { res["Fiellers.Theorem.CI"] = paste0("Confidence interval is the union of two open intervals (-Inf, a], [b, +Inf): a=", r2, ", b=", r1, sep="") }
}
if(D<0|D==0) { res["Fiellers.Theorem.CI"] = paste0("No finite confidence interval exists other than the entire real line.") }
cat("\n# EA -> RE - age + sex (not weighted)\n", file=outfile, append=T)
all_results[[length(all_results) + 1]] = res
write.table(t(data.frame(res)), col.names=F, quote=F, sep=":\t", file=outfile, append=T)

# RE -> EA - age + sex (weighted)
colnames(data)
xmod = lm(ave_MSE_clean ~ myopia_AS_dosage + sex_genetic + age, weights= ipw, data=data)
bx   = xmod$coef[2]
bxse = coeftest(xmod, vcov = sandwich)[2,2]
ymod = glm(eduyearsMore16 ~ myopia_AS_dosage + sex_genetic + age, weights= ipw, data=data, family="binomial")
by   = ymod$coef[2]
byse = coeftest(ymod, vcov = sandwich)[2,2]
beta_ratio = by/bx
N = nrow(data)
res = list()
res["g"] = "myopia_AS_dosage"
res["x"] = "ave_MSE_clean"
res["y"] = "eduyearsMore16"
res["covars"] = "sex_genetic + age"
res["w"] = "ipw"
res["N"] = N
res["bx"] = bx
res["by"] = by
res["ratio"] = beta_ratio
res["bxse"] = bxse
res["byse"] = byse
se_ratio_approx_1 = byse/bx
se_ratio_approx_2 = sqrt(byse^2/bx^2 + by^2*bxse^2/bx^4 - 2*1*by/bx^3)
res["se_ratio_approx_1"] = se_ratio_approx_1
res["se_ratio_approx_2"] = se_ratio_approx_2
f0 = by^2 - qt(0.975, N)^2 * byse^2
f1 = bx^2 - qt(0.975, N)^2 * bxse^2
f2 = by*bx
D = f2^2 - f0*f1
if(D>0) {
  r1 = (f2-sqrt(D))/f1
  r2 = (f2+sqrt(D))/f1
  if(f1>0) { res["Fiellers.Theorem.CI"] = paste0("Confidence interval is a closed interval [a,b]: a=", r1, ", b=", r2, sep="") }
  if(f1<0) { res["Fiellers.Theorem.CI"] = paste0("Confidence interval is the union of two open intervals (-Inf, a], [b, +Inf): a=", r2, ", b=", r1, sep="") }
}
if(D<0|D==0) { res["Fiellers.Theorem.CI"] = paste0("No finite confidence interval exists other than the entire real line.") }
cat("\n# EA -> RE - age + sex (not weighted)\n", file=outfile, append=T)
all_results[[length(all_results) + 1]] = res
write.table(t(data.frame(res)), col.names=F, quote=F, sep=":\t", file=outfile, append=T)

#
# Sensitivity analyses - Other covariates with standard MR
#

covList = list(c("genechip"),
               paste0("PC", 1:10),
               "birthweight",
               "breastfed",
               "northing",
               "easting",
               "tdi_log")

for (i in 1:length(covList)) {
  
  # RE -> EA - age + sex (not weighted)
  cat("\n# RE -> EA\n", file=outfile, append=T)
  res = run.MR(data, "myopia_AS_dosage", "ave_MSE_clean", "eduyears_clean", ccols=c("sex_genetic", "age", covList[[i]]))
  all_results[[length(all_results) + 1]] = res
  write.table(t(data.frame(res)), col.names=F, quote=F, sep=":\t", file=outfile, append=T)
  # RE -> EA - age + sex (weighted)
  cat("\n# RE -> EA\n", file=outfile, append=T)
  res = run.MR(data, "myopia_AS_dosage", "ave_MSE_clean", "eduyears_clean", wcol="ipw", ccols=c("sex_genetic", "age", covList[[i]]))
  all_results[[length(all_results) + 1]] = res
  write.table(t(data.frame(res)), col.names=F, quote=F, sep=":\t", file=outfile, append=T)
  
  # EA -> RE - age + sex (weighted)
  cat("\n# EA -> RE\n", file=outfile, append=T)
  res = run.MR(data, "ea_AS_dosage", "eduyears_clean", "ave_MSE_clean", ccols=c("sex_genetic", "age", covList[[i]]))
  all_results[[length(all_results) + 1]] = res
  write.table(t(data.frame(res)), col.names=F, quote=F, sep=":\t", file=outfile, append=T)
  # EA -> RE (not weighted)
  cat("\n# EA -> RE\n", file=outfile, append=T)
  res = run.MR(data, "ea_AS_dosage", "eduyears_clean", "ave_MSE_clean", wcol="ipw", ccols=c("sex_genetic", "age", covList[[i]]))
  all_results[[length(all_results) + 1]] = res
  write.table(t(data.frame(res)), col.names=F, quote=F, sep=":\t", file=outfile, append=T)

}



#
# Convert list to a dataframe #################################################
#

# Make df
cols = unique(unlist(sapply(all_results, function(x) {names(x)})))
resdf = do.call(rbind, lapply(all_results, function(x) { x[cols] }))
colnames(resdf) = cols

# Write df
outname = "output/MR_result_table.tsv"
write.table(resdf, file=outname, sep="\t", row.name=F)

