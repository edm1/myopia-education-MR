#!/usr/bin/env Rscript
#

setwd("/mnt/seconddrive/data/phd/ukbiobank/myopia_EA_project/8_MRBase/using_split_sample/1_split_sample_weights/")

library("ggplot2")
library("rio")
library("tidyverse")

# Args
indata = "../../../2_derive_variables/output/phenotypes_alleleScores_170724.Rdata"
indosage = "../../../1_extract_instruments/output/EA_UKBB_Okbay.dosage"
split = "A"

# Load sample info
load(indata)

# Only keep data from correct split
data = data %>% filter(coreSampleGeno & !genoExcl & splitSample == split)
data = data.frame(data)
data$eid = as.character(data$eid)
# rownames(data) = data$eid

# Load dosage data
dosagedata = rio::import(indosage, format="tsv", header=T)

# Create data subset in same order as variants
sample_order = gsub("X", "", colnames(dosagedata)[8:ncol(dosagedata)])
data_sub = data[match(sample_order, data$eid), ]



#
# EA instruments ###########################################################
#

#
# Calc betas for EA ~ EA instruments
#

# Define regression function
run_regression = function(row) {
  dosages = as.numeric(row[8:length(row)])
  summ = summary(lm(data_sub$eduyears_clean ~ dosages +
               data_sub$sex_genetic +
               data_sub$age +
               data_sub$PC1 + data_sub$PC2 + data_sub$PC3 + data_sub$PC4 +
               data_sub$PC5 + data_sub$PC6 + data_sub$PC7 + data_sub$PC8 +
               data_sub$PC9 + data_sub$PC10,
               na.action=na.omit))
  data.frame(rsid=row["rsid"],
             b=summ$coefficients[2, 1],
             se=summ$coefficients[2, 2],
             p=summ$coefficients[2, 4],
             n=length(summ$residuals))
}

# Get results
results = as.data.frame(do.call(rbind, apply(dosagedata, 1, run_regression)))
rownames(results) = results[, 1]
results = results[, -c(1)]
colnames(results) = c("b_ea", "se_ea", "p_ea", "n_ea")
eaInstr_ea_results = results


#
# Calc betas for Myopia ~ EA instruments (logistic regression)
#

# Define regression function
run_regression = function(row) {
  dosages = as.numeric(row[8:length(row)])
  summ = summary(glm(data_sub$isMyopic ~ dosages +
                      data_sub$sex_genetic +
                      data_sub$age +
                      data_sub$PC1 + data_sub$PC2 + data_sub$PC3 + data_sub$PC4 +
                      data_sub$PC5 + data_sub$PC6 + data_sub$PC7 + data_sub$PC8 +
                      data_sub$PC9 + data_sub$PC10,
                    na.action=na.omit, family=binomial))
  data.frame(rsid=row["rsid"],
             b=summ$coefficients[2, 1],
             se=summ$coefficients[2, 2],
             p=summ$coefficients[2, 4],
             n=summ$df.residual)
}

# Get results
results = as.data.frame(do.call(rbind, apply(dosagedata, 1, run_regression)))
rownames(results) = results[, 1]
results = results[, -c(1)]
colnames(results) = c("b_myopia", "se_myopia", "p_myopia", "n_myopia")
eaInstr_myopia_results = results

#
# Calc betas for MSE ~ EA instruments
#

# Define regression function
run_regression = function(row) {
  dosages = as.numeric(row[8:length(row)])
  summ = summary(lm(data_sub$ave_MSE_clean ~ dosages +
                      data_sub$sex_genetic +
                      data_sub$age +
                      data_sub$PC1 + data_sub$PC2 + data_sub$PC3 + data_sub$PC4 +
                      data_sub$PC5 + data_sub$PC6 + data_sub$PC7 + data_sub$PC8 +
                      data_sub$PC9 + data_sub$PC10,
                    na.action=na.omit))
  data.frame(rsid=row["rsid"],
             b=summ$coefficients[2, 1],
             se=summ$coefficients[2, 2],
             p=summ$coefficients[2, 4],
             n=length(summ$residuals))
}

# Get results
results = as.data.frame(do.call(rbind, apply(dosagedata, 1, run_regression)))
rownames(results) = results[, 1]
results = results[, -c(1)]
colnames(results) = c("b_mse", "se_mse", "p_mse", "n_mse")
eaInstr_mse_results = results

#
# Calc betas for VA ~ EA instruments
#

# Define regression function
run_regression = function(row) {
  dosages = as.numeric(row[8:length(row)])
  summ = summary(lm(data_sub$ave_logmar_clean_log ~ dosages +
                      data_sub$sex_genetic +
                      data_sub$age +
                      data_sub$ave_MSE_clean +
                      data_sub$PC1 + data_sub$PC2 + data_sub$PC3 + data_sub$PC4 +
                      data_sub$PC5 + data_sub$PC6 + data_sub$PC7 + data_sub$PC8 +
                      data_sub$PC9 + data_sub$PC10,
                    na.action=na.omit))
  data.frame(rsid=row["rsid"],
             b=summ$coefficients[2, 1],
             se=summ$coefficients[2, 2],
             p=summ$coefficients[2, 4],
             n=length(summ$residuals))
}

# Get results
results = as.data.frame(do.call(rbind, apply(dosagedata, 1, run_regression)))
rownames(results) = results[, 1]
results = results[, -c(1)]
colnames(results) = c("b_va", "se_va", "p_va", "n_va")
eaInstr_va_results = results

#
# Merge and plot
#

# Split variant information
var_info = dosagedata[,1:7]
variants = dosagedata[,8:ncol(dosagedata)]
# Calc EAF
var_info$eaf = apply(variants, 1, function(x) {sum(x) / (2*length(x))})

# Combine
eaInstr_res = cbind(var_info,
                    eaInstr_myopia_results,
                    eaInstr_ea_results,
                    eaInstr_mse_results,
                    eaInstr_va_results)

write.table(eaInstr_res, file=paste0("output/eaInstr_stats_split", split, ".tsv"), sep="\t",
            quote=F, col.name=NA)

