#!/usr/bin/env Rscript
#

library("ggplot2")
library("rio")
library("MASS")
library("dplyr")

# Args
indata = "../2_derive_variables/output/phenotypes_alleleScores_170724.Rdata"
indosage = "../1_extract_instruments/output/Myopia_UKBB_Pickrell.dosage"

# Load
load(indata)
dosagedata = import(indosage, format="tsv", header=T)

# Make exclusions
data = data %>% filter(!genoExcl & coreSampleGeno)

# Create data subset in same order as variants
sample_order = gsub("X", "", colnames(dosagedata)[8:ncol(dosagedata)])
data_sub = data[match(sample_order, data$eid), ]

#
# Myopia instruments ###########################################################
#

#
# Calc betas for EA ~ myopia instruments
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
myopiaInstr_ea_results = results


#
# Calc betas for Myopia ~ myopia instruments (logistic regression)
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
myopiaInstr_myopia_results = results

#
# Calc betas for MSE ~ myopia instruments
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
myopiaInstr_mse_results = results


#
# Merge and plot
#

# Split variant information
var_info = dosagedata[,1:7]
variants = dosagedata[,8:ncol(dosagedata)]
# Calc EAF
var_info$eaf = apply(variants, 1, function(x) {sum(x) / (2*length(x))})

# Combine
myopiaInstr_res = cbind(var_info,
                    myopiaInstr_myopia_results,
                    myopiaInstr_ea_results,
                    myopiaInstr_mse_results)
write.table(myopiaInstr_res, file="output/myopiaInstr_stats.tsv", sep="\t",
            quote=F, col.name=NA)

# Plot myopia vs. EA
ggplot(myopiaInstr_res, aes(x=b_myopia, y=b_ea)) +
  geom_hline(yintercept=0, colour="orangered3", linetype="dotted") +
  geom_vline(xintercept=0, colour="orangered3", linetype="dotted") +
  geom_point() +
  geom_errorbar(aes(ymax=b_ea+1.96*se_ea, ymin=b_ea-1.96*se_ea), alpha=0.5) +
  geom_errorbarh(aes(xmax=b_myopia+1.96*se_myopia, xmin=b_myopia-1.96*se_myopia), alpha=0.5) +
  geom_smooth(method='rlm',formula=y~x, se=T, alpha=0.2, fullrange=T) +
  labs(title="Pickrell et al. Myopia instruments in UKBB",
       x="β for Myopia (<-0.75 SphE) in UKBB [w 95% CI]",
       y="β for Eduyears in UKBB [w 95% CI]")
ggsave(file="output/myopiaInstr_myopia_EA.png", height=10, width=12, dpi=300, units="cm")

# Plot MSE vs. EA
ggplot(myopiaInstr_res, aes(x=b_mse, y=b_ea)) +
  geom_hline(yintercept=0, colour="orangered3", linetype="dotted") +
  geom_vline(xintercept=0, colour="orangered3", linetype="dotted") +
  geom_point() +
  geom_errorbar(aes(ymax=b_ea+1.96*se_ea, ymin=b_ea-1.96*se_ea), alpha=0.5) +
  geom_errorbarh(aes(xmax=b_mse+1.96*se_mse, xmin=b_mse-1.96*se_mse), alpha=0.5) +
  geom_smooth(method='rlm',formula=y~x, se=T, alpha=0.2, fullrange=T) +
  labs(title="Pickrell et al. Myopia instruments in UKBB",
       x="β for Refractive Error (MSE) in UKBB [w 95% CI]",
       y="β for Eduyears in UKBB [w 95% CI]")
ggsave(file="output/myopiaInstr_RE_EA.png", height=10, width=12, dpi=300, units="cm")

# Plot myopia vs. MSE
ggplot(myopiaInstr_res, aes(x=b_myopia, y=b_mse)) +
  geom_hline(yintercept=0, colour="orangered3", linetype="dotted") +
  geom_vline(xintercept=0, colour="orangered3", linetype="dotted") +
  geom_point() +
  geom_errorbar(aes(ymax=b_mse+1.96*se_mse, ymin=b_mse-1.96*se_mse), alpha=0.5) +
  geom_errorbarh(aes(xmax=b_myopia+1.96*se_myopia, xmin=b_myopia-1.96*se_myopia), alpha=0.5) +
  geom_smooth(method='rlm',formula=y~x, se=T, alpha=0.2, fullrange=T) +
  labs(title="Pickrell et al. Myopia instruments in UKBB",
       x="β for Myopia (<-0.75 SphE) in UKBB [w 95% CI]",
       y="β for Refractive Error (MSE) in UKBB [w 95% CI]")
ggsave(file="output/myopiaInstr_myopia_RE.png", height=10, width=12, dpi=300, units="cm")
