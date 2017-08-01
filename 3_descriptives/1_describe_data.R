#!/usr/bin/env Rscript
#

library("ggplot2")
library("reshape2")
library("gmodels")
library("dplyr")
library("magrittr")
library("sandwich")
library("lmtest")

setwd("/mnt/seconddrive/data/phd/ukbiobank/myopia_EA_project/3_descriptives/")
#setwd("/Users/ed/Downloads/temp_ukbb/3_descriptives/")

# Args
indata = "../2_derive_variables/output/phenotypes_alleleScores_170728.Rdata"

# Load
load(indata)
data_full = data

# Keep actual sample
#data = data %>% filter(coreSample)
data = data %>% filter(coreSample & coreSampleGeno & !genoExcl)


# Num samples per variable
apply(data, 2, function(x) sum(!is.na(x)))

# Count total MR sample
data_full %>% filter(coreSample) %>% nrow
data_full %>% filter(coreSample & coreSampleGeno) %>% nrow
data_full %>% filter(coreSample & coreSampleGeno & !genoExcl) %>% nrow

prop.table((table(data$eduyears_clean)))

#
# Basic descriptives
#

table(data$splitSample)

# Make myopic, emmotropic, hyperopic variable
data$MSE_cat = NA
data$MSE_cat[data$ave_MSE_clean <= -0.75] = "myopic"
data$MSE_cat[data$ave_MSE_clean >= 0.75] = "hyperopic"
data$MSE_cat[data$ave_MSE_clean > -0.75 & data$ave_MSE_clean < 0.75] = "emmotropic"
data$MSE_cat = factor(data$MSE_cat, levels=c("myopic", "emmotropic", "hyperopic"))

# Summarise data by myopia, emmetropia and hypermetropia
data_melt = data %>% melt(id.vars=c("eid", "MSE_cat")) 
desc_cat = data_melt %>%
              group_by(variable, MSE_cat) %>%
              summarise(mean = round(mean(as.numeric(value), na.rm=T), 3),
                        n = n(),
                        sd = round(sd(value, na.rm=T), 3),
                        lowerCI = mean - 1.96 * sd / sqrt(n),
                        upperCI = mean + 1.96 * sd / sqrt(n),
                        freq = table(value)[1] / n)
# Summarise by all
desc_all = data_melt %>%
  group_by(variable) %>%
  summarise(mean = round(mean(as.numeric(value), na.rm=T), 3),
              n = n(),
              sd = round(sd(value, na.rm=T), 3),
              lowerCI = mean - 1.96 * sd / sqrt(n),
              upperCI = mean + 1.96 * sd / sqrt(n),
              freq = table(value)[1] / n)
# Combine category and all
desc_all$MSE_cat = "all"
desc_all = desc_all[, colnames(desc_cat)]
desc = bind_rows(desc_cat, desc_all)
# desc = desc[order(desc$variable), ]

write.table(desc, file="output/descriptives_v1.tsv", sep="\t", quote=F, row.names=F)

# SD using 23,612 sample
sd(data$ave_MSE_clean)**2
sd(data$eduyears_clean)**2

#
# Observational associations
#



# ModelA with age and sex covars
colnames(data)
modA = lm(ave_MSE_clean ~ eduyears_clean + age + sex, data=data)
coeftest(modA, vcov = sandwich)
coefci(modA, vcov = sandwich)
length(modA$residuals)

# ModelB with age, sex, tdi, ...
colnames(data)
modB = lm(ave_MSE_clean ~ eduyears_clean + age + sex + tdi_log + birthweight + breastfed + northing + easting, data=data)
coeftest(modB, vcov = sandwich)
coefci(modB, vcov = sandwich)
length(modB$residuals)

# refractive error on education with age and sex covars
colnames(data)
modC = lm(eduyears_clean ~ ave_MSE_clean + age + sex, data=data)
coeftest(modC, vcov = sandwich)
coefci(modC, vcov = sandwich)
length(modC$residuals)

# refractive error on education with age, sex, tdi, ...
colnames(data)
modD = lm(eduyears_clean ~ ave_MSE_clean + age + sex + tdi_log + birthweight + breastfed + northing + easting, data=data)
coeftest(modD, vcov = sandwich)
coefci(modD, vcov = sandwich)
length(modD$residuals)

coeftest(lm(eduyears_clean ~ ave_MSE_clean, data=data), vcov = sandwich)
coeftest(lm(ave_MSE_clean ~ eduyears_clean, data=data), vcov = sandwich)

#
# Direct genetic associations
#

library("asbio")

# Remove genetic exclusions
data_clean = data %>% filter(!genoExcl)

# Education allele score with education - Partial R2
lm.without = lm(eduyears_clean ~ age + sex, data=data_clean)
lm.with = lm(eduyears_clean ~ age + sex + ea_AS_dosage, data=data_clean)
partial.R2(lm.without, lm.with) * 100
anova(lm.without, lm.with)
coeftest(lm.with, vcov = sandwich)
coefci(lm.with, vcov = sandwich)
length(lm.with$residuals)

# Education allele score with RE - Partial R2
lm.without = lm(ave_MSE_clean ~ age + sex, data=data_clean)
lm.with = lm(ave_MSE_clean ~ age + sex + ea_AS_dosage, data=data_clean)
partial.R2(lm.without, lm.with) * 100
coeftest(lm.with, vcov = sandwich)
coefci(lm.with, vcov = sandwich)
length(lm.with$residuals)

# Myopia allele score with education - Partial R2
lm.without = lm(eduyears_clean ~ age + sex, data=data_clean)
lm.with = lm(eduyears_clean ~ age + sex + myopia_AS_dosage, data=data_clean)
partial.R2(lm.without, lm.with) * 100
coeftest(lm.with, vcov = sandwich)
coefci(lm.with, vcov = sandwich)
length(lm.with$residuals)

# Myopia allele score with RE - Partial R2
lm.without = lm(ave_MSE_clean ~ age + sex, data=data_clean)
lm.with = lm(ave_MSE_clean ~ age + sex + myopia_AS_dosage, data=data_clean)
partial.R2(lm.without, lm.with) * 100
anova(lm.without, lm.with)
coeftest(lm.with, vcov = sandwich)
coefci(lm.with, vcov = sandwich)
length(lm.with$residuals)

#
# Birthweight and breastfed missing at random?
#

# Birthweight
summary(lm(ave_MSE_clean ~ is.na(data$birthweight), data=data))
summary(lm(ave_MSE_clean ~ birthweight, data=data))
# Breastfed
summary(lm(ave_MSE_clean ~ is.na(data$breastfed), data=data))
summary(lm(ave_MSE_clean ~ breastfed, data=data))
# TDI
summary(lm(ave_MSE_clean ~ is.na(data$tdi_log), data=data))
summary(lm(ave_MSE_clean ~ breastfed, data=data))



#
# Plot refractive error vs educational attainment #############################
#

# Plot mean and CI education against RE
data %>% 
  group_by(eduyears_clean) %>% 
  summarise(mean=mean(ave_MSE_clean),
            se=sd(ave_MSE_clean) / sqrt(n()),
            lowerCI=mean-1.96*se,
            upperCI=mean+1.96*se) %>%
  ggplot(aes(x=eduyears_clean, y=mean)) +
    geom_point() +
    geom_errorbar(aes(ymin=lowerCI, ymax=upperCI), width=.1) +
    xlab("Educational attainment (years)") +
    ylab("Refractive error (D) ± 95%CI") +
    ggtitle("Observational relationship between\neducation and refractive error") +
    scale_x_continuous(breaks=seq(15, 21, 1))
ggsave(file="output/observational_edu_RE.png", h=4, w=5, dpi=150)

data %>% ggplot(aes(x=ave_MSE_clean)) + geom_histogram()
data %>% summarise(mean=mean(ave_MSE_clean))

#
# Test - Num of samples with eye questions #####################################
#
colnames(data)
sum(complete.cases(data[, "q_glassesReason_myopia", drop=F]))
sum(complete.cases(data[, c("q_glassesReason_myopia", "sex_genetic"), drop=F]))
sum(complete.cases(data[, c("ave_MSE_clean", "q_glassesReason_myopia", "sex_genetic"), drop=F]))

# Compare age first wearing glasses to RE
summary(lm(ave_MSE_clean ~ q_glassesAge + q_glassesReason_myopia + q_glassesReason_hypermetropia, data=data))
ggplot(data, aes(x=q_glassesAge, y=ave_MSE_clean)) +
  geom_point(alpha=0.1) +
  geom_smooth(method=lm)

#
# What questions correlate best with visual acuity #############################
#

colnames(data)
cols = c("q_glassesReason_hypermetropia",
         "q_glassesReason_myopia",
         "q_glassesReason_presbyopia",
         "q_glassesReason_astigmatism",
          "q_glassesReason_strabismusm",
          "q_glassesReason_amblyopia",
          "q_wearsGlasses",
          "q_glassesAge")


for (col in cols) {
  form = paste0("worst_logmar_clean_log ~ ", col)
  summ = summary(lm(form, data=data))
  b = summ$coefficients[2,1]
  p = summ$coefficients[2,4]
  print(paste(col, b, p))
}


summary(lm(paste0("worst_logmar_clean_log ~ ", paste0(cols, collapse=" + ")), data=data))



# Plot phenotypes ##############################################################
#

# Test regress
colnames(data)
datacc = data[, c("ave_MSE_clean", "eduyears_clean")]
datacc = datacc[complete.cases(datacc), ]
x = scale(datacc$eduyears_clean)
xres = as.numeric(lm(eduyears_clean ~ ave_MSE_clean, data=datacc)$residuals)
mean(x)
mean(xres)
sd(x)
sd(xres)
scale

#
# Compare phenotypes
#

# Plot histograms of numeric data
datam = melt(data[,sapply(data, class) %in% c("numeric", "integer")])
ggplot(datam, aes(x=value)) + 
  facet_wrap(~ variable, scales="free_x") + 
  geom_histogram()
ggsave(file="output/all_histograms.png", height=12, width=12, dpi=150)

# histogram of refractiver error
colnames(data)
ggplot(data, aes(x=ave_MSE_clean)) +
  geom_histogram() +
  labs(x="Refractive error (average MSE)")
ggsave(file="output/mse_histogram.png", height=8, width=10, units="cm")

# histogram of eduyears
colnames(data)
ggplot(data, aes(x=eduyears_clean)) +
  geom_histogram() +
  labs(x="Eduyears")
ggsave(file="output/eduyears_histogram.png", height=8, width=10, units="cm")

# histogram of eduyears okbay
colnames(data)
ggplot(data, aes(x=eduyearsOkbay_clean)) +
  geom_histogram() +
  labs(x="Eduyears (Okbay)")
ggsave(file="output/eduyears_okbay_histogram.png", height=8, width=10, units="cm")

# Compare Refractive error to eduyears
model = lm(eduyears_clean ~ ave_MSE_clean, data=data)
summary(model)
confint(model)

# Compare eduyears error to refractive error
model = lm(ave_MSE_clean ~ eduyears_clean, data=data)
summary(model)
confint(model)

#
# Compare allele scores to phenotypes ##########################################
#

colnames(data)

#
# RE allele score vs. refractive error
#

ggplot(data, aes(x=re_AS_dosage, y=ave_MSE_clean)) +
  geom_point(alpha=0.1) +
  geom_smooth(method=lm)

summary(lm(ave_MSE_clean ~ re_AS_dosage, data=data))

#
# VA allele score vs. va
#

ggplot(data, aes(x=va_AS_dosage, y=ave_logmar_clean_log)) +
  geom_point(alpha=0.1) +
  geom_smooth(method=lm)

summary(lm(ave_logmar_clean_log ~ va_AS_dosage, data=data))

#
# Myopia allele score vs. refractive error
#

ggplot(data, aes(x=myopia_AS_dosage, y=ave_MSE_clean)) +
  geom_point(alpha=0.03) +
  geom_smooth(method=lm) + 
  labs(x="Myopia allele score", y="Refractive error (Mean Spherical Eqiv)")
ggsave(file="output/myopiaAlleleScore_v_RE.png", h=10, w=10, unit="cm")
model = lm(ave_MSE_clean ~ myopia_AS_dosage, data=data)
summ = summary(model)
confint(model, "myopia_AS_dosage")
length(sum$residuals)

#
# Myopia allele score vs. myopia (<-0.75 SphE)
#

data$isMyopic = factor(data$isMyopic)
table(data$isMyopic)
ggplot(subset(data, !is.na(isMyopic)), aes(x=isMyopic, y=myopia_AS_dosage)) +
  geom_boxplot() + 
  labs(x="Is myopic (<-0.75 Sph)", y="Myopia allele score")
ggsave(file="output/myopiaAlleleScore_v_myopia.png", h=10, w=10, unit="cm")

#summary(glm(isMyopic ~ myopia_AS_dosage, data=data, family=binomial))
sum = summary(lm(myopia_AS_dosage ~ isMyopic, data=data))
length(sum$residuals)

#
# EA allele score vs. EA
#

ggplot(data, aes(x=ea_AS_dosage, y=eduyears_clean)) +
  geom_point(alpha=0.03) +
  geom_smooth(method=lm) + 
  labs(x="Eduyears allele score", y="Eduyears in UKBB")
ggsave(file="output/eduyearsAlleleScore_v_eduyears.png", h=10, w=10, unit="cm")

model = lm(eduyears_clean ~ ea_AS_dosage, data=data)
sum = summary(model)
confint(model)
length(sum$residuals)

#
# EA allele score vs. refractive error
#

ggplot(data, aes(x=ea_AS_dosage, y=ave_MSE_clean)) +
  geom_point(alpha=0.1) +
  geom_smooth(method=lm)

summary(lm(ave_MSE_clean ~ ea_AS_dosage + tdi_log + sex_genetic + age, data=data))

#
# Compare allele scores to potential confounders ###############################
#

# Define regression function
do_regression = function(x, y, data) {

  # Do regression
  mod = lm(data[, y] ~ data[, x])
  summ = summary(mod)
  conf = confint(mod)
  # Extract results
  results = list()
  results["x"] = x
  results["y"] = y
  results["beta"] = summ$coefficients[2,1]
  results["se"] = summ$coefficients[2,2]
  results["beta_lower"] = conf[2,1]
  results["beta_upper"] = conf[2,2]
  results["pval"] = summ$coefficients[2,4]
  results["n"] = length(summ$residuals)
  return(results)

}

# Run for each covar and pheno
x_list = c("birthweight", "breastfed", "northing", "easting", "tdi_log", "age", "sex_genetic", paste0("PC", 1:10))
y_list = c("myopia_AS_dosage", "ave_MSE_clean", "ea_AS_dosage", "eduyears_clean")
# Run loop
outresults = list()
for (x in x_list) {
  for (y in y_list) {
    outresults[[paste0(x, "-", y)]] = do_regression(x, y, data)
  }
}
# Make res data frame table
resdf = do.call(rbind, lapply(outresults, function(x) { x[] }))
# Write
write.table(resdf, file="output/confounders_regression.tsv", sep="\t", row.names=F, quote=F)
