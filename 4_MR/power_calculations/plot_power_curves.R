#!/usr/bin/env Rscript
#
# Plot power curves calculated using http://cnsgenomics.com/shiny/mRnd/
#
# Citation: Calculating statistical power in Mendelian randomization studies Marie-Jo A Brion, Konstantin Shakhbazov, Peter M Visscher International Journal of Epidemiology 2013 42: 1497-1501
#

library("rio")
library("tidyverse")

setwd("/mnt/seconddrive/data/phd/ukbiobank/myopia_EA_project/4_MR/power_calculations/")

#
# RE -> EA
#

# Load
df = rio::import("re_ea_power_data.tsv")
colnames(df) = c("BetaIV", "n23612", "n90000", "n69798")
# Plot
df %>% filter(BetaIV < 0.10) %>% ggplot(aes(x=BetaIV, y=n69798)) +
  geom_point() +
  geom_line(size=1) +
  geom_hline(yintercept=0.8, alpha=0.2, linetype="dashed") +
  geom_vline(xintercept=0.048, alpha=0.2, linetype="dashed") +
  ylab("Power (α=0.05)") +
  xlab("Causal association (years/D)") +
  ggtitle("Refractive error on education")
  #ggtitle("Power to detect causal effect of refractive\nerror on education")
ggsave(file="RE_EA_power_curve.png", units="cm", w=10, h=8, dpi=300)

#
# RE -> EA
#

# Load
df = rio::import("ea_re_power_data.tsv")
colnames(df) = c("BetaIV", "n23612", "n90000", "n69798")
# Plot
df %>% filter(BetaIV < 0.3) %>%
  ggplot(aes(x=BetaIV, y=n69798)) +
  geom_point() +
  geom_line(size=1) +
  geom_hline(yintercept=0.8, alpha=0.2, linetype="dashed") +
  geom_vline(xintercept=0.14, alpha=0.2, linetype="dashed") +
  ylab("Power (α=0.05)") +
  xlab("Causal association (D/year)") +
  scale_color_hue(labels = c("N=23612", "N=69798")) +
  ggtitle("Education on refractive error")
  #ggtitle("Power to detect causal effect of education\non refractive error")
ggsave(file="EA_RE_power_curve.png", units="cm", w=10, h=8, dpi=300)
