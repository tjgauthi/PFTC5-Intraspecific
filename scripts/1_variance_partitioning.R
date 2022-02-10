# Task DA.1 Perform variance analysis of traits

### Call source script----

source(here::here(path = "scripts/0_data_import.R"))

#load in packages
library(tidyverse)
library(lme4)

#check out how many leaves each individual has
leaf_count <- traits_wide %>% 
  group_by(site, plot_uid, individual_uid, taxon) %>% 
  summarize(n = n())

table(leaf_count$n)

#need to decide on how to deal with individuals with 1-2 leaves

#This code all comes from Julie Messier's web site
mod<-lmer(log(leaf_thickness_mm)~1+(1|taxon/individual_uid)+(1|site), data=traits_wide, na.action=na.omit)
variances<-c(unlist(lapply(VarCorr(mod),diag)), attr(VarCorr(mod),"sc")^2) #get variances

#this is the same as: print(VarCorr(mod),comp="Variance")

# raw variance

variances

# % Variance

(var.comp<-variances/sum(variances))
