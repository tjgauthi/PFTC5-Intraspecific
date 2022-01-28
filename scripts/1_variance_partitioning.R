# Task DA.1 Perform variance analysis of traits

### Call source script----

source(here::here(path = "scripts/0_data_import.R"))

#load in packages
library(tidyverse)
library(lme4)

#check out how many leaves each individual has
leaf_count <- traits %>% 
  group_by(site, plot_id, individual_nr, trait, taxon) %>% 
  summarize(n = n())

table(leaf_count$n)

#most individuals have only 1 measurement. This seems like it would make it difficult to assess intraindividual trait variation. Is that a goal of ours?

#test running a model
test_dat <- traits %>% filter(trait == "leaf_thickness_mm")


#This code all comes from Julie Messier's web site
mod<-lmer(log(value)~1+(1|taxon/individual_nr)+(1|site/plot_id), data=test_dat, na.action=na.omit)
variances<-c(unlist(lapply(VarCorr(mod),diag)), attr(VarCorr(mod),"sc")^2) #get variances

#this is the same as: print(VarCorr(mod),comp="Variance")

# raw variance

variances

# % Variance

(var.comp<-variances/sum(variances))
