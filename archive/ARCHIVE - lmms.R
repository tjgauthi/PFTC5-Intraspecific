#trait or traitsCV ~scale(elevation) * taxon + (1|site/plotID)

####Generalized Linear Mixed Models with the value of the trait as responsable variable
library(lme4)

#Provides p-values in type I, II or III anova and summary tables for lmer model fits.
#Call the library once you have choosen the best model to get the summary and then detach the library. 
#library("lmerTest") 

height <- traits %>%
  filter(trait == "plant_height_cm")

#For each of the traits I run a model with and without the log transformation, then I compared the models using AIC
#The AIC of the transformed variable was adjunted following Akaike (1978) in order to enable model comparison

#Plant height
mod_height1 <- lmer(value ~ scale(elevation) * taxon + (1|site/plot_id), data = height,na.action=na.omit)
mod_height2 <- lmer(log(value) ~ scale(elevation) * taxon + (1|site/plot_id), data = height,na.action=na.omit)#Best model

AIC(mod_height1,mod_height2)

AIC(mod_height1,mod_height2) + matrix(ncol=2, c(0,0,0, sum(2*log(height$value)))) #Akaike, H. 1978. "On the likelihood of a time series model," Journal of the Royal Statistical Society, Series D (The Statistician), 27(3/4), pp. 217–235.

summary(mod_height2)

#Dry mass
dry_mass <- traits %>%
  filter(trait == "dry_mass_g")

mod_dry_mass1 <- lmer(value ~ scale(elevation) * taxon + (1|site/plot_id), data = dry_mass,na.action=na.omit)
mod_dry_mass2 <- lmer(log(value) ~ scale(elevation) * taxon + (1|site/plot_id), data = dry_mass,na.action=na.omit)#Best

AIC(mod_dry_mass1,mod_dry_mass2)

AIC(mod_dry_mass1,mod_dry_mass2) + matrix(ncol=2, c(0,0,0, sum(2*log(dry_mass$value)))) 

summary(mod_dry_mass2)

#Leaf area
area <- traits %>%
  filter(trait == "leaf_area_cm2")

mod_area1 <- lmer(value ~ scale(elevation) * taxon + (1|site/plot_id), data = area,na.action=na.omit)
mod_area2 <- lmer(log(value) ~ scale(elevation) * taxon + (1|site/plot_id), data = area,na.action=na.omit)#Best

AIC(mod_area1,mod_area2)

AIC(mod_area1,mod_area2) + matrix(ncol=2, c(0,0,0, sum(2*log(area$value)))) 

summary(mod_area2)

#Specific Leaf Area
sla <- traits %>%
  filter(trait == "sla_cm2_g")

mod_sla1 <- lmer(value ~ scale(elevation) * taxon + (1|site/plot_id), data = sla, na.action=na.omit)##Singular
mod_sla1.1 <- lmer(value ~ scale(elevation) * taxon + (1|site) + (1|plot_id), data = sla, na.action=na.omit)##Singular
mod_sla1.2 <- lmer(value ~ scale(elevation) * taxon + (1|site), data = sla, na.action=na.omit)##Singular
mod_sla1.2 <- lmer(value ~ scale(elevation) * taxon + (1|plot_id), data = sla, na.action=na.omit)

#This model has the structure that we want
mod_sla2 <- lmer(log(value) ~ scale(elevation) * taxon + (1|site/plot_id), data = sla,na.action=na.omit)##Best 

summary(mod_sla2)

#Leaf Dry Matter Content (LDMC)
ldmc <- traits %>%
  filter(trait == "ldmc")

mod_ldmc1 <- lmer(value ~ scale(elevation) * taxon + (1|site/plot_id), data = ldmc,na.action=na.omit)##Singular
mod_ldmc2 <- lmer(log(value) ~ scale(elevation) * taxon + (1|site/plot_id), data = ldmc,na.action=na.omit)##Singular

mod_ldmc1.1 <- lmer(value ~ scale(elevation) * taxon + (1|site) + (1|plot_id), data = ldmc,na.action=na.omit)##Singular
mod_ldmc2.1 <- lmer(log(value) ~ scale(elevation) * taxon + (1|site) + (1|plot_id), data = ldmc,na.action=na.omit)##Singular

mod_ldmc1.2 <- lmer(value ~ scale(elevation) * taxon + (1|site), data = ldmc,na.action=na.omit)##Singular
mod_ldmc2.2 <- lmer(log(value) ~ scale(elevation) * taxon + (1|site), data = ldmc,na.action=na.omit)##Singular

mod_ldmc1.3 <- lmer(value ~ scale(elevation) * taxon + (1|plot_id), data = ldmc,na.action=na.omit)#Best
mod_ldmc2.3 <- lmer(log(value) ~ scale(elevation) * taxon + (1|plot_id), data = ldmc,na.action=na.omit)


AIC(mod_ldmc1.3,mod_ldmc2.3)

AIC(mod_ldmc1.3,mod_ldmc2.3) + matrix(ncol=2, c(0,0,0, sum(2*log(ldmc$value)))) 

summary(mod_ldmc1.3)
summary(mod_ldmc2.3)


#Leaf thickness

thickness <- traits %>%
  filter(trait == "leaf_thickness_mm")

mod_thickness1 <- lmer(value ~ scale(elevation) * taxon + (1|site/plot_id), data = thickness,na.action=na.omit)
mod_thickness2 <- lmer(log(value) ~ scale(elevation) * taxon + (1|site/plot_id), data = thickness,na.action=na.omit)#Best

AIC(mod_thickness1,mod_thickness2)

AIC(mod_thickness1,mod_thickness2) + matrix(ncol=2, c(0,0,0, sum(2*log(thickness$value)))) 

summary(mod_thickness2)

detach("package:lmerTest")

#to test the overall effect of the interactions
anova(mod_height2)
anova(mod_dry_mass2)
anova(mod_area2)
anova(mod_sla2)
anova(mod_ldmc1.3)
anova(mod_ldmc2.3)
anova(mod_thickness2)


#Graphs
library(ggplot2)
library(effects)


ggplot(data.frame(Effect(c("elevation","taxon"),mod_height2)),
       aes(x=elevation,y=fit,color=taxon,group=taxon))+
  geom_line() + 
  ylab("log height")

ggplot(data.frame(Effect(c("elevation","taxon"),mod_dry_mass2)),
       aes(x=elevation,y=fit,color=taxon,group=taxon))+
  geom_line() + 
  ylab("log dry mass")

ggplot(data.frame(Effect(c("elevation","taxon"),mod_area2)),
       aes(x=elevation,y=fit,color=taxon,group=taxon))+
  geom_line() + 
  ylab("log leaf area")

ggplot(data.frame(Effect(c("elevation","taxon"),mod_sla2)),
       aes(x=elevation,y=fit,color=taxon,group=taxon))+
  geom_line() + 
  ylab("log sla")

ggplot(data.frame(Effect(c("elevation","taxon"),mod_ldmc2.3)),
       aes(x=elevation,y=fit,color=taxon,group=taxon))+
  geom_line() + 
  ylab("log ldmc")

ggplot(data.frame(Effect(c("elevation","taxon"),mod_thickness2)),
       aes(x=elevation,y=fit,color=taxon,group=taxon))+
  geom_line() + 
  ylab("log thickness")




#####
# Function for computing coefficient of variation
cv = function(x) { return(sd(x)/mean(x)) } #Josef function

#coef_var <- function(x, na.rm = FALSE) {
  #sd(x, na.rm=na.rm) / mean(x, na.rm=na.rm)
#}


# Compute CV for each trait for each individual
trait_cvs = traits %>% 
  group_by(site, plot_id, individual_nr, functional_group, family, taxon, trait, elevation) %>% 
  summarize(cv_trait = coef_var(value))

#GLMMs with the CV_value of the trait as responsable variable

#CV_Plant Height
height <- trait_cvs %>%
  filter(trait == "plant_height_cm")

mod_height1 <- lmer(cv_trait ~ scale(elevation) * taxon + (1|site/plot_id), data = height,na.action=na.omit)#Singular
mod_height1.1 <- lmer(cv_trait ~ scale(elevation) * taxon + (1|site) + (1|plot_id), data = height,na.action=na.omit)#Best
mod_height2 <- lmer(log(cv_trait) ~ scale(elevation) * taxon + (1|site/plot_id), data = height,na.action=na.omit)#NA/NaN/Inf in 'y'
mod_height2.1 <- lmer(log(cv_trait) ~ scale(elevation) * taxon + (1|site) + (1|plot_id), data = height,na.action=na.omit)#NA/NaN/Inf in 'y'

summary(mod_height1.1)

#CV_dry mass
dry_mass <- trait_cvs %>%
  filter(trait == "dry_mass_g")

mod_dry_mass1 <- lmer(cv_trait ~ scale(elevation) * taxon + (1|site/plot_id), data = dry_mass,na.action=na.omit)#Best
mod_dry_mass2 <- lmer(log(cv_trait) ~ scale(elevation) * taxon + (1|site/plot_id), data = dry_mass,na.action=na.omit)

AIC(mod_dry_mass1,mod_dry_mass2)
AIC(mod_dry_mass1,mod_dry_mass2) + matrix(ncol=2, c(0,0,0, sum(2*log(dry_mass$cv_trait)))) 

summary(mod_dry_mass1)

#CV_leaf area
area <- trait_cvs %>%
  filter(trait == "leaf_area_cm2")

mod_area1 <- lmer(cv_trait ~ scale(elevation) * taxon + (1|site/plot_id), data = area,na.action=na.omit)#Singular
mod_area1.1 <- lmer(cv_trait ~ scale(elevation) * taxon + (1|site) + (1|plot_id), data = area,na.action=na.omit)#Best

#mod_area2 <- lmer(log(cv_trait) ~ scale(elevation) * taxon + (1|site/plot_id), data = area,na.action=na.omit)
#AIC(mod_area1,mod_area2)
#AIC(mod_area1,mod_area2) + matrix(ncol=2, c(0,0,0, sum(2*log(area$value)))) 

summary(mod_area1.1)

#CV_Specific Leaf Area
sla <- trait_cvs %>%
  filter(trait == "sla_cm2_g")

mod_sla1 <- lmer(cv_trait ~ scale(elevation) * taxon + (1|site/plot_id), data = sla, na.action=na.omit)##Singular
mod_sla1.1 <- lmer(cv_trait ~ scale(elevation) * taxon + (1|site) + (1|plot_id), data = sla, na.action=na.omit)##Singular
mod_sla1.2 <- lmer(cv_trait ~ scale(elevation) * taxon + (1|site), data = sla, na.action=na.omit)#Best
mod_sla1.3 <- lmer(cv_trait ~ scale(elevation) * taxon + (1|plot_id), data = sla, na.action=na.omit)##Singular
#mod_sla2 <- lmer(log(value) ~ scale(elevation) * taxon + (1|site/plot_id), data = sla,na.action=na.omit)

summary(mod_sla1.2)

#CV_Leaf Dry Matter Content (LDMC)
ldmc <- trait_cvs %>%
  filter(trait == "ldmc")

mod_ldmc1 <- lmer(cv_trait ~ scale(elevation) * taxon + (1|site/plot_id), data = ldmc,na.action=na.omit)##Singular
mod_ldmc2 <- lmer(log(cv_trait) ~ scale(elevation) * taxon + (1|site/plot_id), data = ldmc,na.action=na.omit)##Singular

mod_ldmc1.1 <- lmer(cv_trait ~ scale(elevation) * taxon + (1|site) + (1|plot_id), data = ldmc,na.action=na.omit)##Singular
mod_ldmc2.2 <- lmer(log(cv_trait) ~ scale(elevation) * taxon + (1|site) + (1|plot_id), data = ldmc,na.action=na.omit)##Singular

mod_ldmc1.3 <- lmer(cv_trait ~ scale(elevation) * taxon + (1|site), data = ldmc,na.action=na.omit)##Singular
mod_ldmc2.3 <- lmer(log(cv_trait) ~ scale(elevation) * taxon + (1|site), data = ldmc,na.action=na.omit)#Best

mod_ldmc1.4 <- lmer(cv_trait ~ scale(elevation) * taxon + (1|plot_id), data = ldmc,na.action=na.omit)#Singular
mod_ldmc2.4 <- lmer(log(cv_trait) ~ scale(elevation) * taxon + (1|plot_id), data = ldmc,na.action=na.omit)#Singular

summary(mod_ldmc2.3)

#CV_leaf thickness

thickness <- trait_cvs %>%
  filter(trait == "leaf_thickness_mm")

mod_thickness1 <- lmer(cv_trait ~ scale(elevation) * taxon + (1|site/plot_id), data = thickness,na.action=na.omit)#Singular
mod_thickness1.1 <- lmer(cv_trait ~ scale(elevation) * taxon + (1|site)+(1|plot_id), data = thickness,na.action=na.omit)#Singular
mod_thickness1.2 <- lmer(cv_trait ~ scale(elevation) * taxon + (1|plot_id), data = thickness,na.action=na.omit)#Singular
mod_thickness1.3 <- lmer(cv_trait ~ scale(elevation) * taxon + (1|site), data = thickness,na.action=na.omit)#Best

summary(mod_thickness1.3)

mod_thickness2 <- lmer(log(cv_trait) ~ scale(elevation) * taxon + (1|site/plot_id), data = thickness,na.action=na.omit)#Singular
mod_thickness2.1 <- lmer(log(cv_trait) ~ scale(elevation) * taxon + (1|site) + (1|plot_id), data = thickness,na.action=na.omit)#Singular
mod_thickness2.2 <- lmer(log(cv_trait) ~ scale(elevation) * taxon + (1|plot_id), data = thickness,na.action=na.omit)#Singular
mod_thickness2.3 <- lmer(log(cv_trait) ~ scale(elevation) * taxon + (1|site), data = thickness,na.action=na.omit)

#AIC(mod_thickness1.3, mod_thickness2.3)
#AIC(mod_thickness1.3,mod_thickness2.3) + matrix(ncol=2, c(0,0,0, sum(2*log(thickness$value)))) 


detach("package:lmerTest")

anova(mod_height1.1)
anova(mod_dry_mass1)
anova(mod_area1.1)
anova(mod_sla1.2)
anova(mod_ldmc2.3)
anova(mod_thickness1.3)


#Graphs
library(ggplot2)
library(effects)


ggplot(data.frame(Effect(c("elevation","taxon"),mod_height1.1)),
       aes(x=elevation,y=fit,color=taxon,group=taxon))+
  geom_line() + 
  ylab("CV height")

ggplot(data.frame(Effect(c("elevation","taxon"),mod_dry_mass1)),
       aes(x=elevation,y=fit,color=taxon,group=taxon))+
  geom_line() + 
  ylab("CV dry mass")

ggplot(data.frame(Effect(c("elevation","taxon"),mod_area1.1)),
       aes(x=elevation,y=fit,color=taxon,group=taxon))+
  geom_line() + 
  ylab("CV leaf area")

ggplot(data.frame(Effect(c("elevation","taxon"),mod_sla1.2)),
       aes(x=elevation,y=fit,color=taxon,group=taxon))+
  geom_line() + 
  ylab("CV sla")

ggplot(data.frame(Effect(c("elevation","taxon"),mod_ldmc2.3)),
       aes(x=elevation,y=fit,color=taxon,group=taxon))+
  geom_line() + 
  ylab("CV ldmc")

ggplot(data.frame(Effect(c("elevation","taxon"),mod_thickness1.3)),
       aes(x=elevation,y=fit,color=taxon,group=taxon))+
  geom_line() + 
  ylab("CV thickness")





  