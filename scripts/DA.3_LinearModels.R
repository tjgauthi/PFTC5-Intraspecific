#trait or traitsCV ~scale(elevation) * taxon + (1|site/individual)
#source(here::here(path = "scripts/0_data_import.R"))

#Generalized Linear Mixed Models with the value of the trait as responsible variable
#

#### Call source script----

source(here::here(path = "scripts/0_data_import.R"))


####SETUP####
library(plyr)
library(Rmisc) #for summarySE function
library(dplyr)
library(lme4)
library(ggplot2)
#library(effects)
library(gridExtra)
library(lmerTest)
library(ggpubr)
library(grid)
library(cowplot)
library(MuMIn)#for r2
library(utils)#load csv
library(tidyr)

# Plotting theme
my_theme = theme_bw() + 
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        text = element_text (size = 10),
        axis.text = element_text(size = 10))

# Color palette
pal_lm <- c("#016392", "#A0CBE8", "#E19825", "#F7C480", "#3E8853", "#9FCD99")
names(pal_lm) <- c("Halenia umbellata", 
                   "Lachemilla orbiculata", 
                   "Paspalum bonplandianum", 
                   "Rhynchospora macrochaeta", 
                   "Gaultheria glomerata", 
                   "Vaccinium floribundum")

####Data Organization####

# Set up data frame for fitting models to log-transformed data
traits_wide_lm = traits_wide %>% 
  mutate(ln_height = log(plant_height_cm),
         ln_dry_mass = log(dry_mass_g),
         ln_area = log(leaf_area_cm2),
         ln_sla = log(sla_cm2_g),
         ln_ldmc = log(ldmc),
         ln_thickness = log(leaf_thickness_mm))

#creating a summarized dataset for plant height
plant_height <- summarySE (data=traits_wide_lm,
                      measurevar = "ln_height",
                      groupvars = c("site", "taxon", "elevation","individual_uid"))

#### Plot 1 - variation in trait values with elevation####

# Fit models to each trait
mod_dry_mass1 <- lmer(ln_dry_mass ~ scale(elevation) * taxon + (1|site/individual_uid), 
                    data = traits_wide_lm,
                    na.action=na.omit)
mod_ldmc1 <- lmer(ln_ldmc ~ scale(elevation) * taxon + (1|site/individual_uid), 
                  data = traits_wide_lm,
                  na.action=na.omit)
mod_area1 <- lmer(ln_area ~ scale(elevation) * taxon + (1|site/individual_uid), 
                      data = traits_wide_lm,
                      na.action=na.omit)
mod_sla1 <- lmer(ln_sla ~ scale(elevation) * taxon + (1|site/individual_uid), 
                  data = traits_wide_lm,
                  na.action=na.omit)
mod_thickness1 <- lmer(ln_thickness ~ scale(elevation) * taxon + (1|site/individual_uid), 
                       data = traits_wide_lm,
                       na.action=na.omit)
mod_height1 <- lmer(ln_height ~ scale(elevation) * taxon + (1|site), #no individual_uid because there is only one height per individual
                    data = plant_height,
                    na.action=na.omit)
# Singular fit LDMC
# Removing the plot_id nested effect does not fix the issue
# It can be resolved by removing random effects altogether
# Not sure if this is a substantial issue since this is going in the sup

#Anova on each regression for P and F values
anova(mod_dry_mass1, type = 3)
r.squaredGLMM(mod_dry_mass1)

anova(mod_ldmc1, type = 3)
r.squaredGLMM(mod_ldmc1)

anova(mod_area1, type = 3)
r.squaredGLMM(mod_area1)

anova(mod_sla1, type = 3)
r.squaredGLMM(mod_sla1)

anova(mod_thickness1, type = 3)
r.squaredGLMM(mod_thickness1)

anova(mod_height1, type = 3)
r.squaredGLMM(mod_height1)

# Generate data that fit the linear models from above 

lm_pred = traits_wide_lm[,1:9] %>% 
  mutate(ln_dry_mass = predict(mod_dry_mass1, re.form=NA, newdata=.),
         ln_area = predict(mod_area1, re.form=NA, newdata=.),
         ln_sla = predict(mod_sla1, re.form=NA, newdata=.),
         ln_ldmc = predict(mod_ldmc1, re.form=NA, newdata=.),
         ln_thickness = predict(mod_thickness1, re.form=NA, newdata=.))

#creating separate data for plant height
lm_pred_plant_height = plant_height[,1:4,6] %>% 
  mutate(ln_height = predict(mod_height1, re.form=NA, newdata=.),)

#Setting factor levels for species so the legend is in the correct order
lm_pred$taxon <- factor(lm_pred$taxon, levels = c("Halenia umbellata", 
                                                  "Lachemilla orbiculata", 
                                                  "Paspalum bonplandianum", 
                                                  "Rhynchospora macrochaeta", 
                                                  "Gaultheria glomerata", 
                                                  "Vaccinium floribundum"))

# Plot the raw data and the linear models
(p_height <- ggplot() +
    geom_smooth(data = plant_height, aes(x = elevation, y = ln_height),method = "lm",se=FALSE, colour = "#4F4F4F")+ #plots the line across all species
    geom_line(data = lm_pred_plant_height, aes(x = elevation, y = ln_height, colour = taxon)) + #plots the predicted line for each species
    geom_point(data = plant_height, aes(x = elevation, y = ln_height, colour=taxon), #plots the original data points
            size = 0.1,alpha = 0.3) +
    ylab("ln Plant Height")+ 
    xlab(NULL)+
    scale_x_continuous(breaks= c(3101,3468,3715))+
    scale_color_manual(values = pal_lm) +
    my_theme)

(p_drymass <- ggplot() +
  geom_smooth(data = traits_wide_lm, aes(x = elevation, y = ln_dry_mass),method = "lm",se=FALSE, colour = "#4F4F4F")+
  geom_line(data = lm_pred, aes(x = elevation, y = ln_dry_mass, color = taxon)) + 
  geom_point(data = traits_wide_lm, aes(x = elevation, y = ln_dry_mass,color = taxon), 
               size = 0.1, alpha = 0.3) +
  ylab("ln Dry Mass")+
  xlab(NULL)+
  scale_x_continuous(breaks= c(3101,3468,3715))+
  scale_color_manual(values = pal_lm) +
  my_theme) 

(p_area <- ggplot() +
  geom_smooth(data = traits_wide_lm, aes(x = elevation, y = ln_area),method = "lm",se=FALSE, colour = "#4F4F4F")+
  geom_line(data = lm_pred, aes(x = elevation, y = ln_area, color = taxon)) +
  geom_point(data = traits_wide_lm, aes(x = elevation, y = ln_area, colour = taxon), 
              size = 0.1,alpha = 0.3) +
  ylab("ln Leaf Area")+
  xlab(NULL)+
  scale_x_continuous(breaks= c(3101,3468,3715))+
  scale_color_manual(values = pal_lm) +
  my_theme) 


(p_sla <- ggplot() +
  geom_smooth(data = traits_wide_lm, aes(x = elevation, y = ln_sla),method = "lm",se=FALSE, colour = "#4F4F4F")+
  geom_line(data = lm_pred, aes(x = elevation, y = ln_sla, color = taxon)) + 
  geom_point(data = traits_wide_lm, aes(x = elevation, y = ln_sla,colour = taxon), 
              size = 0.1, alpha = 0.3) +
  ylab("ln SLA")+
  xlab(NULL) +
  scale_x_continuous(breaks= c(3101,3468,3715))+
  scale_color_manual(values = pal_lm) +
  my_theme)

(p_ldmc <- ggplot() +
  geom_smooth(data = traits_wide_lm, aes(x = elevation, y = ln_ldmc),method = "lm",se=FALSE, colour = "#4F4F4F")+
  geom_line(data = lm_pred, aes(x = elevation, y = ln_ldmc, color = taxon)) + 
  geom_point(data = traits_wide_lm, aes(x = elevation, y = ln_ldmc,colour = taxon), 
              size = 0.1,alpha = 0.3) +
  ylim(-1.8,-0.5) +
  ylab("ln LDMC")+
  xlab(NULL)+
  scale_x_continuous(breaks= c(3101,3468,3715))+
  scale_color_manual(values = pal_lm) +
  my_theme)


(p_thickness <- ggplot() +
  geom_smooth(data = traits_wide_lm, aes(x = elevation, y = ln_thickness),method = "lm",se=FALSE, colour = "#4F4F4F")+
  geom_line(data = lm_pred, aes(x = elevation, y = ln_thickness, color = taxon)) + 
  geom_point(data = traits_wide_lm, aes(x = elevation, y = ln_thickness,colour = taxon), 
              size = 0.1, alpha = 0.3) +
  ylab("ln Leaf Thickness") +
  xlab("Elevation (m)")+
  scale_x_continuous(breaks= c(3101,3468,3715))+
  scale_color_manual(values = pal_lm) +
  my_theme)


# Assemble plot
png("lm_Plot1.png", width = 6,height = 5, units = "in", res=800)

ggarrange(p_drymass, p_ldmc,p_area, p_sla, p_thickness,p_height,
          nrow =2, ncol = 3, heights = c(1,1), labels="AUTO", align = "hv",
          common.legend = TRUE, legend="bottom")
dev.off()

#### Plot 2 - Intra-individual Coefficient of Variation in traits with elevation####

# Coefficient of variation computed at the individual level
# Function for computing coefficient of variation
cv = function(x) { return(sd(x)/mean(x)) }

# Compute CV for each trait for each individual
trait_cvs = traits %>% 
  subset(trait!="plant_height_cm")%>% #removes plant height
  group_by(individual_uid, trait, site, elevation, taxon) %>%
  dplyr::summarize(cv_trait=cv(value))

# plotting coefficient of variation for each trait to see if there is a skew 
trait_cvs %>% 
  ggplot(aes(x = cv_trait)) +
  geom_histogram() +
  facet_wrap(~trait, scales = "free")

#log transforming coefficient of variation for each trait
trait_cvs<- trait_cvs %>% 
  mutate(log_cv_trait = log(cv_trait))%>%
  select(-c(cv_trait))

# plotting log transformed CV for each trait to see if skew is gone
trait_cvs %>% 
  ggplot(aes(x = log_cv_trait)) +
  geom_histogram() +
  facet_wrap(~trait, scales = "free")

# Reorganize as wide data - easier for what follows
trait_cvs_wide = spread(trait_cvs, key = trait, value = log_cv_trait)

# create model for each trait CV
# note that CVs are naturally log transformed already
mod_CVdry_mass1 = lmer(dry_mass_g ~ scale(elevation) * taxon + (1|site),
                       data = trait_cvs_wide,
                       na.action = na.omit)
mod_CVarea1 = lmer(leaf_area_cm2 ~ scale(elevation) * taxon + (1|site),
                   data = trait_cvs_wide,
                   na.action = na.omit)
mod_CVsla1 = lmer(sla_cm2_g ~ scale(elevation) * taxon + (1|site),
                  data = trait_cvs_wide,
                  na.action = na.omit)
mod_CVldmc1 = lmer(ldmc ~ scale(elevation) * taxon + (1|site),
                   data = trait_cvs_wide,
                   na.action = na.omit)
mod_CVthickness1 = lmer(leaf_thickness_mm ~ scale(elevation) * taxon + (1|site),
                        data = trait_cvs_wide,
                        na.action = na.omit)

#calculating p and F values for each model with anova
anova(mod_CVdry_mass1, type = 3)
r.squaredGLMM(mod_CVdry_mass1)

anova(mod_CVldmc1, type = 3)
r.squaredGLMM(mod_CVldmc1)

anova(mod_CVarea1, type = 3)
r.squaredGLMM(mod_CVarea1)

anova(mod_CVsla1, type = 3)
r.squaredGLMM(mod_CVsla1)

anova(mod_CVthickness1, type = 3)
r.squaredGLMM(mod_CVthickness1)


# Generate data that matches linear models above
cv_pred = trait_cvs_wide[,1:4]

cv_pred$CVdry_mass = predict(mod_CVdry_mass1, re.form=NA, newdata=cv_pred)
cv_pred$CVarea = predict(mod_CVarea1, re.form=NA, newdata=cv_pred)
cv_pred$CVsla = predict(mod_CVsla1, re.form=NA, newdata=cv_pred)
cv_pred$CVldmc = predict(mod_CVldmc1, re.form=NA, newdata=cv_pred)
cv_pred$CVthickness = predict(mod_CVthickness1, re.form=NA, newdata=cv_pred)

#Setting factor levels for species so the legend is in the correct order
cv_pred$taxon <- factor(cv_pred$taxon, levels = c("Halenia umbellata", 
                                                  "Lachemilla orbiculata", 
                                                  "Paspalum bonplandianum", 
                                                  "Rhynchospora macrochaeta", 
                                                  "Gaultheria glomerata", 
                                                  "Vaccinium floribundum"))
# Plotting CV linear models
(p_CVdrymass <- ggplot() +
  geom_smooth(data = trait_cvs_wide, aes(x = elevation, y = dry_mass_g),method = "lm",se=FALSE, colour = "#4F4F4F")+ #plots the line across all species
  geom_line(data = cv_pred, aes(x = elevation, y = (CVdry_mass), color = taxon)) + #plots the linear models for each species
  geom_point(data = trait_cvs_wide, aes(x = elevation, y = dry_mass_g,color = taxon) , #plots the raw data points
              size = 0.1,alpha = 0.3) +
  scale_color_manual(values = pal_lm) +
  scale_x_continuous(breaks= c(3101,3468,3715))+
  ylim(-5,0.5) +
  xlab(NULL) +
  ylab("ln CV Dry mass")+
  my_theme)

(p_CVarea <- ggplot() +
  geom_smooth(data = trait_cvs_wide, aes(x = elevation, y = leaf_area_cm2),method = "lm",se=FALSE, colour = "#4F4F4F")+
  geom_line(data = cv_pred, aes(x = elevation, y = (CVarea), color = taxon)) + 
  geom_point(data = trait_cvs_wide, aes(x = elevation, y = leaf_area_cm2,color = taxon), 
             size = 0.1, alpha = 0.3) +
  scale_color_manual(values = pal_lm) +
  scale_x_continuous(breaks= c(3101,3468,3715))+
  ylim(-5,0.5) +
  xlab(NULL) +
  ylab("ln CV Leaf area")+
  my_theme)

(p_CVsla <- ggplot() +
  geom_smooth(data = trait_cvs_wide, aes(x = elevation, y = sla_cm2_g),method = "lm",se=FALSE, colour = "#4F4F4F")+
  geom_line(data = cv_pred, aes(x = elevation, y = (CVsla), color = taxon)) + 
  geom_point(data = trait_cvs_wide, aes(x = elevation, y = sla_cm2_g,color = taxon), 
              size = 0.1,alpha = 0.3) +
  scale_color_manual(values = pal_lm) +
  scale_x_continuous(breaks= c(3101,3468,3715))+
  ylim(-5,0.5) +
  xlab(NULL) +
  ylab("ln CV SLA")+
  my_theme)

(p_CVldmc <- ggplot() +
  geom_smooth(data = trait_cvs_wide, aes(x = elevation, y = ldmc),method = "lm",se=FALSE, colour = "#4F4F4F")+
  geom_line(data = cv_pred, aes(x = elevation, y = (CVldmc), color = taxon)) + 
  geom_point(data = trait_cvs_wide, aes(x = elevation, y = ldmc,color = taxon), 
              size = 0.1,alpha = 0.3) +
  scale_color_manual(values = pal_lm) +
  scale_x_continuous(breaks= c(3101,3468,3715))+
  ylim(-5,0.5) +
  xlab(NULL) +
  ylab("ln CV LDMC")+
  my_theme)

(p_CVthickness <- ggplot() +
  geom_smooth(data = trait_cvs_wide, aes(x = elevation, y = leaf_thickness_mm),method = "lm",se=FALSE, colour = "#4F4F4F")+
  geom_line(data = cv_pred, aes(x = elevation, y = (CVthickness), color = taxon))+
  geom_point(data = trait_cvs_wide, aes(x = elevation, y = leaf_thickness_mm,color = taxon), 
              size = 0.1,alpha = 0.3) +
  scale_color_manual(values = pal_lm) +
  scale_x_continuous(breaks= c(3101,3468,3715))+
  ylim(-5,0.5) +
  xlab("Elevation (m)") +
  ylab("ln CV Leaf thickness")+
  my_theme)

# Assemble plot and save as png

png("lm_Plot2.png", width = 6,height = 5, units = "in", res=800)
ggarrange(p_CVdrymass, p_CVarea, p_CVsla, p_CVldmc,p_CVthickness,
          nrow =2, ncol = 3, heights = c(1,1), labels="AUTO", align = "hv",
          common.legend = TRUE, legend="bottom")
dev.off()


#### ARCHIVE - Test of climate data in model ####
#Climate Data
Climate_raw <- read.csv("C:/Users/tjgau/Desktop/Work/PFTC5/PFTC5-Intraspecific-DA.3_LinearModels/PFTC5-Intraspecific-DA.3_LinearModels/PFTC3_Puna_PFTC5_2019_2020_Climate_clean.csv")

Climate_Data <- Climate_raw %>% filter (site == "ACJ" | site == "WAY"| site == "TRE") %>% filter (treatment == "C") %>% 
  summarySE (measurevar = "value", groupvars = c("site", "variable")) %>% 
  select(-c("N","sd","se","ci")) %>%
  pivot_wider(names_from = "variable", values_from = "value")

traits_wide_lm <-  merge(traits_wide_lm, Climate_Data,
                         by=("site"),
                         all.x = T)

plant_height <-  merge(plant_height, Climate_Data,
                       by=("site"),
                       all.x = T)

# Fit models to each trait
mod_dry_mass2 <- lmer(ln_dry_mass ~ scale(elevation) * taxon + air_temperature + soilmoisture+ (1|site/individual_uid), 
                      data = traits_wide_lm,
                      na.action=na.omit)
mod_ldmc2 <- lmer(ln_ldmc ~ scale(elevation) * taxon + air_temperature + soilmoisture+(1|site/individual_uid), 
                  data = traits_wide_lm,
                  na.action=na.omit)
mod_area2 <- lmer(ln_area ~ scale(elevation) * taxon + air_temperature + soilmoisture+ (1|site/individual_uid), 
                  data = traits_wide_lm,
                  na.action=na.omit)
mod_sla2 <- lmer(ln_sla ~ scale(elevation) * taxon + air_temperature + soilmoisture+(1|site/individual_uid), 
                 data = traits_wide_lm,
                 na.action=na.omit)
mod_thickness2 <- lmer(ln_thickness ~ scale(elevation) * taxon + air_temperature + soilmoisture+(1|site/individual_uid), 
                       data = traits_wide_lm,
                       na.action=na.omit)
mod_height2 <- lmer(ln_height ~ scale(elevation) * taxon + air_temperature + soilmoisture+(1|site), #no individual_uid because there is only one height per individual
                    data = plant_height,
                    na.action=na.omit)
# Singular fit LDMC
# Removing the plot_id nested effect does not fix the issue
# It can be resolved by removing random effects altogether
# Not sure if this is a substantial issue since this is going in the sup

#Anova on each regression for P and F values
anova(mod_dry_mass2)
r.squaredGLMM(mod_dry_mass2)
anova(mod_ldmc2)
r.squaredGLMM(mod_ldmc2)
anova(mod_area2)
r.squaredGLMM(mod_area2)
anova(mod_sla2)
r.squaredGLMM(mod_sla2)
anova(mod_thickness2)
r.squaredGLMM(mod_thickness2)
anova(mod_height2)
r.squaredGLMM(mod_height2)


#### ARCHIVE - Plot 3 - Coefficient of Variation at plot level####
#
# Similar to previous plot, but with CVs computed at the plot level



# Compute CV for each trait for each plot

trait_cvs_plot = traits %>% 

trait_cvs_plot = traits_gathered %>% 
  group_by(site, plot_id, functional_group, family, taxon, trait, elevation) %>% 
  summarize(cv_trait = cv(value)) %>% 
  subset(cv_trait != 0) # chop zeroes - won't play nice later

# Reorganize as wide data - easier for what follows
trait_cvs_wide_plot = spread(trait_cvs_plot, key = trait, value = cv_trait)


mod_CVplot_height1 = lm(log(plant_height_cm) ~ scale(elevation) * taxon,
                                           data = trait_cvs_wide_plot,
                                           na.action = na.omit)
mod_CVplot_dry_mass1 = lm(log(dry_mass_g) ~ scale(elevation) * taxon,
                        data = trait_cvs_wide_plot,
                        na.action = na.omit)
mod_CVplot_area1 = lm(log(leaf_area_cm2) ~ scale(elevation) * taxon,
                        data = trait_cvs_wide_plot,
                        na.action = na.omit)
mod_CVplot_sla1 = lm(log(sla_cm2_g) ~ scale(elevation) * taxon,
                        data = trait_cvs_wide_plot,
                        na.action = na.omit)
mod_CVplot_ldmc1 = lm(log(ldmc) ~ scale(elevation) * taxon,
                        data = trait_cvs_wide_plot,
                        na.action = na.omit)
mod_CVplot_thickness1 = lm(log(leaf_thickness_mm) ~ scale(elevation) * taxon,
                        data = trait_cvs_wide_plot,
                        na.action = na.omit)




