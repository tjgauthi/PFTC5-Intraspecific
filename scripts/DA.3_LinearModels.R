#trait or traitsCV ~scale(elevation) * taxon + (1|site/individual)
#source(here::here(path = "scripts/0_data_import.R"))

#Generalized Linear Mixed Models with the value of the trait as responsible variable
#
# Run 0_data_import.R prior to running this script

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

#remove incomplete rows of data
traits_wide<-na.omit(traits_wide)

#filter out individuals that have 1-2 leaves only
traits_wide <-subset (traits_wide, individual_uid != "TRE_Vaccinium floribundum_5_4" &
                        individual_uid !="TRE_Vaccinium floribundum_5_1" &
                        individual_uid !="TRE_Vaccinium floribundum_3_3" &
                        individual_uid !="TRE_Vaccinium floribundum_1_2" &
                        individual_uid !="TRE_Vaccinium floribundum_1_1" &
                        individual_uid !="TRE_Lachemilla orbiculata_1_12" &
                        individual_uid !="ACJ_Rhynchospora macrochaeta_5_3" &
                        individual_uid !="ACJ_Rhynchospora macrochaeta_1_3" &
                        individual_uid !="WAY_Vaccinium floribundum_5_5" &
                        individual_uid !="WAY_Lachemilla orbiculata_1_2" &
                        individual_uid !="TRE_Vaccinium floribundum_5_2" &
                        individual_uid !="TRE_Paspalum bonplandianum_4_1" &
                        individual_uid !="TRE_Lachemilla orbiculata_5_12" &
                        individual_uid !="TRE_Lachemilla orbiculata_1_1")

#creates a long-format version of traits_wide. 
traits <- traits_wide |> 
  pivot_longer(names_to = "trait", values_to = "value", cols = 10:16)


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
mod_height1 <- lmer(ln_height ~ scale(elevation) * taxon + (1|site), #no individual_uid because there is only one height per individual
                    data = plant_height,
                    na.action=na.omit)
anova(mod_height1)#p-value and Fvalue for table 1

mod_dry_mass1 <- lmer(ln_dry_mass ~ scale(elevation) * taxon + (1|site/individual_uid), 
                    data = traits_wide_lm,
                    na.action=na.omit)
anova(mod_dry_mass1)#p-value and Fvalue for table 1

mod_area1 <- lmer(ln_area ~ scale(elevation) * taxon + (1|site/individual_uid), 
                      data = traits_wide_lm,
                      na.action=na.omit)
anova(mod_area1)#p-value and Fvalue for table 1

mod_sla1 <- lmer(ln_sla ~ scale(elevation) * taxon + (1|site/individual_uid), 
                  data = traits_wide_lm,
                  na.action=na.omit)
anova(mod_sla1)#p-value and Fvalue for table 1

mod_thickness1 <- lmer(ln_thickness ~ scale(elevation) * taxon + (1|site/individual_uid), 
                       data = traits_wide_lm,
                       na.action=na.omit)
anova(mod_thickness1)#p-value and Fvalue for table 1

mod_ldmc1 <- lmer(ln_ldmc ~ scale(elevation) * taxon + (1|site/individual_uid), 
                  data = traits_wide_lm,
                  na.action=na.omit)
anova(mod_ldmc1)#p-value and Fvalue for table 1

# Singular fit LDMC
# Removing the plot_id nested effect does not fix the issue
# It can be resolved by removing random effects altogether
# Not sure if this is a substantial issue since this is going in the sup

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
p_height <- ggplot(data = lm_pred_plant_height, aes(x = elevation, y = ln_height, colour = taxon)) +
  geom_line() + #plots the predicted line
  geom_jitter(data = plant_height,  #plots the original data points
              aes(x = elevation, y = ln_height, colour=taxon), 
              width = 10, 
              height=0, 
              size = 0.1,
              alpha = 0.3) +
  ylab("ln Plant Height")+ 
  xlab(NULL)+
  scale_x_continuous(breaks= c(3101,3468,3715))+
  scale_color_manual(values = pal_lm) +
  my_theme
p_height

p_drymass <- ggplot(data = lm_pred, aes(x = elevation, y = ln_dry_mass, color = taxon)) +
  geom_line() + 
  geom_jitter(data = traits_wide_lm, 
              aes(x = elevation, y = ln_dry_mass), 
              width = 10,
              height=0,
              size = 0.1,
              alpha = 0.3) +
  ylab("ln Dry Mass")+
  xlab(NULL)+
  scale_x_continuous(breaks= c(3101,3468,3715))+
  scale_color_manual(values = pal_lm) +
  my_theme 

p_area <- ggplot(data = lm_pred, aes(x = elevation, y = ln_area, color = taxon)) +
  geom_line() +
  geom_jitter(data = traits_wide_lm, 
              aes(x = elevation, y = ln_area), 
              width = 10, 
              height=0,
              size = 0.1,
              alpha = 0.3) +
  ylab("ln Leaf Area")+
  xlab(NULL)+
  scale_x_continuous(breaks= c(3101,3468,3715))+
  scale_color_manual(values = pal_lm) +
  my_theme 


p_sla <- ggplot(data = lm_pred, aes(x = elevation, y = ln_sla, color = taxon)) +
  geom_line() + 
  geom_jitter(data = traits_wide_lm, 
              aes(x = elevation, y = ln_sla), 
              width = 10,
              height=0,
              size = 0.1,
              alpha = 0.3) +
  ylab("ln SLA")+
  xlab(NULL) +
  scale_x_continuous(breaks= c(3101,3468,3715))+
  scale_color_manual(values = pal_lm) +
  my_theme

p_ldmc <- ggplot(data = lm_pred, aes(x = elevation, y = ln_ldmc, color = taxon)) +
  geom_line() + 
  geom_jitter(data = traits_wide_lm, 
              aes(x = elevation, y = ln_ldmc), 
              width = 10, 
              height=0,
              size = 0.1,
              alpha = 0.3) +
  ylim(-1.8,-0.5) +
  ylab("ln LDMC")+
  xlab(NULL)+
  scale_x_continuous(breaks= c(3101,3468,3715))+
  scale_color_manual(values = pal_lm) +
  my_theme


p_thickness <- ggplot(data = lm_pred, aes(x = elevation, y = ln_thickness, color = taxon)) +
  geom_line() + 
  geom_jitter(data = traits_wide_lm, 
              aes(x = elevation, y = ln_thickness), 
              width = 10, 
              size = 0.1,
              alpha = 0.3) +
  ylab("ln Leaf Thickness") +
  xlab("Elevation (m)")+
  scale_x_continuous(breaks= c(3101,3468,3715))+
  scale_color_manual(values = pal_lm) +
  my_theme 
  theme(legend.position = "bottom") #my theme removes legends. This adds one back in


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
  summarize(cv_trait=cv(value))

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
anova(mod_CVdry_mass1)#P value and F stats for Table lmm1

mod_CVarea1 = lmer(leaf_area_cm2 ~ scale(elevation) * taxon + (1|site),
                   data = trait_cvs_wide,
                   na.action = na.omit)
anova(mod_CVarea1)#P value and F stats for Table lmm1

mod_CVsla1 = lmer(sla_cm2_g ~ scale(elevation) * taxon + (1|site),
                  data = trait_cvs_wide,
                  na.action = na.omit)
anova(mod_CVsla1)

mod_CVldmc1 = lmer(ldmc ~ scale(elevation) * taxon + (1|site),
                   data = trait_cvs_wide,
                   na.action = na.omit)
anova(mod_CVldmc1)

mod_CVthickness1 = lmer(leaf_thickness_mm ~ scale(elevation) * taxon + (1|site),
                        data = trait_cvs_wide,
                        na.action = na.omit)
anova(mod_CVthickness1)


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
p_CVdrymass <- ggplot(data = cv_pred, aes(x = elevation, y = (CVdry_mass), color = taxon)) +
  geom_line() + 
  geom_jitter(data = trait_cvs_wide, 
              aes(x = elevation, y = dry_mass_g), 
              width = 10, 
              height=0,
              size = 0.1,
              alpha = 0.3) +
  scale_color_manual(values = pal_lm) +
  ylim(-5,0.5) +
  xlab(NULL) +
  ylab("ln CV Dry mass")+
  my_theme

p_CVarea <- ggplot(data = cv_pred, aes(x = elevation, y = (CVarea), color = taxon)) +
  geom_line() + 
  geom_jitter(data = trait_cvs_wide, 
              aes(x = elevation, y = leaf_area_cm2), 
              width = 10, 
              size = 0.1,
              alpha = 0.3) +
  scale_color_manual(values = pal_lm) +
  ylim(-5,0.5) +
  xlab(NULL) +
  ylab("ln CV Leaf area")+
  my_theme

p_CVsla <- ggplot(data = cv_pred, aes(x = elevation, y = (CVsla), color = taxon)) +
  geom_line() + 
  geom_jitter(data = trait_cvs_wide, 
              aes(x = elevation, y = sla_cm2_g), 
              width = 10, 
              size = 0.1,
              alpha = 0.3) +
  scale_color_manual(values = pal_lm) +
  ylim(-5,0.5) +
  xlab(NULL) +
  ylab("ln CV SLA")+
  my_theme

p_CVldmc <- ggplot(data = cv_pred, aes(x = elevation, y = (CVldmc), color = taxon)) +
   geom_line() + 
   geom_jitter(data = trait_cvs_wide, 
              aes(x = elevation, y = ldmc), 
              width = 10, 
              size = 0.1,
              alpha = 0.3) +
  scale_color_manual(values = pal_lm) +
  ylim(-5,0.5) +
  xlab(NULL) +
  ylab("ln CV LDMC")+
  my_theme

p_CVthickness <- ggplot(data = cv_pred, aes(x = elevation, y = (CVthickness), color = taxon)) +
  geom_line()+
  geom_jitter(data = trait_cvs_wide, 
              aes(x = elevation, y = leaf_thickness_mm), 
              width = 10, 
              size = 0.1,
              alpha = 0.3) +
  scale_color_manual(values = pal_lm) +
  ylim(-5,0.5) +
  xlab("Elevation (m)") +
  ylab("ln CV Leaf thickness")+
  my_theme

# Assemble plot and save as png

png("lm_Plot2.png", width = 6,height = 5, units = "in", res=800)
ggarrange(p_CVdrymass, p_CVarea, p_CVsla, p_CVldmc,p_CVthickness,
          nrow =2, ncol = 3, heights = c(1,1), labels="AUTO", align = "hv",
          common.legend = TRUE, legend="bottom")
dev.off()


#### ARCHIVE - Plot 3 - Coefficient of Variation at plot level####
#
# Similar to previous plot, but with CVs computed at the plot level
#


# Compute CV for each trait for each plot
trait_cvs_plot = traits %>% 
  group_by(site, plot_id, functional_group, family, taxon, trait, elevation) %>% 
  summarize(cv_trait = cv(value)) %>% 
  subset(cv_trait != 0) # chop zeroes - won't play nice later

# Reorganize as wide data - easier for what follows
trait_cvs_wide_plot = spread(trait_cvs_plot, key = trait, value = cv_trait)

# mod_CVheight1 = lmer(log(plant_height_cm) ~ scale(elevation) * taxon + (1|site),
#                      data = trait_cvs_wide_plot,
#                      na.action = na.omit)
##### SINGULAR - note most of the plant_height CVs are zero b/c measured on the
##### same individuals.

#mod_CVheight2 = lmer(log(plant_height_cm) ~ scale(elevation) * taxon + (1|site),
#                     data = subset(trait_cvs_wide, plant_height_cm != 0),
#                     na.action = na.omit)
##### OK

#mod_CVdry_mass1 = lmer(dry_mass_g ~ scale(elevation) * taxon + (1|site/plot_id),
#                   data = trait_cvs_wide,
#                   na.action = na.omit)
# OK

# mod_CVdry_mass1 = lmer(log(dry_mass_g) ~ scale(elevation) * taxon + (1|site),
#                        data = trait_cvs_wide_plot,
#                        na.action = na.omit)
# # OK
# 
# #mod_CVarea1 = lmer(leaf_area_cm2 ~ scale(elevation) * taxon + (1|site/plot_id),
# #                   data = trait_cvs_wide,
# #                   na.action = na.omit)
# ##### SINGULAR
# 
# mod_CVarea1 = lmer(log(leaf_area_cm2) ~ scale(elevation) * taxon + (1|site),
#                    data = trait_cvs_wide_plot,
#                    na.action = na.omit)
# # OK
# 
# #mod_CVsla1 = lmer(sla_cm2_g ~ scale(elevation) * taxon + (1|site/plot_id),
# #                   data = trait_cvs_wide,
# #                   na.action = na.omit)
# ##### SINGULAR
# mod_CVsla1 = lmer(log(sla_cm2_g) ~ scale(elevation) * taxon + (1|site),
#                   data = trait_cvs_wide_plot,
#                   na.action = na.omit)
# # OK
# 
# #mod_CVldmc1 = lmer(ldmc ~ scale(elevation) * taxon + (1|site/plot_id),
# #                   data = trait_cvs_wide,
# #                   na.action = na.omit)
# ##### SINGULAR
# mod_CVldmc1 = lmer(log(ldmc) ~ scale(elevation) * taxon + (1|site),
#                    data = trait_cvs_wide_plot,
#                    na.action = na.omit)
# 
# #mod_CVthickness1 = lmer(leaf_thickness_mm ~ scale(elevation) * taxon + (1|site/plot_id),
# #                 data = trait_cvs_wide,
# #                 na.action = na.omit)
# ##### SINGULAR
# mod_CVthickness1 = lmer(log(leaf_thickness_mm) ~ scale(elevation) * taxon + (1|site),
#                         data = trait_cvs_wide_plot,
#                         na.action = na.omit)
# OK


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




# Generate model predictions
cv_pred = trait_cvs_wide_plot[,1:7]
cv_pred$CVheight = predict(mod_CVplot_height1, re.form=NA, newdata=cv_pred)
cv_pred$CVdry_mass = predict(mod_CVplot_dry_mass1, re.form=NA, newdata=cv_pred)
cv_pred$CVarea = predict(mod_CVplot_area1, re.form=NA, newdata=cv_pred)
cv_pred$CVsla = predict(mod_CVplot_sla1, re.form=NA, newdata=cv_pred)
cv_pred$CVldmc = predict(mod_CVplot_ldmc1, re.form=NA, newdata=cv_pred)
cv_pred$CVthickness = predict(mod_CVplot_thickness1, re.form=NA, newdata=cv_pred)


# Build panels
p_CVheight <- ggplot(data = cv_pred, aes(x = elevation, y = (CVheight), color = taxon)) +
  geom_jitter(data = (trait_cvs_wide_plot), 
              aes(x = elevation, y = log(plant_height_cm)), 
              width = 10, 
              size = 0.1,
              alpha = 0.3) +
  my_theme +
  scale_color_manual(values = pal_lm) +
  geom_line() + 
  ylim(-4,0) +
  xlab("Elevation (m)") +
  ylab("ln(CV Height)")

p_CVdrymass <- ggplot(data = cv_pred, aes(x = elevation, y = (CVdry_mass), color = taxon)) +
  geom_jitter(data = (trait_cvs_wide_plot), 
              aes(x = elevation, y = log(dry_mass_g)), 
              width = 10, 
              size = 0.1,
              alpha = 0.3) +
  my_theme +
  scale_color_manual(values = pal_lm) +
  geom_line() + 
  ylim(-4,0) +
  xlab("Elevation (m)") +
  ylab("ln(CV Dry mass)")

p_CVarea <- ggplot(data = cv_pred, aes(x = elevation, y = (CVarea), color = taxon)) +
  geom_jitter(data = (trait_cvs_wide_plot), 
              aes(x = elevation, y = log(leaf_area_cm2)), 
              width = 10, 
              size = 0.1,
              alpha = 0.3) +
  my_theme +
  scale_color_manual(values = pal_lm) +
  geom_line() + 
  ylim(-4,0) +
  xlab("Elevation (m)") +
  ylab("ln(CV Leaf area)")

p_CVsla <- ggplot(data = cv_pred, aes(x = elevation, y = (CVsla), color = taxon)) +
  geom_jitter(data = (trait_cvs_wide_plot_plot), 
              aes(x = elevation, y = log(sla_cm2_g)), 
              width = 10, 
              size = 0.1,
              alpha = 0.3) +
  my_theme +
  scale_color_manual(values = pal_lm) +
  geom_line() + 
  ylim(-4,0) +
  xlab("Elevation (m)") +
  ylab("ln(CV SLA)")

p_CVldmc <- ggplot(data = cv_pred, aes(x = elevation, y = (CVldmc), color = taxon)) +
  geom_jitter(data = (trait_cvs_wide_plot), 
              aes(x = elevation, y = log(ldmc)), 
              width = 10, 
              size = 0.1,
              alpha = 0.3) +
  my_theme +
  scale_color_manual(values = pal_lm) +
  geom_line() + 
  ylim(-4,0) +
  xlab("Elevation (m)") +
  ylab("ln(CV LDMC)")

p_CVthickness <- ggplot(data = cv_pred, aes(x = elevation, y = (CVthickness), color = taxon)) +
  geom_jitter(data = (trait_cvs_wide_plot), 
              aes(x = elevation, y = log(leaf_thickness_mm)), 
              width = 10, 
              size = 0.1,
              alpha = 0.3) +
  my_theme +
  scale_color_manual(values = pal_lm) +
  geom_line() + 
  ylim(-4,0) +
  xlab("Elevation (m)") +
  ylab("ln(CV Leaf thickness)")

# Assemble plot and save PDF
pdf("cv_plot_fig.pdf", width=6, height=4)
grid.arrange(p_CVheight, p_CVdrymass, p_CVarea, p_CVsla, p_CVldmc,p_CVthickness, nrow =2)
dev.off()





#
#### ARCHIVE - Plot 4 - Summary figure of slopes for plot 2 & 3 ####
#


# Pull slope values associated with each trait and species
# species trait slope


slopes.cv.ind = data.frame(species = rep(rel_sp[c(1,5,6,2,4,3)],6),
                           trait = c(rep("area", 6),rep("dry_mass", 6),rep("height", 6),rep("ldmc", 6),rep("sla", 6),rep("thickness", 6)),
                           slope = c(coef(summary(mod_CVarea1))[,"Estimate"][c(2,8:12)],
                                     coef(summary(mod_CVdry_mass1))[,"Estimate"][c(2,8:12)],
                                     coef(summary(mod_CVheight1))[,"Estimate"][c(2,8:12)],
                                     coef(summary(mod_CVldmc1))[,"Estimate"][c(2,8:12)],
                                     coef(summary(mod_CVsla1))[,"Estimate"][c(2,8:12)],
                                     coef(summary(mod_CVthickness1))[,"Estimate"][c(2,8:12)]))



slopes.gago = subset(slopes.cv.ind, species == "Gaultheria glomerata")
slopes.others = subset(slopes.cv.ind, species != "Gaultheria glomerata")

slopes.cv.ind = merge(slopes.others, slopes.gago %>% select(-species), by=c( "trait")) %>% 
  mutate(slope = slope.x+slope.y) %>% select(-slope.x,-slope.y) %>% 
  bind_rows(slopes.gago)



slopes.cv.plot = data.frame(species = rep(rel_sp[c(1,5,6,2,4,3)],6),
                            trait = c(rep("area", 6),rep("dry_mass", 6),rep("height", 6),rep("ldmc", 6),rep("sla", 6),rep("thickness", 6)),
                            slope = c(coef(summary(mod_CVplot_area1))[,"Estimate"][c(2,8:12)],
                                      coef(summary(mod_CVplot_dry_mass1))[,"Estimate"][c(2,8:12)],
                                      coef(summary(mod_CVplot_height1))[,"Estimate"][c(2,8:12)],
                                      coef(summary(mod_CVplot_ldmc1))[,"Estimate"][c(2,8:12)],
                                      coef(summary(mod_CVplot_sla1))[,"Estimate"][c(2,8:12)],
                                      coef(summary(mod_CVplot_thickness1))[,"Estimate"][c(2,8:12)]))


slopes.gago = subset(slopes.cv.plot, species == "Gaultheria glomerata")
slopes.others = subset(slopes.cv.plot, species != "Gaultheria glomerata")

slopes.cv.plot = merge(slopes.others, slopes.gago %>% select(-species), by=c( "trait")) %>% 
  mutate(slope = slope.x+slope.y) %>% select(-slope.x,-slope.y) %>% 
  bind_rows(slopes.gago)


p = ggplot(slopes.cv.ind, aes(y = trait, x = slope, color = species)) +
  geom_point(pch = "|",size = 6,position = position_nudge(y = 0.2)) +
  geom_point(data=slopes.cv.plot, pch="|", size=6,position = position_nudge(y=-0.2)) +
  scale_color_manual(values = pal_lm) +
  my_theme +
  geom_vline(xintercept = 0, lty = 2) +
  xlim(-0.3,0.5) +
  geom_hline(yintercept = 1:6) +
  theme(axis.title.y = element_blank())

pdf("plot_overview.pdf", 4,4)
p
dev.off()







#
# Generate data table with p values (do we also need parameter estimates?)
#



pval.traits = data.frame(area = anova(mod_area1)$`Pr(>F)`,
           dry_mass = anova(mod_dry_mass1)$`Pr(>F)`,
           height = anova(mod_height1)$`Pr(>F)`,
           ldmc = anova(mod_ldmc1)$`Pr(>F)`,
           sla = anova(mod_sla1)$`Pr(>F)`,
           thickness = anova(mod_thickness1)$`Pr(>F)`)

rownames(pval.traits) =  c("scale(elevation)", "taxon", "scale(elevation):taxon")
write.csv(pval.traits, "pval.traits.csv")


pval.CV.ind = data.frame(area = anova(mod_CVarea1)$`Pr(>F)`,
                         dry_mass = anova(mod_CVdry_mass1)$`Pr(>F)`,
                         height = anova(mod_CVheight1)$`Pr(>F)`,
                         ldmc = anova(mod_CVldmc1)$`Pr(>F)`,
                         sla = anova(mod_CVsla1)$`Pr(>F)`,
                         thickness = anova(mod_CVthickness1)$`Pr(>F)`)

rownames(pval.CV.ind) =  c("scale(elevation)", "taxon", "scale(elevation):taxon")

write.csv(pval.traits, "pval.CV.ind.csv")



pval.CV.plot = data.frame(area = anova(mod_area1)$`Pr(>F)`,
                         dry_mass = anova(mod_dry_mass1)$`Pr(>F)`,
                         height = anova(mod_height1)$`Pr(>F)`,
                         ldmc = anova(mod_ldmc1)$`Pr(>F)`,
                         sla = anova(mod_sla1)$`Pr(>F)`,
                         thickness = anova(mod_thickness1)$`Pr(>F)`)

rownames(pval.CV.plot) =  c("scale(elevation)", "taxon", "scale(elevation):taxon")

write.csv(pval.traits, "pval.CV.plot.csv")




#
# END
#



