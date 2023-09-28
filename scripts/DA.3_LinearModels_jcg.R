#trait or traitsCV ~scale(elevation) * taxon + (1|site/plotID)
#source(here::here(path = "scripts/0_data_import.R"))



# Josef: Per the recent email from the larger group, we are log-transforming all
# traits and using a consistent random effects structure on all models. I have gone
# through the code to try to make this all consistent



####Generalized Linear Mixed Models with the value of the trait as responsable variable
library(lme4)
library(ggplot2)
#library(effects)
library(gridExtra)
library(lmerTest)

# Plotting theme
my_theme = theme_bw() + 
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Color palette
pal_lm <- c("#016392", "#A0CBE8", "#E19825", "#F7C480", "#3E8853", "#9FCD99")


#
# Run 0_data_import.R prior to running this script
#









#
# Begin code for Plot 1 - variation in trait values with elevation
#


# Set up dataframe for fitting models to log-transformed data
traits_wide_lm = traits_wide %>% 
  mutate(ln_height = log(plant_height_cm),
         ln_dry_mass = log(dry_mass_g),
         ln_area = log(leaf_area_cm2),
         ln_sla = log(sla_cm2_g),
         ln_ldmc = log(ldmc),
         ln_thickness = log(leaf_thickness_mm))

# Fit models to each trait
mod_height1 <- lmer(ln_height ~ scale(elevation) * taxon + (1|site/plot_id), 
                    data = traits_wide_lm,
                    na.action=na.omit)

mod_dry_mass1 <- lmer(ln_dry_mass ~ scale(elevation) * taxon + (1|site/plot_id), 
                    data = traits_wide_lm,
                    na.action=na.omit)

mod_area1 <- lmer(ln_area ~ scale(elevation) * taxon + (1|site/plot_id), 
                      data = traits_wide_lm,
                      na.action=na.omit)

mod_sla1 <- lmer(ln_sla ~ scale(elevation) * taxon + (1|site/plot_id), 
                  data = traits_wide_lm,
                  na.action=na.omit)

mod_ldmc1 <- lmer(ln_ldmc ~ scale(elevation) * taxon + (1|site/plot_id), 
                  data = traits_wide_lm,
                  na.action=na.omit)
###### Singular fit LDMC ######
# Removing the plot_id nested effect does not fix the issue
# It can be resolved by removing random effects altogether
# Not sure if this is a substantial issue since this is going in the sup

mod_thickness1 <- lmer(ln_thickness ~ scale(elevation) * taxon + (1|site/plot_id), 
                  data = traits_wide_lm,
                  na.action=na.omit)


# Generate model predictions
lm_pred = traits_wide_lm[,1:9] %>% 
  mutate(ln_height = predict(mod_height1, re.form=NA, newdata=.),
         ln_dry_mass = predict(mod_dry_mass1, re.form=NA, newdata=.),
         ln_area = predict(mod_area1, re.form=NA, newdata=.),
         ln_sla = predict(mod_sla1, re.form=NA, newdata=.),
         ln_ldmc = predict(mod_ldmc1, re.form=NA, newdata=.),
         ln_thickness = predict(mod_thickness1, re.form=NA, newdata=.),
  )

# Generate panels
p_height <- ggplot(data = lm_pred, aes(x = elevation, y = ln_height, color = taxon)) +
  geom_jitter(data = traits_wide_lm, 
              aes(x = elevation, y = ln_height), 
              width = 10, 
              size = 0.1,
              alpha = 0.3) +
  my_theme +
  scale_color_manual(values = pal_lm) +
  geom_line() + 
  xlab("Elevation (m)") +
  ylab("ln Height")

p_drymass <- ggplot(data = lm_pred, aes(x = elevation, y = ln_dry_mass, color = taxon)) +
  geom_jitter(data = traits_wide_lm, 
              aes(x = elevation, y = ln_dry_mass), 
              width = 10, 
              size = 0.1,
              alpha = 0.3) +
  my_theme +
  scale_color_manual(values = pal_lm) +
  geom_line() + 
  ylim(-6,-0.5) +
  xlab("Elevation (m)") +
  ylab("ln Dry mass")

p_area <- ggplot(data = lm_pred, aes(x = elevation, y = ln_area, color = taxon)) +
  geom_jitter(data = traits_wide_lm, 
              aes(x = elevation, y = ln_area), 
              width = 10, 
              size = 0.1,
              alpha = 0.3) +
  my_theme +
  scale_color_manual(values = pal_lm) +
  geom_line() +
  ylim(-2,4) +
  xlab("Elevation (m)") +
  ylab("ln Leaf area")

p_sla <- ggplot(data = lm_pred, aes(x = elevation, y = ln_sla, color = taxon)) +
  geom_jitter(data = traits_wide_lm, 
              aes(x = elevation, y = ln_sla), 
              width = 10, 
              size = 0.1,
              alpha = 0.3) +
  my_theme +
  scale_color_manual(values = pal_lm) +
  geom_line() + 
  ylim(3.5,6) +
  xlab("Elevation (m)") +
  ylab("ln SLA")

p_ldmc <- ggplot(data = lm_pred, aes(x = elevation, y = ln_ldmc, color = taxon)) +
  geom_jitter(data = traits_wide_lm, 
              aes(x = elevation, y = ln_ldmc), 
              width = 10, 
              size = 0.1,
              alpha = 0.3) +
  my_theme +
  scale_color_manual(values = pal_lm) +
  geom_line() + 
  ylim(-1.8,-0.5) +
  xlab("Elevation (m)") +
  ylab("ln LDMC")

p_thickness <- ggplot(data = lm_pred, aes(x = elevation, y = ln_thickness, color = taxon)) +
  geom_jitter(data = traits_wide_lm, 
              aes(x = elevation, y = ln_thickness), 
              width = 10, 
              size = 0.1,
              alpha = 0.3) +
  my_theme +
  scale_color_manual(values = pal_lm) +
  geom_line() + 
  ylim(-3,-0.5) +
  xlab("Elevation (m)") +
  ylab("ln Leaf thickness")

# Assemble plot, save to PDF
pdf("lm_fig.pdf", width = 6, height = 4)
grid.arrange(p_height, p_drymass, p_area, p_sla, p_ldmc, p_thickness, nrow =2)
dev.off()










#
# Begin code for Plot 2 - Intra-individual variation in traits with elevation
#
# CVs computed at the individual level
#


#####
# Function for computing coefficient of variation
cv = function(x) { return(sd(x)/mean(x)) } #Josef function


# Compute CV for each trait for each individual
trait_cvs = traits %>% 
  group_by(site, plot_id, individual_nr, functional_group, family, taxon, trait, elevation) %>% 
  summarize(cv_trait = cv(value)) %>% 
  subset(cv_trait != 0) # Snip zeroes - these cause problems later

# Reorganize as wide data - easier for what follows
trait_cvs_wide = spread(trait_cvs, key = trait, value = cv_trait)

#mod_CVheight1 = lmer(log(plant_height_cm) ~ scale(elevation) * taxon + (1|site/plot_id),
#                   data = trait_cvs_wide,
#                   na.action = na.omit)
##### SINGULAR - note most of the plant_height CVs are zero b/c measured on the
##### same individuals, but these were removed above

mod_CVheight1 = lmer(log(plant_height_cm) ~ scale(elevation) * taxon + (1|site),
                     data = subset(trait_cvs_wide, plant_height_cm != 0),
                     na.action = na.omit)
##### OK

#mod_CVdry_mass1 = lmer(dry_mass_g ~ scale(elevation) * taxon + (1|site/plot_id),
#                   data = trait_cvs_wide,
#                   na.action = na.omit)
# OK

mod_CVdry_mass1 = lmer(log(dry_mass_g) ~ scale(elevation) * taxon + (1|site),
                       data = trait_cvs_wide,
                       na.action = na.omit)
# OK

#mod_CVarea1 = lmer(leaf_area_cm2 ~ scale(elevation) * taxon + (1|site/plot_id),
#                   data = trait_cvs_wide,
#                   na.action = na.omit)
##### SINGULAR

mod_CVarea1 = lmer(log(leaf_area_cm2) ~ scale(elevation) * taxon + (1|site),
                   data = trait_cvs_wide,
                   na.action = na.omit)
# OK

#mod_CVsla1 = lmer(sla_cm2_g ~ scale(elevation) * taxon + (1|site/plot_id),
#                   data = trait_cvs_wide,
#                   na.action = na.omit)
##### SINGULAR
mod_CVsla1 = lmer(log(sla_cm2_g) ~ scale(elevation) * taxon + (1|site),
                  data = trait_cvs_wide,
                  na.action = na.omit)
# OK

#mod_CVldmc1 = lmer(ldmc ~ scale(elevation) * taxon + (1|site/plot_id),
#                   data = trait_cvs_wide,
#                   na.action = na.omit)
##### SINGULAR
mod_CVldmc1 = lmer(log(ldmc) ~ scale(elevation) * taxon + (1|site),
                   data = trait_cvs_wide,
                   na.action = na.omit)

#mod_CVthickness1 = lmer(leaf_thickness_mm ~ scale(elevation) * taxon + (1|site/plot_id),
#                 data = trait_cvs_wide,
#                 na.action = na.omit)
##### SINGULAR
mod_CVthickness1 = lmer(log(leaf_thickness_mm) ~ scale(elevation) * taxon + (1|site),
                        data = trait_cvs_wide,
                        na.action = na.omit)
# OK


# Generate model predictions
cv_pred = trait_cvs_wide[,1:7]
cv_pred$CVheight = predict(mod_CVheight1, re.form=NA, newdata=cv_pred)
cv_pred$CVdry_mass = predict(mod_CVdry_mass1, re.form=NA, newdata=cv_pred)
cv_pred$CVarea = predict(mod_CVarea1, re.form=NA, newdata=cv_pred)
cv_pred$CVsla = predict(mod_CVsla1, re.form=NA, newdata=cv_pred)
cv_pred$CVldmc = predict(mod_CVldmc1, re.form=NA, newdata=cv_pred)
cv_pred$CVthickness = predict(mod_CVthickness1, re.form=NA, newdata=cv_pred)


# Build panels
p_CVheight <- ggplot(data = cv_pred, aes(x = elevation, y = (CVheight), color = taxon)) +
  geom_jitter(data = (trait_cvs_wide), 
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
  geom_jitter(data = (trait_cvs_wide), 
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
  geom_jitter(data = (trait_cvs_wide), 
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
  geom_jitter(data = (trait_cvs_wide), 
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
  geom_jitter(data = (trait_cvs_wide), 
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
  geom_jitter(data = (trait_cvs_wide), 
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
pdf("cv_individual_fig.pdf", width=6, height=4)
grid.arrange(p_CVheight, p_CVdrymass, p_CVarea, p_CVsla, p_CVldmc,p_CVthickness, nrow =2)
dev.off()












#
# Begin code for Plot 3 - Intraspecific variation in traits with elevation
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
# Draft of an overview figure
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



