library(lme4)
library(lmerTest)
library(nlme)
library(gridExtra)

source(here::here(path = "scripts/DC.2_add_scaffold_chemical_traits.R"))

# Basic plotting theme
my_theme = theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Add new column with uid for each individual plant
traits_chem_final$ind_uid_new = paste(traits_chem_final$site,
                                      traits_chem_final$plot_id,
                                      traits_chem_final$taxon,
                                      traits_chem_final$individual_nr) %>% as.factor()

# Reorder taxon levels for correct colour assignment in plots
traits_chem_final$taxon <- fct_reorder(traits_chem_final$taxon, 
                                       as.numeric(as.factor(traits_chem_final$functional_group)))

#
#
#
#
#

# Idea 1 - trait-trait relationships at the overall, species, and individual level
# Ordinary least squares approach

# Scale and centre trait values
traits_scaled = traits_chem_final %>% 
  mutate(ldmc_sc = scale(ldmc), sla_sc = scale(sla_cm2_g), 
         leaf_thickness_sc = scale(leaf_thickness_mm))

# Make new uniqueID
traits_scaled$ind_uid_new = paste(traits_scaled$site,
                                  traits_scaled$plot_id,
                                  traits_scaled$taxon,
                                  traits_scaled$individual_nr) %>% as.factor()

# Generate plot showing overall, species, and individual regressions
p1 = ggplot(traits_scaled, aes(x=ldmc_sc, y=leaf_thickness_sc, color=taxon)) +
  my_theme +
  scale_color_manual(values = pal_lm) +
  ylim(-2,3) +
  geom_point(size=0.9) +
  geom_smooth(method="lm", se=F,lwd=0.4, color="grey80", aes(group=ind_uid_new)) +
  geom_smooth(method="lm", se=F, lwd=0.8, aes(color=taxon)) +
  geom_smooth(method="lm", se=F, color="black", lty=2) 

# Make inset plot showing distribution of slopes for individuals
trait.lm.list = nlme::lmList(leaf_thickness_sc ~ ldmc_sc | ind_uid_new, data = traits_scaled, na.action = na.exclude)
model_coef = subset(coef(trait.lm.list), !is.na(ldmc_sc))

z = lm(leaf_thickness_sc ~ ldmc_sc, traits_scaled)
species.lm.list = nlme::lmList(leaf_thickness_sc ~ ldmc_sc | taxon, data = traits_scaled, na.action = na.exclude)
species.slopes = data.frame(coef(species.lm.list))
species.slopes$taxon = rownames(species.slopes)

inset_plot = ggplotGrob(
  ggplot(data = model_coef, aes(x = ldmc_sc)) + 
    geom_density() + 
    theme_classic() +
    geom_vline(xintercept = 0, lty = 3) +
    geom_vline(xintercept = coef(z)[2], lty = 2) +
    geom_vline(data = species.slopes, aes(xintercept = ldmc_sc, color=taxon)) +
    xlim(c(-3,3)) +
    xlab("Slope") +
    #ylim(c(0,22)) +
    scale_color_manual(values = pal_lm) +
    theme(legend.position = "none") +
    theme(axis.title = element_text( size=8),
          rect = element_rect(fill = "transparent"),
          plot.background = element_rect(colour = "transparent") 
    )
)

p1 + annotation_custom(grob = inset_plot, xmin = -3.4, xmax = -0.2, ymin = 1.6, ymax = 3)

#
#
#
#
#
# Same plot, but using linear mixed model
# Build model with random slopes for individul nested in taxon
z1 = lmer(leaf_thickness_sc ~ ldmc_sc +  (ldmc_sc|taxon/ind_uid_new),
          REML = F,
          data = traits_scaled)

# Loop over species, grab species name and min and max ldmc value
pred.dat = c()
for (i in unique(traits_scaled$taxon)) {
  cur.dat = subset(traits_scaled, taxon == i & !is.na(ldmc_sc))
  new.dat = data.frame(taxon = i, ldmc_sc = c(min(cur.dat$ldmc_sc), max(cur.dat$ldmc_sc)))
  pred.dat = bind_rows(pred.dat,new.dat)
}

# Generate predicted lines at the taxon level
pred.dat$leaf_thickness_sc = predict(z1, pred.dat, re.form =   (~ldmc_sc|taxon))

# Generate oveall relationship prediction
pred.overall = data.frame(ldmc_sc = c(min(traits_scaled$ldmc_sc, na.rm=T),
                                      max(traits_scaled$ldmc_sc, na.rm=T)))
pred.overall$leaf_thickness_sc = predict(z1, pred.overall, re.form=NA)

# Generate individual level predictions
pred.ind = c()
for (i in unique(traits_scaled$ind_uid_new)) {
  cur.dat = subset(traits_scaled, ind_uid_new == i & !is.na(ldmc_sc))
  new.dat = data.frame(ind_uid_new = i, taxon = unique(cur.dat$taxon), ldmc_sc = c(min(cur.dat$ldmc_sc), max(cur.dat$ldmc_sc)))
  pred.ind = bind_rows(pred.ind,new.dat)
}
pred.ind$leaf_thickness_sc = predict(z1, pred.ind)


# Generate plot showing overall, species, and individual regressions
p2 = ggplot(traits_scaled, aes(x=ldmc_sc, y=leaf_thickness_sc, color=taxon)) +
  my_theme +
  scale_color_manual(values = pal_lm) +
  ylim(-2,3) +
  geom_point(size=0.9) +
  geom_line(data = pred.ind, color="grey80", aes(group=ind_uid_new)) +
  geom_line(data = pred.dat, lwd=0.8, aes(color=taxon)) +
  geom_line(data = pred.overall, color = "black", lty = 2, aes(color = NA))


model_coef = c()
model_coef$ldmc = ranef(z1)$`ind_uid_new:taxon`$ldmc_sc
model_coef = as.data.frame(model_coef)

species.slopes = ranef(z1)$taxon
species.slopes$taxon = rownames(species.slopes)
as.data.frame(species.slopes)

inset_plot = ggplotGrob(
  ggplot(data = model_coef, aes(x = ldmc)) + 
    geom_density() + 
    theme_classic() +
    geom_vline(xintercept = 0, lty = 3) +
    geom_vline(xintercept = fixef(z1)[2], lty = 2) +
    geom_vline(data = species.slopes, aes(xintercept = ldmc_sc, color=taxon)) +
    xlim(c(-0.5,0.5)) +
    xlab("Slope") +
    #ylim(c(0,22)) +
    scale_color_manual(values = pal_lm) +
    theme(legend.position = "none") +
    theme(axis.title = element_text( size=8),
          rect = element_rect(fill = "transparent"),
          plot.background = element_rect(colour = "transparent") 
    )
)

p2 + annotation_custom(grob = inset_plot, xmin = -3.4, xmax = -0.2, ymin = 1.6, ymax = 3)






# 
# 
# 
# 
# 
# 
# # Generate plot showing overall, species, and individual regressions
# p1 = ggplot(traits_chem_final, aes(x=ldmc, y=leaf_thickness_mm, color=taxon)) +
#   #geom_smooth(method="lm") + 
#   #geom_smooth(method="lm", aes(color=NULL)) +
#   my_theme +
#   scale_color_manual(values = pal_lm) +
#   ylim(0,0.625) +
#   geom_point(size=0.9) +
#   geom_smooth(method="lm", se=F,lwd=0.4, color="grey80", aes(group=ind_uid_new)) +
#   geom_smooth(method="lm", se=F, lwd=0.8, aes(color=taxon)) +
#   geom_smooth(method="lm", se=F, color="black", lty=2) #+
#   #theme(legend.position = "none")
# 
# # Make inset plot showing distribution of slopes for individuals
# trait.lm.list = nlme::lmList(leaf_thickness_mm ~ ldmc | ind_uid_new, data = traits_chem_final, na.action = na.exclude)
# model_coef = subset(coef(trait.lm.list), !is.na(ldmc))
# 
# z = lm(leaf_thickness_mm ~ ldmc, traits_chem_final)
# species.lm.list = nlme::lmList(leaf_thickness_mm ~ ldmc | taxon, data = traits_chem_final, na.action = na.exclude)
# species.slopes = data.frame(coef(species.lm.list))
# species.slopes$taxon = rownames(species.slopes)
# 
# inset_plot = ggplotGrob(
#   ggplot(data = model_coef, aes(x = ldmc)) + 
#     geom_density() + 
#     theme_classic() +
#     geom_vline(xintercept = 0, lty = 3) +
#     geom_vline(xintercept = coef(z)[2], lty = 2) +
#     geom_vline(data = species.slopes, aes(xintercept = ldmc, color=taxon)) +
#     xlim(c(-3,3)) +
#     xlab("Slope") +
#     #ylim(c(0,22)) +
#     scale_color_manual(values = pal_lm) +
#     theme(legend.position = "none") +
#     theme(axis.title = element_text( size=8),
#           rect = element_rect(fill = "transparent"),
#           plot.background = element_rect(colour = "transparent") 
#           )
# )
# 
# p1 + annotation_custom(grob = inset_plot, xmin = 0, xmax = 0.3, ymin = 0.45, ymax = 0.63)
# 
# # Another example - ldmc vs. sla
# 
# # Generate plot showing overall, species, and individual regressions
# p1 = ggplot(traits_chem_final, aes(x=ldmc, y=sla_cm2_g, color=taxon)) +
#   #geom_smooth(method="lm") + 
#   #geom_smooth(method="lm", aes(color=NULL)) +
#   my_theme +
#   scale_color_manual(values = pal_lm) +
#   #ylim(0,0.625) +
#   geom_point(size=0.9) +
#   geom_smooth(method="lm", se=F,lwd=0.4, color="grey80", aes(group=ind_uid_new)) +
#   geom_smooth(method="lm", se=F, lwd=0.8, aes(color=taxon)) +
#   geom_smooth(method="lm", se=F, color="black", lty=2) #+
# #theme(legend.position = "none")
# 
# # Make inset plot showing distribution of slopes for individuals
# trait.lm.list = nlme::lmList(sla_cm2_g ~ ldmc | ind_uid_new, data = traits_chem_final, na.action = na.exclude)
# model_coef = subset(coef(trait.lm.list), !is.na(ldmc))
# 
# z = lm(sla_cm2_g ~ ldmc, traits_chem_final)
# species.lm.list = nlme::lmList(sla_cm2_g ~ ldmc | taxon, data = traits_chem_final, na.action = na.exclude)
# species.slopes = data.frame(coef(species.lm.list))
# species.slopes$taxon = rownames(species.slopes)
# 
# inset_plot = ggplotGrob(
#   ggplot(data = model_coef, aes(x = ldmc)) + 
#     geom_density() + 
#     theme_classic() +
#     geom_vline(xintercept = 0, lty = 3) +
#     geom_vline(xintercept = coef(z)[2], lty = 2) +
#     geom_vline(data = species.slopes, aes(xintercept = ldmc, color=taxon)) +
#     xlim(c(-1200,1200)) +
#     xlab("Slope") +
#     #ylim(c(0,22)) +
#     scale_color_manual(values = pal_lm) +
#     theme(legend.position = "none") +
#     theme(axis.title = element_text( size=8),
#           rect = element_rect(fill = "transparent"),
#           plot.background = element_rect(colour = "transparent") 
#     )
# )
# 
# p1 + annotation_custom(grob = inset_plot, xmin = 0.38, xmax = 0.65, ymin = 300, ymax = 450)


#
#
#
#
#

# Idea 2 - trait-elevation relationships broken out by species
p1 = ggplot(traits_chem_final, aes(x = elevation, y = ldmc, color = taxon)) +
  geom_point() +
  geom_smooth(method = "lm") +
  my_theme +
  theme(legend.position = "none")

p2 = ggplot(traits_chem_final, aes(x = elevation, y = leaf_thickness_mm, color = taxon)) +
  geom_point() +
  geom_smooth(method = "lm") +
  my_theme +
  theme(legend.position = "none")

p3 = ggplot(traits_chem_final, aes(x = elevation, y = sla_cm2_g, color = taxon)) +
  geom_point() +
  geom_smooth(method = "lm") +
  my_theme +
  theme(legend.position = "none")

p4 = ggplot(traits_chem_final, aes(x = elevation, y = plant_height_cm, color = taxon)) +
  geom_point() +
  geom_smooth(method = "lm") +
  my_theme +
  theme(legend.position = "none")

grid.arrange(p1,p2,p3,p4,ncol=2)



#
#
#
#
#

# Idea 3 - traitCV-elevation relationships broken out by species
# Function to compute coefficient of variation
cv = function(x) { return(sd(x)/mean(x)) }

# Select traits of interest
which_traits = c("plant_height_cm", "wet_mass_g", "dry_mass_g", "leaf_area_cm2", 
                 "sla_cm2_g", "ldmc", "leaf_thickness_mm", "carbon_derived",
                 "c_n_derived", "nitrogen_derived", "phosphorus_derived", 
                 "c_p_derived", "delta_15N", "delta_13C")

# Compute cv for each trait
traits_long_with_cv = traits_chem_final %>% 
  gather(trait_name, trait_value, which_traits) %>% 
  group_by(site, plot_id, individual_nr, ind_uid_new, functional_group, family, 
           taxon, genus, species, elevation, plot_uid, individual_uid, trait_name) %>% 
  summarize(traitCV = cv(trait_value))

traits_data_cut = subset(traits_long_with_cv, traitCV != 0)

# Broken out by species
ggplot(traits_data_cut, aes(x = elevation, y = traitCV, color = taxon)) +
  geom_point() +
  geom_smooth(method = "lm") +
  my_theme #+
  #theme(legend.position = "none")

# Broken out by trait
ggplot(traits_data_cut, aes(x = elevation, y = traitCV, color = trait_name)) +
  geom_jitter(width = 20) +
  geom_smooth(method = "lm") +
  my_theme +
  ylim(0,0.6)#+
  #theme(legend.position = "none")



# Attempt to plot both trait and cv trait together?
ggplot(subset(traits, trait == "sla_cm2_g"), aes(x = elevation, y = value)) +
  geom_jitter(width = 20) +
  geom_smooth(method = "lm") +
  geom_smooth(method = "lm", data = subset(traits_data_cut, trait_name == "sla_cm2_g"), aes(x = elevation, y = traitCV), lty = 2, color = "black")
  my_theme #+
  #ylim(0,0.6)#+

#
#
#
#
#

# Hypothesis 2 - structural leaf traits show less intraspecific variation than 
# physiological and chemical traits
# Hypothesis 3 - the functional groups will show different amounts of intraindividual
# and intraspecific variability with more variation in the forbs, less in graminoids, 
# and least in woody species; and
# Hypothesis 4 -  there is relatively more intraindividual and intraspecific 
# variability at the lower elevation populations than at the higher elevation.

# I believe these can be evaluated in a common model

# Set as factors
traits_long_with_cv$taxon = as.factor(traits_long_with_cv$taxon)
traits_long_with_cv$functional_group = as.factor(traits_long_with_cv$functional_group)

# Classify traits as structural or chemical
traits_long_with_cv = merge(traits_long_with_cv, 
      data.frame(trait_name = as.factor(which_traits), 
                 trait_type = as.factor(c(rep("structural", 7), rep("chemical", 7)))),
      by = "trait_name")


# Start with a model like this, determine best model by progressive elimination of terms
# Can't do this yet; need full chemical traits to determine best model
z1 = lmer(traitCV ~ trait_type*elevation*functional_group + (1|trait_name) + (1|taxon) + (1|plot_uid/ind_uid_new),
          REML = F,
          data = traits_long_with_cv)

# Note that here I have computed CV for intra-individual, not intraspecific trait
# variation - not sure how we would want to compute intraspecific CV? We will
# likely need to build a separate model for that

# Could we use stepAIC or similar for model selection?

