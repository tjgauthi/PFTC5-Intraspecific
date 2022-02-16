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



#
#
#
#
#

# Idea 1 - trait-trait relationships at the overall, species, and individual level

# Generate plot showing overall, species, and individual regressions
p1 = ggplot(traits_chem_final, aes(x=ldmc, y=leaf_thickness_mm, color=taxon)) +
  #geom_smooth(method="lm") + 
  #geom_smooth(method="lm", aes(color=NULL)) +
  my_theme +
  ylim(0,0.625) +
  geom_point(size=0.9) +
  geom_smooth(method="lm", se=F,lwd=0.4, color="darkgrey", aes(group=ind_uid_new)) +
  geom_smooth(method="lm", se=F, lwd=0.6, aes(color=taxon)) +
  geom_smooth(method="lm", se=F, color="black", lty=2) #+
  #theme(legend.position = "none")

# Make inset plot showing distribution of slopes for individuals
trait.lm.list = nlme::lmList(leaf_thickness_mm ~ ldmc | ind_uid_new, data = traits_chem_final, na.action = na.exclude)
model_coef = subset(coef(trait.lm.list), !is.na(ldmc))

z = lm(leaf_thickness_mm ~ ldmc, traits_chem_final)
species.lm.list = nlme::lmList(leaf_thickness_mm ~ ldmc | taxon, data = traits_chem_final, na.action = na.exclude)
species.slopes = data.frame(coef(species.lm.list))
species.slopes$taxon = rownames(species.slopes)

inset_plot = ggplotGrob(
  ggplot(data = model_coef, aes(x = ldmc)) + 
    geom_density() + 
    theme_classic() +
    geom_vline(xintercept = 0, lty = 3) +
    geom_vline(xintercept = coef(z)[2], lty = 2) +
    geom_vline(data = species.slopes, aes(xintercept = ldmc, color=taxon)) +
    xlim(c(-3,3)) +
    xlab("Slope") +
    #ylim(c(0,22)) +
    theme(legend.position = "none") +
    theme(axis.title = element_text( size=8),
          rect = element_rect(fill = "transparent"),
          plot.background = element_rect(colour = "transparent") 
          )
)

p1 + annotation_custom(grob = inset_plot, xmin = 0, xmax = 0.3, ymin = 0.45, ymax = 0.63)



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
  geom_point() +
  geom_smooth(method = "lm") +
  my_theme #+
  #theme(legend.position = "none")



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

