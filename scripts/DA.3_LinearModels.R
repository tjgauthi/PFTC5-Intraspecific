library(lme4)
library(nlme)
library(gridExtra)

# look at means, how they change across all sites, how they change w/in species, etc.
# Think abt.that anderegg paper

# I did that quickly, it doesn't look very interesting

# Another idea is to look at trait-trait relationships w/in individuals, w/in species,
# and overall

my_theme = theme_bw() + 
  #theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

traits_chem_final$ind_uid_new = paste(traits_chem_final$site,
                                      traits_chem_final$plot_id,
                                      traits_chem_final$taxon,
                                      traits_chem_final$individual_nr) %>% as.factor()

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









# Something like this - this shows overall and species lines
ggplot(traits_chem_final, aes(x=ldmc, y=leaf_thickness_mm,color=ind_uid_new)) +
  geom_point() +
  geom_smooth(method="lm") + 
  geom_smooth(method="lm", aes(color=NULL))

# Function for computing coefficient of variation
cv = function(x) { return(sd(x)/mean(x)) }

# Compute CV for each trait for each individual
trait_cvs = traits %>% 
  group_by(site, plot_id, individual_nr, functional_group, family, taxon, trait, elevation) %>% 
  summarize(cv_trait = cv(value))

# Before we can use "individual" as an effect in our model, we need to assign a 
# unique individual number to each individual, otherwise all the "1"s will be
# treated as the same individual, etc.

# This assigns a unique number "unique_id" to each individual in the dataset
trait_cvs = trait_cvs %>% 
  ungroup() %>% 
  group_by(site, plot_id, individual_nr, taxon) %>% 
  mutate(unique_id = group_indices())
  
trait_cvs$taxon = as.factor(trait_cvs$taxon)
trait_cvs$functional_group = as.factor(trait_cvs$functional_group)


#z = lmer(traitCV ~ Species + FunctionalGroup + Site/Elevation + Individuals + (1|Species/Individuals) + (1|Site/Species))

z = lm(cv_trait ~ taxon + functional_group,
         data = subset(trait_cvs, trait == "ldmc"))
# I think we basically can't get an independent effect of functional group here because
# taxon is fully nested within functional group? Not sure.

# Basically I think this model is not going to work and we need to figure out what
# kinds of questions we are trying to ask with the model and then build a model
# which will answer that question. From the MS draft, it looks like we are trying to 
# test hypotheses 2, 3, and 4 with the LME models. These are:

# 2. structural leaf traits show less intraspecific variation than physiological and chemical traits
# For this, perhaps we need to classify traits as structural, physiological, and chemical
# We also need to make sure that whatever metric we use to compare intraspecific variation is 
# comparable between different traits, i.e. things with different units. I think the CV satisfies this.
#
# All the traits we have so far seem like structural traits to me - no chemical yet
# Do we even have any physiological traits? We could get some for paspalum and gaultheria

# Here, trait_cvs is measuring within-individual variation, but what we actually want is
# within-species. So we could take the CV of all the leaves of a species ... but just doing
# that washes out any differences between elevations, individuals, plots, etc.
# Maybe what you want here is an ANOVA/ANCOVA? bc what we want is to ask whether
# the groups differ from each other. If we just computed intraspecific variation (e.g.
# trait CV at the individual level), then we could do an ANOVA on that with the "trait
# type" as the factor.

# To do:
# Compute "intraspecific variation" for each trait (CV of all leaves, or of individual means)
# Classify traits as chemical, structural, and physiological
# --Add dummy traits as needed
# --Do we have any physiological traits?
# Then we basically want to construct a model:
# --Response is ITV
# --Fixed effect = trait_type (struct,physio,chem)
# --Random effects = species, (fn group?), site/elev (or elev as cov.), plot, trait? 


# 3. the functional groups will show different amounts of intraindividual and 
# intraspecific variability with more variation in the forbs, less in graminoids, 
# and least in woody species;

# So again, it seems like an ANOVA is appropriate here. Maybe we can use an ANCOVA
# with the elevation as a covariate to make a model which would also test hyp. 4.
# Another thing to consider here is a Tukey test, which would give us pairwise
# comparisons between group means.

# To do:
# Compute WITV and ITV - need two models here for two response variables
# Model:
# --Response = WITV or ITV
# --Fixed = Fn_group + elevation
# --Random = 

library(nlme)
z = lme(cv_trait ~ functional_group,
    data = trait_cvs,
    random = ~ 1|plot_id/site/unique_id + 1|trait,
    na.action = na.exclude)

z1 = lmer(cv_trait ~ functional_group + (1|trait) + (1|taxon) + (1|siteplot/unique_id),
         REML = F,
         data = trait_cvs)

# z2 is preferred
z2 = lmer(cv_trait ~ functional_group + (1|trait) +  (1|siteplot/unique_id),
          REML = F,
          data = trait_cvs)
# z2 is much preferred
z3 = lmer(cv_trait ~ functional_group + (1|siteplot/unique_id),
          REML = F,
          data = trait_cvs)

# I guess z2 is preferred?
z4 = lmer(cv_trait ~ functional_group + (1|trait) +  (1|plot_id/site/unique_id),
          REML = F,
          data = trait_cvs)
# z2 still
z5 = lmer(cv_trait ~ functional_group + (1|trait) +  (1|plot_id/site),
          REML = F,
          data = trait_cvs)

# z2 still; possibly z4?
z6 = lmer(cv_trait ~ functional_group + (1|trait) +  (1|unique_id),
          REML = F,
          data = trait_cvs)

# z7 wins I guess?
z7 = lmer(cv_trait ~ functional_group + elevation +(1|trait) +  (1|siteplot/unique_id),
          REML = F,
          data = trait_cvs)

# z8 wins I guess, or no diff so take the simpler one.
z8 = lmer(cv_trait ~ functional_group + elevation +(1|trait) +  (1|unique_id),
          REML = F,
          data = trait_cvs)

# Apparently z9 is better by both aic and bic, but no different according to p value
z9 = lmer(cv_trait ~ elevation +(1|trait) +  (1|unique_id),
          REML = F,
          data = trait_cvs)

z0 = lmer(cv_trait ~ functional_group + elevation +(1|trait) +  (1|unique_id),
          REML = T,
          data = trait_cvs)


# Could it just be as simple as this?
z = lm(cv_trait ~ functional_group, data = trait_cvs)
anova(z)
# I'm not sure what benefit we would gain from random effects here
# Based on this test, there is no difference by fn group in our data, so no
# need to go on to do a Tukey test.


# 4. there is relatively more intraindividual and intraspecific variability at 
# the lower elevation populations than at the higher elevation.

# We could evaluate this by computing traitCV (within-individual) and regressing
# it against elevation - this would get at the intRAspecific question. We could
# add appropriate random effects for species, and plot. Individual is 
# nested within species and plot, but plot and species are crossed I think, bc 
# each plot has all species. Something like that.

# This will also give us an estimate of within-species variation due to the random
# effects structure, but I think it won't let us determine how that variation is
# affected by elevation. I think. We need something like an interaction term for that
# but I don't think you can do that for random effects?








