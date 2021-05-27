# Task DC.2 Add scaffolded chemical traits
# Retrieve chemical traits from online database as placeholders

# We currently don't have any trait data for: P, B, C, and 15N, 13C..
# So to scaffold all the models, we'll take data from related species from BIEN

### Call source script----

source(here::here(path = "scripts/0_data_import.R"))

### Packages ----
library(BIEN)
library(tidyverse)

### Prepping dataset ----

# Species and genus level only retrieval doesn't work, as our species don't have data

# Higher level instead
traits_chem <- traits_wide %>% 
  separate(taxon,
           into = c("genus", "species"),
           sep = "\\s", 2,
           remove = FALSE,
           # keep authority and species together
           extra = "merge")


# Get family
family_list <- traits_chem %>% 
  distinct(family) %>% 
  pull(family)


# Retrieve relevant records from BIEN (NB. Not all traits there, i.e. delta 13C/delta N15, C:P ratio)
BIEN_chem_traits_family <- BIEN_trait_traitbyfamily(family_list,
                    c("leaf nitrogen content per leaf dry mass",
                      "leaf phosphorus content per leaf dry mass", 
                      "leaf carbon content per leaf nitrogen content",
                      "leaf carbon content per leaf dry mass"))

# Pulls families we have matches for
family_matches = 
  inner_join(traits_chem %>%
               distinct(family),
             BIEN_chem_traits_family,
             by = c('family' = 'scrubbed_family')) %>%
  select(-c(url_source,project_pi, project_pi_contact, access, id, unit, method)) %>%
  mutate(match_level = rep("family", nrow(.))) %>%
  #rename to make joining easier later
  rename(taxon = scrubbed_species_binomial,
         genus = scrubbed_genus)


# Summarise trait values for each of the six families
traits_chem_family_summary <- family_matches %>% 
  group_by(trait_name, family) %>% 
  summarise(mean_trait = mean(as.numeric(trait_value)))

# Relabel attributes
traits_chem_family_summary$trait_name[traits_chem_family_summary$trait_name == "leaf carbon content per leaf nitrogen content"] <- "c_n_derived"
traits_chem_family_summary$trait_name[traits_chem_family_summary$trait_name == "leaf nitrogen content per leaf dry mass"] <- "nitrogen_derived"
traits_chem_family_summary$trait_name[traits_chem_family_summary$trait_name == "leaf phosphorus content per leaf dry mass"] <- "phosphorus_derived"
traits_chem_family_summary$trait_name[traits_chem_family_summary$trait_name == "leaf carbon content per leaf dry mass"] <- "carbon_derived"

# Pivot wider to match traits_chem format                                        
traits_chem_family_summary_wide <- traits_chem_family_summary %>% 
  pivot_wider(names_from = trait_name,
            values_from = c(mean_trait))

# Add dummy variables for those that aren't on BIEN
# C:P, Delta_N15 and Delta_13C
traits_chem_family_summary_wide <- traits_chem_family_summary_wide %>% 
  mutate(c_p_derived = rnorm(n(), 10,3),
         delta_15N = rnorm(n(),20,3),
         delta_13C = rnorm(n(),30,3))


# Merge Datasets
traits_chem_final <- left_join(traits_chem, traits_chem_family_summary_wide,
                       by = "family")            

# To be integrated into the 0 cleaning script until we get our own data
# END ----