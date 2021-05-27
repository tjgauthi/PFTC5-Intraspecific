# Task DC.1 Add scaffolded chemical traits
# Retrieve chemical traits from online database as placeholders

### Call source script----

source(here::here(path = "scripts/0_data_import.R"))

### Packages ----
library(BIEN)
library(tidyverse)

### Prepping dataset ----

# Species and genus level only retrieval doesn't work, as our species don't have data
BIEN_species_traits_species <- BIEN_trait_species(c("Halenia umbellata",
                                                 "Lachemilla orbiculata",
                                                 "Gaultheria glomerata",
                                                 "Paspalum bonplandianum",
                                                 "Rhynchospora macrochaeta",
                                                 "Vaccinium floribundum"))
# Higher level instead
traits_chem <- traits_wide %>% 
  separate(taxon,
           into = c("genus", "species"),
           sep = "\\s", 2,
           remove = FALSE,
           # keep authority and species together
           extra = "merge")

genus_list <- traits_chem %>% 
  distinct(genus) %>% 
  pull(genus)


# Retrieve records from BIEN and filter to relevant
BIEN_chem_traits_genus <- BIEN_trait_genus(genus_list)

BIEN_chem_traits_genus <- BIEN_chem_traits_genus %>% 
  filter(trait_name == c("leaf nitrogen content per leaf dry mass",
           "leaf phosphorus content per leaf dry mass", 
           "leaf carbon content per leaf nitrogen content",
           "leaf carbon content per leaf dry mass"))
# Not many records for each genus other than Vaccinium, so shift to family level

# Get family
family_list <- traits_chem %>% 
  distinct(family) %>% 
  pull(family)


# Retrieve records from BIEN
BIEN_chem_traits_family <- BIEN_trait_traitbyfamily(family_list,
                    c("leaf nitrogen content per leaf dry mass",
                      "leaf phosphorus content per leaf dry mass", 
                      "leaf carbon content per leaf nitrogen content",
                      "leaf carbon content per leaf dry mass"))

test <- anti_join(traits_chem %>% 
            distinct(family),
          BIEN_chem_traits_family,
          by = c("family" = "scrubbed_family"))

#pulls species we have matches for
spp_matches =
  inner_join(traits_chem %>%
               ungroup() %>%
               distinct(taxon),
             BIEN_chem_traits_family,
             by = c('taxon' = 'scrubbed_species_binomial')) %>%
  select(-c(url_source,project_pi, project_pi_contact, access, id, unit, method)) %>%
  mutate(match_level = rep("species", nrow(.))) %>%
  #rename to make joining easier later
  rename(family = scrubbed_family,
         genus = scrubbed_genus)

#pulls genera we have matches for
genus_matches = 
  inner_join(traits_chem %>%
               distinct(genus),
             BIEN_chem_traits_family,
             by = c('genus' = 'scrubbed_genus')) %>%
  select(-c(url_source,project_pi, project_pi_contact, access, id, unit, method)) %>%
  mutate(match_level = rep("genus", nrow(.))) %>%
  #removes traits that we have at spp level
  anti_join(.,
            spp_matches,
            by = c("scrubbed_family" = "family", "trait_name")) %>%
  #rename to make joining easier later
  rename(family = scrubbed_family,
         taxon = scrubbed_species_binomial)

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

# Merge back to traits_chem with pivot_longer first
traits_chem_long <- traits_chem %>% 
  pivot_longer(cols = plant_height_cm:leaf_thickness_mm,
               names_to = "trait",
               values_to = "value")

testing <- rbind(
  left_join(traits_chem_long %>% 
              ungroup() %>% 
              select(-c(trait,value)),
            family_matches %>% 
              group_by(trait_name, family) %>% 
              summarise(mean_trait = mean(as.numeric(trait_value))) %>% 
              ungroup() %>% 
              pivot_wider(id_cols = -c(trait_name, mean_trait),
                          names_from = trait_name,
                          values_from = mean_trait,
                          values_fill = NA),
            by = 'family')) %>% 
  rename(c_n = 'leaf carbon content per leaf nitrogen content',
         nitrogen = 'leaf nitrogen content per leaf dry mass',
         phosphorus = 'leaf phosphorus content per leaf dry mass',
         carbon = 'leaf carbon content per leaf dry mass') %>% 
  pivot_longer(cols = c(c_n, nitrogen, phosphorus, carbon),
               names_to = 'trait',
               values_to = 'value')

# Summarise trait values for each of the six families
traits_chem_family_summary <- family_matches %>% 
  group_by(trait_name, family) %>% 
  summarise(mean_trait = mean(as.numeric(trait_value))) %>% 

# Relabel attributes
traits_chem_family_summary$trait_name[traits_chem_family_summary$trait_name == "leaf carbon content per leaf nitrogen content"] <- "c_n"
traits_chem_family_summary$trait_name[traits_chem_family_summary$trait_name == "leaf nitrogen content per leaf dry mass"] <- "nitrogen"
traits_chem_family_summary$trait_name[traits_chem_family_summary$trait_name == "leaf phosphorus content per leaf dry mass"] <- "phosphorus"
traits_chem_family_summary$trait_name[traits_chem_family_summary$trait_name == "leaf carbon content per leaf dry mass"] <- "carbon"

# Pivot wider to match traits_chem format                                        
traits_chem_family_summary_wide <- traits_chem_family_summary %>% 
  pivot_wider(names_from = trait_name,
            values_from = c(mean_trait))


# Merge Datasets
traits_chem_final <- left_join(traits_chem, traits_chem_family_summary_wide,
                       by = "family")            
