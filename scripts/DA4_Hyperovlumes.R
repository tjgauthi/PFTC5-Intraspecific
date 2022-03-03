#' ####################################################################### #
#' PROJECT: [DA.4 - Hypervolumes] 
#' CONTENTS: 
#'  - Calculating hypervolumes and comparing their sizes
#' AUTHOR: [Erik Kusch]
#' ####################################################################### #

# PREAMBLE ================================================================
rm(list=ls())

## Directories ------------------------------------------------------------


## Packages ---------------------------------------------------------------
install.load.package <- function(x) {
  if (!require(x, character.only = TRUE))
    install.packages(x, repos='http://cran.us.r-project.org')
  require(x, character.only = TRUE)
}
package_vec <- c(
  "hypervolume"# names of the packages required placed here as character objects
)
sapply(package_vec, install.load.package)

## Functionality ----------------------------------------------------------
`%nin%` <- Negate(`%in%`) # a function for negation of %in% function 

# DATA LOADING ============================================================
source("scripts/0_data_import.R") # sourcing data import script
traits_wide$ID <- paste(traits_wide$individual_uid, traits_wide$taxon, sep="_")
traits_df <- na.omit(traits_wide[ , c(-2:-4, -16:-18)])



# ANALYSES ================================================================
FUN.Hypervolumes <- function(data = traits_df, # the data frame with traits and the grouping column
                             TraitCols = 6:12, # which columns of the data frame contain the trait values
                             Grouping = "family" # by which column to establish groups
){
  data_scale <- scale(data[,TraitCols])
  data_ls <- split(data_scale, data[,Grouping])
  hv_ls <- lapply(data_ls, hypervolume)
  return(hv_ls)
  # volumes_ls <- lapply(hv_ls, get_volume)
  # return(volumes_ls)
}


## Creating Hypervolumes --------------------------------------------------
if(!file.exists(file.path("./data", "Hypervolumes.RData"))){
  
  ### Individual Hypervolumes ####
  ID_hv <- FUN.Hypervolumes(Grouping = "ID")
  ID_vols <- lapply(ID_hv, get_volume)
  
  ### Species Hypervolumes ####
  taxon_hv <- FUN.Hypervolumes(Grouping = "taxon")
  taxon_vols <- lapply(taxon_hv, get_volume)
  
  ### Functional Group Hypervolumes ####
  functional_group_hv <- FUN.Hypervolumes(Grouping = "functional_group")
  functional_group_vols <- lapply(functional_group_hv, get_volume)
  
  ### Saving Hypervolume Lists ####
  save(ID_hv, ID_vols, 
       taxon_hv, taxon_vols, 
       functional_group_hv, functional_group_vols,
       file = file.path("./data", "Hypervolumes.RData"))
}else{
  load(file.path("./data", "Hypervolumes.RData"))
}

## Comparison of Volumes --------------------------------------------------
Vols_df <- data.frame(values = unlist(ID_vols),
                      ID = names(unlist(ID_vols))
)
Vols_df$ID <- gsub(Vols_df$ID, pattern = ".untitled", replacement = "")
plot_df <- plyr::join(x = Vols_df, y = traits_df, by = "ID")

ggplot(plot_df, aes(y = values, x = taxon)) + 
  geom_boxplot() + 
  theme_bw()

ggplot(plot_df, aes(y = values, x = functional_group)) + 
  geom_boxplot() + 
  theme_bw()
