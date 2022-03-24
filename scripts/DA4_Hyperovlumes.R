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
  "hypervolume",
  "car",
  "multcomp"
)
sapply(package_vec, install.load.package)

## Functionality ----------------------------------------------------------
`%nin%` <- Negate(`%in%`) # a function for negation of %in% function 

# DATA LOADING ============================================================
source("scripts/0_data_import.R") # sourcing data import script
traits_wide$ID <- paste(traits_wide$individual_uid, traits_wide$taxon, sep="_")

# POOL ELEVATION WITHIN SITE. WAY: 3101; ACJ: 3468; TRE: 3715
traits_df <- na.omit(traits_wide[ , c(-2:-4, -16:-18)])
traits_df$elevation <- ifelse(traits_df$site == 'ACJ', 3468,
					   ifelse(traits_df$site == 'WAY', 3101,
					   ifelse(traits_df$site == 'TRE', 3715, NA)))

Grouping_df <- as.data.frame(table(traits_df[,c("taxon", "functional_group")]))
Grouping_df <- Grouping_df[Grouping_df$Freq != 0, -3]
Grouping_df$taxon <- as.character(Grouping_df$taxon)

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

FUN.Overlap <- function(data, # the list of hypervolumes from which to pull
                        names # the names of the hypervolumes which to compare
                        ){
  overlap <- hypervolume_set(data[[names[1]]], data[[names[2]]], 
                             check.memory = FALSE)
  overlap <- hypervolume::hypervolume_overlap_statistics(overlap)[1] # jaccard distance
  return(overlap)
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

# Data prep
Vols_df <- data.frame(values = unlist(ID_vols),
                      ID = names(unlist(ID_vols))
)
Vols_df$ID <- gsub(Vols_df$ID, pattern = ".untitled", replacement = "")
plot_df <- plyr::join(x = Vols_df, y = traits_df, by = "ID") %>%
			distinct(values, ID, .keep_all = TRUE) %>%
			filter(values != max(.$values)) # This filters out the outlier which is TRE_1_NA_Rhynchospora macrochaeta

plot_df$taxon <- factor(plot_df$taxon, levels = c('Halenia umbellata', 'Lachemilla orbiculata', 'Paspalum bonplandianum', 'Rhynchospora macrochaeta', 'Gaultheria glomerata', 'Vaccinium floribundum'))
plot_df$functional_group <- as.factor(plot_df$functional_group)

# Function to get post-hoc test stats
get_posthoc_stats <- function(model, group_var) {

	cmp <- do.call(mcp, setNames(list("Tukey"), group_var))

	stats <- glht(model, linfct = cmp) %>%
				summary(test = adjusted("bonferroni"))
				
	stats.df <- cbind('Estimate' = stats$test$coefficients,
					  'Std_Error' = stats$test$sigma,
					  'z_value' = stats$test$tstat,
					  'p_value' = stats$test$pvalues
					  ) %>%
					data.frame() %>%
					tibble::rownames_to_column(var = 'Comparison')
		
	return(stats.df)
	
}

# Function to get compact letter display of post-hoc groups
get_letters <- function(model, group_var) {

	cmp <- do.call(mcp, setNames(list("Tukey"), group_var))

	letters <- glht(model, linfct = cmp) %>%
					cld() %>%
					.[[10]] %>%
					.[[1]] %>%
					toupper() %>%
					as.list()
					
	return(letters)
	
}

# Taxon analysis and boxplot
taxon.lm <- lm(values ~ taxon, data = plot_df) # fit isn't great, but I've seen worse
taxon.results <- car::Anova(taxon.lm, type = 2) # Not sure what type we should be using

taxon.posthoc.stats <- get_posthoc_stats(model = taxon.lm, group_var = 'taxon')
taxon.posthoc.letters <- get_letters(model = taxon.lm, group_var = 'taxon')

taxon.boxplot <- ggplot(plot_df, aes(y = values, x = taxon, color = functional_group)) +
	geom_boxplot() + 
	stat_summary(geom = 'text',
		label = taxon.posthoc.letters,
		fun = max,
		vjust = -1
		) +
  theme_bw()
print(taxon.boxplot)

# Functional group analysis and boxplot
fg.lm <- lm(values ~ functional_group, data = plot_df) # fit isn't great, but I've seen worse
fg.results <- car::Anova(fg.lm, type = 2) # Not sure what type we should be using

fg.posthoc.stats <- get_posthoc_stats(model = fg.lm, group_var = 'functional_group')
fg.posthoc.letters <- get_letters(model = fg.lm, group_var = 'functional_group')

fg.boxplot <- ggplot(plot_df, aes(y = values, x = functional_group)) + 
  geom_boxplot() + 
  stat_summary(geom = 'text',
  	label = fg.posthoc.letters,
	fun = max, 
	vjust = -1 
	) +
  theme_bw()
print(fg.boxplot)

# Elevation analysis and scatterplot
elevation.lm <- lm(values ~ elevation, data = plot_df) # fit isn't great, but I've seen worse
elevation.results <- car::Anova(elevation.lm, type = 2) # Not sure what type we should be using

elevation.scatterplot <- ggplot(plot_df, aes(y = values, x = elevation, color = elevation)) + 
  geom_point() + 
  theme_bw() +
  annotate("text", label = "womp womp :(", x = 3200, y = 8)
print(elevation.scatterplot)

### bar plots ----------
plot_df <- data.frame(sp = c(Grouping_df$taxon[Grouping_df$functional_group == "Forb"], "Overlap",
                             Grouping_df$taxon[Grouping_df$functional_group == "Graminoid"], "Overlap",
                             Grouping_df$taxon[Grouping_df$functional_group == "Woody"], "Overlap"),
                      fg = rep(c("Forb", "Graminoid", "Woody"), each = 3),
                      value = as.numeric(
                        c(taxon_vols[Grouping_df$taxon[Grouping_df$functional_group == "Forb"]],
                          FUN.Overlap(data = taxon_hv,
                                      names = Grouping_df$taxon[Grouping_df$functional_group == "Forb"]),
                          taxon_vols[Grouping_df$taxon[Grouping_df$functional_group == "Graminoid"]],
                          FUN.Overlap(data = taxon_hv,
                                      names = Grouping_df$taxon[Grouping_df$functional_group == "Graminoid"]),
                          taxon_vols[Grouping_df$taxon[Grouping_df$functional_group == "Woody"]],
                          FUN.Overlap(data = taxon_hv,
                                      names = Grouping_df$taxon[Grouping_df$functional_group == "Woody"]))
                      ),
                      sp2 = rep(c("sp1", "sp2", "Overlap"), 3)
)
plot_df$value[1:2] <- plot_df$value[1:2]-plot_df$value[3]/2
plot_df$value[4:5] <- plot_df$value[4:5]-plot_df$value[6]/2
plot_df$value[7:8] <- plot_df$value[7:8]-plot_df$value[9]/2


plot_df$sp2 <- as.factor(plot_df$sp2)
plot_df$sp2 <- relevel(plot_df$sp2, 'sp2')

gplot <- ggplot(plot_df, aes(x = fg, y = value, fill = sp2, label = sp)) + 
  geom_bar(position="stack", stat="identity") + 
  geom_text(size = 3, position = position_stack(vjust = 0.6)) + 
  theme_bw() + theme(legend.position = "none") + labs(x = "Functional Group", y = "Hypervolume Size")
print(gplot)

### venn diagram ----
hv_list = hypervolume_join(functional_group_hv)

plot_df <- data.frame(Area = c("Forb", "Graminoid", "Woody",
                               "Forb & Graminoid",
                               "Forb & Woody",
                               "Graminoid & Woody",
                               "All"),
                      Value = as.numeric(c(
                        functional_group_vols,
                        FUN.Overlap(data = functional_group_hv, names = c("Forb", "Graminoid")),
                        FUN.Overlap(data = functional_group_hv, names = c("Forb", "Woody")),
                        FUN.Overlap(data = functional_group_hv, names = c("Graminoid", "Woody")),
                        get_volume(hypervolume_set_n_intersection(hv_list))
                      ))
)

# draw.triple.venn(plot_df$Value[1], plot_df$Value[2], plot_df$Value[3], 
#                  plot_df$Value[4], plot_df$Value[5], plot_df$Value[6], 0.001,
#                  c("Forb", "Graminoid", "Woody"), fill = c("red", "green", " blue"))

# myV3 <- createVennObj(nSets = nrow(plot_df), sSizes = plot_df$Value)
# myV3 <- plotVenn(nVennObj = myV3, nCycles = 10000, outFile = "test.png")

message("#############################################")
print("We now need to run some simple tests like Mann-Whitney U to assess significance of difference in volume size")