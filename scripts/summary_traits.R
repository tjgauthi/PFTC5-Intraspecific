library(dplyr)
library(ggplot2)

# Upload data
data <- read_csv("data/raw/PFTC3-Puna-PFTC5_Peru_2018-2020_LeafTraits_clean.csv")

# select data needed
data <- data %>%
  filter(site %in% c("WAY", "ACJ", "TRE") &
           year == 2020 & treatment %in% "C" &
           trait %in%c("plant_height_cm", "leaf_area_cm2", "sla_cm2_g", "ldmc") &
           functional_group %in% c("Forb", "Graminoid", "Woody"))

# order the site acccording to altitude
data$site <- factor(data$site, levels=c("WAY", "ACJ", "TRE"))

## summary traits plots
# 1
data %>% ggplot(aes(y = value, x = site, fill = site)) + 
  geom_boxplot() + 
  facet_wrap(~ trait, scales = "free") +
   scale_fill_grey()+ theme_classic() +
   theme(legend.position="bottom")

# 2
data %>% ggplot(aes(y = value, x = trait, fill = site)) + 
  geom_boxplot() + 
  facet_wrap(~ site, scales = "free") +
  scale_fill_grey()+ theme_classic() +
  theme(legend.position="bottom")

# 3   
data %>% ggplot(aes(y = value, x = site, fill = site)) + 
  geom_violin() + 
  facet_wrap(~ trait, scales = "free") +
  scale_fill_grey()+ 
  theme_classic() +
  theme(legend.position="bottom")

# 4 
data %>% ggplot(aes(y = value, x = functional_group, fill = functional_group)) + 
  geom_boxplot()+
  facet_wrap(~ trait, scales = "free")+
  scale_fill_grey()+ 
  theme_classic()+
  theme(legend.position="bottom")

# 5
data %>% ggplot(aes(y = value, x = functional_group, fill = functional_group))+ 
  geom_violin()+ 
  facet_wrap(~ trait, scales = "free")+
  scale_fill_grey()+ 
  theme_classic()+
  theme(legend.position="bottom")


data %>% ggplot(aes(y = value, x = family, fill = family)) + 
  geom_boxplot()+
  facet_wrap(~ trait, scales = "free")+
  scale_fill_grey()+ theme_classic()+
  theme(legend.position="bottom")



data %>% ggplot(aes(y = value, x = family, fill = family))+ geom_violin()+ facet_wrap(~ trait, scales = "free")+
  scale_fill_grey()+ theme_classic()+
  theme(legend.position="bottom")+theme(axis.text.x = element_text(angle = 90))

