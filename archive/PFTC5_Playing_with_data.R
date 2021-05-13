#Data analysis- site ACJ-------
#1. Site ACJ plant height 
data_final %>% 
  filter(Site == "ACJ") %>%
  ggplot(aes(x = ID, y = Plant_Height_cm, col = Genus)) +
  geom_point()

#plot thikness color by Genus 
data_final %>% 
  filter(Site == "ACJ") %>%
  ggplot(aes(x = Leaf_Thickness_1_mm, y = Leaf_Thickness_2_mm, col = Genus)) +
  geom_point()+
  geom_smooth(method = "lm", se = FALSE)

data_final %>% 
  filter(Site == "ACJ") %>%
  ggplot(aes(x = Leaf_Thickness_1_mm, y = Leaf_Thickness_2_mm)) +
  geom_point()+
  geom_smooth(method = "lm", se = FALSE)

#Plot together- mean_Plant_Height_cm and mean_Plant_Length_cm 
data_final %>% 
  filter(Site == "ACJ") %>%
  group_by(Genus) %>%
  summarise(n = n(), 
            mean_Plant_Height_cm = mean(Plant_Height_cm, na.rm = TRUE),
            mean_Plant_Length_cm = mean(Plant_Length_cm, na.rm = TRUE)) %>% 
  pivot_longer(cols = c(n, mean_Plant_Height_cm, mean_Plant_Length_cm),
               names_to = c("measures"), values_to = "values") %>% 
  ggplot(aes(x = measures, y = values, fill = Genus)) + 
  geom_col(position = "dodge")

