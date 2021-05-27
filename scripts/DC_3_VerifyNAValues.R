#pulls out any rows with NA values in it
df_non_na <- traits_wide[rowSums(is.na(traits_wide)) > 0, ]