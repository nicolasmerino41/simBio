library(ggplot2)
library(magrittr)
library(ggcorrplot)
arch = read.csv("drago_results_with_herbivory_new.csv")
arch1 = subset(arch, complete.cases(arch))
# arch1$richness_similarity = 100 - arch1$richness_similarity

# Assuming the data is stored in 'arch1'
# Scatter plot of 'sigma' vs 'avg_shannon'
arch1 %>%
  dplyr::filter(alfa != 0.001 & alfa != 0.01) %>%
  ggplot(aes(x = as.factor(sigma), y = avg_shannon, color = as.factor(epsilon))) +
  geom_point() +
  facet_wrap(~as.factor(epsilon)) +
  labs(
    title = "Sigma vs Avg Shannon Index",
    x = "Sigma", 
    y = "Avg Shannon Index"
  ) +
  rcartocolor::scale_color_carto_d(palette = "Earth") +
  theme_minimal()

# Scatter plot of 'sigma' vs 'avg_shannon', colored by 'Suitability Type'
ggplot(arch1, aes(x = sigma, y = avg_shannon, color = sigma_comp)) +
  geom_point() +
  facet_wrap(~sigma_comp) +
  labs(title = "Average Shannon Vs Sigma by Sigma Comp",
       x = "Sigma", y = "Average Shannon",
       color = "Sigma comp") +  # This changes the legend title
  theme_minimal()

# Scatter plot of 'epsilon' vs 'Average Shannon' faceted by 'Suitability Type'
ggplot(arch1, aes(x = epsilon, y = avg_shannon)) +
  geom_point() +
  facet_wrap(~sigma) +
  labs(title = "Average Shannon Vs Epsilon by Interaction Strength",
       x = "Epsilon", y = "Average Shannon") +
  theme_classic()

# Box plot of 'avg_bbp' by 'sigma'
ggplot(arch1, aes(x = as.character(sigma), y = avg_bbp)) +
  geom_point() +
  labs(title = "Avg BBP by sigma",
       x = "Sigma", y = "Avg BBP") +
  theme_minimal()

# Boxplot of 'richness_similarity' vs 'sigma' faceted by 'Suitability Type' (k_DA_name)
ggplot(arch1, aes(x = as.factor(sigma), y = richness_similarity)) +
  geom_boxplot() +
  facet_wrap(~ epsilon) +  # Facet by Suitability Type
  labs(title = "Richness Similarity vs Sigma by Epsilon",
       x = "Sigma", y = "Richness Similarity", 
       facet = "Epsilon") +
  theme_minimal()
# Boxplot of 'alive_predator' vs 'sigma' faceted by 'Suitability Type' (k_DA_name)

arch1 %>%
  dplyr::filter(alfa != 0.001 & alfa != 0.01) %>%
  ggplot(aes(x = as.factor(sigma), y = alive_predators, color = alfa)) +
  geom_point() +
  facet_wrap(~ epsilon) +  # Facet by Suitability Type
  labs(title = "Alive predators vs Sigma by Epsilon",
       x = "Sigma", y = "Alive predators", 
       facet = "Epsilon") +
  scale_color_carto_c(palette = "Earth") +
  theme_classic()


# CORRELATION
numeric_cols <- arch1[, c("sigma", "epsilon", "alfa", "avg_shannon", "avg_bbp", "richness_similarity", "alive_predators", "mean_tl")]
corr_matrix <- cor(numeric_cols, use = "complete.obs")

# Plot the lower triangle of the correlation matrix
ggcorrplot(corr_matrix, lab = TRUE, type = "upper", title = "Correlation Matrix")












