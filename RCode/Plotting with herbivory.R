library(ggplot2)
library(magrittr)
library(ggcorrplot)
arch = read.csv("drago_result_with_herbivory_new.csv")
arch1 = subset(arch, complete.cases(arch))
# arch1$richness_similarity = 100 - arch1$richness_similarity

# Assuming the data is stored in 'arch1'
# Scatter plot of 'sigma' vs 'avg_shannon'
arch1 %>%
  # dplyr::filter(alfa != 0.001 & alfa != 0.01) %>%
  ggplot(aes(x = as.factor(sigma_comp), y = avg_shannon, color = as.factor(alfa))) +
  geom_point() +
  facet_wrap(~as.factor(epsilon)) +
  labs(
    title = "Avg Shannon Index vs Sigma by Epsilon",
    x = "Sigma", 
    y = "Avg Shannon Index"
  ) +
  rcartocolor::scale_color_carto_d(palette = "Earth") +
  theme_minimal()

# Scatter plot of 'Average Shannon' VS 'epsilon' Faceted by 'Sigma'
ggplot(arch1, aes(x = as.factor(epsilon), y = avg_shannon)) +
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


model = lm(avg_shannon ~ sigma + epsilon + alfa + sigma_comp + assymetry, data = arch1)
model_with_interaction = lm(avg_shannon ~ (sigma + epsilon + alfa + sigma_comp + assymetry)^2, data = arch1)

summary(model)
summary(model_with_interaction)

# Plot residuals vs fitted values
plot(model_with_interaction, which = 1)

# Normal Q-Q plot
plot(model_with_interaction, which = 2)

# Scale-Location plot
plot(model_with_interaction, which = 3)

# Residuals vs Leverage
plot(model_with_interaction, which = 5)

# Check multicollinearity using Variance Inflation Factor (VIF)
library(car)
vif(model)

# Backward elimination based on AIC
model_step <- step(model_with_interaction, direction = "backward")

model_step



