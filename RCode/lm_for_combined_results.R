library(lm)
library(readr)

df = read.csv("C:/Users/MM-1/OneDrive/PhD/GitHub/simBio/combined_results.csv")
abs_df = subset(df, df$community_type != "Real")
model <- lm(log(resilience*-1) ~ conn * avg_degree + mu + eps + m_alpha +
              reactivity + prop_omn + total_species + avg_conn,
            data = abs_df)

# View summary of the model
summary(model)

