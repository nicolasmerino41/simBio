library(dplyr)
df = read.csv("C:/Users/MM-1/Downloads/nico_spp.csv")

df_ibe_pen = subset(df, X > -9.5 & X < 3.5)

df_ibe_pen = subset(df, Y > 36.0 & Y < 44)

average_densities = df_ibe_pen %>% group_by(names) %>% summarise(mean = mean(dens_km))
