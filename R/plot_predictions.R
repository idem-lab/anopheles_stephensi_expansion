# plot rasters for presentation
library(terra)
library(ggplot2)
library(tidyterra)

potential_distribution <- rast("output/rasters/derived/predicted_potential_distribution.tif")

# and for An. gambiae
ggplot() +
  geom_spatraster(
    data = potential_distribution,
  ) +
  scale_fill_distiller(
    palette = "YlGnBu",
    direction = 1,
    na.value = grey(0.9)
  ) +
  labs(fill = "Probability of suitability") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  ggtitle("Predicted potential distribution of An. stephensi")

ggsave("figures/An_stephensi_potential_distribution.png",
       bg = "white",
       width = 7,
       height = 5)
