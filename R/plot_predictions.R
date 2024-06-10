# plot rasters for presentation
library(terra)
library(ggplot2)
library(tidyterra)

potential_distribution <- rast("output/rasters/derived/predicted_potential_distribution.tif")
realised_distributions_all <- rast("output/rasters/derived/predicted_occupancy.tif")
realised_distributions_mode_1 <- rast("output/rasters/derived/predicted_occupancy_mode_1.tif")
realised_distributions_mode_2 <- rast("output/rasters/derived/predicted_occupancy_mode_2.tif")
larval_covs <- rast("output/rasters/derived/larval_habitat_covariates.tif")

# load detection data for plotting

first_detection <- readRDS("output/tabular/first_detection.RDS") %>%
  filter(ever_detected == 1)

ggplot() +
  geom_spatraster(
    data = potential_distribution,
    maxcell = Inf
  ) +
  scale_fill_distiller(
    palette = "YlGnBu",
    direction = 1,
    na.value = "white",
    limits = c(0, 1)
  ) +
  labs(fill = "Suitability") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  ggtitle("Predicted potential distribution of An. stephensi") +
  geom_point(
    aes(
      x = x,
      y = y
    ),
    data = first_detection,
    pch = 21,
    colour = "black",
    fill = "red",
    # stroke = 0.8
  ) +
  xlab("") +
  ylab("")

ggsave("figures/An_stephensi_potential_distribution.png",
       bg = "white",
       width = 9,
       height = 6,
       dpi = 600)

ggplot() +
  geom_spatraster(
    data = potential_distribution,
    maxcell = Inf
  ) +
  scale_fill_distiller(
    palette = "YlGnBu",
    direction = 1,
    na.value = "white",
    limits = c(0, 1)
  ) +
  labs(fill = "Suitability") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  ggtitle("Predicted potential distribution of An. stephensi") +
  xlab("") +
  ylab("")

ggsave("figures/An_stephensi_potential_distribution_nopoints.png",
       bg = "white",
       width = 9,
       height = 6,
       dpi = 600)


years_keep <- c("2000", "2010", "2022")

# plot all realised distributions

ggplot() +
  geom_spatraster(
    data = realised_distributions_all[[years_keep]],
  ) +
  scale_fill_distiller(
    palette = "YlGnBu",
    direction = 1,
    na.value = "white"
  ) +
  labs(fill = "Probability of presence") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  facet_wrap(~lyr) +
  ggtitle("Predicted realised distribution of An. stephensi")

ggsave("figures/An_stephensi_realised_distribution.png",
       bg = "white",
       width = 12,
       height = 4,
       dpi = 600)

# now for the two modes

ggplot() +
  geom_spatraster(
    data = realised_distributions_mode_1[[years_keep]],
  ) +
  scale_fill_distiller(
    palette = "YlGnBu",
    direction = 1,
    na.value = "white"
  ) +
  labs(fill = "Probability of presence") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  facet_wrap(~lyr) +
  ggtitle("Predicted realised distribution of An. stephensi (mode 1)")

ggsave("figures/An_stephensi_realised_distribution_mode_1.png",
       bg = "white",
       width = 12,
       height = 4,
       dpi = 600)


ggplot() +
  geom_spatraster(
    data = realised_distributions_mode_2[[years_keep]],
  ) +
  scale_fill_distiller(
    palette = "YlGnBu",
    direction = 1,
    na.value = "white"
  ) +
  labs(fill = "Probability of presence") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  facet_wrap(~lyr) +
  ggtitle("Predicted realised distribution of An. stephensi (mode 1)")

ggsave("figures/An_stephensi_realised_distribution_mode_2.png",
       bg = "white",
       width = 12,
       height = 4,
       dpi = 600)

# plot detection covariates
as_detection_density <- rast("output/rasters/derived/an_stephensi_detection_density.tif")

ggplot() +
  geom_spatraster(
    data = as_detection_density[[c("2000", "2010", "2022")]],
  ) +
  scale_fill_binned(
    # guide = guide_coloursteps(
    #   direction = -1
    # ),
    # type = "gradient",
    low = "light blue",
    high = "#132B43",
    breaks = c(1e-4, 0.01, 0.025, 0.05),
    na.value = "white"
  ) +
  labs(fill = "Increased detection probability") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  facet_wrap(~lyr) +
  ggtitle("Increase in likelihood of detection")

ggsave("figures/An_stephensi_detection_effort.png",
       bg = "white",
       width = 12,
       height = 4,
       dpi = 600)

# plot larval habitat rasters too
larval_covs_keep <- c("built_volume", "tcb", "tcw")

ggplot() +
  geom_spatraster(
    data = larval_covs[[larval_covs_keep]],
  ) +
  scale_fill_distiller(
    palette = "PRGn",
    direction = 1,
    na.value = "white"
  ) +
  labs(fill = "Scaled covariate value") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  facet_wrap(~lyr) +
  ggtitle("Larval habitat covariates")

ggsave("figures/An_stephensi_larval_covariates.png",
       bg = "white",
       width = 12,
       height = 4,
       dpi = 600)
