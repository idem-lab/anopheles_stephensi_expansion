# plot rasters for presentation
library(terra)
library(ggplot2)
library(tidyterra)
library(ggtext)

climatic_rel_abund <- rast("output/rasters/derived/climatic_rel_abund.tif")
potential_distribution <- rast("output/rasters/derived/predicted_potential_distribution.tif")
realised_distributions_all <- rast("output/rasters/derived/predicted_occupancy.tif")
larval_covs <- rast("output/rasters/derived/larval_habitat_covariates.tif")

# climatic rel_abund doesn't have the extra region beyond the target region
mask <- climatic_rel_abund * 0

# resample the monthly predictions of climatic suitability, for masking and
# plotting in the same way as the others
climatic_rel_abund_monthly <- rast("output/rasters/derived/An_stephensi_mechanistic_abundance.tif")
climatic_rel_abund_monthly <- resample(climatic_rel_abund_monthly,
                                       mask,
                                       method = "bilinear")

# mask and crop for plotting
climatic_rel_abund <- mask(climatic_rel_abund, mask)
climatic_rel_abund_monthly <- mask(climatic_rel_abund_monthly, mask)
potential_distribution <- mask(potential_distribution, mask)
realised_distributions_all <- mask(realised_distributions_all, mask)
larval_covs <- mask(larval_covs, mask)

plot_ext <- ext(-19.8, 103.7, -38.2, 41)
climatic_rel_abund <- crop(climatic_rel_abund, plot_ext)
climatic_rel_abund_monthly <- crop(climatic_rel_abund_monthly, plot_ext)
potential_distribution <- crop(potential_distribution, plot_ext)
realised_distributions_all <- crop(realised_distributions_all, plot_ext)
larval_covs <- crop(larval_covs, plot_ext)

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
    panel.grid.minor = element_blank(),
    plot.title = element_markdown(),
    plot.subtitle = element_markdown()
  ) +
  ggtitle("Predicted potential distribution of *An. stephensi*") +
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
       height = 5.5,
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
    panel.grid.minor = element_blank(),
    plot.title = element_markdown(),
    plot.subtitle = element_markdown()
  ) +
  ggtitle("Predicted potential distribution of *An. stephensi*") +
  xlab("") +
  ylab("")

ggsave("figures/An_stephensi_potential_distribution_nopoints.png",
       bg = "white",
       width = 9,
       height = 5.5,
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
    panel.grid.minor = element_blank(),
    plot.title = element_markdown(),
    plot.subtitle = element_markdown()
  ) +
  facet_wrap(~lyr) +
  ggtitle("Predicted realised distribution of *An. stephensi*")

ggsave("figures/An_stephensi_realised_distribution.png",
       bg = "white",
       width = 12,
       height = 3,
       dpi = 600)

# plot detection covariates
as_detection_density <- rast("output/rasters/derived/an_stephensi_detection_density.tif")

ggplot() +
  geom_spatraster(
    data = as_detection_density[[c("2000", "2010", "2022")]],
  ) +
  scale_fill_binned(
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
    panel.grid.minor = element_blank(),
    plot.title = element_markdown(),
    plot.subtitle = element_markdown()
  ) +
  facet_wrap(~lyr) +
  ggtitle("Increase in likelihood of detection")

ggsave("figures/An_stephensi_detection_effort.png",
       bg = "white",
       width = 12,
       height = 3,
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
    panel.grid.minor = element_blank(),
    plot.title = element_markdown(),
    plot.subtitle = element_markdown()
  ) +
  facet_wrap(~lyr) +
  ggtitle("Larval habitat covariates")

ggsave("figures/An_stephensi_larval_covariates.png",
       bg = "white",
       width = 12,
       height = 3,
       dpi = 600)

# climatic relative abundance
ggplot() +
  geom_spatraster(
    data = climatic_rel_abund,
    maxcell = Inf
  ) +
  scale_fill_distiller(
    palette = "YlGnBu",
    direction = 1,
    na.value = "white",
    limits = c(0, 1)
  ) +
  labs(fill = "Relative carrying capacity") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_markdown(),
    plot.subtitle = element_markdown()
  ) +
  ggtitle("Predicted annual climatic suitability for *An. stephensi*",
          "Relative carrying capacity per larval habitat, if present") +
  xlab("") +
  ylab("")

ggsave("figures/An_stephensi_mechanistic_suitability.png",
       bg = "white",
       width = 9,
       height = 5.5,
       dpi = 600)


ggplot() +
  geom_spatraster(
    data = climatic_rel_abund_monthly,
    maxcell = Inf
  ) +
  scale_fill_distiller(
    palette = "YlGnBu",
    direction = 1,
    na.value = "white",
    limits = c(0, 1)
  ) +
  labs(fill = "Relative\ncarrying\ncapacity") +
  facet_wrap(~lyr, ncol = 3) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_markdown(),
    plot.subtitle = element_markdown()
  ) +
  ggtitle("Predicted monthly climatic suitability for *An. stephensi*",
          "Relative carrying capacity per larval habitat, if present") +
  xlab("") +
  ylab("")

ggsave("figures/An_stephensi_monthly_mechanistic_suitability.png",
       bg = "white",
       width = 12,
       height = 12,
       dpi = 600)
