library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)

# some notes on this seasonality data
# - some entries will be NA because of inability to extract data from paper, but
# record of paper remains in dataset
# - some entries may appear repeated if multiple sites were mentioned in text
# however combined data were presented in the paper figs/text
# - definitions of counts will vary among studies even for the same metric
# - conversion to monthly will vary also by study, in that some studies report
# a monthly metric which is derived as the reported average for the month,
# others might be a reported daily rate which is multiplied to convert to
# monthly

# conclusion:
# it is probable that the best use of these data will be to scale them from
# zero such that max = 1 and all other values = x/max - see graphs at bottom



seasonality_raw <- read_xlsx(
  path = "data/seasonality/seasonality.xlsx",
  sheet = "seasonality"
)


mons <- toupper(month.abb)

seasonality <- seasonality_raw %>%
  filter(
    grepl(
      pattern = "vector density",
      x = METRIC
    )
  ) %>%
  select(
    id = ID,
    year_start = STUDY_YEAR_START,
    year_end = STUDY_YEAR_END,
    country = ADMIN0,
    iso = ISO,
    admin0 = ADMIN0,
    admin1 = ADMIN1,
    admin2 = ADMIN2,
    admin3 = ADMIN3,
    lat = LAT,
    lon = LONG,
    sp = MOSQUITO_SPECIES,
    metric = METRIC,
    all_of(mons)
  ) %>%
  mutate(
    lat = as.numeric(lat),
    lon = as.numeric(lon),
    rowid = row_number()
  ) %>%
  pivot_longer(
    cols = all_of(mons),
    names_to = "month"
  ) %>%
  mutate(
    metric = sub(
      pattern = " vector density",
      replacement = "",
      x = .$metric
    ) %>%
      sub(
        pattern = "vector density ",
        replacement = "",
        x = .
      ) %>%
      gsub(
        pattern = "\\(",
        replacement = "",
        x = .
      ) %>%
      gsub(
        pattern = "\\)",
        replacement = "",
        x = .
      )
  ) %>%
  filter(value != "na") %>%
  mutate(
    value = as.numeric(value),
    month = factor(
      month,
      levels = mons
    ),
    monthnum = as.numeric(month)
  ) %>%
  filter(
    metric %in% c(
      "no. mosquitoes/trap/month",
      "no. mosquitoes/trap/month - indoor",
      #"aquatic habitats containing larvae",
      #"average no. larvae/dip/month",
      "average no. mosquitoes/house/month",
      "no. mosquitoes/house/month",
      "no. mosquitoes/man/month",
      "no. mosquitoes/month",
      "no. mosquitoes/night",
      "no. mosquitoes/trap/month - outdoor",
      #"aquatic habitats containing larvae %",
      "average no. mosquitoes/month",
      "average no. mosquitoes/man/month"
    )
  ) %>%
  filter(value >= 0) %>% # fixing dodgy value
  filter(country != "equatorial guinea") %>% # all zeroes
  group_by(rowid) %>%
  mutate(
    scaled_val = value / max(value)
  )


ggplot(seasonality) +
  geom_point(
    aes(
      x = month,
      y = value,
      colour = country
    )
  ) +
  scale_y_sqrt()

ggplot(
  seasonality
) +
  geom_point(
    aes(
      x = monthnum,
      y = value,
      colour = country
    )
  ) +
  scale_y_sqrt() +
  geom_smooth(
    aes(
      x = monthnum,
      y = value,
      colour = country
    )
  ) +
  facet_wrap(
    ~ country,
    scales = "free"
  )




ggplot(
  seasonality
) +
  geom_point(
    aes(
      x = monthnum,
      y = value
    )
  ) +
  scale_y_sqrt() +
  geom_smooth(
    aes(
      x = monthnum,
      y = value
    ),
    colour = "red",
    method = "loess"
  )  +
  facet_wrap(
    ~ country,
    scales = "free"
  ) #+
  # geom_smooth(
  #   data = seasonality %>%
  #     filter(country != "equatorial guinea", value > 0),
  #   aes(
  #     x = monthnum,
  #     y = value
  #   ),
  #   colour = "orange",
  #   method = "loess"
  # )


# tidy up data set, include trap type, explore all-zero instances


ggplot(seasonality) +
  geom_point(
    aes(
      x = month,
      y = scaled_val,
      colour = country
    )
  )

ggplot(
  seasonality
) +
  geom_point(
    aes(
      x = monthnum,
      y = scaled_val,
      colour = country
    )
  ) +
  geom_smooth(
    aes(
      x = monthnum,
      y = scaled_val,
      colour = country
    ),
    method = "loess"
  ) +
  facet_wrap(
    ~ country,
    scales = "free"
  )




ggplot(
  seasonality
) +
  geom_point(
    aes(
      x = monthnum,
      y = scaled_val,
      colour = metric
    )
  ) +
  geom_smooth(
    aes(
      x = monthnum,
      y = scaled_val
    ),
    colour = "grey40",
    method = "loess"
  )  +
  facet_wrap(
    ~ country,
    scales = "free"
  )
