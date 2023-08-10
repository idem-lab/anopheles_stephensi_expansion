library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)

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
  filter(value >= 0) # fixing dodgy value


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
  seasonality  #%>%
    #filter(country != "equatorial guinea") # all zeroes screws with 
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
