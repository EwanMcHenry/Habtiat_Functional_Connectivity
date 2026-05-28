# what are the road widths
library(U.utilities)
library(tidyverse)

# took random points along lines in qgis and measured the road there or there abouts
## avoided roundabouts and junctions - they were over represeted in the random points (one per line)
## where dueled measured to middle between



rds <- read.csv(paste0(gis.wd, "\\Data\\Roads\\road_widths_sample.csv"))

summary(glm(rds$road.width ~ rds$highway))

# mean of each
rds %>%
  group_by(highway) %>%
  filter(!is.na(road.width)) %>%
  summarise(mean.width = mean(road.width, na.rm = TRUE),
            sd.width = sd(road.width, na.rm = TRUE),
            buffer = 0.5 * mean.width,
            n = n())
