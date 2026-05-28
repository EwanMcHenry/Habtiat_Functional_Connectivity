### select and mask lcm ----
## select years and corresponding index ----
lcm <- load_lcm_year(
    lcm.directs = lcm.directs,
    year = this.year,
    country = this.country,
    resolution = constants$lcm.res
) 

lcm.landscape <- mask_lcm_landscape_year(
  lcm = lcm,
  landscape = ts.buff,
  country = this.country
)

# add roads
roads.rast <- terra::rasterize(
  terra::vect(roads.landscape),
  lcm.landscape,   # template raster
  field = 1,
  background = 0,
  touches = TRUE
)

lcm.landscape[roads.rast == 1] <- 22 # assign roads a unique value in lcm to be able to identify them later

# write roads raster
# writeRaster(roads.rast, paste0(gis.wd, "\\Data\\Roads\\road_widths_sample.tif"), overwrite = TRUE)

rm(lcm)
gc()


