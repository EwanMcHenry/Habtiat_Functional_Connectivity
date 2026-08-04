### select and mask lcm ----

load(paste0(func.conect.path, 
            "\\analysis outputs\\", this.tss[this.ts.num], "\\r_curated data_.RData"))
load(paste0(func.conect.path, 
            "\\analysis outputs\\", this.tss[this.ts.num], "\\r_curated_data.awi_nwss_roads.RData"))


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
names(lcm.landscape) <- paste0("lcm")
trouble_plot(lcm.landscape, 
             name = paste0("lcm_with_roads_", this.tss[this.ts.num], "_",this.year ), 
             do_troubleshooting  = T,
             out_path = paste0(func.conect.path,
                               "\\analysis outputs\\", 
                               this.tss[this.ts.num],
                               "\\", this.year, "\\lcmres", constants$lcm.res)
)


# write roads raster
# writeRaster(roads.rast, paste0(gis.wd, "\\Data\\Roads\\road_widths_sample.tif"), overwrite = TRUE)

rm(lcm,
   roads.landscape, roads.rast)
gc()


print("LCM curation (script03.1) done")

# sort objects in memory by size

# sort(
#   sapply(ls(envir = .GlobalEnv), function(x)
#     object.size(get(x, envir = .GlobalEnv))),
#   decreasing = TRUE
# )