# define Uk wide hexgrid
hex.grid0 = st_make_grid(countries, c(constants$hexdist.h, constants$hexdist.v), what = "polygons", square = F)
hex.grid = st_sf(hex.grid0) %>%
  # add grid ID
  mutate(grid_id = 1:length(lengths(hex.grid0)))
rm(hex.grid0)
gc()


save(hex.grid, 
     file = paste0(func.conect.path, 
                   "\\analysis outputs\\r_curated_global_data_.RData")
)

print("script02.01 done")
