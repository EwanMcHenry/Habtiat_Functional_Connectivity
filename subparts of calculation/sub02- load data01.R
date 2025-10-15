# Ewan McHenry
# functional connectivity metric dev
# script 02 - loading data

# data ----
# hab cost and edge etc
countries = st_read(paste0(gis.wd, "\\Data\\administrative boundaries\\Countries\\R.5countries.simp100m.shp"))%>% st_transform( 27700) %>% arrange(name)

# full lcms ----
## select years and corresponding index ----
t_year <- min(years.considered)
tplus1_year <- max(years.considered)

t_idx <- which(lcm.directs$year == t_year)
tplus1_idx <- which(lcm.directs$year == tplus1_year)

## load rasters ----
# Load GB rasters
lcm_t.rast25.gb <- raster(lcm.directs$gb.25[t_idx])
lcm_tplus1.rast25.gb <- raster(lcm.directs$gb.25[tplus1_idx])

# Load NI rasters
lcm_t.rast25.ni <- raster(lcm.directs$ni.25[t_idx])
lcm_tplus1.rast25.ni <- raster(lcm.directs$ni.25[tplus1_idx])

## replace 0 values with 13 for sea ----
lcm_tplus1.rast25.gb[lcm_tplus1.rast25.gb == 0] = 13 
lcm_tplus1.rast25.ni[lcm_tplus1.rast25.ni == 0] = 13 
lcm_t.rast25.gb[lcm_t.rast25.gb == 0] = 13 
lcm_t.rast25.ni[lcm_t.rast25.ni == 0] = 13 

# woodland info -----
nwss = st_read(nwss.dir) %>% st_transform( 27700)

awi = st_read(awi.dir) %>% st_transform( 27700)

### save ----
save(countries,
     lcm_tplus1.rast25.gb,
     lcm_t.rast25.gb,
     lcm_tplus1.rast25.ni,
     lcm_t.rast25.ni,
     lcm_t.rast25.ni,
     nwss,
     awi, 
     file = paste0(func.conect.path, 
                   "\\analysis outputs\\r_global_data_.RData")
)

print("script02 done")
