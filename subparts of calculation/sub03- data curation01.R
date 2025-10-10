# Ewan McHenry
##------ Fri Feb 25 09:19:50 2022 ------##
# functional connectivity metric dev
# script 03 - data curation

# curation ----


## CUT DATA TO LANDSCAPE AND SAVE ----
### select buffered landscape to avoid edge effects  ----
this.ts = Focal_landscape[this.ts.num,]
ts.buff = this.ts %>% 
  st_simplify( preserveTopology = T, dTolerance = constants$landscape.buffer.simplification.tolerance[1]) %>% # first simplify hack to reduce run time. this is a rough buffer to negate edge effects, so can be v rough
  st_buffer( dist = constants$buffer.roundLandscape) %>% # buffer
  st_simplify( preserveTopology = FALSE, dTolerance = constants$landscape.buffer.simplification.tolerance[2]) # simplify again
### name of country ----
this.country = st_intersection( st_point_on_surface(this.ts) , countries)$name.1 

### select hex grid ----

# intersect hexgrid with landscape - note that grid id is from original UK-wide grid, allows easy cross-ID
tsbuff.hexgrid <- st_intersection(hex.grid, ts.buff) %>% 
  st_make_valid() %>%  st_cast("MULTIPOLYGON") %>% #st_cast("POLYGON") %>% 
  dplyr::select(grid_id)

ts.hexgrid <- st_intersection(hex.grid, this.ts) %>% 
  st_make_valid() %>%  st_cast("MULTIPOLYGON") %>% #st_cast("POLYGON") %>% 
  dplyr::select(grid_id) 
ts.hexgrid$hex.ha = st_area(ts.hexgrid) %>% 
  set_units(value = "ha") %>%
  as.numeric()

### select lcm and transform crs of other data if needed----
# different lcm for NI/GB
if(this.country == "N.Ireland"){
  
  lcm_tplus1.rast25 = lcm_tplus1.rast25.ni #%>% projectRaster(., crs = crs(lcm_tplus1.rast25.gb), method = "ngb")
  lcm_t.rast25 = lcm_t.rast25.ni #%>% projectRaster(., crs = crs(lcm_t.rast25.gb), method = "ngb")

  ## cut LCM to landscape
  tscrop.lcm_tplus1.rast25 <- crop(lcm_tplus1.rast25, extent(ts.buff %>% st_transform(29903) ))
  tsbuff.lcm_tplus1.rast25 <- fast_mask(tscrop.lcm_tplus1.rast25, ts.buff %>% st_transform(29903))
  
  tscrop.lcm_t.rast25 <- crop(lcm_t.rast25, extent(ts.buff %>% st_transform(29903)))
  tsbuff.lcm_t.rast25 <- fast_mask(tscrop.lcm_t.rast25, ts.buff %>% st_transform(29903))
  
  } else{
  lcm_tplus1.rast25 = lcm_tplus1.rast25.gb
  lcm_t.rast25 = lcm_t.rast25.gb

  ## cut LCM to landscape
  tscrop.lcm_tplus1.rast25 <- crop(lcm_tplus1.rast25, extent(ts.buff))
  tsbuff.lcm_tplus1.rast25 <- fast_mask(tscrop.lcm_tplus1.rast25, ts.buff)
  
  tscrop.lcm_t.rast25 <- crop(lcm_t.rast25, extent(ts.buff))
  tsbuff.lcm_t.rast25 <- fast_mask(tscrop.lcm_t.rast25, ts.buff)
  
}

# Create a named list of rasters, keyed by year

raster.list <- setNames(
  list(
    tsbuff.lcm_t.rast25,
    tsbuff.lcm_tplus1.rast25
  ),
  sort(years.considered)
)


# Create lcm.landscape by explicitly matching years
lcm.landscape <- lapply(years.considered, function(yr) {
  list(year = yr, raster = raster.list[[as.character(yr)]])
})

# tscrop.lcm_tplus1.rast25.unpro <- crop(lcm_tplus1.rast25.unpro, extent(ts.buff))
# tsbuff.lcm_tplus1.rast25.unpro <- fast_mask(tscrop.lcm_tplus1.rast25.unpro, ts.buff)
# 
# tscrop.lcm_t.rast25.unpro <- crop(lcm_t.rast25.unpro, extent(ts.buff))
# tsbuff.lcm_t.rast25.unpro <- fast_mask(tscrop.lcm_t.rast25.unpro, ts.buff)


### select awi ----
tsbuff.awi <- st_intersection(awi, ts.buff) %>% 
  st_make_valid() %>%  st_cast("MULTIPOLYGON") %>% 
  st_cast("POLYGON")

### select nwss with more than 50% native canopy and rasterise with ts ----

if(this.country == "Scotland"){
  tsbuff.nwss <- st_intersection(nwss[nwss$NATIVE_PCT>=50,] %>% st_cast("POLYGON"), ts.buff) %>% 
    st_make_valid() %>%  st_cast("MULTIPOLYGON") %>% 
    st_cast("POLYGON")
} else{ tsbuff.nwss = NA}

### rasterise layers with same config as lcm
tsbuff.rast = fasterize(ts.buff, tsbuff.lcm_tplus1.rast25)
tsbuff.awi.raster = fasterize(tsbuff.awi %>% st_as_sf(), tsbuff.lcm_tplus1.rast25)
if(this.country == "Scotland"){
  tsbuff.nwss.raster = fasterize(tsbuff.nwss %>% st_as_sf(), tsbuff.lcm_tplus1.rast25) %>% 
  buffr (., distance = 25, units = "geographic", target_value = 1)
} else{ tsbuff.nwss.raster = NA}

awi.landscape = tsbuff.awi#eg.awi
nwss.landscape = tsbuff.nwss.raster



### save ----
save(awi.landscape, 
     nwss.landscape,
  #tsbuff.lcm_tplus1.rast25,
     #tsbuff.lcm_t.rast25,
     # tsbuff.lcm_tplus1.rast25.unpro, 
     # tsbuff.lcm_t.rast25.unpro, 
     tsbuff.awi, 
     #tsbuff.nwss,
     tsbuff.nwss.raster,
     # tsbuff.rast,
    # tsbuff.awi.raster,
     tsbuff.hexgrid,
     ts.hexgrid,
     this.ts,
     this.country,
    lcm.landscape,
     file = paste0(func.conect.path, 
                   "\\analysis outputs\\", this.tss[this.ts.num], "\\r_curated data_.RData")
    )

# save some space in RAM by removing some objects ----
rm(awi, lcm_tplus1.rast25.ni,hex.grid0, hex.grid, nwss, 
   tsbuff.nwss, 
   tsbuff.nwss.raster,
   tsbuff.awi.raster, lcm_tplus1.rast25, lcm_t.rast25, tsbuff.rast, 
   tsbuff.lcm_tplus1.rast25, tsbuff.lcm_t.rast25,
   tscrop.lcm_tplus1.rast25, tscrop.lcm_t.rast25, 
   lcm_tplus1.rast25.gb,
   lcm_t.rast25.ni,
   lcm_t.rast25.gb#, 
   #tsbuff.lcm_tplus1.rast25.unpro, tsbuff.lcm_t.rast25.unpro
   )
gc()
print("Data curation (script03) done")
