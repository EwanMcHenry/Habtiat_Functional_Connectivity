# Ewan McHenry
##------ Fri Feb 25 09:19:50 2022 ------##
# functional connectivity metric dev
# script 03 - data curation

# curation ----

## make folder for this TS's objects ----
dir.create(paste0(func.conect.path, 
                  "\\analysis outputs\\", this.tss[this.ts.num]))
dir.create(paste0(func.conect.path, 
                  "\\analysis outputs\\", this.tss[this.ts.num],"\\", this.year))

## costs and edge effects of habitat types ----

## CUT DATA TO LANDSCAPE AND SAVE ----
### select buffered landscape to avoid edge effects  ----
this.ts = Focal_landscape[this.ts.num,]
ts.buff = this.ts %>% st_simplify( preserveTopology = T, dTolerance = constants$landscape.buffer.simplification.tolerance[1]) %>% # first simplify hack to reduce run time. this is a rough buffer to negate edge effects, so can be v rough
  st_buffer( dist = constants$buffer.roundLandscape) %>% 
  st_simplify( preserveTopology = FALSE, dTolerance = constants$landscape.buffer.simplification.tolerance[2])
### name country
this.country = st_intersection( st_point_on_surface(this.ts) , countries)$name.1 

### select hex grid ----
# define Uk wide hexgrid
hex.grid0 = st_make_grid(countries, c(constants$hexdist.h, constants$hexdist.v), what = "polygons", square = F)
hex.grid = st_sf(hex.grid0) %>%
  # add grid ID
  mutate(grid_id = 1:length(lengths(hex.grid0)))

# intersect with landscape - note that grid id is from original UK-wide grid, allows easy cross-ID
tsbuff.hexgrid <- st_intersection(hex.grid, ts.buff) %>% 
  st_make_valid() %>%  st_cast("MULTIPOLYGON") %>% st_cast("POLYGON") %>% 
  dplyr::select(grid_id)

ts.hexgrid <- st_intersection(hex.grid, this.ts) %>% 
  st_make_valid() %>%  st_cast("MULTIPOLYGON") %>% st_cast("POLYGON") %>% 
  dplyr::select(grid_id) 
ts.hexgrid$hex.ha = st_area(ts.hexgrid) %>% 
  set_units(value = "ha") %>%
  as.numeric()

### select lcm and transform crs of other data if needed----
# different lcm for NI/GB
if(this.country == "N.Ireland"){
  lcm19.rast25 = lcm19.rast25.ni #%>% projectRaster(., crs = crs(lcm19.rast25.gb), method = "ngb")
  lcm90.rast25 = lcm90.rast25.ni #%>% projectRaster(., crs = crs(lcm90.rast25.gb), method = "ngb")

  ## cut LCM to landscape
  tscrop.lcm19.rast25 <- crop(lcm19.rast25, extent(ts.buff %>% st_transform(29903) ))
  tsbuff.lcm19.rast25 <- fast_mask(tscrop.lcm19.rast25, ts.buff %>% st_transform(29903))
  
  tscrop.lcm90.rast25 <- crop(lcm90.rast25, extent(ts.buff %>% st_transform(29903)))
  tsbuff.lcm90.rast25 <- fast_mask(tscrop.lcm90.rast25, ts.buff %>% st_transform(29903))
  
  } else{
  lcm19.rast25 = lcm19.rast25.gb
  lcm90.rast25 = lcm90.rast25.gb

  ## cut LCM to landscape
  tscrop.lcm19.rast25 <- crop(lcm19.rast25, extent(ts.buff))
  tsbuff.lcm19.rast25 <- fast_mask(tscrop.lcm19.rast25, ts.buff)
  
  tscrop.lcm90.rast25 <- crop(lcm90.rast25, extent(ts.buff))
  tsbuff.lcm90.rast25 <- fast_mask(tscrop.lcm90.rast25, ts.buff)
  
}

if(this.year == 2019){
  lcm.landscape = tsbuff.lcm19.rast25 #eg.lcm19.rast25
}
if(this.year == 1990){
  lcm.landscape = tsbuff.lcm90.rast25 #eg.lcm19.rast25
}



# tscrop.lcm19.rast25.unpro <- crop(lcm19.rast25.unpro, extent(ts.buff))
# tsbuff.lcm19.rast25.unpro <- fast_mask(tscrop.lcm19.rast25.unpro, ts.buff)
# 
# tscrop.lcm90.rast25.unpro <- crop(lcm90.rast25.unpro, extent(ts.buff))
# tsbuff.lcm90.rast25.unpro <- fast_mask(tscrop.lcm90.rast25.unpro, ts.buff)


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
tsbuff.rast = fasterize(ts.buff, tsbuff.lcm19.rast25)
tsbuff.awi.raster = fasterize(tsbuff.awi %>% st_as_sf(), tsbuff.lcm19.rast25)
if(this.country == "Scotland"){
  tsbuff.nwss.raster = fasterize(tsbuff.nwss %>% st_as_sf(), tsbuff.lcm19.rast25) %>% 
  buffr (., distance = 25, units = "geographic", target_value = 1)
} else{ tsbuff.nwss.raster = NA}

awi.landscape = tsbuff.awi#eg.awi
nwss.landscape = tsbuff.nwss.raster



### save ----
save(awi.landscape, 
     nwss.landscape,
  #tsbuff.lcm19.rast25,
     #tsbuff.lcm90.rast25,
     # tsbuff.lcm19.rast25.unpro, 
     # tsbuff.lcm90.rast25.unpro, 
     tsbuff.awi, 
     #tsbuff.nwss,
     tsbuff.nwss.raster,
     # tsbuff.rast,
    # tsbuff.awi.raster,
     tsbuff.hexgrid,
     ts.hexgrid,
     this.ts,
     this.country, 
     file = paste0(func.conect.path, 
                   "\\analysis outputs\\", this.tss[this.ts.num], "\\", this.year, "\\r_curated data_.RData")
    )

# save some space in RAM by removing some objects ----
rm(awi, lcm19.rast25.ni,hex.grid0, hex.grid, 
   tsbuff.nwss, 
   tsbuff.nwss.raster,
   tsbuff.awi.raster, lcm19.rast25, lcm90.rast25, tsbuff.rast, 
   tsbuff.lcm19.rast25, tsbuff.lcm90.rast25,
   tscrop.lcm19.rast25, tscrop.lcm90.rast25, 
   lcm19.rast25.gb,
   lcm90.rast25.ni,
   lcm90.rast25.gb#, 
   #tsbuff.lcm19.rast25.unpro, tsbuff.lcm90.rast25.unpro
   )
gc()
print("Data curation (script03) done")
