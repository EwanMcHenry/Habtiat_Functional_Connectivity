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
  st_buffer( dist = constants$max.dispersal.considered) %>% # buffer
  st_simplify( preserveTopology = FALSE, dTolerance = constants$landscape.buffer.simplification.tolerance[2]) # simplify again
### name of country ----
this.country = st_intersection( st_point_on_surface(this.ts) , countries)$nation 

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


# transform landscape and data to crs ----

crs_use <- if (this.country == "Northern Ireland") {
  29903
} else {
  27700
}

awi <- st_transform(awi, crs_use)
nwss <- st_transform(nwss, crs_use)

this.ts <- st_transform(this.ts, crs_use)
ts.buff <- st_transform(ts.buff, crs_use)
tsbuff.hexgrid <- st_transform(tsbuff.hexgrid, crs_use)
ts.hexgrid <- st_transform(ts.hexgrid, crs_use)


### select awi ----
awi.landscape <- st_intersection(awi, ts.buff) %>% 
  st_make_valid() %>%  st_cast("MULTIPOLYGON") %>% 
  st_cast("POLYGON")

### select nwss with more than 50% native canopy and rasterise with ts ----

if(this.country == "Scotland"){
  nwss.native.conifer.landscape <- st_intersection(nwss[nwss$NATIVE_PCT>=50,] %>% st_cast("POLYGON"), ts.buff) %>% 
    st_make_valid() %>%  st_cast("MULTIPOLYGON") %>% 
    st_cast("POLYGON")
} else{ nwss.native.conifer.landscape = NA}

# ### rasterise layers with same config as lcm
# tsbuff.rast = fasterize(ts.buff, tsbuff.lcm_tplus1.rast25)
# tsbuff.awi.raster = fasterize(tsbuff.awi %>% st_as_sf(), tsbuff.lcm_tplus1.rast25)
# if(this.country == "Scotland"){
#   tsbuff.nwss.raster = fasterize(tsbuff.nwss %>% st_as_sf(), tsbuff.lcm_tplus1.rast25) %>% 
#   buffr (., distance = 25, units = "geographic", target_value = 1)
# } else{ tsbuff.nwss.raster = NA}
# 
# awi.landscape = tsbuff.awi#eg.awi
# nwss.native.conifer.landscape = tsbuff.nwss.raster
# 

### select roads ----
roads.landscape <- st_intersection(roads_uk_major, ts.buff) %>% 
  st_make_valid() 
roads.landscape <- roads.landscape %>%
  st_buffer(dist = .$buffer_size) %>% 
  st_union() %>% st_make_valid() %>% st_cast("POLYGON") #%>% st_simplify(preserveTopology = TRUE, dTolerance = 5) # buffer roads by width, union to avoid overlaps, simplify to reduce complexity


# troubleshooting saves ----
if(troubleshooting ==T){
  st_write(tsbuff.hexgrid, 
           paste0(func.conect.path, "\\troubleshooting_saves\\tsbuff_hexgrid.gpkg"), delete_dsn = TRUE)
  st_write(ts.hexgrid, paste0(func.conect.path, "\\troubleshooting_saves\\ts_hexgrid.gpkg"), delete_dsn = TRUE)
  st_write(this.ts, paste0(func.conect.path, "\\troubleshooting_saves\\ts.gpkg"), delete_dsn = TRUE)
  st_write(ts.buff, paste0(func.conect.path, "\\troubleshooting_saves\\tsbuff.gpkg"), delete_dsn = TRUE)
  st_write(awi.landscape, paste0(func.conect.path, "\\troubleshooting_saves\\awi_landscape.gpkg"), delete_dsn = TRUE)
  st_write(roads.landscape, paste0(func.conect.path, "\\troubleshooting_saves\\roads_landscape.gpkg"), delete_dsn = TRUE)
  
  if(this.country == "Scotland"){
    st_write(nwss.native.conifer.landscape, paste0(func.conect.path, "\\troubleshooting_saves\\nwss_landscape.gpkg"), delete_dsn = TRUE)
  }
  
}

### save ----
save(awi.landscape, 
     nwss.native.conifer.landscape,
     roads.landscape,

     file = paste0(func.conect.path, 
                   "\\analysis outputs\\", this.tss[this.ts.num], "\\r_curated_data.awi_nwss_roads.RData")
)

save(tsbuff.hexgrid,
     awi.landscape,
     ts.hexgrid,
     this.ts,
     this.country,
     roads.landscape,
     countries,
     file = paste0(func.conect.path, 
                   "\\analysis outputs\\", this.tss[this.ts.num], "\\r_curated data_.RData")
    )

# save some space in RAM by removing some objects ----
rm(roads.landscape, roads_uk_major,
   countries,
   awi, lcm_tplus1.rast25.ni,hex.grid0, hex.grid, nwss, 
   tsbuff.nwss, 
   tsbuff.nwss.raster,
   tsbuff.awi.raster, lcm_tplus1.rast25, lcm_t.rast25, tsbuff.rast, 
   tsbuff.lcm_tplus1.rast25, tsbuff.lcm_t.rast25,
   tscrop.lcm_tplus1.rast25, tscrop.lcm_t.rast25, 
   lcm_tplus1.rast25.gb,
   lcm_t.rast25.ni,
   lcm_t.rast25.gb,awi.landscape, 
     nwss.native.conifer.landscape,
     
     tsbuff.hexgrid,
   ts.hexgrid,
   this.ts,
   
   this.country,
   
   lcm.landscape
   #, 
   #tsbuff.lcm_tplus1.rast25.unpro, tsbuff.lcm_t.rast25.unpro
)
gc()
print("Data curation (script03) done")
