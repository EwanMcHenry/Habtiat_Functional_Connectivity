# Ewan McHenry
##------ Fri Feb 25 09:19:50 2022 ------##
# functional connectivity metric dev

# script 05 - matrix and patch isolation


# ## load curated data ----
load(paste0(func.conect.path, "\\analysis outputs\\", this.tss[this.ts.num], "\\r_curated data_.RData"))
## LOAD PATCH DATA ----
load(paste0(func.conect.path, "\\analysis outputs\\", this.tss[this.ts.num], "\\", this.year, "\\r_funcconnect_patchwork.RData"))

r_stack <- terra::rast(
  paste0(func.conect.path, "\\analysis outputs\\", 
         this.tss[this.ts.num], "\\", this.year, "\\r_stack.tif"))


# THE MATRIX ----

# make cost layer ----
hab.cost.lcm = r_stack$lcm  
values(hab.cost.lcm) = dispers.costs$scaled.ecolog.cost[values(r_stack$lcm)]
trouble_plot(hab.cost.lcm, "habitat cost layer")
#replace gaps with high-cost landscape - these are normally sea, but can be beyond edge of landscape, this shouldnt be a problem, becasue landscape is buffered
values(hab.cost.lcm)[is.na(values(hab.cost.lcm))] = max(dispers.costs$scaled.ecolog.cost[values(r_stack$lcm)],na.rm=T)
# aggragate cost raster by mean to reduce resolution and computational power ----
hab.cost.agg.rast = aggregate(hab.cost.lcm, constants$cost_agg_n, mean, na.rm = T) # 4 x 4 mean aggregation

# LANDSCAPE cost INFO ----

# Mean landscape cost
# mean cost traveling through "not sea" in landscape
landscape.mean.ecolog.cost.not.sea = mean(values(hab.cost.lcm)[!(values(r_stack$lcm  ) %in% coastal_water_lcm) & !is.na(values(r_stack$lcm  ))])*constants$cost.scale.factor
# landscape.median.scaled.ecolog.cost.not.sea = median(values(hab.cost.lcm)[!(values(r_stack$lcm) %in%coastal_water_lcm) & !is.na(values(r_stack$lcm))])*constants$cost.scale.factor

# mean hex cost
# mean costs of hexes "not sea" in land scape
not.sea.cost = hab.cost.lcm
values(not.sea.cost)[(values(r_stack$lcm) %in% coastal_water_lcm) | is.na(values(r_stack$lcm))] = NA
ts.hexgrid$mean.ecolog.cost.not.sea = exact_extract(not.sea.cost, ts.hexgrid, "mean" )*constants$cost.scale.factor
# ts.hexgrid$median.ecolog.cost.not.sea = exact_extract(not.sea.cost, ts.hexgrid, "median" )*constants$cost.scale.factor


# euclidian distances between all patches ----
patch.euc.dists <- as.matrix(dist(st_coordinates(patch_centroid_info), diag = T))


# cost distance between patches ----
#warning that indexing might not work
if(sum(patch_centroid_info$row_id !=1:length(patch_centroid_info$uid))>0){
  print("HEY, LOOK OUT!\n patch ID doesnt equal row number -- this will casue problems in allocating cost distances!! ")
}

time.now = Sys.time()
print("cost distance calculation")


big.cost.dist = matrix (NA, nrow = dim(patch_centroid_info[patch_centroid_info$focal_landscape==1,])[1], ncol = dim(patch_centroid_info)[1]) # dataframe to save cost distances
patches <- patch_centroid_info %>% select(patches, grid_id, focal_landscape )


# build candidate pairs - global neighbourhood graph
nn <- sf::st_is_within_distance(
  patches[patch_centroid_info$focal_landscape==1,], # only interested in within landscape connections - had taken this out before.. though not sure why
  patches,
  dist = constants$max.dispersal.considered # cutoff distance beyond which pair not considered
)

edges <- data.frame(
  from = rep(seq_along(nn), lengths(nn)),
  to   = unlist(nn)) %>%
  dplyr::filter(from < to)

# attach coordinates to edges 
coords <- sf::st_coordinates(patch_centroid_info)

edges$x_from <- coords[edges$from, 1]
edges$y_from <- coords[edges$from, 2]
edges$x_to   <- coords[edges$to, 1]
edges$y_to   <- coords[edges$to, 2]


pairs <- data.frame(
  from = rep(seq_along(nn), lengths(nn)),
  to   = unlist(nn)
)


for (hthhex in 1:length(unique(patch_centroid_info$grid_id[patch_centroid_info$focal_landscape ==1]))){
  # for each unique hexgrid ID with patches in the lansdacep 
  this.grid_id = unique(patch_centroid_info$grid_id[patch_centroid_info$in.landscape==T])[hthhex]
  # focal patches (dispersal sources) in this grid cell, also within the landscape
  grid.focal_centroids = patch_centroid_info[patch_centroid_info$grid_id == this.grid_id &
                                                    patch_centroid_info$in.landscape == T,]
  
  # if there are any patches in this cell, in this landscape ()
#  if(length(grid.focal_centroids)>0){ # ive taken this out, I dont think it nessessary any more now loop only considering focal patches that are in landscape

    #buffer focal grid cell
  if(this.country == "N.Ireland"){ # transform crs
    this.grid_buffered <- ts.hexgrid[ts.hexgrid$grid_id == this.grid_id,] %>%
      st_simplify( preserveTopology = T, dTolerance = 100) %>% # first simplify hack to reduce run time. this is a rough buffer to negate edge effects, so can be v rough
      st_buffer( dist = constants$max.dispersal.considered, joinStyle = "MITRE" ) %>% 
      st_union() %>% # in case different island within hex create multiple polygons
      st_as_sf() %>% 
      st_simplify( preserveTopology = FALSE, dTolerance = constants$max.dispersal.considered/50) %>% # v rough, but contribution to connectivity max dispers km away is sooo tiny
      st_transform(29903)    
  } else{
    this.grid_buffered <- ts.hexgrid[ts.hexgrid$grid_id == this.grid_id,] %>%
      st_simplify( preserveTopology = T, dTolerance = 100) %>% # first simplify hack to reduce run time. this is a rough buffer to negate edge effects, so can be v rough
      st_buffer( dist = constants$max.dispersal.considered, joinStyle = "MITRE" ) %>% 
      st_union() %>% # in case different island within hex create multiple polygons
      st_as_sf() %>% 
      st_simplify( preserveTopology = FALSE, dTolerance = constants$max.dispersal.considered/50) # v rough, but contribution to connectivity max dispers km away is sooo tiny
  }
    # id patch centroids within buffered and focal hex as sp objects -- needed for fast  costdistance()
    grid.focal_centroids.sp = as( grid.focal_centroids %>% st_transform( 27700) , Class = "Spatial")
    
    grid.buffered_centroids = st_intersection(patch_centroid_info, this.grid_buffered %>% st_transform( 27700))
    grid.buffered_centroids.sp = as( grid.buffered_centroids , Class = "Spatial")
    
    # mask cost layer to buffered hex
    gridbuff_cost <- crop(hab.cost.agg.rast, this.grid_buffered ) %>%  
      fast_mask(., this.grid_buffered ) %>%  
      projectRaster(., crs = crs(grid.buffered_centroids.sp), res = res(hab.cost.agg.rast), method = "ngb" )
    #gridbuff_cost[is.na(gridbuff_cost)] = max(values(gridbuff_cost), na.rm = T)
    
    # transition matrix and correct for diagonal/ortho -- note scale = F, to get real cost distance measure, rather than scaled
    gridbuff_Tcost <- transition(gridbuff_cost,function(x) 1/mean(x),8)
    gridbuff_Tcost.C <- geoCorrection(gridbuff_Tcost, type="c", multpl=FALSE, scl=F)
    # cost distances between patches in focal hex and buffer around 
    big.cost.dist[grid.focal_centroids$hex_splitpatch_ID , grid.buffered_centroids$hex_splitpatch_ID  ] = # potential trouble here if ID doesnt == row number in OG patch centroids data
      costDistance(gridbuff_Tcost.C,grid.focal_centroids.sp, grid.buffered_centroids.sp)
    
 # }
  svMisc::progress(hthhex)
}

time.now - Sys.time()
# rescale cost distance to actual effective distance (before costs scaled to be close to 1 for computational efficiency)
effective.distance = big.cost.dist*constants$cost.scale.factor
rm(big.cost.dist)
# and make all that are beyond maximum considered distance NA. the efficiecny hack emplyed here (doing cost.distance by hexes buffered by max considered distance) will have calculated for some it didnt need to.
effective.distance [patch.euc.dists>constants$max.dispersal.considered] = NA
print("cost distance done")


# SAVE PATCH DATA ----
save(effective.distance,
   #  big.cost.dist,
     landscape.mean.scaled.ecolog.cost.not.sea,
     landscape.median.scaled.ecolog.cost.not.sea,
     ts.hexgrid,
     patch.euc.dists,
     patch_centroid_info, 
     file = 
     paste0(func.conect.path, 
            "\\analysis outputs\\", this.ts.for.loop[this.ts.num], "\\", this.year, "\\r_funcconnect_MatrixCostDists.RData")
)

if(grepl("Illustrative", this.ts.for.loop[this.ts.num]) ){
  save(gridbuff_Tcost.C,
       file = 
         paste0(func.conect.path, 
                "\\analysis outputs\\", 
                     this.ts.for.loop[this.ts.num], "\\", this.year, "\\cost.objects.RData")
  )
  
}

rm()
gc()
     
print("Landscape matrix and patch linkage (script05) done")





