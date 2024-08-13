# Ewan McHenry
##------ Fri Feb 25 09:19:50 2022 ------##
# functional connectivity metric dev
# script 05 - matrix and patch isolation

# THE MATRIX ----

# Ncells of landtypes in each hex

ts.hexgrid$lcm.ncells = exact_extract(lcm.landscape, ts.hexgrid, "count" )
bl.cells = lcm.landscape
values(bl.cells)[values(lcm.landscape) != 1] = 0
ts.hexgrid$bl.ncells = exact_extract(bl.cells, ts.hexgrid, "sum" )
landnotcoastal.cells = lcm.landscape
values(landnotcoastal.cells) = 0
values(landnotcoastal.cells)[values(lcm.landscape) %in% c(1:12,14,20:21 )] = 1
ts.hexgrid$landnotcoastal.ncells = exact_extract(landnotcoastal.cells, ts.hexgrid, "sum" )
rm(bl.cells, landnotcoastal.cells)

# make cost layer 
hab.cost.lcm = lcm.landscape  
values(hab.cost.lcm) = dispers.costs$scaled.ecolog.cost[values(lcm.landscape)]
#replace gaps with high-cost landscape - these are normally sea, but can be beyond edge of landscape, this shouldnt be a problem, becasue landscape is buffered
values(hab.cost.lcm)[is.na(values(hab.cost.lcm))] = max(values(hab.cost.lcm),na.rm=T)
# aggragate cost raster by mean to reduce resolution and computational power ----
hab.cost.lcm.100mres.mean = aggregate(hab.cost.lcm, constants$cost.res, mean, na.rm = T) # 4 x 4 mean aggregation

# mean cost traveling through "not sea" in landscape
landscape.mean.scaled.ecolog.cost.not.sea = mean(values(hab.cost.lcm)[!(values(lcm.landscape) %in% c(13, 15:19 )) & !is.na(values(lcm.landscape))])*constants$cost.res
landscape.median.scaled.ecolog.cost.not.sea = median(values(hab.cost.lcm)[!(values(lcm.landscape) %in% c(13, 15:19 )) & !is.na(values(lcm.landscape))])*constants$cost.res
# mean costs of hexes "not sea" in land scape
not.sea.cost = hab.cost.lcm
values(not.sea.cost)[(values(lcm.landscape) %in% c(13, 15:19 )) | is.na(values(lcm.landscape))] = NA
ts.hexgrid$mean.scaled.ecolog.cost.not.sea = exact_extract(not.sea.cost, ts.hexgrid, "mean" )*constants$cost.res
ts.hexgrid$median.scaled.ecolog.cost.not.sea = exact_extract(not.sea.cost, ts.hexgrid, "median" )*constants$cost.res

# euclidian distances between all patches ----
patch.euc.dists <- as.matrix(dist(st_coordinates(bl.patch.hexid.centroids), diag = T))

#warning that indexing might not work
if(sum(bl.patch.hexid.centroids$hex_splitpatch_ID !=1:length(bl.patch.hexid.centroids$hex_splitpatch_ID))>0){
  print("HEY, LOOK OUT!\n patch ID doesnt equal row number -- this will casue problems in allocating cost distances!! ")
}

time.now = Sys.time()
big.cost.dist = matrix (NA, nrow = dim(bl.patch.hexid.centroids)[1], ncol = dim(bl.patch.hexid.centroids)[1])
print("cost distance calculation")
for (hthhex in 1:length(unique(bl.patch.hexid.centroids$grid_id[bl.patch.hexid.centroids$in.landscape==T]))){
  # for each unique hexgrid ID with patches in the lansdacep 
  this.grid_id = unique(bl.patch.hexid.centroids$grid_id[bl.patch.hexid.centroids$in.landscape==T])[hthhex]
  # focal patches (dispersal sources) in this grid cell, also within the landscape
  grid.focal_centroids = bl.patch.hexid.centroids[bl.patch.hexid.centroids$grid_id == this.grid_id &
                                                    bl.patch.hexid.centroids$in.landscape == T,]
  
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
    
    grid.buffered_centroids = st_intersection(bl.patch.hexid.centroids, this.grid_buffered %>% st_transform( 27700))
    grid.buffered_centroids.sp = as( grid.buffered_centroids , Class = "Spatial")
    
    # mask cost layer to buffered hex
    gridbuff_cost <- crop(hab.cost.lcm.100mres.mean, this.grid_buffered ) %>%  
      fast_mask(., this.grid_buffered ) %>%  
      projectRaster(., crs = crs(grid.buffered_centroids.sp), res = res(hab.cost.lcm.100mres.mean), method = "ngb" )
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
effective.distance = big.cost.dist*constants$cost.res
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
     bl.patch.hexid.centroids, 
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





