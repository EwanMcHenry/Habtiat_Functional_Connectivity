
# Ewan McHenry
##------ Fri Feb 25 09:19:50 2022 ------##
# functional connectivity metric dev
# script 04 - patch work

## id broadleaf patches from raster of land cover 
## find patch area
## intersect with grid
## calculate area and awi area within patches
## find area of hostile edge within patch and within-patch-awi hostile edge
## patch-grid centroids

# PATCH WORK ----


## id habitat ----

## edit native conifer in scotland
if(this.country == "Scotland" & constants$focal.hab.num.lcm ==1){ # if scotland and focal habitat is broadleaf
  # change conifer cells with >50% native canopy in NWSS to focal habitat in lcm landscape
  ## (for habitat id, intensive landuse edge effect and movement calcs)
  values(lcm.landscape)[
    values(lcm.landscape) == constants$alt.hab.scot.nwss & 
      values(nwss.native.conifer.landscape) == 1] = constants$focal.hab.num.lcm
}

bl.landscape = lcm.landscape
bl.landscape[lcm.landscape != constants$focal.hab.num.lcm] = NA
bl.landscape[lcm.landscape == constants$focal.hab.num.lcm] = 1
trouble_plot(bl.landscape, "broadleaf (and conifer in scotland) habitat raster")

## id patches ----
## first buffer to join patches that are super close
bl.buff <- buffr (bl.landscape, distance = constants$buffer.for.patchid, units = "geographic", target_value = 1) # buffer neighbouring patches
trouble_plot(bl.buff, "buffered habitat raster, before cutting by intensive landuse")

## cut through buffered joins with any intensive landuse types
bl.buff[ lcm.landscape %in% dispers.costs$hab.num[dispers.costs$patch_breaking]] = NA
trouble_plot(bl.buff, "buffered habitat raster, cut by intensive landuse")

## id initial patches from cut habitat buffer
bl.patch.id <- bl.buff %>% 
  terra::patches(directions = 8) %>% 
  terra::mask(bl.landscape)
trouble_plot(bl.patch.id, "patch ID raster")


## id non-habitat for edge effect calcs ----
intensive.landscape = lcm.landscape 
#remove all semi-natural landcover types
intensive.landscape[lcm.landscape %in% dispers.costs$hab.num[dispers.costs$semi.natural]] = NA
trouble_plot(intensive.landscape, "intensive landuse raster")
#remove all beyond possible affecting range  from focal hab, to speed up edge effect calcs
edge.max.effecting <- buffr(bl.landscape, distance = max(dispers.costs$edge.extent, na.rm = T)*1.1, units = "geographic", target_value = 1) # 
edge.effecting.landscape =   intensive.landscape %>% 
  terra::mask(edge.max.effecting)
trouble_plot(edge.effecting.landscape, "edge effecting landscape raster")
# buffer edge.effecting.landscape by area of effect
## individual buffering by type then combine
classes <- dispers.costs$hab.num[!dispers.costs$semi.natural]
edge.effecting_list <- lapply(classes, function(cl) {
  r <- edge.effecting.landscape == cl
  r[is.na(r)] <- 0
  r
})
trouble_plot(edge.effecting_list[[2]], "edge effecting landscape - e.g arable")

grow_class <- function(r, dist) {
  d <- terra::distance(r) 
  r_grown <- d <= (dist*2) # is the raster cell within the distance of effect of this class? multiplied by 2 becasue the distance given is the maxis from the edge of the class (e.g. 25 for the neighbouring 25m cell with res 25m)... so here its like, is the max buffer distance closer to this cell or the next
  trouble_plot(d, "distance to edge effecting class")
  
  r_grown
}
grown_edge.effecting_list <- lapply(classes, function(cl) {
  r <- 1*(edge.effecting_list[[which(classes==cl)]])
  r[r==0] <- NA
  dist <-  dispers.costs$edge.extent[dispers.costs$hab.num == cl] 
  g  <-  grow_class(r, dist)
  g
})
trouble_plot(grown_edge.effecting_list[[2]], "grown edge effecting landscape - e.g arable")
# combine all edge effecting classes into one raster
grown_edge.effecting.rast <- Reduce(`|`, grown_edge.effecting_list)
trouble_plot(grown_edge.effecting.rast, "grown edge effecting landscape - combined")
sub_optimal_edge_habitat <-  grown_edge.effecting.rast %>% 
  terra::mask(bl.landscape)
trouble_plot(bl.landscape, "bl.landscape")
trouble_plot(sub_optimal_edge_habitat, "suboptimal edge habitat")
core_habitat <- !grown_edge.effecting.rast %>% 
  terra::mask(bl.landscape)
core_habitat[core_habitat == 0] = NA
trouble_plot(core_habitat, "core habitat")



## id area of suboptimal edge habiat ----

## split initial patches by hexgrid
### hex raster
hex.r <- terra::rasterize(
  terra::vect(tsbuff.hexgrid),
  bl.patch.id,
  field = "grid_id"
) %>% 
  terra::resample(., bl.patch.id, method = "near")
trouble_plot(hex.r, "hex ID raster")
## df of patch and hexid
cells <- data.frame(
  patch = terra::values(bl.patch.id),
  hex   = terra::values(hex.r)
)
cells <- cells[complete.cases(cells), ]
cells$uid <- interaction(cells$patches, cells$grid_id, drop = TRUE)

bl.hexid.rast <- terra::rasterize(tsbuff.hexgrid %>% 
                                    st_transform(crs(bl.patch.id)), 
                                  bl.patch.id, field = "grid_id") %>% 
  mask(., bl.patch.id)
trouble_plot(bl.hexid.rast, "habitat hex ID raster")


### subpatches polygonise and split by grid cell ----
# allows reporting on that scale

# polygonise patches
# if(this.country == "N.Ireland"){
#   bl.patch.id.poly = stars::st_as_stars(bl.patch.id) %>% 
#     st_as_sf() %>% group_by(clumps ) %>% summarize(geometry = st_union(geometry)) %>%  st_make_valid() %>%  st_cast("MULTIPOLYGON") %>%   st_cast("POLYGON") %>% # to make work for NI -- overcoming some error to do with projection triggered by st_as_sf(merge=T)
#     sf::as_Spatial() %>% 
#     # sf::as_Spatial(sf::st_as_sf(stars::st_as_stars(bl.patch.id), 
#     #                                              as_points = FALSE, merge = TRUE)) %>% 
#     rgeos::gBuffer( byid = TRUE, width = 0) %>% st_as_sf()%>% 
#     st_transform(27700) %>% 
#     st_simplify(dTolerance = 5)  %>% 
#     rowid_to_column( "subpatch_ID")
#   
# } else{
  bl.patch.id.poly = stars::st_as_stars(bl.patch.id) %>% 
    st_as_sf(merge = T, connect8 = T) %>% 
    sf::as_Spatial() %>% 
    # sf::as_Spatial(sf::st_as_sf(stars::st_as_stars(bl.patch.id), 
    #                                              as_points = FALSE, merge = TRUE)) %>% 
    rgeos::gBuffer( byid = TRUE, width = 0) %>% st_as_sf()%>% 
    st_transform(27700) %>% 
    st_simplify(dTolerance = 5)  %>% 
    rowid_to_column( "subpatch_ID")
# }


# Intersect with grid
bl.patch.id.poly.hexid = st_intersection(bl.patch.id.poly, tsbuff.hexgrid) %>% 
  st_make_valid() %>%  
  # st_cast("MULTIPOLYGON") %>% st_cast("POLYGON") %>%  # for consistant geometries across features
  # group_by(paste(clumps,grid_id )) %>% summarize(geometry = st_union(geometry)) %>%  # stop diagonal joins within patch breaking off
  rowid_to_column( "hex_splitpatch_ID")

## Patch areas ---- 
# sub patch.hex, sub patch and patch clump (not including clumping buffer)
bl.patch.id.poly.hexid$subpatch.hex.ha = st_area(bl.patch.id.poly.hexid) %>% 
  set_units(value = "ha") %>%
  as.numeric()
bl.patch.id.poly.hexid$subpatch.ha = full_join(bl.patch.id.poly.hexid, aggregate(bl.patch.id.poly.hexid$subpatch.hex.ha,list(bl.patch.id.poly.hexid$subpatch_ID) , sum) %>% 
                                                 rename(subpatch_ID = Group.1, subpatch.area = x))$subpatch.area
bl.patch.id.poly.hexid$clump.ha = full_join(bl.patch.id.poly.hexid, aggregate(bl.patch.id.poly.hexid$subpatch.hex.ha,list(bl.patch.id.poly.hexid$clumps) , sum) %>% 
                                              rename(clumps = Group.1, clump.area = x))$clump.area

## AWI area within patches ----
awi.bl.patch.hexid0 = st_intersection(bl.patch.id.poly.hexid, awi.landscape[,]) 
awi.bl.patch.hexid = awi.bl.patch.hexid0[npts(awi.bl.patch.hexid0, by_feature = T) >= 4,] %>%  # hot fix to stop an error, only 2 in SSR landscape get dropped out
  st_make_valid() %>% 
  # st_cast("MULTIPOLYGON") %>% st_cast("POLYGON") %>%  # for consistant geometries across features
  st_simplify(dTolerance = 5) 
rm(awi.bl.patch.hexid0)

awi.bl.patch.hexid$awi.subpatch.hex.ha <- st_area(awi.bl.patch.hexid) %>% 
  set_units(value = "ha") %>%
  as.numeric()

awiAREA.subpatch.hexid <- aggregate(awi.bl.patch.hexid$awi.subpatch.hex.ha,list(awi.bl.patch.hexid$hex_splitpatch_ID) , sum) %>%
  rename(hex_splitpatch_ID = Group.1, awi.subpatch.hex.ha = x)

bl.patch.id.poly.hexid$awi.subpatch.hex.ha = full_join(bl.patch.id.poly.hexid, awiAREA.subpatch.hexid %>% as.data.frame(), 
                                                       by = "hex_splitpatch_ID"
)$awi.subpatch.hex.ha %>% 
  replace_na(0)

## Negative edge effects ----
### hostile edge polygon ----
### polygonise lcm without woodland and join edge info

  lcm.poly = stars::st_as_stars(lcm.landscape) %>% 
    st_as_sf(as_points = FALSE, merge = T, connect8 = T ) %>% #group_by(layer) %>% summarize(geometry = st_union(geometry)) %>%  st_make_valid() %>%  st_cast("MULTIPOLYGON") %>%   st_cast("POLYGON") %>% # to make work for NI -- overcoming some error to do with projection triggered by st_as_sf(merge=T)
    sf::as_Spatial() %>% 
    rgeos::gBuffer( byid = TRUE, width = 0) %>% st_as_sf()%>% 
    st_transform(27700) %>% 
    rename(hab.num = layer) %>% 
    mutate(hab.num = factor(hab.num)) %>% 
    st_simplify(dTolerance = 5) 
  


# add info on extent of negative edges
lcm.poly.no.patch = lcm.poly[lcm.poly$hab.num != 1,] %>% 
  left_join(dispers.costs, by = "hab.num")
# create hostile area
edge <- lcm.poly.no.patch %>% 
  st_buffer(dist = lcm.poly.no.patch$edge.extent, endCapStyle = "SQUARE", joinStyle = "MITRE" ) %>% 
  st_union()%>% 
  st_simplify(dTolerance = 5) 

### subpatch edge area ----
patch.edge = st_intersection(bl.patch.id.poly.hexid, edge) #%>% 
  # st_make_valid() %>%   st_cast("MULTIPOLYGON") %>% st_cast("POLYGON") # for consistant geometries across features
patch.edge$edge.subpatch.hex.ha <- st_area(patch.edge) %>% 
  set_units(value = "ha") %>%
  as.numeric()
edge.subpatch.hexid <- aggregate(patch.edge$edge.subpatch.hex.ha,list(patch.edge$hex_splitpatch_ID) , sum) %>%
  rename(hex_splitpatch_ID = Group.1, edge.subpatch.hex.ha = x)

bl.patch.id.poly.hexid$edge.subpatch.hex.ha = full_join(bl.patch.id.poly.hexid, edge.subpatch.hexid %>% as.data.frame(), 
                                                        by = "hex_splitpatch_ID")$edge.subpatch.hex.ha %>% 
  replace_na(0)

### subpatch awi edge area ----
awi.edge = st_intersection(awi.bl.patch.hexid, edge) #%>% 
  # st_make_valid() %>% st_cast("MULTIPOLYGON") %>% st_cast("POLYGON") # for consistant geometries across features

awi.edge$edge.awi.hex.ha <- st_area(awi.edge) %>% 
  set_units(value = "ha") %>%
  as.numeric()# this wont work -- need to aggregate by sub hex patch

edge.awi.subpatch.hexid <- aggregate(awi.edge$edge.awi.hex.ha,list(awi.edge$hex_splitpatch_ID) , sum) %>%
  rename(hex_splitpatch_ID = Group.1, edge.awi.subpatch.hex.ha = x)


bl.patch.id.poly.hexid$edge.awi.subpatch.hex.ha = full_join(bl.patch.id.poly.hexid, edge.awi.subpatch.hexid %>% as.data.frame(), 
                                                            by = "hex_splitpatch_ID")$edge.awi.subpatch.hex.ha %>% 
  replace_na(0)

## patch centroids ----
bl.patch.hexid.centroids = st_point_on_surface(bl.patch.id.poly.hexid) %>% st_transform( 27700) # point on surface to make sure centroid is actually in teh patch


bl.patch.hexid.centroids$in.landscape = lengths(st_intersects(bl.patch.hexid.centroids, this.ts)) > 0
bl.patch.hexid.centroids.sp = as( bl.patch.hexid.centroids , Class = "Spatial")


# SAVE PATCH DATA ----
save(lcm.landscape,
     bl.landscape,
     bl.patch.hexid.centroids,
     bl.patch.hexid.centroids.sp,
     bl.patch.id.poly.hexid,
     file = 
       paste0(func.conect.path, "\\analysis outputs\\", 
              this.tss[this.ts.num], "\\", this.year, "\\r_funcconnect_patchwork.RData")
)

if(grepl("Illustrative", this.tss[this.ts.num]) ){
  save(awi.bl.patch.hexid,
       patch.edge,
       awi.edge,
       awi.landscape,
       edge,
       intensive.landscape,
       file = 
         paste0(func.conect.path, "\\analysis outputs\\", 
                this.tss[this.ts.num], "\\", this.year, "\\edge_awi.polys.RData")
  )
  
  }

rm(lcm.landscape,
   bl.landscape,
   awi.bl.patch.hexid,
   awi.edge, awi.landscape, tsbuff.awi, awiAREA.subpatch.hexid, 
   bl.buff,  
   edge, edge.awi.subpatch.hexid, edge.subpatch.hexid, patch.edge,
   lcm.poly, lcm.poly.no.patch, 
   intensive.landscape)
gc()
print("Patch definition (script04) done")
