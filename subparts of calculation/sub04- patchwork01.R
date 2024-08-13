
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
## ID patches ----
### based on nehbouring clumps ---- 
bl.lcm = lcm.landscape#eg.lcm19.rast25 
bl.lcm[lcm.landscape != constants$focal.hab.num.lcm] = NA

nobl.lcm = lcm.landscape 
nobl.lcm[lcm.landscape == constants$focal.hab.num.lcm] = NA

if(this.country == "Scotland" & constants$focal.hab.num.lcm %in% 1:2){ # if scotland conifer withini native woodland (NWSS>50%) == patch  
  values(bl.lcm)[values(lcm.landscape) == 2 & values(nwss.landscape) == 1] = constants$focal.hab.num.lcm # conifer is focal habitat in scotland
  values(lcm.landscape)[values(lcm.landscape) == 2 & values(nwss.landscape) == 1] = constants$focal.hab.num.lcm
}

bl.buff <- buffr (bl.lcm, distance = constants$buffer.for.patchid, units = "geographic", target_value = 1) # buffer neighbouring patches
bl.patch.id <- bl.buff %>%  clump() %>% mask(.,bl.lcm)# patch ids to buffered clumps, then remove buffered area so patch ID = clump ID

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
       file = 
         paste0(func.conect.path, "\\analysis outputs\\", 
                this.tss[this.ts.num], "\\", this.year, "\\edge_awi.polys.RData")
  )
  
  }

rm(awi.bl.patch.hexid,
   awi.edge, awi.landscape, tsbuff.awi, awiAREA.subpatch.hexid, 
   bl.buff, bl.lcm, 
   edge, edge.awi.subpatch.hexid, edge.subpatch.hexid, patch.edge,
   lcm.poly, lcm.poly.no.patch, 
   nobl.lcm)
gc()
print("Patch definition (script04) done")
