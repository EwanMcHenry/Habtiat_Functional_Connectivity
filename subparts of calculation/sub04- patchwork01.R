
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


## id habitat - bl.landscape ----

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

## id patches - bl.patch.id.rast ----
## first buffer to join patches that are super close
bl.buff <- buffr (bl.landscape, distance = constants$buffer.for.patchid, units = "geographic", target_value = 1) # buffer neighbouring patches
trouble_plot(bl.buff, "buffered habitat raster, before cutting by intensive landuse")

## cut through buffered joins with any intensive landuse types
bl.buff[ lcm.landscape %in% dispers.costs$hab.num[dispers.costs$patch_breaking]] = NA
trouble_plot(bl.buff, "buffered habitat raster, cut by intensive landuse")

## id initial patches from cut habitat buffer
bl.patch.id.rast <- bl.buff %>% 
  terra::patches(directions = 8) %>% 
  terra::mask(bl.landscape)
trouble_plot(bl.patch.id.rast, "patch ID raster")


## id suboptimal edge and core habitat - sub_optimal_edge_habitat.rast - core_habitat.rast ----
# sub_optimal_edge_habitat.rast
# core_habitat.rast

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
sub_optimal_edge_habitat.rast <-  grown_edge.effecting.rast %>% 
  terra::mask(bl.landscape)
sub_optimal_edge_habitat.rast[sub_optimal_edge_habitat.rast == 0] = NA
names(sub_optimal_edge_habitat.rast) <- "suboptimal_edge_habitat"
trouble_plot(bl.landscape, "bl.landscape")
trouble_plot(sub_optimal_edge_habitat.rast, "suboptimal edge habitat")

core_habitat.rast <- !grown_edge.effecting.rast %>% 
  terra::mask(bl.landscape)
core_habitat.rast[core_habitat.rast == 0] = NA
names(core_habitat.rast) <- "core_habitat"
trouble_plot(core_habitat.rast, "core habitat")



# raster hex - hex_id.r ----

## split initial patches by hexgrid
### hex raster
hex_id.r <- terra::rasterize(
  terra::vect(tsbuff.hexgrid),
  bl.patch.id.rast,
  field = "grid_id"
) 
trouble_plot(hex_id.r, "hex ID raster")


# raster awi - awi_hab.r ----
awi_hab.r <- terra::rasterize(
  terra::vect(awi.landscape),
  bl.patch.id.rast) %>%
  resample(bl.patch.id.rast) %>%
  mask(bl.patch.id.rast) # mask to patch raster - only interested in awi within patches
names(awi_hab.r) <- "awi_hab"
awi_hab.r[is.na(awi_hab.r)] <- 0 # make NA = 0, so can sum with edge habitat to find awi edge habitat
trouble_plot(awi_hab.r, "awi_hab.r")

# raster landscape - focal_landscape.rast ----

focal_landscape.rast <- terra::rasterize(
  terra::vect(this.ts),
  bl.landscape
)
#0 nas
focal_landscape.rast[is.na(focal_landscape.rast)] <- 0
names(focal_landscape.rast) <- "focal_landscape"


#  combining ----

## stack rasters and combine in df 
r_stack <- c(
  bl.patch.id.rast,
  lcm.landscape,
  hex_id.r,
  sub_optimal_edge_habitat.rast,
  core_habitat.rast,
  awi_hab.r,
  focal_landscape.rast
)
cells <- as.data.frame(r_stack, xy = TRUE, na.rm = FALSE)

cells2 <- cells %>%
  dplyr::filter(
    !is.na(patches),
    !is.na(grid_id)
  ) %>%
  dplyr::mutate(
    edge = suboptimal_edge_habitat == 1,
    core = core_habitat == 1,
    awi  = awi_hab == 1,
    focal.landscape = focal_landscape == 1
  ) %>%
  dplyr::mutate(
    class = dplyr::case_when(
      edge & awi  ~ "awi_edge",
      edge & !awi ~ "nonawi_edge",
      core & awi  ~ "awi_core",
      core & !awi ~ "nonawi_core",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(class))

## aggregate area by patch, hex and class
cell_area <- prod(terra::res(bl.patch.id.rast))

agg <- cells2 %>%
  dplyr::group_by(patches, grid_id, focal_landscape, class) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::mutate(area = n * cell_area)
agg_wide <- agg %>%
  select(!n) %>%
  tidyr::pivot_wider(
    names_from = class,
    values_from = area,
    values_fill = 0
  )

## find patch centroids and join with area data
centroids <- cells2 %>%
  dplyr::group_by(patches, grid_id, focal_landscape) %>%
  dplyr::summarise(
    x = mean(x),
    y = mean(y),
    .groups = "drop"
  )  %>%
  dplyr::left_join(agg_wide, by = c("patches", "grid_id", "focal_landscape"))

## spatialise centroids
patch_centroid_info <- sf::st_as_sf(centroids, coords = c("x", "y"), crs = terra::crs(bl.patch.id.rast))
patch_centroid_info <- patch_centroid_info %>% 
  mutate(edge_tot = awi_edge + nonawi_edge,
         core_tot = awi_core + nonawi_core,
         patch_area = edge_tot + core_tot,
         uid = paste0(patches, "_", grid_id),
         row_id = 1:nrow(patch_centroid_info),)
trouble_plot(patch_centroid_info, "patch_centroid_info_points")

# hex summary info ----
## want to find for heach hexgrid cell: 
### area of lcm cells - total
### area of each lcm type
### area of habitat
### are of non-coastal lcm types
coastal_water_lcm <- c(13, 15:19)
### area of awi-core, awi-edge, nonawi-core, nonawi-edge 

hex_summary <- cells %>%
  dplyr::filter(
    !is.na(grid_id))  %>%
  dplyr::group_by(grid_id) %>%
  dplyr::summarise(
    total_area = dplyr::n() * cell_area,
    habitat_area = sum(!is.na(patches)) * cell_area,
    noncoastal_water_area = sum(!(lcm %in% coastal_water_lcm)) * cell_area,
    awi_edge_area = sum(suboptimal_edge_habitat == TRUE &
                     awi_hab == 1, na.rm = TRUE) * cell_area,
    nonawi_edge_area = sum(suboptimal_edge_habitat == TRUE &
                        awi_hab == 0, na.rm = TRUE) * cell_area,
    awi_core_area = sum(core_habitat == TRUE & awi_hab == 1,
                   na.rm = TRUE) * cell_area,
    nonawi_core_area =sum(core_habitat == TRUE &
                       awi_hab == 0, na.rm = TRUE) * cell_area,
    .groups = "drop"
  )

lcm_summary <- cells %>%
  dplyr::filter(
    !is.na(grid_id))  %>%
  dplyr::group_by(grid_id, lcm) %>%
  dplyr::summarise(
    area = dplyr::n() * cell_area,
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    names_from = lcm,
    values_from = area,
    names_prefix = "lcm_",
    values_fill = 0
  )
hex_summary <- hex_summary %>%
  dplyr::left_join(lcm_summary, by = "grid_id")



# SAVE PATCH DATA ----
save(patch_centroid_info,
     hex_summary,
     file = 
       paste0(func.conect.path, "\\analysis outputs\\", 
              this.tss[this.ts.num], "\\", this.year, "\\r_funcconnect_patchwork.RData")
)
terra::writeRaster(
  r_stack,
  file = paste0(func.conect.path, "\\analysis outputs\\", 
                               this.tss[this.ts.num], "\\", this.year, "\\r_stack.tif"),
  overwrite = TRUE
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


# 

rm(lcm_summary,
   hex_summary,
   patch_centroid_info,
   countries,
   tsbuff.hexgrid, ts.hexgrid,
   cells, cells2,
   agg_wide, agg,
   centroids,
   lcm.landscape,
   bl.landscape,
   # awi.bl.patch.hexid,awi.edge, tsbuff.awi, awiAREA.subpatch.hexid, 
   # edge, edge.awi.subpatch.hexid, edge.subpatch.hexid, patch.edge,
   # lcm.poly, lcm.poly.no.patch, 
   awi.landscape, 
   bl.buff,  r_stack,
   intensive.landscape)
gc()
print("Patch definition (script04) done")
# sort(sapply(ls(), function(x) object.size(get(x))), decreasing = TRUE)
