# Ewan McHenry
##------ Fri Mar 04 11:59:09 2022 ------##
# functional connectivity metric dev
# script 06 - equivelent connected area

## LOAD PATCH DATA ----
load(paste0(func.conect.path, "\\analysis outputs\\", this.tss[this.ts.num], "\\", this.year, "\\lcmres", constants$lcm.res, 
            "\\04_r_funcconnect_patchwork.RData"))
## Load cost distance data
load(paste0(func.conect.path, "\\analysis outputs\\", this.tss[this.ts.num], "\\", this.year, "\\lcmres", constants$lcm.res, 
            "\\05_r_funcconnect_MatrixCostDists.RData")
)

landscape_stats[[this.ts.num]]$year_stats$n.patch = sum(patch_centroid_info$focal_landscape==1)

# convert areas to ha -----
patch_centroid_info <- patch_centroid_info %>%
  mutate(
    patch_area_ha = patch_area / 10000,
    edge_tot_ha   = edge_tot / 10000,
    core_tot_ha   = core_tot / 10000,
    awi_tot_ha    = awi_tot / 10000,
    nonawi_tot_ha = nonawi_tot / 10000
  )

#EFFECTIVE PATCH AREA (EPA) ----

patch_centroid_info <- patch_centroid_info %>% 
  mutate(
    effective_patch_area_ha = core_tot_ha * constants$non.awi.qual.eff + # core non awi value
      awi_tot_ha * constants$awi.qual.eff + # core awi value
      edge_tot_ha * constants$non.awi.qual.eff * constants$relative.edge.quality + # edge non-awi value
      edge_tot_ha * constants$awi.qual.eff * constants$relative.edge.quality # edge awi value
  )

# patch contributions to ECA ----
## EPA to pairs----
candidate_pairs = candidate_pairs %>% 
  left_join(patch_centroid_info %>% as.data.frame() %>% dplyr::select(uid, effective_patch_area_ha, -geometry ), by = c("uid_from" = "uid")) %>%
  rename(effective_patch_area_ha_from = effective_patch_area_ha) %>%
  left_join(patch_centroid_info %>% as.data.frame() %>% dplyr::select(uid, effective_patch_area_ha, -geometry ), by = c("uid_to" = "uid")) %>%
  rename(effective_patch_area_ha_to = effective_patch_area_ha)

## probability of connection ----
candidate_pairs$prob_connect_lcd <-  exp(-constants$alpha * candidate_pairs$lcd)
candidate_pairs$prob_connect_euclid <-  exp(-constants$alpha * candidate_pairs$euclid_dist)

## ECA value of each pair ----
candidate_pairs$eca_lcd_contribution = candidate_pairs$effective_patch_area_ha_from * candidate_pairs$effective_patch_area_ha_to * candidate_pairs$prob_connect_lcd
candidate_pairs$eca_euclidean_contribution = candidate_pairs$effective_patch_area_ha_from * candidate_pairs$effective_patch_area_ha_to * candidate_pairs$prob_connect_euclid
candidate_pairs$eca_value_ha_lcd = sqrt(candidate_pairs$eca_lcd_contribution)
candidate_pairs$eca_value_ha_euclidean = sqrt(candidate_pairs$eca_euclidean_contribution)


## sum patch contributions ----
### only the "from" which considers the contribution from patches outside the focal area, within the buffer

source_patch_centroid_info <- patch_centroid_info %>%
  filter(focal_landscape == 1) %>%
  # left_join(
  #   candidate_pairs %>% group_by(uid_to) %>% 
  #     summarise(eca_lcd_contribution_to = sum(eca_lcd_contribution, na.rm = TRUE),
  #               eca_euclidean_contribution_to = sum(eca_euclidean_contribution, na.rm = TRUE)),
  #   by = c("uid" = "uid_to")
  # ) %>%
  left_join(
    candidate_pairs %>% group_by(uid_from) %>% 
      summarise(eca_lcd_contribution_from = sum(eca_lcd_contribution, na.rm = TRUE),
                eca_euclidean_contribution_from = sum(eca_euclidean_contribution, na.rm = TRUE)),
    by = c("uid" = "uid_from")
  )

# source_patch_centroid_info$eca_value_ha_lcd_to <- sqrt(source_patch_centroid_info$eca_lcd_contribution_to)
source_patch_centroid_info$eca_value_ha_lcd_from <- sqrt(source_patch_centroid_info$eca_lcd_contribution_from)
# source_patch_centroid_info$eca_value_ha_euclidean_to <- sqrt(source_patch_centroid_info$eca_euclidean_contribution_to)
source_patch_centroid_info$eca_value_ha_euclidean_from <- sqrt(source_patch_centroid_info$eca_euclidean_contribution_from)

# LANDSCAPE SCALE METRICS ----

landscape_stats[[this.ts.num]][[this.year]]$landscape.metrics = 
  data.frame(name = this.tss[this.ts.num],
             year = this.year,
             lcd.ECA = sqrt(sum(source_patch_centroid_info$eca_lcd_contribution_from)),
             euclid.ECA = sqrt(sum(source_patch_centroid_info$eca_euclidean_contribution_from)),
             n.clumps = length(unique(source_patch_centroid_info$uid)),
             tot.patch.area.effective.ha = sum(source_patch_centroid_info$effective_patch_area_ha),
             tot.patch.ha = sum(source_patch_centroid_info$patch_area_ha),
             tot.aw.patch.ha = sum(source_patch_centroid_info$awi_tot_ha),
             tot.edge.patch.ha = sum(source_patch_centroid_info$edge_tot_ha),
             tot.awi_edge.ha = sum(source_patch_centroid_info$awi_edge),
             med.patch.ha = median(source_patch_centroid_info$patch_area_ha),
             mean.patch.ha = mean(source_patch_centroid_info$patch_area_ha),
             largest_patch_ha = max(source_patch_centroid_info$patch_area_ha),
             largest_patch_quality_weighted_ha = max(source_patch_centroid_info$effective_patch_area_ha)
             
             
             )

# sensitive.landscape.lcd.ECA = sqrt(sum(bl.patch.hexid.centroids$sensitive.lcd.contrib.to[bl.patch.hexid.centroids$in.landscape == T]))
# sensitive.landscape.scaled.lcd.ECA = sqrt(sum(bl.patch.hexid.centroids$sensitive.scaled.lcd.contrib.to[bl.patch.hexid.centroids$in.landscape == T]))
# sensitive.landscape.euclid.ECA = sqrt(sum(bl.patch.hexid.centroids$sensitive.euclid.contrib.to[bl.patch.hexid.centroids$in.landscape == T]))

# stats per grid cell ----

# summarise patch metrics to hex level
hex_summary <- source_patch_centroid_info %>%
  st_drop_geometry() %>%
  group_by(grid_id) %>%
  summarise(
    hex.lcd.eca     = sqrt(sum(eca_lcd_contribution_from, na.rm = TRUE)),
    hex.euclid.eca        = sqrt(sum(eca_euclidean_contribution_from , na.rm = TRUE)),
    
    hex.lcd.eca.contrib     = sum(eca_lcd_contribution_from, na.rm = TRUE),
    hex.euclid.eca.contrib        = sum(eca_euclidean_contribution_from , na.rm = TRUE),
    
    
    n.clumps              = n_distinct(uid),
    
    tot.patch.ha          = sum(patch_area_ha, na.rm = TRUE),
    tot.aw.patch.ha       = sum(awi_tot_ha, na.rm = TRUE),
    tot.edge.patch.ha     = sum(edge_tot_ha, na.rm = TRUE),
    tot.awedge.patch.ha   = sum(awi_edge, na.rm = TRUE),
    tot.effective.patch.ha  = sum(effective_patch_area_ha, na.rm = TRUE),
    
    .groups = "drop"
  )

ts.hexgrid <- ts.hexgrid %>%
  left_join(hex_summary, by = "grid_id") %>%
  mutate(
    across(
      c(
        hex.lcd.eca,
        hex.euclid.eca,
        hex.lcd.eca.contrib,
        hex.euclid.eca.contrib,
        n.clumps,
        tot.patch.ha,
        tot.aw.patch.ha,
        tot.edge.patch.ha,
        tot.awedge.patch.ha,
        tot.effective.patch.ha
      ),
      ~replace_na(.x, 0)
    )
  )

# RETIRED

# # sqrt of sum of all summed contributions of patches in the defined landscape by hex
# hex.eca.lcd =  aggregate(source_patch_centroid_info$eca_value_ha_lcd_from, # contrib to each patch
#                                           by = list(source_patch_centroid_info$grid_id),  # by hex
#                                           function(x){sqrt(sum(x))}) %>% rename(grid_id = Group.1, hex.eca = x)
# hex.eca.euclid =  aggregate(source_patch_centroid_info$eca_value_ha_euclidean_from, # contrib to each patch
#                                           by = list(source_patch_centroid_info$grid_id),  # by hex
#                                           function(x){sqrt(sum(x))}) %>% rename(grid_id = Group.1, hex.eca = x)
# hex.n.clumps =  aggregate(source_patch_centroid_info$uid,  # scaled (so mean landscaep cost ==1) contribution to patch
#                           by = list(source_patch_centroid_info$grid_id), function(x){ # over each grid
#                             length(unique(x))}) %>% rename(grid_id = Group.1, n.clumps = x)
# hex.tot.patch.ha = aggregate(source_patch_centroid_info$patch_area_ha,  # scaled (so mean landscaep cost ==1) contribution to patch
#                           by = list(source_patch_centroid_info$grid_id), function(x){ # over each grid
#                             sum(x)}) %>% rename(grid_id = Group.1, tot.patch.ha = x)  
# hex.tot.aw.patch.ha =  aggregate(source_patch_centroid_info$awi_tot_ha,  # scaled (so mean landscaep cost ==1) contribution to patch
#                               by = list(source_patch_centroid_info$grid_id), function(x){ # over each grid
#                                 sum(x)}) %>% rename(grid_id = Group.1, tot.aw.patch.ha = x)
# hex.tot.edge.patch.ha =  aggregate(source_patch_centroid_info$edge_tot_ha,  # scaled (so mean landscaep cost ==1) contribution to patch
#                               by = list(source_patch_centroid_info$grid_id), function(x){ # over each grid
#                                 sum(x)}) %>% rename(grid_id = Group.1, tot.edge.patch.ha = x)
# 
# hex.tot.awedge.patch.ha =  aggregate(source_patch_centroid_info$awi_edge,  # scaled (so mean landscaep cost ==1) contribution to patch
#                               by = list(source_patch_centroid_info$grid_id), function(x){ # over each grid
#                                 sum(x)}) %>% rename(grid_id = Group.1, tot.awedge.patch.ha = x)
# 
# # add summed contributions to hex shape object
# ts.hexgrid$hex.lcd.eca = left_join(ts.hexgrid, hex.sum.lcd.contrib.to)$hex.eca
# ts.hexgrid$hex.lcd.eca[is.na(ts.hexgrid$hex.lcd.eca)] = 0 # give all NAs (hexes with no contributions) 0 eca
# ts.hexgrid$hex.scaled.lcd.eca = left_join(ts.hexgrid, hex.sum.scaled.lcd.contrib.to)$hex.eca
# ts.hexgrid$hex.scaled.lcd.eca[is.na(ts.hexgrid$hex.scaled.lcd.eca)] = 0 # give all NAs (hexes with no contributions) 0 eca
# ts.hexgrid$hex.euclid.eca = left_join(ts.hexgrid, hex.sum.euclid.contrib.to)$hex.eca
# ts.hexgrid$hex.euclid.eca[is.na(ts.hexgrid$hex.euclid.eca)] = 0 # give all NAs (hexes with no contributions) 0 eca
# 
# ts.hexgrid$n.clumps = left_join(ts.hexgrid, hex.n.clumps)$n.clumps
# ts.hexgrid$n.clumps[is.na(ts.hexgrid$n.clumps)] = 0 # give all NAs (hexes with no contributions) 0 eca
# ts.hexgrid$tot.patch.ha = left_join(ts.hexgrid, hex.tot.patch.ha)$tot.patch.ha
# ts.hexgrid$tot.patch.ha[is.na(ts.hexgrid$tot.patch.ha)] = 0 # give all NAs (hexes with no contributions) 0 eca
# ts.hexgrid$tot.aw.patch.ha = left_join(ts.hexgrid, hex.tot.aw.patch.ha)$tot.aw.patch.ha
# ts.hexgrid$tot.aw.patch.ha[is.na(ts.hexgrid$tot.aw.patch.ha)] = 0 # give all NAs (hexes with no contributions) 0 eca
# ts.hexgrid$tot.edge.patch.ha = left_join(ts.hexgrid, hex.tot.edge.patch.ha)$tot.edge.patch.ha
# ts.hexgrid$tot.edge.patch.ha[is.na(ts.hexgrid$tot.edge.patch.ha)] = 0 # give all NAs (hexes with no contributions) 0 eca
# ts.hexgrid$tot.awedge.patch.ha = left_join(ts.hexgrid, hex.tot.awedge.patch.ha)$tot.awedge.patch.ha
# ts.hexgrid$tot.awedge.patch.ha[is.na(ts.hexgrid$tot.awedge.patch.ha)] = 0 # give all NAs (hexes with no contributions) 0 eca

# ts.hexgrid$sensitive.hex.lcd.eca = left_join(ts.hexgrid, sensitive.hex.sum.lcd.contrib.to)$hex.eca
# ts.hexgrid$sensitive.hex.lcd.eca[is.na(ts.hexgrid$sensitive.hex.lcd.eca)] = 0 # give all NAs (hexes with no contributions) 0 eca
# ts.hexgrid$sensitive.hex.scaled.lcd.eca = left_join(ts.hexgrid, sensitive.hex.sum.scaled.lcd.contrib.to)$hex.eca
# ts.hexgrid$sensitive.hex.scaled.lcd.eca[is.na(ts.hexgrid$sensitive.hex.scaled.lcd.eca)] = 0 # give all NAs (hexes with no contributions) 0 eca
# ts.hexgrid$sensitive.hex.euclid.eca = left_join(ts.hexgrid, sensitive.hex.sum.euclid.contrib.to)$hex.eca
# ts.hexgrid$sensitive.hex.euclid.eca[is.na(ts.hexgrid$sensitive.hex.euclid.eca)] = 0 # give all NAs (hexes with no contributions) 0 eca

# make uniform areas for hexes with different split landmasses

if(length(unique(ts.hexgrid$grid_id[duplicated(ts.hexgrid$grid_id)])) > 0){
  print(paste0("WARNING: ", length(unique(ts.hexgrid$grid_id[duplicated(ts.hexgrid$grid_id)])), " hexes have been split into multiple polygons. Summing areas and patch counts for these hexes"))
}

for(i in 1: length(unique(ts.hexgrid$grid_id[duplicated(ts.hexgrid$grid_id)]))){
  ts.hexgrid$hex.ha[ts.hexgrid$grid_id == unique(ts.hexgrid$grid_id[duplicated(ts.hexgrid$grid_id)])[i] ] = ts.hexgrid$hex.ha[ts.hexgrid$grid_id == unique(ts.hexgrid$grid_id[duplicated(ts.hexgrid$grid_id)])[i] ] %>% sum()
  ts.hexgrid$lcm.ncells[ts.hexgrid$grid_id == unique(ts.hexgrid$grid_id[duplicated(ts.hexgrid$grid_id)])[i] ] = ts.hexgrid$lcm.ncells[ts.hexgrid$grid_id == unique(ts.hexgrid$grid_id[duplicated(ts.hexgrid$grid_id)])[i] ] %>% sum()
  ts.hexgrid$bl.ncells[ts.hexgrid$grid_id == unique(ts.hexgrid$grid_id[duplicated(ts.hexgrid$grid_id)])[i] ] = ts.hexgrid$bl.ncells[ts.hexgrid$grid_id == unique(ts.hexgrid$grid_id[duplicated(ts.hexgrid$grid_id)])[i] ] %>% sum()
  ts.hexgrid$landnotcoastal.ncells[ts.hexgrid$grid_id == unique(ts.hexgrid$grid_id[duplicated(ts.hexgrid$grid_id)])[i] ] = ts.hexgrid$landnotcoastal.ncells[ts.hexgrid$grid_id == unique(ts.hexgrid$grid_id[duplicated(ts.hexgrid$grid_id)])[i] ] %>% sum()
  
}

#standardise by terrestrial area of hex - i.e. this is eca if hex was complete  - max number of cells in any hex/ number of non-coastal land cells in this hex
ts.hexgrid$hex.standardised.lcd.eca = ts.hexgrid$hex.lcd.eca *(max(ts.hexgrid$non_coastal_land_area)/ts.hexgrid$non_coastal_land_area)
ts.hexgrid$hex.standardised.lcd.eca[is.na(ts.hexgrid$hex.standardised.lcd.eca)] = 0

ts.hexgrid$hex.standardised.euclid.eca = ts.hexgrid$hex.euclid.eca *(max(ts.hexgrid$non_coastal_land_area)/ts.hexgrid$non_coastal_land_area)
ts.hexgrid$hex.standardised.euclid.eca[is.na(ts.hexgrid$hex.standardised.euclid.eca)] = 0




# ts.hexgrid$sensitive.hex.standardised.lcd.eca = ts.hexgrid$sensitive.hex.lcd.eca *(max(ts.hexgrid$lcm.ncells)/ts.hexgrid$landnotcoastal.ncells)
# ts.hexgrid$sensitive.hex.standardised.lcd.eca[is.na(ts.hexgrid$sensitive.hex.standardised.lcd.eca)] = 0
# ts.hexgrid$sensitive.hex.standardised.scaled.lcd.eca = ts.hexgrid$sensitive.hex.scaled.lcd.eca *(max(ts.hexgrid$lcm.ncells)/ts.hexgrid$landnotcoastal.ncells)
# ts.hexgrid$sensitive.hex.standardised.scaled.lcd.eca[is.na(ts.hexgrid$sensitive.hex.standardised.scaled.lcd.eca)] = 0
# ts.hexgrid$sensitive.hex.standardised.euclid.eca = ts.hexgrid$sensitive.hex.euclid.eca *(max(ts.hexgrid$lcm.ncells)/ts.hexgrid$landnotcoastal.ncells)
# ts.hexgrid$sensitive.hex.standardised.euclid.eca[is.na(ts.hexgrid$sensitive.hex.standardised.euclid.eca)] = 0

st_write(ts.hexgrid, paste0(func.conect.path, "\\analysis outputs\\", this.tss[this.ts.num], "\\", this.year, "\\lcmres", constants$lcm.res, "\\r_funcconnect_hexgrid.gpkg"), delete_layer = TRUE)
st_write(source_patch_centroid_info, paste0(func.conect.path, "\\analysis outputs\\", this.tss[this.ts.num], "\\", this.year, "\\lcmres", constants$lcm.res, "\\r_funcconnect_patch_centroids.gpkg"), delete_layer = TRUE)
# save ----
save(ts.hexgrid,
     landscape_stats, 
    # sensitive.landscape.lcd.ECA, sensitive.landscape.scaled.lcd.ECA, sensitive.landscape.euclid.ECA,
    source_patch_centroid_info,
     file = 
      paste0(func.conect.path, 
             "\\analysis outputs\\", this.tss[this.ts.num], "\\", this.year, "\\lcmres", constants$lcm.res, 
             "\\06_r_funcconnect_EffectiveAreas_ECAobs_.RData")
)

print("Equivelent Connected Area calculations (script06) done")



# testing

