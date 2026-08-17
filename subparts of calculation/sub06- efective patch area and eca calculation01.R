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
candidate_pairs$prob_connect_lcd <-  exp(-constants$alpha * (candidate_pairs$lcd/constants$magic.lcd.refactoror ))
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

landscape_stats[[this.tss[this.ts.num]]][[this.year]]$landscape.summary = 
  data.frame(name = this.tss[this.ts.num],
             year = this.year,
             
             lcd.ECA = sqrt(sum(source_patch_centroid_info$eca_lcd_contribution_from)),
             euclid.ECA = sqrt(sum(source_patch_centroid_info$eca_euclidean_contribution_from)),

             n.clumps = length(unique(source_patch_centroid_info$uid)),
             
             sum.patch.area.effective.ha = sum(source_patch_centroid_info$effective_patch_area_ha),
             sum.patch.ha = sum(source_patch_centroid_info$patch_area_ha),
             sum.aw.patch.ha = sum(source_patch_centroid_info$awi_tot_ha),
             sum.edge.patch.ha = sum(source_patch_centroid_info$edge_tot_ha),
             sum.awi_edge.ha = sum(source_patch_centroid_info$awi_edge),
             
             med.patch.ha = median(source_patch_centroid_info$patch_area_ha),
             mean.patch.ha = mean(source_patch_centroid_info$patch_area_ha),
             
             largest_patch_ha = max(source_patch_centroid_info$patch_area_ha),
             largest_patch_quality_weighted_ha = max(source_patch_centroid_info$effective_patch_area_ha),
             
             landscape.ecolog.cost.not.sea_mean = landscape_stats[[this.tss[this.ts.num]]][[this.year]]$landscape.ecolog.cost.not.sea_mean,
             landscape.ecolog.cost.not.sea_median = landscape_stats[[this.tss[this.ts.num]]][[this.year]]$landscape.ecolog.cost.not.sea_median,
             landscape.ecolog.cost.not.sea_p70 = landscape_stats[[this.tss[this.ts.num]]][[this.year]]$landscape.ecolog.cost.not.sea_p70
             
             )

landscape_stats[[this.tss[this.ts.num]]][[this.year]]$model_params = list(
  
  
  lcm_res = constants$lcm.res,
  buffer.for.patchid = constants$buffer.for.patchid,
  
  relative.edge.quality = constants$relative.edge.quality,
  non.awi.qual.eff = constants$non.awi.qual.eff,
  awi.qual.eff = constants$awi.qual.eff,
  
  alpha = constants$alpha,
  costs = dispers.costs,
  focalhab.cost = constants$focalhab.cost,
  cost.scale.factor = constants$cost.scale.factor,
  agg_cost_surface = constants$agg_cost_surface,
  cost.res = constants$cost.res
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
    
    sum.patch.area.effective.ha  = sum(effective_patch_area_ha, na.rm = TRUE),
    sum.patch.ha = sum(patch_area_ha, na.rm = TRUE),
    sum.aw.patch.ha = sum(awi_tot_ha, na.rm = TRUE),
    sum.edge.patch.ha = sum(edge_tot_ha, na.rm = TRUE),
    sum.awi_edge.ha = sum(awi_edge, na.rm = TRUE),
    
    med.patch.ha = median(patch_area_ha, na.rm = TRUE),
    mean.patch.ha = mean(patch_area_ha, na.rm = TRUE),
    
    largest_patch_ha = max(patch_area_ha, na.rm = TRUE),
    largest_patch_quality_weighted_ha = max(effective_patch_area_ha, na.rm = TRUE),

    .groups = "drop"
  )

### deal with split-hexes with duplicated grid_ids ----

dup_ids <- unique(ts.hexgrid$grid_id[duplicated(ts.hexgrid$grid_id)])
lcm_cols <- grep("^lcm_", names(ts.hexgrid), value = TRUE)

if(length(dup_ids) > 0){
  print(paste0("WARNING: ", length(dup_ids), " hexes have been split into multiple polygons. Summing areas and patch counts for these hexes"))
}

for(i in seq_along(dup_ids)) {
  rows <- ts.hexgrid$grid_id == dup_ids[i]
  
  ts.hexgrid$hex.ha[rows] <- sum(ts.hexgrid$hex.ha[rows], na.rm = TRUE)
  ts.hexgrid$total_area[rows] <- sum(ts.hexgrid$total_area[rows], na.rm = TRUE)
  
  ts.hexgrid[rows, lcm_cols] <-
    matrix(colSums(ts.hexgrid[rows, lcm_cols], na.rm = TRUE),
           nrow = sum(rows), 
           ncol = length(lcm_cols), byrow = TRUE)
}

### deal with NA entries


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
        sum.patch.area.effective.ha,
        sum.patch.ha,
        sum.aw.patch.ha,
        sum.edge.patch.ha,
        sum.awi_edge.ha,
        med.patch.ha,
        mean.patch.ha,
        largest_patch_ha,
        largest_patch_quality_weighted_ha      ),
      ~replace_na(.x, 0)
    )
  )


#standardise by terrestrial area of hex - i.e. this is eca if hex was complete  - max number of cells in any hex/ number of non-coastal land cells in this hex
ts.hexgrid$hex.standardised.lcd.eca = ts.hexgrid$hex.lcd.eca *(max(ts.hexgrid$noncoastal_water_area)/ts.hexgrid$noncoastal_water_area)
ts.hexgrid$hex.standardised.lcd.eca[is.na(ts.hexgrid$hex.standardised.lcd.eca)] = 0

ts.hexgrid$hex.standardised.euclid.eca = ts.hexgrid$hex.euclid.eca *(max(ts.hexgrid$noncoastal_water_area)/ts.hexgrid$noncoastal_water_area)
ts.hexgrid$hex.standardised.euclid.eca[is.na(ts.hexgrid$hex.standardised.euclid.eca)] = 0



landscape_stats[[this.tss[this.ts.num]]][[this.year]]$ts.hexgrid = ts.hexgrid
landscape_stats[[this.tss[this.ts.num]]][[this.year]]$source_patch_centroid_info = source_patch_centroid_info
landscape_stats[[this.tss[this.ts.num]]][[this.year]]$candidate_pairs = candidate_pairs


st_write(ts.hexgrid, paste0(func.conect.path, "\\analysis outputs\\", this.tss[this.ts.num], "\\", this.year, "\\lcmres", constants$lcm.res, "\\r_funcconnect_hexgrid.gpkg"), delete_layer = TRUE)
st_write(source_patch_centroid_info, paste0(func.conect.path, "\\analysis outputs\\", this.tss[this.ts.num], "\\", this.year, "\\lcmres", constants$lcm.res, "\\r_funcconnect_patch_centroids.gpkg"), delete_layer = TRUE)
# save ----
save(landscape_stats, 
     file = 
      paste0(func.conect.path, 
             "\\analysis outputs\\", this.tss[this.ts.num], "\\", this.year, "\\lcmres", constants$lcm.res, 
             "\\06_r_funcconnect_EffectiveAreas_ECAobs_.RData")
)

print("Equivelent Connected Area calculations (script06) done")



# testing

