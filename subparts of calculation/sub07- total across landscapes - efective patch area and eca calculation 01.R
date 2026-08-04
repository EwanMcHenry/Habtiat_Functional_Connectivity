##------ Wed Jul 20 12:32:04 2022 ------##
# Ewan McHenry
# calculation of ECA for entirity of all landscapes considered, then combine hexgrids for each ts considered
# based heavily on subscript 06
# assums no overlap in landscapes and no contribution of one landscape to antoher
# adds a connectivity metric for all landscapes considered I think. I dont really like this and am removing it from the master as of ##------ Thu Aug 22 16:09:20 2024 ------##
# consider it unnesseary


# LOAD AND COMBINE PATCH INFO FROM ALL individual LANDSCAPES ----
# load data needed for landscapes in loop
individ.ts.outputs = list (NA)
  for(y in seq_along(years.considered)){
    this.year = years.considered[y]
    load(paste0(func.conect.path, 
                "\\analysis outputs\\", this.tss[this.ts.num], "\\", this.year, "\\lcmres", constants$lcm.res, 
                "\\06_r_funcconnect_EffectiveAreas_ECAobs_.RData"))
    # give an identifier for the landscape stuff came from
    ts.hexgrid$ts.name = this.tss
    source_patch_centroid_info$ts.name = this.tss
    # hexgrid with stuff 
    individ.ts.outputs[[this.tss]][[this.year]] = list(name = this.tss, 
                                                       year = this.year,
                                                       alpha = constants$alpha,
                                                       costs = dispers.costs,
                                                       cost.scale.factor = constants$cost.scale.factor,
                                                       lcm_res = constants$lcm.res,
                                                       hex_grid = ts.hexgrid,
                                                       source_patch_centroid_info = source_patch_centroid_info,
                                                       landscape_stats = c(
                                                         name = this.tss,
                                                         year = this.year,
                                                         lcd.ECA = landscape_stats[[this.ts.num]][[this.year]]$lcd.ECA,
                                                         euclid.ECA = landscape_stats[[this.ts.num]][[this.year]]$euclid.ECA,
                                                         

                                                         n.patch = landscape_stats[[this.ts.num]][[this.year]]$n.clumps,
                                                         quality_weighted_habitat_ha_sum = landscape_stats[[this.ts.num]][[this.year]]$tot.patch.area.effective.ha,
                                                         habitat_ha_sum = landscape_stats[[this.ts.num]][[this.year]]$tot.patch.ha,
                                                         aw_habitat_ha_sum = landscape_stats[[this.ts.num]][[this.year]]$tot.aw.patch.ha,
                                                         edge_habitat_ha_sum = landscape_stats[[this.ts.num]][[this.year]]$tot.edge.patch.ha,
                                                         awi_edge_ha_sum = landscape_stats[[this.ts.num]][[this.year]]$tot.awi_edge.ha,
                                                         
                                                         med.patch.ha = landscape_stats[[this.ts.num]][[this.year]]$med.patch.ha,
                                                         mean.patch.ha = landscape_stats[[this.ts.num]][[this.year]]$mean.patch.ha,
                                                         largest_patch_ha = landscape_stats[[this.ts.num]][[this.year]]$largest_patch_ha,
                                                         largest_patch_quality_weighted_ha =landscape_stats[[this.ts.num]][[this.year]]$largest_patch_quality_weighted_ha,
                                                         
                                                         
                                                         ecolog.cost.not.sea_mean = landscape_stats[[this.ts.num]][[this.year]]$landscape.ecolog.cost.not.sea_mean,
                                                         ecolog.cost.not.sea_median = landscape_stats[[this.ts.num]][[this.year]]$landscape.ecolog.cost.not.sea_median,
                                                         ecolog.cost.not.sea_p70 = landscape_stats[[this.ts.num]][[this.year]]$landscape.ecolog.cost.not.sea_p70
                                                         )
                                                       )
  }

# SET UP FOR LOOP FOR YEARS
for (this.year in years.considered) {
  # pull together patch info from  all landscapes from that year
  source_patch_centroid_info = do.call(rbind,
                                     map(individ.ts.outputs, "source_patch_centroid_info")[map(individ.ts.outputs, "year") %>% unlist() == this.year])
  ts.hexgrid = do.call(rbind,
                       map(individ.ts.outputs, "ts.hexgrid")[map(individ.ts.outputs, "year") %>% unlist() == this.year])
  # combine hexgrid info from same hexgrid from different treescapes ----
  # ## arrg i give up (weighted median), its an idea for future version!
  # ts.hexgrid = ts.hexgrid[,names(ts.hexgrid) %in% c("grid_id", "hex.ha", "lcm.ncells", "bl.ncells", "landnotcoastal.ncells",
  #                                                   "mean.scaled.ecolog.cost.not.sea", "median.scaled.ecolog.cost.not.sea", 
  #                                                   "ts.name")]
  # boardering.to.combine =  ts.hexgrid$grid_id %in% ts.hexgrid$grid_id[duplicated(ts.hexgrid$grid_id) & !duplicated(paste(ts.hexgrid$grid_id, ts.hexgrid$ts.name))] # grid ids that have matching same grid id, but in different ts 
  # temp.hexgrid1 = ts.hexgrid[!boardering.to.combine,]
  # temp.hexgrid2 = ts.hexgrid[0,]
  # temp.hexgrid3 = ts.hexgrid[boardering.to.combine,] %>%
  #   group_by(grid_id) %>% 
  #   summarize()
  # 
  # for ( ii seq_along(unique(ts.hexgrid$grid_id[boardering.to.combine])))){
  #   these.polys =  ts.hexgrid$grid_id %in% unique(ts.hexgrid$grid_id[boardering.to.combine])[ii]
  #   temp.hexgrid3$hex.ha = sum(ts.hexgrid$hex.ha[these.polys] )
  #   temp.hexgrid3$lcm.ncells = sum(ts.hexgrid$lcm.ncells[these.polys] )
  #   temp.hexgrid3$bl.ncells = sum(ts.hexgrid$bl.ncells[these.polys] )
  #   temp.hexgrid3$landnotcoastal.ncells = sum(ts.hexgrid$landnotcoastal.ncells[these.polys] )
  #   temp.hexgrid3$mean.scaled.ecolog.cost.not.sea = weighted.mean(ts.hexgrid$mean.scaled.ecolog.cost.not.sea[these.polys] , ts.hexgrid$landnotcoastal.ncells[these.polys])
  #   temp.hexgrid3$median.scaled.ecolog.cost.not.sea = weighted.median(ts.hexgrid$mean.scaled.ecolog.cost.not.sea[these.polys] , ts.hexgrid$landnotcoastal.ncells[these.polys])
  #   temp.hexgrid3$landnotcoastal.ncells = sum(ts.hexgrid$landnotcoastal.ncells[these.polys] )
  #   
  #                              
  #   
  # } ----
  

  # LANDSCAPE SCALE METRICS ----
  # square root of sum of all summed contributions of patches in the defined landscape
  
  landscape.metrics = data.frame(name = ts.andAll.lcm.names[length(ts.andAll.lcm.names)],
                                 year = this.year,
                                 leastcost.ECA =         sqrt(sum(source_patch_centroid_info$leastcost.contrib.to[source_patch_centroid_info$in.landscape == T])),
                                 # scaled.leastcost.ECA =  sqrt(sum(source_patch_centroid_info$scaled.leastcost.contrib.to[source_patch_centroid_info$in.landscape == T])),
                                 euclid.ECA = sqrt(sum(source_patch_centroid_info$euclid.contrib.to[source_patch_centroid_info$in.landscape == T])),
                                 n.clumps = length(unique(source_patch_centroid_info$clumps[source_patch_centroid_info$in.landscape == T])),
                                 tot.patch.ha = sum(source_patch_centroid_info$subpatch.hex.ha [source_patch_centroid_info$in.landscape == T]),
                                 tot.aw.patch.ha = sum(source_patch_centroid_info$awi.subpatch.hex.ha [source_patch_centroid_info$in.landscape == T]),
                                 tot.edge.patch.ha = sum(source_patch_centroid_info$edge.subpatch.hex.ha [source_patch_centroid_info$in.landscape == T]),
                                 tot.awedge.patch.ha = sum(source_patch_centroid_info$edge.awi.subpatch.hex.ha [source_patch_centroid_info$in.landscape == T]),
                                 med.patch.ha = median(source_patch_centroid_info$clump.ha[!duplicated(source_patch_centroid_info$clumps)])
  )
  landscape.metrics$mean.patch.ha = landscape.metrics$tot.patch.ha/ landscape.metrics$n.clumps
  
  # sensitive.landscape.leastcost.ECA = sqrt(sum(source_patch_centroid_info$sensitive.leastcost.contrib.to[source_patch_centroid_info$in.landscape == T]))
  # sensitive.landscape.scaled.leastcost.ECA = sqrt(sum(source_patch_centroid_info$sensitive.scaled.leastcost.contrib.to[source_patch_centroid_info$in.landscape == T]))
  # sensitive.landscape.euclid.ECA = sqrt(sum(source_patch_centroid_info$sensitive.euclid.contrib.to[source_patch_centroid_info$in.landscape == T]))
  
  # # ECA per grid cell ----
  # # sqrt of sum of all summed contributions of patches in the defined landscape by hex
  # hex.sum.leastcost.contrib.to =  aggregate(source_patch_centroid_info$leastcost.contrib.to[source_patch_centroid_info$in.landscape == T], # contrib to each patch
  #                                           by = list(source_patch_centroid_info$grid_id[source_patch_centroid_info$in.landscape == T]),  # by hex
  #                                           function(x){sqrt(sum(x))}) %>% rename(grid_id = Group.1, hex.eca = x)
  # hex.sum.scaled.leastcost.contrib.to =  aggregate(source_patch_centroid_info$scaled.leastcost.contrib.to[source_patch_centroid_info$in.landscape == T],  # scaled (so mean landscaep cost ==1) contribution to patch
  #                                                  by = list(source_patch_centroid_info$grid_id[source_patch_centroid_info$in.landscape == T]), function(x){ # over each grid
  #                                                    sqrt(sum(x))}) %>% rename(grid_id = Group.1, hex.eca = x)
  # hex.sum.euclid.contrib.to =  aggregate(source_patch_centroid_info$euclid.contrib.to[source_patch_centroid_info$in.landscape == T], 
  #                                        by = list(source_patch_centroid_info$grid_id[source_patch_centroid_info$in.landscape == T]), function(x){
  #                                          sqrt(sum(x))}) %>% rename(grid_id = Group.1, hex.eca = x)
  # hex.n.clumps =  aggregate(source_patch_centroid_info$clumps[source_patch_centroid_info$in.landscape == T],  # scaled (so mean landscaep cost ==1) contribution to patch
  #                           by = list(source_patch_centroid_info$grid_id[source_patch_centroid_info$in.landscape == T]), function(x){ # over each grid
  #                             length(unique(x))}) %>% rename(grid_id = Group.1, n.clumps = x)
  # hex.tot.patch.ha =  aggregate(source_patch_centroid_info$subpatch.hex.ha[source_patch_centroid_info$in.landscape == T],  # scaled (so mean landscaep cost ==1) contribution to patch
  #                               by = list(source_patch_centroid_info$grid_id[source_patch_centroid_info$in.landscape == T]), function(x){ # over each grid
  #                                 sum(x)}) %>% rename(grid_id = Group.1, tot.patch.ha = x)
  # hex.tot.aw.patch.ha =  aggregate(source_patch_centroid_info$awi.subpatch.hex.ha[source_patch_centroid_info$in.landscape == T],  # scaled (so mean landscaep cost ==1) contribution to patch
  #                                  by = list(source_patch_centroid_info$grid_id[source_patch_centroid_info$in.landscape == T]), function(x){ # over each grid
  #                                    sum(x)}) %>% rename(grid_id = Group.1, tot.aw.patch.ha = x)
  # hex.tot.edge.patch.ha =  aggregate(source_patch_centroid_info$edge.subpatch.hex.ha[source_patch_centroid_info$in.landscape == T],  # scaled (so mean landscaep cost ==1) contribution to patch
  #                                    by = list(source_patch_centroid_info$grid_id[source_patch_centroid_info$in.landscape == T]), function(x){ # over each grid
  #                                      sum(x)}) %>% rename(grid_id = Group.1, tot.edge.patch.ha = x)
  # hex.tot.awedge.patch.ha =  aggregate(source_patch_centroid_info$edge.awi.subpatch.hex.ha[source_patch_centroid_info$in.landscape == T],  # scaled (so mean landscaep cost ==1) contribution to patch
  #                                      by = list(source_patch_centroid_info$grid_id[source_patch_centroid_info$in.landscape == T]), function(x){ # over each grid
  #                                        sum(x)}) %>% rename(grid_id = Group.1, tot.awedge.patch.ha = x)
  # 
  # 
  # 
  # # sensitive.hex.sum.leastcost.contrib.to =  aggregate(source_patch_centroid_info$sensitive.leastcost.contrib.to[source_patch_centroid_info$in.landscape == T], 
  # #                                           by = list(source_patch_centroid_info$grid_id[source_patch_centroid_info$in.landscape == T]), function(x){
  # #                                             sqrt(sum(x))}) %>% 
  # #   rename(grid_id = Group.1, hex.eca = x)
  # # sensitive.hex.sum.scaled.leastcost.contrib.to =  aggregate(source_patch_centroid_info$sensitive.leastcost.contrib.to[source_patch_centroid_info$in.landscape == T], 
  # #                                                  by = list(source_patch_centroid_info$grid_id[source_patch_centroid_info$in.landscape == T]), function(x){
  # #                                                    sqrt(sum(x))}) %>% 
  # #   rename(grid_id = Group.1, hex.eca = x)
  # # sensitive.hex.sum.euclid.contrib.to =  aggregate(source_patch_centroid_info$sensitive.euclid.contrib.to[source_patch_centroid_info$in.landscape == T], 
  # #                                        by = list(source_patch_centroid_info$grid_id[source_patch_centroid_info$in.landscape == T]), function(x){
  # #                                          sqrt(sum(x))}) %>% 
  # #   rename(grid_id = Group.1, hex.eca = x)
  # 
  # # add summed contributions to hex shape object
  # ts.hexgrid$hex.leastcost.eca = left_join(ts.hexgrid, hex.sum.leastcost.contrib.to)$hex.eca
  # ts.hexgrid$hex.leastcost.eca[is.na(ts.hexgrid$hex.leastcost.eca)] = 0 # give all NAs (hexes with no contributions) 0 eca
  # ts.hexgrid$hex.scaled.leastcost.eca = left_join(ts.hexgrid, hex.sum.scaled.leastcost.contrib.to)$hex.eca
  # ts.hexgrid$hex.scaled.leastcost.eca[is.na(ts.hexgrid$hex.scaled.leastcost.eca)] = 0 # give all NAs (hexes with no contributions) 0 eca
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
  # 
  # # ts.hexgrid$sensitive.hex.leastcost.eca = left_join(ts.hexgrid, sensitive.hex.sum.leastcost.contrib.to)$hex.eca
  # # ts.hexgrid$sensitive.hex.leastcost.eca[is.na(ts.hexgrid$sensitive.hex.leastcost.eca)] = 0 # give all NAs (hexes with no contributions) 0 eca
  # # ts.hexgrid$sensitive.hex.scaled.leastcost.eca = left_join(ts.hexgrid, sensitive.hex.sum.scaled.leastcost.contrib.to)$hex.eca
  # # ts.hexgrid$sensitive.hex.scaled.leastcost.eca[is.na(ts.hexgrid$sensitive.hex.scaled.leastcost.eca)] = 0 # give all NAs (hexes with no contributions) 0 eca
  # # ts.hexgrid$sensitive.hex.euclid.eca = left_join(ts.hexgrid, sensitive.hex.sum.euclid.contrib.to)$hex.eca
  # # ts.hexgrid$sensitive.hex.euclid.eca[is.na(ts.hexgrid$sensitive.hex.euclid.eca)] = 0 # give all NAs (hexes with no contributions) 0 eca
  # 
  # # make uniform areas for hexes with different split landmasses
  # for(i in 1: length(unique(ts.hexgrid$grid_id[duplicated(ts.hexgrid$grid_id)]))){
  #   ts.hexgrid$hex.ha[ts.hexgrid$grid_id == unique(ts.hexgrid$grid_id[duplicated(ts.hexgrid$grid_id)])[i] ] = ts.hexgrid$hex.ha[ts.hexgrid$grid_id == unique(ts.hexgrid$grid_id[duplicated(ts.hexgrid$grid_id)])[i] ] %>% sum()
  #   ts.hexgrid$lcm.ncells[ts.hexgrid$grid_id == unique(ts.hexgrid$grid_id[duplicated(ts.hexgrid$grid_id)])[i] ] = ts.hexgrid$lcm.ncells[ts.hexgrid$grid_id == unique(ts.hexgrid$grid_id[duplicated(ts.hexgrid$grid_id)])[i] ] %>% sum()
  #   ts.hexgrid$bl.ncells[ts.hexgrid$grid_id == unique(ts.hexgrid$grid_id[duplicated(ts.hexgrid$grid_id)])[i] ] = ts.hexgrid$bl.ncells[ts.hexgrid$grid_id == unique(ts.hexgrid$grid_id[duplicated(ts.hexgrid$grid_id)])[i] ] %>% sum()
  #   ts.hexgrid$landnotcoastal.ncells[ts.hexgrid$grid_id == unique(ts.hexgrid$grid_id[duplicated(ts.hexgrid$grid_id)])[i] ] = ts.hexgrid$landnotcoastal.ncells[ts.hexgrid$grid_id == unique(ts.hexgrid$grid_id[duplicated(ts.hexgrid$grid_id)])[i] ] %>% sum()
  #   
  # }
  # 
  # #standardise by terrestrial area of hex - i.e. this is eca if hex was complete  - max number of cells in any hex/ number of non-coastal land cells in this hex
  # ts.hexgrid$hex.standardised.leastcost.eca = ts.hexgrid$hex.leastcost.eca *(max(ts.hexgrid$lcm.ncells)/ts.hexgrid$landnotcoastal.ncells)
  # ts.hexgrid$hex.standardised.leastcost.eca[is.na(ts.hexgrid$hex.standardised.leastcost.eca)] = 0
  # 
  # ts.hexgrid$hex.standardised.scaled.leastcost.eca = ts.hexgrid$hex.scaled.leastcost.eca *(max(ts.hexgrid$lcm.ncells)/ts.hexgrid$landnotcoastal.ncells)
  # ts.hexgrid$hex.standardised.scaled.leastcost.eca[is.na(ts.hexgrid$hex.standardised.scaled.leastcost.eca)] = 0
  # 
  # ts.hexgrid$hex.standardised.euclid.eca = ts.hexgrid$hex.euclid.eca *(max(ts.hexgrid$lcm.ncells)/ts.hexgrid$landnotcoastal.ncells)
  # ts.hexgrid$hex.standardised.euclid.eca[is.na(ts.hexgrid$hex.standardised.euclid.eca)] = 0
  # 
  # 
  # 
  
  # ts.hexgrid$sensitive.hex.standardised.leastcost.eca = ts.hexgrid$sensitive.hex.leastcost.eca *(max(ts.hexgrid$lcm.ncells)/ts.hexgrid$landnotcoastal.ncells)
  # ts.hexgrid$sensitive.hex.standardised.leastcost.eca[is.na(ts.hexgrid$sensitive.hex.standardised.leastcost.eca)] = 0
  # ts.hexgrid$sensitive.hex.standardised.scaled.leastcost.eca = ts.hexgrid$sensitive.hex.scaled.leastcost.eca *(max(ts.hexgrid$lcm.ncells)/ts.hexgrid$landnotcoastal.ncells)
  # ts.hexgrid$sensitive.hex.standardised.scaled.leastcost.eca[is.na(ts.hexgrid$sensitive.hex.standardised.scaled.leastcost.eca)] = 0
  # ts.hexgrid$sensitive.hex.standardised.euclid.eca = ts.hexgrid$sensitive.hex.euclid.eca *(max(ts.hexgrid$lcm.ncells)/ts.hexgrid$landnotcoastal.ncells)
  # ts.hexgrid$sensitive.hex.standardised.euclid.eca[is.na(ts.hexgrid$sensitive.hex.standardised.euclid.eca)] = 0
  
  # save ----
  save(ts.hexgrid,
       landscape.metrics, 
       # sensitive.landscape.leastcost.ECA, sensitive.landscape.scaled.leastcost.ECA, sensitive.landscape.euclid.ECA,
       source_patch_centroid_info,
       file = 
         paste0(func.conect.path, 
                "\\analysis outputs\\All treescapes\\", this.year, "\\r_funcconnect_EffectiveAreas_ECAobs_.RData")
  )
}
print("Total Lanscapes Equivelent Connected Area calculations (script07) done")
