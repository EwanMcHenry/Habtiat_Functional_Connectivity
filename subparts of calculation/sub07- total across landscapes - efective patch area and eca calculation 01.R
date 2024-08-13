##------ Wed Jul 20 12:32:04 2022 ------##
# Ewan McHenry
# calculation of ECA for entirity of all landscapes considered, then combine hexgrids for each ts considered
# based heavily on subscript 06
# assums no overlap in landscapes and no contribution of one landscape to antoher

# make directory to save objects in
dir.create(paste0(func.conect.path, 
                  "\\analysis outputs\\", this.tss)
)

# LOAD AND COMBINE PATCH INFO FROM ALL individual LANDSCAPES ----
load(paste0(func.conect.path, 
            "\\analysis outputs\\", this.tss[1], "\\", years.considered[1], "\\r_funcconnect_EffectiveAreas_ECAobs_.RData"))
# load data needed for landscapes in loop
individ.ts.outputs = list (NA)
for(i in seq_along(this.tss)){
  for(y in seq_along(years.considered)){
    load(paste0(func.conect.path, 
                "\\analysis outputs\\", this.tss[i], "\\", years.considered[y], "\\r_funcconnect_EffectiveAreas_ECAobs_.RData"))
    # give an identifier for the landscape stuff came from
    ts.hexgrid$ts.name = this.tss[i]
    bl.patch.hexid.centroids$ts.name = this.tss[i]
    # hexgrid with stuff 
    individ.ts.outputs[[(i-1)*length(years.considered)+y]] = list(name = this.tss[i], 
                                                            year = years.considered[y],
                                                            ts.hexgrid = ts.hexgrid, 
                                                            bl.patch.hexid.centroids = bl.patch.hexid.centroids)
  }}

# SET UP FOR LOOP FOR YEARS
for (this.year in years.considered) {
  # pull together patch info from  all landscapes from that year
  bl.patch.hexid.centroids = do.call(rbind,
                                     map(individ.ts.outputs, "bl.patch.hexid.centroids")[map(individ.ts.outputs, "year") %>% unlist() == this.year])
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
  
  # make directory to save full landscape objects in
  dir.create(paste0(func.conect.path, 
                    "\\analysis outputs\\All treescapes\\", this.year)
  )
  
  # LANDSCAPE SCALE METRICS ----
  # square root of sum of all summed contributions of patches in the defined landscape
  
  landscape.metrics = data.frame(name = ts.andAll.lcm.names[length(ts.andAll.lcm.names)],
                                 year = this.year,
                                 leastcost.ECA =         sqrt(sum(bl.patch.hexid.centroids$leastcost.contrib.to[bl.patch.hexid.centroids$in.landscape == T])),
                                 scaled.leastcost.ECA =  sqrt(sum(bl.patch.hexid.centroids$scaled.leastcost.contrib.to[bl.patch.hexid.centroids$in.landscape == T])),
                                 euclid.ECA = sqrt(sum(bl.patch.hexid.centroids$euclid.contrib.to[bl.patch.hexid.centroids$in.landscape == T])),
                                 n.clumps = length(unique(bl.patch.hexid.centroids$clumps[bl.patch.hexid.centroids$in.landscape == T])),
                                 tot.patch.ha = sum(bl.patch.hexid.centroids$subpatch.hex.ha [bl.patch.hexid.centroids$in.landscape == T]),
                                 tot.aw.patch.ha = sum(bl.patch.hexid.centroids$awi.subpatch.hex.ha [bl.patch.hexid.centroids$in.landscape == T]),
                                 tot.edge.patch.ha = sum(bl.patch.hexid.centroids$edge.subpatch.hex.ha [bl.patch.hexid.centroids$in.landscape == T]),
                                 tot.awedge.patch.ha = sum(bl.patch.hexid.centroids$edge.awi.subpatch.hex.ha [bl.patch.hexid.centroids$in.landscape == T]),
                                 med.patch.ha = median(bl.patch.hexid.centroids$clump.ha[!duplicated(bl.patch.hexid.centroids$clumps)])
  )
  landscape.metrics$mean.patch.ha = landscape.metrics$tot.patch.ha/ landscape.metrics$n.clumps
  
  # sensitive.landscape.leastcost.ECA = sqrt(sum(bl.patch.hexid.centroids$sensitive.leastcost.contrib.to[bl.patch.hexid.centroids$in.landscape == T]))
  # sensitive.landscape.scaled.leastcost.ECA = sqrt(sum(bl.patch.hexid.centroids$sensitive.scaled.leastcost.contrib.to[bl.patch.hexid.centroids$in.landscape == T]))
  # sensitive.landscape.euclid.ECA = sqrt(sum(bl.patch.hexid.centroids$sensitive.euclid.contrib.to[bl.patch.hexid.centroids$in.landscape == T]))
  
  # # ECA per grid cell ----
  # # sqrt of sum of all summed contributions of patches in the defined landscape by hex
  # hex.sum.leastcost.contrib.to =  aggregate(bl.patch.hexid.centroids$leastcost.contrib.to[bl.patch.hexid.centroids$in.landscape == T], # contrib to each patch
  #                                           by = list(bl.patch.hexid.centroids$grid_id[bl.patch.hexid.centroids$in.landscape == T]),  # by hex
  #                                           function(x){sqrt(sum(x))}) %>% rename(grid_id = Group.1, hex.eca = x)
  # hex.sum.scaled.leastcost.contrib.to =  aggregate(bl.patch.hexid.centroids$scaled.leastcost.contrib.to[bl.patch.hexid.centroids$in.landscape == T],  # scaled (so mean landscaep cost ==1) contribution to patch
  #                                                  by = list(bl.patch.hexid.centroids$grid_id[bl.patch.hexid.centroids$in.landscape == T]), function(x){ # over each grid
  #                                                    sqrt(sum(x))}) %>% rename(grid_id = Group.1, hex.eca = x)
  # hex.sum.euclid.contrib.to =  aggregate(bl.patch.hexid.centroids$euclid.contrib.to[bl.patch.hexid.centroids$in.landscape == T], 
  #                                        by = list(bl.patch.hexid.centroids$grid_id[bl.patch.hexid.centroids$in.landscape == T]), function(x){
  #                                          sqrt(sum(x))}) %>% rename(grid_id = Group.1, hex.eca = x)
  # hex.n.clumps =  aggregate(bl.patch.hexid.centroids$clumps[bl.patch.hexid.centroids$in.landscape == T],  # scaled (so mean landscaep cost ==1) contribution to patch
  #                           by = list(bl.patch.hexid.centroids$grid_id[bl.patch.hexid.centroids$in.landscape == T]), function(x){ # over each grid
  #                             length(unique(x))}) %>% rename(grid_id = Group.1, n.clumps = x)
  # hex.tot.patch.ha =  aggregate(bl.patch.hexid.centroids$subpatch.hex.ha[bl.patch.hexid.centroids$in.landscape == T],  # scaled (so mean landscaep cost ==1) contribution to patch
  #                               by = list(bl.patch.hexid.centroids$grid_id[bl.patch.hexid.centroids$in.landscape == T]), function(x){ # over each grid
  #                                 sum(x)}) %>% rename(grid_id = Group.1, tot.patch.ha = x)
  # hex.tot.aw.patch.ha =  aggregate(bl.patch.hexid.centroids$awi.subpatch.hex.ha[bl.patch.hexid.centroids$in.landscape == T],  # scaled (so mean landscaep cost ==1) contribution to patch
  #                                  by = list(bl.patch.hexid.centroids$grid_id[bl.patch.hexid.centroids$in.landscape == T]), function(x){ # over each grid
  #                                    sum(x)}) %>% rename(grid_id = Group.1, tot.aw.patch.ha = x)
  # hex.tot.edge.patch.ha =  aggregate(bl.patch.hexid.centroids$edge.subpatch.hex.ha[bl.patch.hexid.centroids$in.landscape == T],  # scaled (so mean landscaep cost ==1) contribution to patch
  #                                    by = list(bl.patch.hexid.centroids$grid_id[bl.patch.hexid.centroids$in.landscape == T]), function(x){ # over each grid
  #                                      sum(x)}) %>% rename(grid_id = Group.1, tot.edge.patch.ha = x)
  # hex.tot.awedge.patch.ha =  aggregate(bl.patch.hexid.centroids$edge.awi.subpatch.hex.ha[bl.patch.hexid.centroids$in.landscape == T],  # scaled (so mean landscaep cost ==1) contribution to patch
  #                                      by = list(bl.patch.hexid.centroids$grid_id[bl.patch.hexid.centroids$in.landscape == T]), function(x){ # over each grid
  #                                        sum(x)}) %>% rename(grid_id = Group.1, tot.awedge.patch.ha = x)
  # 
  # 
  # 
  # # sensitive.hex.sum.leastcost.contrib.to =  aggregate(bl.patch.hexid.centroids$sensitive.leastcost.contrib.to[bl.patch.hexid.centroids$in.landscape == T], 
  # #                                           by = list(bl.patch.hexid.centroids$grid_id[bl.patch.hexid.centroids$in.landscape == T]), function(x){
  # #                                             sqrt(sum(x))}) %>% 
  # #   rename(grid_id = Group.1, hex.eca = x)
  # # sensitive.hex.sum.scaled.leastcost.contrib.to =  aggregate(bl.patch.hexid.centroids$sensitive.leastcost.contrib.to[bl.patch.hexid.centroids$in.landscape == T], 
  # #                                                  by = list(bl.patch.hexid.centroids$grid_id[bl.patch.hexid.centroids$in.landscape == T]), function(x){
  # #                                                    sqrt(sum(x))}) %>% 
  # #   rename(grid_id = Group.1, hex.eca = x)
  # # sensitive.hex.sum.euclid.contrib.to =  aggregate(bl.patch.hexid.centroids$sensitive.euclid.contrib.to[bl.patch.hexid.centroids$in.landscape == T], 
  # #                                        by = list(bl.patch.hexid.centroids$grid_id[bl.patch.hexid.centroids$in.landscape == T]), function(x){
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
       bl.patch.hexid.centroids,
       file = 
         paste0(func.conect.path, 
                "\\analysis outputs\\All treescapes\\", this.year, "\\r_funcconnect_EffectiveAreas_ECAobs_.RData")
  )
}
print("Total Lanscapes Equivelent Connected Area calculations (script07) done")
