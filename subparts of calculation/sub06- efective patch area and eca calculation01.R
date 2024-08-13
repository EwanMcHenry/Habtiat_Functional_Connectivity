# Ewan McHenry
##------ Fri Mar 04 11:59:09 2022 ------##
# functional connectivity metric dev
# script 06 - equivelent connected area

#EFFECTIVE PATCH AREA ----

bl.patch.hexid.centroids$edge.nonawi.subpatch.hex.ha = bl.patch.hexid.centroids$edge.subpatch.hex.ha - bl.patch.hexid.centroids$edge.awi.subpatch.hex.ha 
bl.patch.hexid.centroids$core.subpatch.hex.ha = bl.patch.hexid.centroids$subpatch.hex.ha - bl.patch.hexid.centroids$edge.subpatch.hex.ha
bl.patch.hexid.centroids$core.awi.subpatch.hex.ha =   bl.patch.hexid.centroids$awi.subpatch.hex.ha - bl.patch.hexid.centroids$edge.awi.subpatch.hex.ha
bl.patch.hexid.centroids$core.nonawi.subpatch.hex.ha = bl.patch.hexid.centroids$core.subpatch.hex.ha - bl.patch.hexid.centroids$core.awi.subpatch.hex.ha

bl.patch.hexid.centroids$effective.ha = 
  bl.patch.hexid.centroids$core.nonawi.subpatch.hex.ha * constants$non.awi.qual.eff + # core non awi value
  bl.patch.hexid.centroids$core.awi.subpatch.hex.ha * constants$awi.qual.eff + # core awi value
  bl.patch.hexid.centroids$edge.nonawi.subpatch.hex.ha * constants$non.awi.qual.eff * constants$relative.edge.quality + # edge non-awi value
  bl.patch.hexid.centroids$edge.awi.subpatch.hex.ha * constants$awi.qual.eff * constants$relative.edge.quality # edge awi value
  

# patch contributions to ECA ----
n.patch = dim(bl.patch.hexid.centroids)[1]
least.cost.contrib.from.patch.hexid = matrix(NA, nrow = n.patch, ncol = n.patch ) 
scaled.least.cost.contrib.from.patch.hexid = matrix(NA, nrow = n.patch, ncol = n.patch ) 
euclid.contrib.from.patch.hexid = matrix(NA, nrow = n.patch, ncol = n.patch ) 

# sensitive.least.cost.contrib.from.patch.hexid = matrix(NA, nrow = n.patch, ncol = n.patch ) 
# sensitive.scaled.least.cost.contrib.from.patch.hexid = matrix(NA, nrow = n.patch, ncol = n.patch ) 
# sensitive.euclid.contrib.from.patch.hexid = matrix(NA, nrow = n.patch, ncol = n.patch ) 

print("ECA Calculations")
for ( i in 1:n.patch){
  least.cost.contrib.from.patch.hexid[i,] =         bl.patch.hexid.centroids$effective.ha [i] * bl.patch.hexid.centroids$effective.ha * exp(-constants$alpha * effective.distance[,i])
  scaled.least.cost.contrib.from.patch.hexid[i,] =  bl.patch.hexid.centroids$effective.ha [i] * bl.patch.hexid.centroids$effective.ha * exp(-constants$alpha * effective.distance[,i]* (1/landscape.mean.scaled.ecolog.cost.not.sea)) # scaled so that average cell has cost of 1
  euclid.contrib.from.patch.hexid[i,] = bl.patch.hexid.centroids$effective.ha [i] * bl.patch.hexid.centroids$effective.ha * exp(-constants$alpha * patch.euc.dists[,i])
  
  # sensitive.least.cost.contrib.from.patch.hexid[i,] = bl.patch.hexid.centroids$effective.ha [i] * bl.patch.hexid.centroids$effective.ha * exp(-constants$alpha.sensitive * effective.distance[,i])
  # sensitive.scaled.least.cost.contrib.from.patch.hexid[i,] = bl.patch.hexid.centroids$effective.ha [i] * bl.patch.hexid.centroids$effective.ha * exp(-constants$alpha.sensitive * effective.distance[,i]* (1/landscape.mean.scaled.ecolog.cost.not.sea))# scaled so that average cell has cost of 1
  # sensitive.euclid.contrib.from.patch.hexid[i,] = bl.patch.hexid.centroids$effective.ha [i] * bl.patch.hexid.centroids$effective.ha * exp(-constants$alpha.sensitive * patch.euc.dists[,i])
  
  svMisc::progress(i)
  
}
# sum of all contributions to/from every patch
bl.patch.hexid.centroids$leastcost.contrib.to = rowSums(least.cost.contrib.from.patch.hexid, na.rm = T)
bl.patch.hexid.centroids$scaled.leastcost.contrib.to = rowSums(scaled.least.cost.contrib.from.patch.hexid, na.rm = T)
bl.patch.hexid.centroids$euclid.contrib.to = rowSums(euclid.contrib.from.patch.hexid, na.rm = T)

# bl.patch.hexid.centroids$sensitive.leastcost.contrib.to = rowSums(sensitive.least.cost.contrib.from.patch.hexid, na.rm = T)
# bl.patch.hexid.centroids$sensitive.scaled.leastcost.contrib.to = rowSums(sensitive.scaled.least.cost.contrib.from.patch.hexid, na.rm = T)
# bl.patch.hexid.centroids$sensitive.euclid.contrib.to = rowSums(sensitive.euclid.contrib.from.patch.hexid, na.rm = T)

# LANDSCAPE SCALE METRICS ----
# square root of sum of all summed contributions of patches in the defined landscape

landscape.metrics = data.frame(name = this.tss[this.ts.num],
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

# ECA per grid cell ----
# sqrt of sum of all summed contributions of patches in the defined landscape by hex
hex.sum.leastcost.contrib.to =  aggregate(bl.patch.hexid.centroids$leastcost.contrib.to[bl.patch.hexid.centroids$in.landscape == T], # contrib to each patch
                                          by = list(bl.patch.hexid.centroids$grid_id[bl.patch.hexid.centroids$in.landscape == T]),  # by hex
                                          function(x){sqrt(sum(x))}) %>% rename(grid_id = Group.1, hex.eca = x)
hex.sum.scaled.leastcost.contrib.to =  aggregate(bl.patch.hexid.centroids$scaled.leastcost.contrib.to[bl.patch.hexid.centroids$in.landscape == T],  # scaled (so mean landscaep cost ==1) contribution to patch
                                          by = list(bl.patch.hexid.centroids$grid_id[bl.patch.hexid.centroids$in.landscape == T]), function(x){ # over each grid
                                            sqrt(sum(x))}) %>% rename(grid_id = Group.1, hex.eca = x)
hex.sum.euclid.contrib.to =  aggregate(bl.patch.hexid.centroids$euclid.contrib.to[bl.patch.hexid.centroids$in.landscape == T], 
                                       by = list(bl.patch.hexid.centroids$grid_id[bl.patch.hexid.centroids$in.landscape == T]), function(x){
                                         sqrt(sum(x))}) %>% rename(grid_id = Group.1, hex.eca = x)
hex.n.clumps =  aggregate(bl.patch.hexid.centroids$clumps[bl.patch.hexid.centroids$in.landscape == T],  # scaled (so mean landscaep cost ==1) contribution to patch
                          by = list(bl.patch.hexid.centroids$grid_id[bl.patch.hexid.centroids$in.landscape == T]), function(x){ # over each grid
                            length(unique(x))}) %>% rename(grid_id = Group.1, n.clumps = x)
hex.tot.patch.ha =  aggregate(bl.patch.hexid.centroids$subpatch.hex.ha[bl.patch.hexid.centroids$in.landscape == T],  # scaled (so mean landscaep cost ==1) contribution to patch
                          by = list(bl.patch.hexid.centroids$grid_id[bl.patch.hexid.centroids$in.landscape == T]), function(x){ # over each grid
                            sum(x)}) %>% rename(grid_id = Group.1, tot.patch.ha = x)
hex.tot.aw.patch.ha =  aggregate(bl.patch.hexid.centroids$awi.subpatch.hex.ha[bl.patch.hexid.centroids$in.landscape == T],  # scaled (so mean landscaep cost ==1) contribution to patch
                              by = list(bl.patch.hexid.centroids$grid_id[bl.patch.hexid.centroids$in.landscape == T]), function(x){ # over each grid
                                sum(x)}) %>% rename(grid_id = Group.1, tot.aw.patch.ha = x)
hex.tot.edge.patch.ha =  aggregate(bl.patch.hexid.centroids$edge.subpatch.hex.ha[bl.patch.hexid.centroids$in.landscape == T],  # scaled (so mean landscaep cost ==1) contribution to patch
                                 by = list(bl.patch.hexid.centroids$grid_id[bl.patch.hexid.centroids$in.landscape == T]), function(x){ # over each grid
                                   sum(x)}) %>% rename(grid_id = Group.1, tot.edge.patch.ha = x)
hex.tot.awedge.patch.ha =  aggregate(bl.patch.hexid.centroids$edge.awi.subpatch.hex.ha[bl.patch.hexid.centroids$in.landscape == T],  # scaled (so mean landscaep cost ==1) contribution to patch
                                   by = list(bl.patch.hexid.centroids$grid_id[bl.patch.hexid.centroids$in.landscape == T]), function(x){ # over each grid
                                     sum(x)}) %>% rename(grid_id = Group.1, tot.awedge.patch.ha = x)



# sensitive.hex.sum.leastcost.contrib.to =  aggregate(bl.patch.hexid.centroids$sensitive.leastcost.contrib.to[bl.patch.hexid.centroids$in.landscape == T], 
#                                           by = list(bl.patch.hexid.centroids$grid_id[bl.patch.hexid.centroids$in.landscape == T]), function(x){
#                                             sqrt(sum(x))}) %>% 
#   rename(grid_id = Group.1, hex.eca = x)
# sensitive.hex.sum.scaled.leastcost.contrib.to =  aggregate(bl.patch.hexid.centroids$sensitive.leastcost.contrib.to[bl.patch.hexid.centroids$in.landscape == T], 
#                                                  by = list(bl.patch.hexid.centroids$grid_id[bl.patch.hexid.centroids$in.landscape == T]), function(x){
#                                                    sqrt(sum(x))}) %>% 
#   rename(grid_id = Group.1, hex.eca = x)
# sensitive.hex.sum.euclid.contrib.to =  aggregate(bl.patch.hexid.centroids$sensitive.euclid.contrib.to[bl.patch.hexid.centroids$in.landscape == T], 
#                                        by = list(bl.patch.hexid.centroids$grid_id[bl.patch.hexid.centroids$in.landscape == T]), function(x){
#                                          sqrt(sum(x))}) %>% 
#   rename(grid_id = Group.1, hex.eca = x)

# add summed contributions to hex shape object
ts.hexgrid$hex.leastcost.eca = left_join(ts.hexgrid, hex.sum.leastcost.contrib.to)$hex.eca
ts.hexgrid$hex.leastcost.eca[is.na(ts.hexgrid$hex.leastcost.eca)] = 0 # give all NAs (hexes with no contributions) 0 eca
ts.hexgrid$hex.scaled.leastcost.eca = left_join(ts.hexgrid, hex.sum.scaled.leastcost.contrib.to)$hex.eca
ts.hexgrid$hex.scaled.leastcost.eca[is.na(ts.hexgrid$hex.scaled.leastcost.eca)] = 0 # give all NAs (hexes with no contributions) 0 eca
ts.hexgrid$hex.euclid.eca = left_join(ts.hexgrid, hex.sum.euclid.contrib.to)$hex.eca
ts.hexgrid$hex.euclid.eca[is.na(ts.hexgrid$hex.euclid.eca)] = 0 # give all NAs (hexes with no contributions) 0 eca

ts.hexgrid$n.clumps = left_join(ts.hexgrid, hex.n.clumps)$n.clumps
ts.hexgrid$n.clumps[is.na(ts.hexgrid$n.clumps)] = 0 # give all NAs (hexes with no contributions) 0 eca
ts.hexgrid$tot.patch.ha = left_join(ts.hexgrid, hex.tot.patch.ha)$tot.patch.ha
ts.hexgrid$tot.patch.ha[is.na(ts.hexgrid$tot.patch.ha)] = 0 # give all NAs (hexes with no contributions) 0 eca
ts.hexgrid$tot.aw.patch.ha = left_join(ts.hexgrid, hex.tot.aw.patch.ha)$tot.aw.patch.ha
ts.hexgrid$tot.aw.patch.ha[is.na(ts.hexgrid$tot.aw.patch.ha)] = 0 # give all NAs (hexes with no contributions) 0 eca
ts.hexgrid$tot.edge.patch.ha = left_join(ts.hexgrid, hex.tot.edge.patch.ha)$tot.edge.patch.ha
ts.hexgrid$tot.edge.patch.ha[is.na(ts.hexgrid$tot.edge.patch.ha)] = 0 # give all NAs (hexes with no contributions) 0 eca
ts.hexgrid$tot.awedge.patch.ha = left_join(ts.hexgrid, hex.tot.awedge.patch.ha)$tot.awedge.patch.ha
ts.hexgrid$tot.awedge.patch.ha[is.na(ts.hexgrid$tot.awedge.patch.ha)] = 0 # give all NAs (hexes with no contributions) 0 eca

# ts.hexgrid$sensitive.hex.leastcost.eca = left_join(ts.hexgrid, sensitive.hex.sum.leastcost.contrib.to)$hex.eca
# ts.hexgrid$sensitive.hex.leastcost.eca[is.na(ts.hexgrid$sensitive.hex.leastcost.eca)] = 0 # give all NAs (hexes with no contributions) 0 eca
# ts.hexgrid$sensitive.hex.scaled.leastcost.eca = left_join(ts.hexgrid, sensitive.hex.sum.scaled.leastcost.contrib.to)$hex.eca
# ts.hexgrid$sensitive.hex.scaled.leastcost.eca[is.na(ts.hexgrid$sensitive.hex.scaled.leastcost.eca)] = 0 # give all NAs (hexes with no contributions) 0 eca
# ts.hexgrid$sensitive.hex.euclid.eca = left_join(ts.hexgrid, sensitive.hex.sum.euclid.contrib.to)$hex.eca
# ts.hexgrid$sensitive.hex.euclid.eca[is.na(ts.hexgrid$sensitive.hex.euclid.eca)] = 0 # give all NAs (hexes with no contributions) 0 eca

# make uniform areas for hexes with different split landmasses
for(i in 1: length(unique(ts.hexgrid$grid_id[duplicated(ts.hexgrid$grid_id)]))){
  ts.hexgrid$hex.ha[ts.hexgrid$grid_id == unique(ts.hexgrid$grid_id[duplicated(ts.hexgrid$grid_id)])[i] ] = ts.hexgrid$hex.ha[ts.hexgrid$grid_id == unique(ts.hexgrid$grid_id[duplicated(ts.hexgrid$grid_id)])[i] ] %>% sum()
  ts.hexgrid$lcm.ncells[ts.hexgrid$grid_id == unique(ts.hexgrid$grid_id[duplicated(ts.hexgrid$grid_id)])[i] ] = ts.hexgrid$lcm.ncells[ts.hexgrid$grid_id == unique(ts.hexgrid$grid_id[duplicated(ts.hexgrid$grid_id)])[i] ] %>% sum()
  ts.hexgrid$bl.ncells[ts.hexgrid$grid_id == unique(ts.hexgrid$grid_id[duplicated(ts.hexgrid$grid_id)])[i] ] = ts.hexgrid$bl.ncells[ts.hexgrid$grid_id == unique(ts.hexgrid$grid_id[duplicated(ts.hexgrid$grid_id)])[i] ] %>% sum()
  ts.hexgrid$landnotcoastal.ncells[ts.hexgrid$grid_id == unique(ts.hexgrid$grid_id[duplicated(ts.hexgrid$grid_id)])[i] ] = ts.hexgrid$landnotcoastal.ncells[ts.hexgrid$grid_id == unique(ts.hexgrid$grid_id[duplicated(ts.hexgrid$grid_id)])[i] ] %>% sum()
  
}

#standardise by terrestrial area of hex - i.e. this is eca if hex was complete  - max number of cells in any hex/ number of non-coastal land cells in this hex
ts.hexgrid$hex.standardised.leastcost.eca = ts.hexgrid$hex.leastcost.eca *(max(ts.hexgrid$lcm.ncells)/ts.hexgrid$landnotcoastal.ncells)
ts.hexgrid$hex.standardised.leastcost.eca[is.na(ts.hexgrid$hex.standardised.leastcost.eca)] = 0

ts.hexgrid$hex.standardised.scaled.leastcost.eca = ts.hexgrid$hex.scaled.leastcost.eca *(max(ts.hexgrid$lcm.ncells)/ts.hexgrid$landnotcoastal.ncells)
ts.hexgrid$hex.standardised.scaled.leastcost.eca[is.na(ts.hexgrid$hex.standardised.scaled.leastcost.eca)] = 0

ts.hexgrid$hex.standardised.euclid.eca = ts.hexgrid$hex.euclid.eca *(max(ts.hexgrid$lcm.ncells)/ts.hexgrid$landnotcoastal.ncells)
ts.hexgrid$hex.standardised.euclid.eca[is.na(ts.hexgrid$hex.standardised.euclid.eca)] = 0




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
             "\\analysis outputs\\", this.tss[this.ts.num], "\\", this.year, "\\r_funcconnect_EffectiveAreas_ECAobs_.RData")
)

print("Equivelent Connected Area calculations (script06) done")
