##------ Wed Aug 31 14:12:23 2022 ------##
# functional connectivity metric dev
# sub 01.2 - setting constants




# SET MODEL CONSTANTS ----
# patch identification
# often single real life "patch" might be split by 1-2 cells in teh data, this is the distance to join such split-patches into one
buffer.for.patchid = 25 # distance to buffer around LCM patchest to define patch ID, this is 1/2 the max distance separating clumps within the same patch

# patch quality modifiers  - taken from delphi expert opinion
read.csv(paste0(func.conect.path, "\\delphi point estimate summary.csv"))
relative.edge.quality = 0.6 # relative value of edge habitat compared to core(1)
non.awi.qual.eff = 1 # relative value of non awi
awi.qual.eff =  1/0.45 # relative value of awi within patches

# dispersal parameters - just change the dispersal.dist.95
dispersal.dist.95  = 5000 # "max" dispersal distance (95th percentile)
alpha  = -log(0.05)/dispersal.dist.95 # dispersal kernal scaleing parameter
max.dispersal.considered  = dispersal.dist.95  * 2 # to limit max dispersal distance between patches for computationall efficiency
max.p.considered = exp(-alpha * max.dispersal.considered)
mean.effective.dispersal = -log(0.5)/alpha

#another set of dispersal distances 
dispersal.dist.95.sensitive = 1000
alpha.sensitive =  -log(0.05)/dispersal.dist.95.sensitive
max.dispersal.considered.sensitive  = dispersal.dist.95.sensitive  * 2 # to limit max dispersal distance between patches for computationall efficiency
max.p.considered.sensitive = exp(-alpha * max.dispersal.considered.sensitive)
mean.effective.dispersal.sensitive = -log(0.5)/alpha.sensitive

buffer.roundLandscape = max.dispersal.considered # buffer around landscape, to account for contribution of neighboring habitat

# dispersal cost
focalhab.cost = 0.05 # cost of moving through focal (bl woodland) habitat - given nominal amount to protect against imapct of super thin corridors
cost.res = 4 # n cells to be mean-aggregated (vert and horiz) for dispersal cost. higher to reduce computing time

# define which cost set is being used here
eycott = read.csv(paste0(gis.wd, "\\Connectivity\\Functional connectivity\\functional conectivity metric dev\\hab costs and edge effects Eycott 2011.csv"))

eycott$guy.cost[eycott$guy.cost == 1000] = 50 # hot fix - saltwater guy cost too high (1000), messes with scaling of costs (need to be scaled to not have range of 1000s to make algorithm run nice), here I make it more reasonable
dispers.costs <- data.frame(hab = eycott$hab,
                            hab.num = eycott$hab.num %>% as.factor(),
                            ecolog.cost = eycott$guy.cost,
                            edge.extent = eycott$eycott.edge.extent
)

# hex grid size  - NOTE - hex size is smaller than max dispersal considered. Im happy with this, but if it is bigger it might cause problems later in cost dist calc subscript 05.. might not ... think about it
hexdist.v = 5000
hexdist.h = 5000

landscape.buffer.simplification.tolerance = c(100, 1000) # used in sub03 - curation to cut data by a very rough buffer beyond focal landscape, to consider outside impacts
