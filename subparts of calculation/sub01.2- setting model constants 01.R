##------ Wed Aug 31 14:12:23 2022 ------##
# functional connectivity metric dev
# sub 01.2 - setting constants
read.csv(paste0(func.conect.path, "\\delphi point estimate summary.csv"))


# SET MODEL CONSTANTS ----
constants <- list(
  ## patch identification ----
  # often single real life "patch" might be split by 1-2 cells in teh data, this is the distance to join such split-patches into one
  buffer.for.patchid = 25, # distance to buffer around LCM patchest to define patch ID, this is 1/2 the max distance separating clumps within the same patch

  ## patch quality modifiers  - taken from delphi expert opinion
  relative.edge.quality = 0.6, # relative value of edge habitat compared to core(1)
  non.awi.qual.eff = 1, # relative value of non awi
  awi.qual.eff =  1/0.45, # relative value of awi within patches
  
  # dispersal parameters - just change the dispersal.dist.set
  dispersal.dist.set = c(5000) # distance where dispersal success probably is set
  dispersal.dist.set = 1000 #another set of dispersal distances
  
  prob.dispersal.at.set = 0.5, # probability of successful dispersal between patches at dispersal.dist.set

  multiple.max.considered = 2, # factor to limit patch searching to an absolute max dispersal distance for computationall efficiency
  
  
  # derivative dispersal parameters
  alpha  = -log(0.05)/dispersal.dist.set # dispersal kernal scaleing parameter
  max.dispersal.considered  = dispersal.dist.set  * multiple.max.considered # to limit max dispersal distance between patches for computationall efficiency
  max.p.considered = exp(-alpha * max.dispersal.considered)
  mean.effective.dispersal = -log(0.5)/alpha
  
#  add dispersal dist at 95 to constants
  alpha.sensitive =  -log(0.05)/dispersal.dist.set.sensitive
  max.dispersal.considered.sensitive  = dispersal.dist.set.sensitive  * multiple.max.considered # to limit max dispersal distance between patches for computationall efficiency
  max.p.considered.sensitive = exp(-alpha * max.dispersal.considered.sensitive)
  mean.effective.dispersal.sensitive = -log(0.5)/alpha.sensitive
  
  buffer.roundLandscape = max.dispersal.considered # buffer around landscape, to account for contribution of neighboring habitat
  
  # dispersal cost
  focalhab.cost = 0.05 # cost of moving through focal (bl woodland) habitat - given nominal amount to protect against imapct of super thin corridors
  cost.res = 4 # n cells to be mean-aggregated (vert and horiz) for dispersal cost. higher to reduce computing time

# define which cost set is being used here
eycott = read.csv("hab costs and edge effects Eycott 2011.csv")

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
