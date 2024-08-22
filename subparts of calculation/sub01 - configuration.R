##------ Wed Aug 31 14:12:23 2022 ------##
# functional connectivity metric dev
# sub 01.2 - setting constants


# Note - also can config "Data\\hab costs and edge effects Eycott 2011.csv"

# read.csv(paste0(func.conect.path, "\\Data\\delphi point estimate summary.csv"))


## working directories ----
# maindrive = "D:\\Users\\Ewan McHenry\\OneDrive - the Woodland Trust"
# #maindrive = "C:\\Users\\emc2\\OneDrive - The Woodland Trust"
# ts.wd = paste0(maindrive , "\\Treescapes analysis")
func.conect.path = paste0(gis.wd, "\\Connectivity\\Habtiat_Functional_Connectivity") # the project directory for code, outputs etc 

## DEFINE LANDSCAPE(S) ----
# landscape must be a single polygon... obviously, st_union is to make sure.
#Focal_landscape = st_read(paste0(gis.wd, "Data\\Treescape boundaries\\Ewan TS_priority_v1.01gbgrid01.shp")) %>% st_transform( 27700) %>% arrange(name) # sf of landscapes for whcih connectivitty is to be calcualted
Focal_landscape = st_read(paste0(gis.wd, "\\Data\\Landscapes\\Usk Catchments\\Usk Catchments.shp")) %>% st_transform( 27700) %>% st_union() %>% st_as_sf()
Focal_landscape$name = "Usk Catchments"

## Define year ----
years.considered = c(2019, 1990) # vector of years to be calcualted over -- must be LCM data availible and comparible for these years
# years.considered = c( 2019) # vector of years to be calcualted over -- must be LCM data availible and comparible for these years


# SET MODEL CONSTANTS ----
constants <- list(
  # hex grid size  
  hexdist.v = 1000,
  hexdist.h = 1000,
  # - NOTE - hex size is smaller than max dispersal considered. Im happy with this, but if it is bigger it might cause problems later in cost dist calc subscript 05.. might not ... think about it
  
  ## patch identification ----
  focal.hab.num.lcm = 1, # 1 is broadleaf
  # patch clumping buffer -   # often single real life "patch" might be split by 1-2 cells in teh data, this is the distance to join such split-patches into one
  buffer.for.patchid = 25, # distance to buffer around LCM patchest to define patch ID, this is 1/2 the max distance separating clumps within the same patch

  ## patch quality modifiers  - taken from delphi expert opinion
  relative.edge.quality = 0.6, # relative value of edge habitat compared to core(1)
  non.awi.qual.eff = 1, # relative value of non awi
  awi.qual.eff =  1/0.45, # relative value of awi within patches
  
  # dispersal distance parameters - just change the dispersal.dist.set
  dispersal.dist.set = 1000, # distance where dispersal success probably is set
  prob.dispersal.at.set = 0.5, # probability of successful dispersal at dispersal.dist.set
  multiple.max.considered = 2, # factor to limit patch searching to an absolute max dispersal distance for computationall efficiency

  # dispersal cost parameters
  focal.hab.cost.num = 1, # lcm code, 1 is broadleaf
  focalhab.cost = 0.05, # cost of moving through focal (bl woodland) habitat - given nominal amount to protect against imapct of super thin corridors
                # could consider to be 1
  cost.res = 4, # n cells to be mean-aggregated (vert and horiz) for dispersal cost. higher to reduce computing time
  
  #landscape buffer
  # used in sub03 - 
  # define resolution of the landscape and the very rough buffer beyond focal landscape, to consider outside impacts
  landscape.buffer.simplification.tolerance = c(100, 1000) 
  
  )  
  
# define dispersal costs set ----
dispers.costs = read.csv("Data\\hab costs and edge effects Eycott 2011.csv")  %>% 
  mutate(guy.cost = ifelse(guy.cost == 1000, 50, guy.cost)) %>% 
         #replace all 1000 costs with 50, as they are too high and mess with scaling
  mutate(hab = hab %>% as.factor(),
         hab.num = hab.num %>% as.factor(),
         ecolog.cost = guy.cost,
         edge.extent = eycott.edge.extent) %>% 
  # assign custom focal hab type cost
  mutate(ecolog.cost = ifelse (hab.num == constants$focal.hab.cost.num, 
                               constants$focalhab.cost, ecolog.cost))

# DERIVED CONSTANTS -----

sub.code.path = paste0(func.conect.path, "\\subparts of calculation") # subparts to code

this.tss = Focal_landscape$name # vector of names of all landscapes to be calculated over

  # derivative dispersal parameters
constants$alpha  = -log(constants$prob.dispersal.at.set)/constants$dispersal.dist.set # dispersal kernal scaleing parameter
constants$dispers_95 = -log(0.05)/constants$alpha # distance where dispersal success is 0.05
constants$max.dispersal.considered  = constants$dispers_95  * constants$multiple.max.considered # to limit max dispersal distance between patches for computationall efficiency
constants$max.p.considered = exp(-constants$alpha * constants$max.dispersal.considered)
constants$mean.effective.dispersal = -log(0.5)/constants$alpha
  
# derived landscape buffer
constants$buffer.roundLandscape = constants$max.dispersal.considered # buffer around landscape, to account for contribution of neighboring habitat

# derive scaling for dispersal costs to optimise cost distance calculations
# scale cost by fraction of maximum cost - cost clacs work best when costs close to 1, only use multiples becasue easy to scale final distances
constants$cost.scale.factor = max(dispers.costs$ecolog.cost)/5
dispers.costs$scaled.ecolog.cost <-  dispers.costs$ecolog.cost/constants$cost.scale.factor


