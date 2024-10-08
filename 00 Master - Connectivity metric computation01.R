# Ewan McHenry
##------ Fri Jan 21 15:50:05 2022 ------##

# development of functional connectivity metric
#  note that this wont work across GB & NI - need to be done seperate if bioth within same landscape
# 
##------ Fri Feb 25 09:19:50 2022 ------##
timestart = Sys.time()

# load libraries and functions - read subscript 01 ----
source("subparts of calculation\\sub00- loading libraries and functions.R")
# source("D:\\Users\\Ewan McHenry\\OneDrive - the Woodland Trust\\GIS\\Ewans gis specifications.R")

# Configuration ----
## SET MODEL CONSTANTS ----
source("subparts of calculation\\sub01 - configuration.R") # configureation file

# read subscript 02 - load uncurated global data ----
# contains data used to define patches, and the landscape
source("subparts of calculation\\sub02- load data01.R")

# set up loop for multiple years and landscapes ----
this.ts.num = 1
this.year = years.considered[1]
this.ts.for.loop = this.tss[this.ts.num]
for(this.ts.num in 1: length(this.tss)){
  ## make folder for this focal landscape ----
  dir.create(paste0(func.conect.path, 
                    "\\analysis outputs\\", this.tss[this.ts.num]))
  # load global data, uncourated
  load(paste0(func.conect.path, 
              "\\analysis outputs\\r_global_data_.RData"))  
  
  # read subscript 03 - data curation ----
  ## make folder for this landscape's rdata
  ## curate costs and edge effects, scaling etc
  ## cut data to buffered treescape
  ### hexgrid
  ### lcm from 90 and 19
  ### awi
  source(paste0(sub.code.path, "\\sub03- data curation01.R"))
  
  for(this.year in years.considered){
    
    # folder for this year for this landscape
    dir.create(paste0(func.conect.path, 
                      "\\analysis outputs\\", this.tss[this.ts.num],"\\", this.year))
    


## load curated data ----
load(paste0(func.conect.path, "\\analysis outputs\\", this.tss[this.ts.num], "\\r_curated data_.RData"))

# read subscript 04 - defining patch attributes ----
## define patches - Broadleaf (or conifer in >50% native NWSS) cells contiguous or within small buffer (buffer.for.patchid x 2) 
### polygonise and intersect with grid
### area and awi area of clumps, patches and patch fragments within each hex
## negative edge effects
### polygonise lcm, buffer by edge effect, dissolve, and find area of edge and awi edge in each hex-patch fragment
## patch centroids

source(paste0(sub.code.path, "\\sub04- patchwork01.R"))

  }}

# savepoint/checkpoint - next loop takes a bit, so break for a checkpoint.

for(this.ts.for.loop in this.tss){#: length(this.tss)]){
  for(this.year in years.considered){
    this.ts.num = which(this.tss == this.ts.for.loop) # this bit should maybe automatically chosen in for loop in future
    # ## load curated data ----
    load(paste0(func.conect.path, "\\analysis outputs\\", this.tss[this.ts.num], "\\r_curated data_.RData"))
    
## LOAD PATCH DATA ----
load(paste0(func.conect.path, "\\analysis outputs\\", this.tss[this.ts.num], "\\", this.year, "\\r_funcconnect_patchwork.RData"))

# read subscript 05 - dispersal cost between patches ----
## count ncells of within each hex of: broadleaf, land not coastal and all cells
## make cost layer using lcm and dispersal dataframe
## mean and median landscape cost per hex
## euclidian distance between patch centroids
## warning if row names dont equal ID, will mess up indexing
## least cost distance and rescale
source(paste0(sub.code.path, "\\sub05- matric and patch isolation01.R"))

## Load cost distance data ----
load(paste0(func.conect.path, "\\analysis outputs\\", this.tss[this.ts.num], "\\", this.year, "\\r_funcconnect_MatrixCostDists.RData")
            )

# read subscript 06 - effective area and eca calculation ----
## effective patch area, usign edge and awi area
## ECA calculation per patch, per hex and per landscape
### hex == contribution of all patches in hex to all patches in buffered landscape
#### standardised to patch area of non-coastal land
### patch = contribution to all in buffered landscape
### landscape = to all in buffered
### euclidian, least cost and least cost scaled (so that mean cost == 1)
source(paste0(sub.code.path, "\\sub06- efective patch area and eca calculation01.R"))

# end landscape-years loop
print(paste(this.ts.for.loop, this.year, "done"))
  } # end year loop
  } # end landscape loop
# time taken ----

timedone = Sys.time()
timetaken = timedone - timestart


###############################################################################
###############################################################################


