# Ewan McHenry
##------ Fri Jan 21 15:50:05 2022 ------##

# development of functional connectivity metric, first focusign on saving scotland Rainforest Treescape
# original script becoming unwhealdy, so splitting up to make mor manageable
##------ Fri Feb 25 09:19:50 2022 ------##
timestart = Sys.time()

# load libraries and functions - read subscript 01 ----
source("subparts of calculation\\sub01- loading libraries and functions 01.R")
source("D:\\Users\\Ewan McHenry\\OneDrive - the Woodland Trust\\GIS\\Ewans functions.R")
source("D:\\Users\\Ewan McHenry\\OneDrive - the Woodland Trust\\GIS\\Ewans gis specifications.R")

# Configuration ----
## working directories ----
maindrive = "D:\\Users\\Ewan McHenry\\OneDrive - the Woodland Trust"
#maindrive = "C:\\Users\\emc2\\OneDrive - The Woodland Trust"
ts.wd = paste0(maindrive , "\\Treescapes analysis")
gis.wd = paste0( maindrive, "\\GIS")
func.conect.path = paste0(gis.wd, "\\Connectivity\\Functional connectivity\\functional conectivity metric dev")
sub.code.path = paste0(func.conect.path, "\\code\\subparts of calculation")
## DEFINE LANDSCAPE(S) ----
#Focal_landscape = st_read(paste0(gis.wd, "\\Data\\Treescape boundaries\\Ewan TS_priority_v1.01gbgrid01.shp")) %>% st_transform( 27700) %>% arrange(name) # sf of landscapes for whcih connectivitty is to be calcualted
Focal_landscape = st_read(paste0(gis.wd, "\\Data\\Rainforest\\welsh rainforest\\2oceaniczone.shp")) %>% st_transform( 27700) 
Focal_landscape$name = "Welsh Rainforest"
this.tss = Focal_landscape$name # vector of names of all landscapes to be calculated over
## Define year ----
#this.years = c( 2019, 1990) # vector of years to be calcualted over -- must be LCM data availible and comparible for these years
this.years = c( 2019) # vector of years to be calcualted over -- must be LCM data availible and comparible for these years
## SET MODEL CONSTANTS ----
source("subparts of calculation\\sub01.2- setting model constants 01.R")


# set up loop for multiple years and landscapes ----

for(this.ts.for.loop in this.tss[1: length(this.tss)]){
  for(this.year in this.years){
    year = this.year
    this.ts.num = which(ts.lcm.names == this.ts.for.loop) # this bit should maybe automatically chosen in for loop in future
    
# read subscript 02 - load uncurated data ----
source(paste0(sub.code.path, "\\sub02- load data01.R"))
# read subscript 03 - data curation ----
## make folder for this landscape's rdata
## curate costs and edge effects, scaling etc
## cut data to buffered treescape
### hexgrid
### lcm from 90 and 19
### awi

source(paste0(sub.code.path, "\\sub03- data curation01.R"))

## load curated data ----
load(paste0(func.conect.path, "\\analysis outputs\\", ts.lcm.names[this.ts.num], "\\", this.year, "\\r_curated data_.RData"))
# optional - define smaller landscape subset for testing ----
# e.size = 30000
# e = as.vector(extent(tsbuff.lcm19.rast25))
# cent = c(mean(e[c(1,2)]),mean(e[c(3,4)]))
# cent[2] = cent[2] + 20200
# cent[1] = cent[1] + 6000
# 
# a = as(extent(cent[1] - e.size, cent[1] + e.size,  cent[2] - e.size, cent[2] + e.size), 'SpatialPolygons')
# crs(a) = crs(tsbuff.awi)
# 
# a.st = st_as_sf(a) %>% st_transform(27700)
# 
# eg.lcm19.rast25 = crop(tsbuff.lcm19.rast25, a)
# eg.awi = st_intersection(tsbuff.awi, a.st) %>% 
#   st_make_valid() %>%  st_cast("MULTIPOLYGON") %>% 
#   st_cast("POLYGON")
# if(this.country == "Scotland"){
#   eg.nwss.raster = crop(tsbuff.nwss.raster, a)
# 
# } else{ eg.nwss.raster = NA}
# 
# 
# eg.hex = st_intersection(tsbuff.hexgrid, a.st) %>% 
#   st_make_valid() %>%  st_cast("MULTIPOLYGON") %>% 
#   st_cast("POLYGON")
# # eg.awi.rast <- crop(tsbuff.awi.raster, a)

## define landscape data ----

if(year == 2019){
  lcm.landscape = tsbuff.lcm19.rast25 #eg.lcm19.rast25
}
if(year == 1990){
  lcm.landscape = tsbuff.lcm90.rast25 #eg.lcm19.rast25
}


awi.landscape = tsbuff.awi#eg.awi
bl.lcm = lcm.landscape#eg.lcm19.rast25 
nwss.landscape = tsbuff.nwss.raster

# awi.rast.landscape = tsbuff.awi.raster#eg.awi.rast
# hex.landscape = tsbuff.hexgrid#eg.hex

# bl.lcm = eg.lcm19.rast25
# awi.landscape = eg.awi
# lcm.landscape = eg.lcm19.rast25
# nwss.landscape = eg.nwss.raster
# hex.landscape = eg.hex
# # awi.rast.landscape = eg.awi.rast

# read subscript 04 - defining patch attributes ----
## define patches - Broadleaf (or conifer in >50% native NWSS) cells contiguous or within small buffer (buffer.for.patchid x 2) 
### polygonise and intersect with grid
### area and awi area of clumps, patches and patch fragments within each hex
## negative edge effects
### polygonise lcm, buffer by edge effect, dissolve, and find area of edge and awi edge in each hex-patch fragment
## patch centroids

source(paste0(sub.code.path, "\\sub04- patchwork01.R"))

  }}
this.ts.for.loop = this.tss[1]
this.year = this.years[1]

for(this.ts.for.loop in this.tss){#: length(this.tss)]){
  for(this.year in this.years){
    year = this.year
    this.ts.num = which(ts.lcm.names == this.ts.for.loop) # this bit should maybe automatically chosen in for loop in future
    # ## load curated data ----
    load(paste0(func.conect.path, "\\analysis outputs\\", ts.lcm.names[this.ts.num], "\\", this.year, "\\r_curated data_.RData"))
    
## LOAD PATCH DATA ----
load(paste0(func.conect.path, "\\analysis outputs\\", ts.lcm.names[this.ts.num], "\\", this.year, "\\r_funcconnect_patchwork.RData")
     )



# read subscript 05 - dispersal cost between patches ----
## count ncells of within each hex of: broadleaf, land not coastal and all cells
## make cost layer using lcm and dispersal dataframe
## mean and median landscape cost per hex
## euclidian distance between patch centroids
## warning if row names dont equal ID, will mess up indexing
## least cost distance and rescale
source(paste0(sub.code.path, "\\sub05- matric and patch isolation01.R"))


## Load cost distance data ----
load(paste0(func.conect.path, "\\analysis outputs\\", ts.lcm.names[this.ts.num], "\\", this.year, "\\r_funcconnect_MatrixCostDists.RData")
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
}}
# read subscript 07 - total combined landscapes effective area and eca calculation ----
# based heavily on sub06, just loading and runign for combined
source(paste0(sub.code.path, "\\sub07- total across landscapes - efective patch area and eca calculation 01.R"))

# time taken ----

timedone = Sys.time()
timetaken = timedone - timestart




# ## Load eca and effective area data ----
load(paste0(func.conect.path, "\\analysis outputs\\", ts.lcm.names[this.ts.num], "\\", this.year, "\\r_funcconnect_EffectiveAreas_ECAobs_.RData")
)
# # subscript 07 - change in eca, landscape stats and spatial hex plots ----
# source(paste0(sub.code.path, "\\sub07- stats, change and map plots01.R"))
# ## load and store all landscape and hexgrid eca info: 
# ### landscape.metrics.all 
# ### all.hexgrids
# ## create change in landscape and hexgrid objects 
# ### landscape.metrics.all.change
# ### change.hexgrids 
# ## ggplots and plotly individual maps -
# ### save pdfs and pngs of individual plots
# ### hexgird eca eca.hexmap
# ### eca.hexmap.plotly
# ## comparison ggplots and save 
# ### comparison.eca.hexmap
# ## change in ECA ggplot and plotly
# ### change.eca.hexmap
# ### change.eca.hexmap.plotly
# 
# # load plots and stats  ----
# load(paste0(func.conect.path, 
#             "\\analysis outputs\\", "\\r_plots_stats_change.RData"))