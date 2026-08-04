

# quick load and show of stuff from Northern forest


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
# read sub 0.2.1 - curation of global data
source("subparts of calculation\\sub02.1 - curate patch_agnostic.R")



load(paste0(func.conect.path, 
            "\\analysis outputs\\", this.tss[this.ts.num], "\\", years.considered[1], "\\lcmres", constants$lcm.res, 
            "\\06_r_funcconnect_EffectiveAreas_ECAobs_.RData"))


bigstats <-  landscape_stats[[1]]$year_stats$landscape.metrics

landscape_stats[[1]]$year_stats

# $landscape.metrics
# name year  lcd.ECA euclid.ECA n.clumps tot.patch.area.effective.ha tot.patch.ha tot.aw.patch.ha tot.edge.patch.ha
# 1 Northern Forest 2020 5103.541   7982.794    34011                    252591.3     159146.5        30087.81          110386.1
# source_patch_centroid_info med.patch.ha mean.patch.ha
# 1                  153981875       1.3125      4.679266

load(paste0(func.conect.path, 
            "\\analysis outputs\\", this.tss[this.ts.num], "\\", years.considered[2], "\\lcmres", constants$lcm.res, 
            "\\06_r_funcconnect_EffectiveAreas_ECAobs_.RData"))

landscape_stats[[1]]$year_stats
bigstats[2,] <-  landscape_stats[[1]]$year_stats$landscape.metrics


# $landscape.metrics
# name year leastcost.ECA euclid.ECA n.clumps tot.patch.area.effective.ha tot.patch.ha tot.aw.patch.ha tot.edge.patch.ha
# 1 Northern Forest 2024       5005.22   8218.719    37500                    271386.3       172229        30183.06          122584.6
# source_patch_centroid_info med.patch.ha mean.patch.ha
# 1                  158894375       1.3125      4.592773

