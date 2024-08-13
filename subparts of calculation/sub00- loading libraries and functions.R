# Ewan McHenry
# functional connectivity metric dev
# script 01 - loading libraries, functions and data



# libraries ----
library(tidyverse)
library(sf) # for gis
library(raster)
library(rgdal)
library(units) # for set_units()
library(stars)
library(sp)
library(rgeos)
library(rgis) # remotes::install_github("Pakillo/rgis") # for fast_mask # -- other of same name... devtools::install_github('jgcri/rgis')
# remotes::install_github("hunzikp/velox") # dependancy needs to be installed before rgis # may need to run R with admin privaleges to install

library(fasterize)
library(gdistance) # for costdistance - https://rdrr.io/rforge/gdistance/f/inst/doc/gdistance1.pdf 
# remotes::install_github("arielfri/buffr")# for fast raster buffering
library(buffr) # remotes::install_github("arielfri/buffr")# for fast raster buffering
library(compiler)
library(foreach)
library(doParallel)
library(svMisc)
library(exactextractr)
library(ggpubr)
library(mapview) # npts()
library(htmltools)
library(plotly)
library(scales)#prety_breaks()
library(leaflet)
library(leaflet.providers)
library(gridExtra) # for saving ggplots etc
library(reldist)# for weighted quantile function
library(ggnewscale)
library(ggpattern)# remotes::install_github("coolbutuseless/ggpattern")
library(ragg) # scaling ggplots - https://www.tidyverse.org/blog/2020/08/taking-control-of-plot-scaling/
library(ggspatial)
library(ggstar)
library(cowplot)
library(RColorBrewer)
library(GGally)
library(network)
library(sna)
library(ggnetwork)
library(ggExtra)
