# Ewan McHenry
# functional connectivity metric dev
# script 01 - loading libraries, functions and data

# libraries ----
library(tidyverse)
library(sf) # for gis
library(raster)
# library(rgdal) # replaced with stars
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
library(giscoR)
library(osmextract) # for roads

library(U.utilities) # Ewan custom functions git_install("EwanMcHenry/U.utilities")

# mask_lcm_landscape <- function(lcm, landscape, country) {
#   
#   lcm.r <- if (country == "Northern Ireland") lcm$ni else lcm$gb
#   
#   crs_use <- if (country == "Northern Ireland") 29903 else 27700
#   
#   ts.vect <- terra::vect(st_transform(landscape, crs_use))
#   
#   template <- terra::crop(terra::rast(lcm.r[[1]]), ts.vect)
#   
#   fast_mask <- function(r) {
#     terra::mask(terra::crop(r, template), ts.vect)
#   }
#   
#   lapply(lcm.r, fast_mask)
# }


trouble_plot <- function(x,name) {
  if(troubleshooting == T){
    #plot raster
    if(class(x) == "SpatRaster"){
      writeRaster(x, 
                  paste0(func.conect.path, 
                         "\\troubleshooting_saves\\", name, ".tif"), overwrite = TRUE)
    }else if(class(x) == "sf"){
      st_write(x, 
               paste0(func.conect.path, 
                      "\\troubleshooting_saves\\", name, ".gpkg"), delete_dsn = TRUE)
    } else if(class(x) == "data.frame"){
      write.csv(x, 
                paste0(func.conect.path, 
                       "\\troubleshooting_saves\\", name, ".csv"), row.names = FALSE)
    }
  }
}
