# Ewan McHenry
##------ Wed Jun  1 15:22:52 2022 ------##
# leaflet maps
# should run independently if "Presenting stats, change..." code run

# working directories ----
maindrive = "D:\\Users\\Ewan McHenry\\OneDrive - the Woodland Trust"
#maindrive = "C:\\Users\\emc2\\OneDrive - The Woodland Trust"
ts.wd = paste0(maindrive , "\\Treescapes analysis")
gis.wd = paste0( maindrive, "\\GIS")
func.conect.path = paste0(gis.wd, "\\Connectivity\\Functional connectivity\\functional conectivity metric dev")
leaflet.path = paste0(func.conect.path, 
                      "\\analysis outputs\\.maps\\leaflet maps")
# libraries ----

library(tidyverse)
library(sf) # for gis
library(htmltools)
library(leaflet)
library(leaflet.providers)
library(leaflet.extras)
library(reldist)# for weighted quantile function
library(scales) # colour_ramp() returns hexcode for colour interpolation 
library(colourvalues) # for  colour_values(temp), default pallete is viridis
library(htmlwidgets)
library(leaflegend)

source("D:\\Users\\Ewan McHenry\\OneDrive - the Woodland Trust\\GIS\\Ewans functions.R")
source("D:\\Users\\Ewan McHenry\\OneDrive - the Woodland Trust\\GIS\\Ewans gis specifications.R")

# LOAD DATA
load(paste0(func.conect.path, "\\analysis outputs\\", "\\r_plots_stats_change.RData"))

# hot fix for dispersal cost scaling ------
hot.fix = 0 # logical - had full loop been run since 6/6/2022 when hotfix of bug in hex median cost put in
# hotfix cost.scale.factor- 
eycott = read.csv(paste0(gis.wd, "\\Connectivity\\Functional connectivity\\functional conectivity metric dev\\hab costs and edge effects Eycott 2011.csv"))
eycott$guy.cost[eycott$guy.cost == 1000] = 50 # hot fix - saltwater guy cost too high (1000), messes with scaling of costs (need to be scaled to not have range of 1000s to make algorithm run nice), here I make it more reasonable

dispers.costs <- data.frame(hab = eycott$hab,
                            hab.num = eycott$hab.num %>% as.factor(),
                            ecolog.cost = eycott$guy.cost,
                            edge.extent = eycott$eycott.edge.extent
                            )
cost.scale.factor = max(dispers.costs$ecolog.cost)/5

# hot fix 
if (hot.fix==1){
  for( ii in 1: length(map(all.hexgrids, "hexgrid"))){
    all.hexgrids[[ii]]$hexgrid$mean.scaled.ecolog.cost.not.sea = all.hexgrids[[ii]]$hexgrid$mean.scaled.ecolog.cost.not.sea * cost.scale.factor
    all.hexgrids[[ii]]$hexgrid$median.scaled.ecolog.cost.not.sea = all.hexgrids[[ii]]$hexgrid$median.scaled.ecolog.cost.not.sea * cost.scale.factor
  }
}

# delete above when hotfix not needed
# create directory for saving leaflet maps as .html ----
dir.create(leaflet.path)

# DEFINE LANDSCAPE(S) and constants ----
this.tss = ts.andAll.lcm.names # vector of names of all landscapes to be calculated over
this.years = c( 2019, 1990) # vector of years to be calcualted over -- must be LCM data availible and comparible for these years
nice.names = ts.andAll.nice.names

# LEAFLET MAP of connectivity in each year and change ----
rr <- tags$div(
  HTML(paste0('<strong>Native woodland Functional conecivity  </strong><br>Dr. Ewan McHenry, ', format(Sys.time(), '%d %B, %Y')))
)  

treescape.leaflet <- lapply(seq_along(this.tss), FUN = function(i) {
  # sf objects for mapping
  grid.change = change.hexgrids[[which(map(change.hexgrids, "name") == this.tss[i])]]$hexgrid %>% st_simplify(dTolerance = 100) %>% st_transform(4326)
  # list of grids for each year
  grid.year =  vector(mode = "list", length = length(this.years))
  names(grid.year) = paste0("year",this.years )
  for( ii in seq_along(this.years)){
    grid.year[[ii]] = all.hexgrids[map(all.hexgrids, "name") == this.tss[i] & map(all.hexgrids, "year") == this.years[ii]][[1]]$hexgrid %>% st_simplify(dTolerance = 100) %>% st_transform(4326)
    
  }
  
  ## colour/fill pallete formatting ----
  ### snapshot absolute ECA and permiabilty----
  #### yearly connectivity ----
  var.name = "hex.standardised.leastcost.eca" # small hex area scales for
  col.lim.var.name = "hex.leastcost.eca" # this solves issue where small portion-hexes with high cover where skewing colour scale. If they are high they get high colour, but dont mess with the scale
  conective.years.considered = this.years[1:2]
  var.landscapes.for.colscale = map(all.hexgrids, "name") == this.tss[i] & map(all.hexgrids, "year") %in% conective.years.considered # ID all grids in considered years for this landscape
  col.lim.var.same.landscape      = as.data.frame(bind_rows(map(all.hexgrids, "hexgrid")[var.landscapes.for.colscale ]))[,col.lim.var.name] # all of col.lim.variable from considered years in this landscape
  col.lim.weights.same.landscape  = as.data.frame(bind_rows(map(all.hexgrids, "hexgrid")[var.landscapes.for.colscale ]))[, "hex.ha"] # all weights from cells in the considered years in this landscape
  big.enough.to.consider = col.lim.weights.same.landscape> (max(col.lim.weights.same.landscape)*0.9) # only consider those with area > 50% land
  colour.limits = c(0, find.lims(var = col.lim.var.same.landscape, quant.weights = col.lim.weights.same.landscape, 
                                 consider = big.enough.to.consider, quant.prob = 0.98, sd.mult = 3 ))
  # squish variable
  connective.scale.vals = lapply(seq_along(conective.years.considered), FUN = function(ii) {
    temp = grid.year[[ii]] %>% as.data.frame() %>% select(all_of(var.name))
    temp[temp[,1]>colour.limits[2],] = colour.limits[2]
    temp[temp[,1]<colour.limits[1],] = colour.limits[1]
    temp %>% as.matrix() %>% c()
  })
  # make pallete function
  connect.year.pal <- colorNumeric(
    palette = "viridis",
    domain = unlist(connective.scale.vals) )
  
  #### yearly permiability ---- 
  var.name = "mean.scaled.ecolog.cost.not.sea" # small hex area scales for
  col.lim.var.name = "mean.scaled.ecolog.cost.not.sea" # this solves issue where small portion-hexes with high cover where skewing colour scale. If they are high they get high colour, but dont mess with the scale
  perm.years.considered = this.years[1]
  var.landscapes.for.colscale = map(all.hexgrids, "name") == this.tss[i] & map(all.hexgrids, "year") %in% perm.years.considered
  col.lim.var.same.landscape      = as.data.frame(bind_rows(map(all.hexgrids, "hexgrid")[var.landscapes.for.colscale ]))[,col.lim.var.name] # all of col.lim.variable from considered years in this landscape
  col.lim.weights.same.landscape  = as.data.frame(bind_rows(map(all.hexgrids, "hexgrid")[var.landscapes.for.colscale ]))[, "hex.ha"] # all weights from cells in the considered years in this landscape
  big.enough.to.consider = col.lim.weights.same.landscape> (max(col.lim.weights.same.landscape)*0.9) # only consider those with area > 50% land
  colour.limits = c(min(col.lim.var.same.landscape), # note this is not 0, 
                    find.lims(var = col.lim.var.same.landscape, quant.weights = col.lim.weights.same.landscape, 
                                 consider = big.enough.to.consider, quant.prob = 0.98, sd.mult = 3 ))
  
  # squish variable
  perm.scale.vals = lapply(seq_along(perm.years.considered), FUN = function(ii) {
    temp = grid.year[[ii]] %>% as.data.frame() %>% select(var.name)%>%as.matrix() %>% c()
    temp[temp > colour.limits[2]] = colour.limits[2]
    temp[temp < colour.limits[1]] = colour.limits[1]
    temp %>% as.matrix() %>% c()
  })
  # make pallete function
  perm.year.pal <- colorNumeric(
    palette = "inferno",
    domain = unlist(perm.scale.vals), reverse = TRUE )

  ###  change in connectivity ----
  #limits to squish colour ramp to
  col.lim.weights.same.landscape = as.data.frame(grid.change)[,names(grid.change) == "hex.ha"] # hex size to weight quantilelimit 
  abs.col.lim = find.lims(var = abs(grid.change$hex.standardised.leastcost.eca), quant.weights = col.lim.weights.same.landscape, 
                            consider = rep(T, length(var)), quant.prob = 0.98, sd.mult = 3 ) # upper limit finder, from Ewans functions
  colour.limits = c(-abs.col.lim, abs.col.lim) 

  # squish variable
  temp.var = grid.change$hex.standardised.leastcost.eca 
  temp.var[temp.var < colour.limits[1]] = colour.limits[1]
  temp.var[temp.var > colour.limits[2]] = colour.limits[2]
  grid.change$scale.vals = temp.var

  # make palletes
  col.to.ramp = c(E.cols$connectiv.low, "#FFFFFF", E.cols$connectiv.high) # select colours for ramp
  cust.pallet = colour_ramp(col.to.ramp)
  pal.connect.change.continuous = colorNumeric(cust.pallet, # make continuous pallette for plot
                                               domain = c(grid.change$scale.vals, colour.limits))
  ## pretty break pallete for legend
  ### breaks
  mybins <- pretty_breaks(n = 6)(grid.change$scale.vals)
  mybins[mybins %in% range(mybins)] = c(-Inf, Inf)
  ### create pallette& colours
  pal.connect.change.bin <- colorBin( palette= cust.pallet, domain = grid.change$scale.vals, na.color="transparent", bins=mybins)
  pal.connect.change.bin.rev <- colorBin( palette= cust.pallet, domain = grid.change$scale.vals, na.color="transparent", bins=mybins, reverse = T)
  connect.change.bincols = pal.connect.change.bin.rev(seq(min(grid.change$scale.vals),max(grid.change$scale.vals), length = length(mybins)+10)) %>% unique()
  ### labels
  mybin.labs = rep("", length(mybins)-1)
  for (ii in 1: sum(mybins<0)){
    if(ii %% 2){ # if remainder from x/2 -- this makes gaps in lables
      mybin.labs[ii] = paste(mybins[ii+1],"to" ,mybins[ii])
      mybin.labs[ii+sum(mybins<0)] = paste(mybins[ii+sum(mybins<0)],"to" ,mybins[ii+1+sum(mybins<0)])
    } else{
      mybin.labs[ii] = ""
      mybin.labs[ii+sum(mybins<0)] = ""
    }
  }
  ## even break pallete for legend - alternative
  even.breaks = seq(min(grid.change$scale.vals), max(grid.change$scale.vals), length = (length(mybins)+ (length(mybins)%%2 -1)))
  ### create pallette & colours
  pal.connect.change.evenbreak <- colorBin( palette= cust.pallet, domain = grid.change$scale.vals, na.color="transparent", bins=length(even.breaks)+1)
  connect.change.even.cols = pal.connect.change.evenbreak(even.breaks ) 
  connect.change.even.cols[ceiling(length(even.breaks)/2)] = col.to.ramp[2] # fix middle colour to middle of diverging pallete
  # labels
  even.break.lab = rep("", length(even.breaks)-1) 
  even.break.lab[1] = "Increased"
  even.break.lab[length(even.breaks)] = "Decreased"
  even.break.lab[ceiling(length(even.breaks)/2)] = "No change"
  
  ## text function ----
  html.conectivity.text.for.fig = function(this.grid, this.year = NULL, change = F){ # function to make text for popup
    if(change == T){title = "<b>Change in:</b> <br>"
    text.round.function = function(x, digits = 0){ # funciton to make nicely formatted numbers for popup
      gsub(" ", "", paste0(c("", "", "+")[sign(x)+2] ,format(round(x,digits), big.mark=",")), fixed = TRUE)
      }} else{title = paste0("<b>",this.year,":</b> <br>")
    text.round.function = function(x, digits = 0){ # funciton to make nicely formatted numbers for popup
      gsub(" ", "", format(round(x,digits), big.mark=","), fixed = TRUE) # positive sign if positive, rounded to 0 digits and , every thousand
    }} 
    paste0( title,
           "Functional connectivity: ECA(PC) = ", text.round.function(this.grid$hex.leastcost.eca), " ha <br>",
           "N patches = ", text.round.function(this.grid$n.clumps), " <br>",
           "Total habitat = ", text.round.function(this.grid$tot.patch.ha), " ha (",
           text.round.function(this.grid$tot.aw.patch.ha), " ha on ancient woodland ) <br>",
           "Core habitat = ", text.round.function(this.grid$tot.patch.ha-this.grid$tot.edge.patch.ha), " ha (",
           text.round.function(this.grid$tot.aw.patch.ha-this.grid$tot.awedge.patch.ha), " ha on ancient woodland ) <br>",
           "Poor quality edge habitat = ", text.round.function(this.grid$tot.edge.patch.ha), " ha (", 
           text.round.function(this.grid$tot.awedge.patch.ha), " ha on ancient woodland ) <br>",
           "Mean dispersal cost = ", text.round.function((this.grid$mean.scaled.ecolog.cost.not.sea), digits = 2), " <br>"
    )
  }
  
  # map params e.g. zoom default ----
  bounds <- grid.change %>% st_transform(4326) %>% st_bbox() %>% as.character()
  
  # run leaflet ----

  leaflet.map = leaflet() %>%
    # addTiles() %>%
    addProviderTiles(providers$CartoDB.Voyager) %>%
    ####### change over time ----
    addPolygons(data = grid.change , stroke = T, color = "grey" ,
                fillColor =  ~pal.connect.change.continuous(grid.change$scale.vals), weight = 0.5, smoothFactor = 0.5,
                opacity = 0.3, fillOpacity = 0.9, 
                popup = html.conectivity.text.for.fig(grid.change, change = T), 
                group = "Connectivity change 1990 to 2019") %>%
    addLegend(position = "bottomright", 
              colors = connect.change.even.cols %>% rev(), 
              # values = seq(min(grid.change$scale.vals), max(grid.change$scale.vals), length = 2) %>% rev(),#grid.change$scale.vals,
              title = "Functional connectivity<br> change: &#916ECA (ha)",
              labels = even.break.lab ,
              labFormat = labelFormat(transform = function(x) sort(x, decreasing = F)),
              opacity = 1, 
              group = paste0("Connectivity change ", this.years[2], " to ", this.years[1]),
              className = paste0("info legend ", keep.only.letters(paste0("Connectivity change ", this.years[2], " to ", this.years[1])))) %>%
    ##### yearly connectivity -----
      ###### 2019 ----
      addPolygons(data = grid.year[[1]] , stroke = F, color = "grey" ,
                  fillColor =  ~connect.year.pal(connective.scale.vals[[1]]) , 
                  weight = 0.5, smoothFactor = 0.5,
                  opacity = 0.3, fillOpacity = 0.8, 
                  popup = html.conectivity.text.for.fig(grid.year[[1]], this.year = conective.years.considered[1], change = F), 
                  group = paste0("Functional connectivity ", conective.years.considered[1])) %>%
      ###### 1990 ----
    addPolygons(data = grid.year[[2]] , stroke = F, color = "grey" ,
                fillColor =  ~connect.year.pal(connective.scale.vals[[2]]) , 
                weight = 0.5, smoothFactor = 0.5,
                opacity = 0.3, fillOpacity = 0.8, 
                popup = html.conectivity.text.for.fig(grid.year[[2]], this.year = conective.years.considered[2], change = F), 
                group = paste0("Functional connectivity ", conective.years.considered[2])) %>%
      ###### shared connectivity legend ----
    addLegend(position = "bottomright", 
              pal = connect.year.pal,
              values = connective.scale.vals[[1]] ,
              title = paste("Equivelent<br> Connected<br> Area (ha)"),
              opacity = 1, 
              group = paste0("Functional connectivity ", conective.years.considered[1]),
              className = paste0("info legend ", keep.only.letters(paste0("Functional connectivity ", conective.years.considered[1])))) %>%
      
    ##### yearly landscape dispersal cost ----
    addPolygons(data = grid.year[[1]] , stroke = F, color = "grey" ,
              fillColor =  ~perm.year.pal(perm.scale.vals[[1]]) , 
              weight = 0.5, smoothFactor = 0.5,
              opacity = 0.3, fillOpacity = 0.8, 
              popup = html.conectivity.text.for.fig(grid.year[[1]], this.year = perm.years.considered[1], change = F), 
              group = paste0("Dispersal permiability ", perm.years.considered[1])) %>%
    addLegend(position = "bottomright", 
              pal = perm.year.pal,
              values = perm.scale.vals[[1]] ,
              title = paste(perm.years.considered[1] ,"<br>Mean dispersal cost"),
              opacity = 1, 
              group = paste0("Dispersal permiability ", perm.years.considered[1]),
              className = paste0("info legend ", keep.only.letters(paste0("Dispersal permiability ", perm.years.considered[1])))) %>%
    ### leaflet options ----
    addLayersControl(baseGroups = c( paste0("Functional connectivity ", conective.years.considered),
                                     "Connectivity change 1990 to 2019",
                                     paste0("Dispersal permiability ", perm.years.considered), "Basemap"),
    options = layersControlOptions(collapsed = F)) %>% 
    addSearchOSM( options = searchOptions(position = "topleft", autoCollapse = TRUE, minLength = 2)) %>%
    # from https://github.com/rstudio/leaflet/issues/477#issuecomment-776802927 -- i dont really understand it...it allows spaces in group names...
    htmlwidgets::onRender("      
    function(el, x) {
      var updateLegend = function () {
          var selectedGroup = document.querySelectorAll('input:checked')[0].nextSibling.innerText.substr(1).replace(/[^a-zA-Z]+/g, ''); 
          document.querySelectorAll('.legend').forEach( a => a.hidden=true );
          document.querySelectorAll('.legend').forEach( l => { if (l.classList.contains(selectedGroup)) l.hidden=false; } );
      };
         updateLegend();
         this.on('baselayerchange', el => updateLegend());
      }"
    ) %>% 
    addControl(rr, position = "bottomleft") 
    
    
  # save ----
    saveWidget(leaflet.map, file=paste0(leaflet.path, "\\",this.tss[i] ,"leaflet.html"))
  
  list(name = this.tss[i],
       leaflet = leaflet.map) 
})

save(treescape.leaflet, 
     file = paste0(func.conect.path, 
                   "\\analysis outputs\\", "\\r_leaflets.RData")
)

  




