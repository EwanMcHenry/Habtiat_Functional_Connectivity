# Ewan McHenry
##------ Fri Mar 04 11:59:09 2022 ------##
# functional connectivity metric dev
# determing change over time, makign maps and tables and storign improtant stats.
# this code should opperate stand alone, so long as the master computation code has already been run
# libraries and ----
source("subparts of calculation\\sub00- loading libraries and functions.R")
source("subparts of calculation\\sub01 - configuration.R") # configureation file
library(ggmap)
library("rnaturalearth")
library("rnaturalearthdata")
library("rnaturalearthhires") # install.packages("rnaturalearthhires", repos = "http://packages.ropensci.org", type = "source")
world <- ne_countries(scale = "large", returnclass = "sf")

# configuration ----

nice.names = Focal_landscape$name

stand.plot.height = 7 # plot height ratio
stand.plot.width = 6 # plot width ratio
hex.line.size = 0.05 # line size for hexgrid
round_to_for_legend_state = 25 # number of significant figures for legend
round_to_for_legend_change = 10
# number of significant figures for legend

n_breaks_in_legend = 4 # number of breaks in legend
n_breaks_in_legend_change = 5

# Configuration ----
## SET MODEL CONSTANTS ----

# LOAD and store all landscape and hexgrid eca info: landscape.metrics.all all.hexgrids ----
#  make objects to store in
load(paste0(func.conect.path, "\\analysis outputs\\", 
            this.tss[1], "\\", years.considered[1], "\\r_funcconnect_EffectiveAreas_ECAobs_.RData"))
landscape.metrics.all = landscape.metrics[0,]
# load up data needed for landscapes in loop

hex.colnames = names(ts.hexgrid)
all.hexgrids = list (NA)
for(i in seq_along(this.tss)){
  for(y in seq_along(years.considered)){
    load(paste0(func.conect.path, "\\analysis outputs\\", 
                this.tss[i], "\\", years.considered[y], "\\r_funcconnect_EffectiveAreas_ECAobs_.RData"))
    
    landscape.metrics.all[(i-1)*length(years.considered)+y,] = landscape.metrics
    
    # hexgrid with stuff 
    all.hexgrids[[(i-1)*length(years.considered)+y]] = list(name = this.tss[i], 
                                                      year = years.considered[y],
                                                      hexgrid = ts.hexgrid[,hex.colnames],
                                                      height.width.ratio = as.numeric((st_bbox(ts.hexgrid)$ymax- st_bbox(ts.hexgrid)$ymin)/(st_bbox(ts.hexgrid)$xmax- st_bbox(ts.hexgrid)$xmin)),
                                                      bl.patch.hexid.centroids = bl.patch.hexid.centroids)
  }}

# CURATE - create change in 
## change in landscape landscape.metrics.all.change ----
landscape.metrics.all.change = slice(landscape.metrics[, -2], 0)
landscape.metrics.all.change[1: length(this.tss),] = NA

for(i in seq_along(this.tss)){
    landscape.metrics.all.change$name[i] = this.tss[i]
    landscape.metrics.all.change[i,2:dim(landscape.metrics.all.change)[2] ] = 
      landscape.metrics.all[landscape.metrics.all$name == this.tss[i] & landscape.metrics.all$year == max(years.considered),3:dim(landscape.metrics)[2]] - 
      landscape.metrics.all[landscape.metrics.all$name == this.tss[i] & landscape.metrics.all$year == min(years.considered),3:dim(landscape.metrics)[2]] 
    }

## change in and hexgrid objects change.hexgrids
change.hexgrids <- lapply(seq_along(this.tss), FUN = function(i) {
  dont.diff = c("grid_id", "hex.ha")
  col.to.diff = c("lcm.ncells", "bl.ncells", "landnotcoastal.ncells",
                  "mean.scaled.ecolog.cost.not.sea", "median.scaled.ecolog.cost.not.sea", "hex.leastcost.eca",
                  "hex.scaled.leastcost.eca", "hex.euclid.eca", "n.clumps", 
                  "tot.patch.ha", "tot.aw.patch.ha", "tot.edge.patch.ha",
                  "tot.awedge.patch.ha", "hex.standardised.leastcost.eca", "hex.standardised.scaled.leastcost.eca",
                  "hex.standardised.euclid.eca") # one of these minus the next years is the difference
  
  
  temp.grid = list(name = map_chr(all.hexgrids, "name")[map_dbl(all.hexgrids, "year") == max(years.considered) & 
                                                          map_chr(all.hexgrids, "name") ==  this.tss[i]], # map names, index this one
                   years = years.considered, 
                   hexgrid = all.hexgrids[[which(map_dbl(all.hexgrids, "year") == max(years.considered) & 
                                                   map_chr(all.hexgrids, "name") ==  this.tss[i])]]$hexgrid[,dont.diff]) # map hexgrids sf, index id, geometry and area of this named ts, most recent ( same for all years)
  
  temp.grid$hexgrid[,(length(dont.diff)+2):(length(dont.diff) + length(col.to.diff)+1) ] = # rest cols == change from 1st to last year of this ts
    # mx year - min, as.df, indexing out ID and geometry
    as.data.frame(all.hexgrids[[which(map_dbl(all.hexgrids, "year") == max(years.considered) & 
                                        map_chr(all.hexgrids, "name") ==  this.tss[i])]]$hexgrid )[,col.to.diff ] -
    as.data.frame(all.hexgrids[[which(map_dbl(all.hexgrids, "year") == min(years.considered) & 
                                        map_chr(all.hexgrids, "name") ==  this.tss[i])]]$hexgrid )[,col.to.diff ]
  temp.grid$height.width.ratio = as.numeric((st_bbox(temp.grid$hexgrid)$ymax- st_bbox(temp.grid$hexgrid)$ymin)/(st_bbox(temp.grid$hexgrid)$xmax- st_bbox(temp.grid$hexgrid)$xmin))

  
  temp.grid
})


## MERGE YEARLY DATA INTO ONE OBJECT FOR EACH TS ----
combo.hexgrid = list()
for(i in seq_along(this.tss)){
  # start with hexgrid id, geometry and hex area, taken from the first year of this ts
    combo.hexgrid[[i]] = all.hexgrids[[(i-1)*length(years.considered)+1]]$hexgrid[,c("grid_id", "hex.grid0", "hex.ha")] 
  for(y in seq_along(years.considered)){
    # this year hexgrid
    temp.hexgrid = all.hexgrids[[(i-1)*length(years.considered)+y]]$hexgrid
    # include year in name of cols
    names(temp.hexgrid)[4:dim(temp.hexgrid)[2]] <-  paste0(names(temp.hexgrid)[4:dim(temp.hexgrid)[2]], "_", all.hexgrids[[(i-1)*length(years.considered)+y]]$year) # rename cols to include year
    # match temp.hexgrid cols 4: length(names()) to combo.hexgrid by gridid and 
    combo.hexgrid[[i]] = combo.hexgrid[[i]] %>%
      left_join(as.data.frame(temp.hexgrid) %>%
                  select(-hex.grid0, -hex.ha),
                by = "grid_id") # join by grid_id, keep geometry and area, then rest cols
    
  }
}


#### ADD CHANGE between years for each important.output.cols ----
for (i in seq_along(combo.hexgrid)) {
  for (j in seq_along(important.output.cols)) {
    # create change col name
    change.col.name <- paste(important.output.cols[j], "change", change_years[1],change_years[2], sep = "_")

    # safely calculate the difference using [[ to get numeric vectors
    combo.hexgrid[[i]][[change.col.name]] <-
      combo.hexgrid[[i]][[paste0(important.output.cols[j], "_", change_years[2])]] -
      combo.hexgrid[[i]][[paste0(important.output.cols[j], "_", change_years[1])]]
  }
}
  
#### ORDER COLUMNS BY INDICATOR AND YEAR, WITH CHANGE ----
  cols.order <- c("grid_id", "hex.grid0", "hex.ha")
  for (ind in important.output.cols) {
    cols.order <- c(
      cols.order,
      paste(ind, sort(years.considered), sep = "_"),
      paste(ind, "change", change_years[1], change_years[2], sep = "_")
    )
  }
  combo.hexgrid[[i]] = combo.hexgrid[[i]][,cols.order] # order cols by indicator, then year
  
  
for(i in seq_along(this.tss)){
  write_sf(combo.hexgrid[[i]], 
           paste0(func.conect.path, "\\analysis outputs\\", 
                  this.tss[i], "\\", 
                  this.tss[i], "_", change_years[2],"-",change_years[1] , "_ECA_output_EMcH.gpkg"))
}

# INDIVIDAUL PLOTS - ggplots list - hexgird eca eca.hexmap  ---- 

  eca.hexmap <- lapply(seq_along(all.hexgrids), FUN = function(i) {
    ### variables ----
    var.name = "hex.standardised.leastcost.eca" # small hex area scales for
    col.lim.var.name = "hex.leastcost.eca" # hex.leastcost.eca solves issue where small portion-hexes with high cover where skewing colour scale. If they are high they get high colour, but dont mess with the scale

    main.title = paste("Bigger, better, more joined up: functional connectivity of native woodland")
    sub.title = paste0(this.tss[this.tss == all.hexgrids[[i]]$name], " ", all.hexgrids[[i]]$year, ". Landscape ECA(PC) = ",landscape.metrics.all$leastcost.ECA[i]  %>%  round( digits = 0) %>% format( big.mark = ","), " ha"  )
    fill.scale.title = "ECA(PC) (ha)"
    
    grid = all.hexgrids[[i]]$hexgrid %>% st_simplify(dTolerance = 100)
    var = as.data.frame(grid)[,names(grid) == var.name]
    
    ##### colour limits & plot ----
    # set variable and weights to make colour scale
    col.lim.var.same.landscape = as.data.frame(bind_rows(map(all.hexgrids, "hexgrid")[# all the hexgrids
      map(all.hexgrids, "name") == all.hexgrids[[i]]$name ]))[,col.lim.var.name] # that have same name as this one (n = n.years), extract this variable and concat into single vector
    col.lim.weights.same.landscape = as.data.frame(bind_rows(map(all.hexgrids, "hexgrid")[# all the hexgrids
      map(all.hexgrids, "name") == all.hexgrids[[i]]$name ]))[, "hex.ha"] # that have same name as this one (n = n.years), extract the hex area and concat into single vector
    col.lim.weights.same.landscape = col.lim.weights.same.landscape[col.lim.var.same.landscape !=0]
    col.lim.var.same.landscape = col.lim.var.same.landscape[col.lim.var.same.landscape !=0]
    
    colour.limits = c(0,min(mean(col.lim.var.same.landscape) + 4* sd(col.lim.var.same.landscape),
                             wtd.quantile (col.lim.var.same.landscape, q = 0.99, 
                                           na.rm = FALSE, weight= col.lim.weights.same.landscape))) 
    dividor = 1
    ### create directory for maps ----
    dir.create(paste0(func.conect.path, "\\analysis outputs\\.maps\\", all.hexgrids[[i]]$name), showWarnings = FALSE)
    
    ### map plot ----

    plot = ggplot(data=grid) +
      geom_sf(data = world %>%  st_transform(27700), size = 0.1) +
      geom_sf(mapping = aes(fill = hex.standardised.leastcost.eca, 
                            text =  map(paste0("ECA(PC) = ", format(round(hex.leastcost.eca,0), big.mark=","), " ha <br>", 
                                               "N patches = ", format(round(n.clumps,0), big.mark=","), " <br>",
                                               "Total habitat = ", format(round(tot.patch.ha,0), big.mark=","), " ha (including ", format(round(tot.aw.patch.ha,0), big.mark=","), " ha ancient woodland ) <br>",
                                               "Core habitat = ", format(round(tot.patch.ha-tot.edge.patch.ha,0), big.mark=","), " ha (including ", format(round(tot.aw.patch.ha-tot.awedge.patch.ha,0), big.mark=","), " ha ancient woodland ) <br>",
                                               "Negative edge habitat = ", format(round(tot.edge.patch.ha,0), big.mark=","), " ha (including ", format(round(tot.awedge.patch.ha,0), big.mark=","), " ha ancient woodland ) <br>"
                            ), HTML)
      ), colour = "black", size = hex.line.size) +
      scale_fill_viridis_c(name = fill.scale.title,
                           limits = range(colour.limits, colour.brks(lims = colour.limits, n = n_breaks_in_legend, round_to = round_to_for_legend_change)),
                           oob = scales::squish, 
                           breaks = colour.brks(lims = colour.limits, n = n_breaks_in_legend, round_to = round_to_for_legend_state),
                           labels = colour.lable(x = col.lim.var.same.landscape ,
                                                 lims = colour.limits , n = n_breaks_in_legend, 
                                                 round_to = round_to_for_legend_state),
                           #option = "magma",direction = -1 
                           guide = guide_colorbar(
                             direction = "horizontal", barheight = unit(2, units = "mm"),
                             barwidth = unit(50, units = "mm"), draw.ulim = F,
                             title.position = 'top', title.hjust = 0.5, label.hjust = 0.5))+
      labs(x = NULL, y = NULL , title = main.title, subtitle = sub.title#, caption = ""
      )+
      coord_sf(xlim = pad.lim(st_bbox(grid) [c(1,3)]), ylim = pad.lim(st_bbox(grid) [c(2,4)]), expand = FALSE) +
      theme_map() +
      theme(legend.position = "bottom"
      ) 
    

    
    ##### save plot ----
    
    # ggsave(plot, filename = paste0(func.conect.path, 
    #                                "\\analysis outputs\\", all.hexgrids[[i]]$name, "\\", all.hexgrids[[i]]$year, 
    #                                "\\hex.ECA.pdf"),
    #        height = 5  , width = stand.plot.width *  max(1,1/(all.hexgrids[[i]]$height.width.ratio)))
    
    ggsave(plot, filename = paste0(func.conect.path, 
                                   "\\analysis outputs\\.maps\\", all.hexgrids[[i]]$name, "\\", 
                                   all.hexgrids[[i]]$name, all.hexgrids[[i]]$year, "_hex.ECA.pdf"),
           height = stand.plot.height  , width = stand.plot.width *  max(1,1/(all.hexgrids[[i]]$height.width.ratio)))
    
    ggsave(plot, filename = paste0(func.conect.path, 
                                   "\\analysis outputs\\.maps\\", all.hexgrids[[i]]$name, "\\", 
                                   all.hexgrids[[i]]$name, all.hexgrids[[i]]$year, "_hex.ECA.png"),
           height = stand.plot.height  , width = stand.plot.width *  max(1,1/(all.hexgrids[[i]]$height.width.ratio)), bg = "white", dpi = 900)
    
    # POINT PLOT ----
    max_pointplot = max(colour.limits, colour.brks(lims = colour.limits, n = n_breaks_in_legend, round_to = round_to_for_legend_change))
    point_data <-  data.frame(var = sort(var), order= 1:length(var))
    target_hex <-  data.frame (var = quantile(point_data$var, probs = quantile_target),
                                order = round(length(var)*quantile_target) )
    point_data$var[point_data$var > max_pointplot] <- max_pointplot # just hold at max
    
    point_plot = ggplot(point_data, aes(x = order, y = var, colour = var)) +
      geom_point(stat = "identity", size = 0.1) +
      scale_y_continuous(name = "ECA(PC) (ha)", 
                         breaks = pretty_breaks(n = n_breaks_in_legend_change, round_to = round_to_for_legend_change)) +
      scale_x_discrete(name = "Order") +
      scale_colour_viridis_c(name = fill.scale.title,
                           limits = range(colour.limits, colour.brks(lims = colour.limits, n = n_breaks_in_legend, round_to = round_to_for_legend_change)),
                           oob = scales::squish,
                           guide = F) +
      labs(title = paste("ECA(PC) for hexes in", all.hexgrids[[i]]$name, all.hexgrids[[i]]$year)) +
      theme_pubr() +
      theme(plot.title = element_text(hjust = 0.5, size = 10),
            axis.text.y = element_text(size = 8),
            axis.text.x = element_text(size = 8))
    
    #pointplot with lines highlighting target quantile ECA
    point_plot_with_target <- point_plot +
      geom_hline(data = target_hex, aes(yintercept = var), linetype = "dashed", color = "red") +
      geom_segment(data = target_hex, aes(x = order, xend = order, y = 0, yend = var), color = "red") +
      geom_text(data = target_hex, aes(x = order/2, y = var, 
                                       label = paste0((1 - quantile_target) * 100, "% threshold: ECA = ", round(var, 2), " ha")), 
                vjust = -0.5, hjust = 1, color = "red", size = 3.5) +
      theme(legend.position = "none")

    #### save point plots ----
    dir.create(paste0(func.conect.path, "\\analysis outputs\\.maps\\", all.hexgrids[[i]]$name, "\\point plots"), showWarnings = FALSE)
    ggsave(point_plot, filename = paste0(func.conect.path, 
                                         "\\analysis outputs\\.maps\\", all.hexgrids[[i]]$name, "\\point plots\\", 
                                         all.hexgrids[[i]]$name, all.hexgrids[[i]]$year, "_hex.ECA.point.pdf"),
           height = stand.plot.height/3  , width = stand.plot.width *  max(1,1/(all.hexgrids[[i]]$height.width.ratio)))
    ggsave(point_plot, filename = paste0(func.conect.path, 
                                         "\\analysis outputs\\.maps\\", all.hexgrids[[i]]$name, "\\point plots\\", 
                                         all.hexgrids[[i]]$name, all.hexgrids[[i]]$year, "_hex.ECA.point.png"),
           height = stand.plot.height/3  , width = stand.plot.width *  max(1,1/(all.hexgrids[[i]]$height.width.ratio)), bg = "white", dpi = 900)
    # save point plot with target
    ggsave(point_plot_with_target, filename = paste0(func.conect.path, 
                                                     "\\analysis outputs\\.maps\\", all.hexgrids[[i]]$name, "\\point plots\\", 
                                                     all.hexgrids[[i]]$name, all.hexgrids[[i]]$year, "_hex.ECA.point.with_target.pdf"),
           height = stand.plot.height/3  , width = stand.plot.width *  max(1,1/(all.hexgrids[[i]]$height.width.ratio)))
    ggsave(point_plot_with_target, filename = paste0(func.conect.path, 
                                                     "\\analysis outputs\\.maps\\", all.hexgrids[[i]]$name, "\\point plots\\", 
                                                     all.hexgrids[[i]]$name, all.hexgrids[[i]]$year, "_hex.ECA.point.with_target.png"),
           height = stand.plot.height/3  , width = stand.plot.width *  max(1,1/(all.hexgrids[[i]]$height.width.ratio)), bg = "white", dpi = 900)
    
    
    ### output ----
    list(map = plot,
         points = list(point_plot,
                       point_plot_with_target),
         target_hex = target_hex)
    
  })

# source("archive\\save hex plots as one pdf.R")

# COMPARISON PLOTS -  ggplot list comparison.eca.hexmap ----
comparison.eca.hexmap <- lapply(seq_along(this.tss), FUN = function(i) {
  joint.title = paste("Bigger, better, more joined up: Functional connectivity of native woodland")#, nice.names[i])
  sub.title01 = paste0("Landscape ECA(PC) = ",landscape.metrics.all$leastcost.ECA[landscape.metrics.all$name == this.tss[i] & landscape.metrics.all$year == min(years.considered)]  %>%  round( digits = 0) %>% format( big.mark = ","), " ha")
  sub.title02 = paste0("Landscape ECA(PC) = ",landscape.metrics.all$leastcost.ECA[landscape.metrics.all$name == this.tss[i] & landscape.metrics.all$year == max(years.considered)]  %>%  round( digits = 0) %>% format( big.mark = ","), " ha")

  comparison.plot = ggarrange( eca.hexmap[map(all.hexgrids, "name") == this.tss[i] & map(all.hexgrids, "year") == min(years.considered)][[1]]$map+
                                 labs(title = min(years.considered), subtitle = sub.title01),
                               eca.hexmap[map(all.hexgrids, "name") == this.tss[i] & map(all.hexgrids, "year") == max(years.considered)][[1]]$map+
                                 labs(title = max(years.considered), subtitle = sub.title02),
                               ncol = 2, common.legend = TRUE, legend = "bottom") %>%
    annotate_figure(top = text_grob(joint.title, face = "bold", size = 11))
  
  # create directory for map plots
  # dir.create(paste0(func.conect.path, "\\analysis outputs\\.maps\\treescape comparison\\", all.hexgrids[[i]]$name), showWarnings = FALSE)

  # ggsave(comparison.plot, filename = paste0(gis.wd, 
  #                                           "\\Connectivity\\Functional connectivity\\functional conectivity metric dev\\analysis outputs\\", this.tss[i], "\\comparison_", min(years.considered), "_", max(years.considered), "_hex.ECA.pdf"),
  #        height = 5  , width = 2 * stand.plot.width *  max(1,1/(change.hexgrids[[i]]$height.width.ratio)))
  
  # hotfix 12.07.22 - to make these more accessible.. hope it works, if no revert to above
  

  ggsave(comparison.plot, filename = paste0(func.conect.path, 
                                            "\\analysis outputs\\.maps\\", all.hexgrids[[i]]$name, "\\comparison - ",
                                            min(years.considered), "_", max(years.considered), "_", this.tss[i], "_hex.ECA.pdf"),
         height = stand.plot.height  , width = stand.plot.width * 2 *  max(1,1/(change.hexgrids[[i]]$height.width.ratio)))

  ggsave(comparison.plot, filename = paste0(func.conect.path, 
                                            "\\analysis outputs\\.maps\\", all.hexgrids[[i]]$name, "\\comparison - ",
                                            min(years.considered), "_", max(years.considered), "_", this.tss[i], "_hex.ECA.png"),
         height = stand.plot.height  , width = stand.plot.width * 2 *  max(1,1/(change.hexgrids[[i]]$height.width.ratio)), dpi = 900, bg = "white")
  
  comparison.plot
})
 
  
# TODO - add point plots for comparison plots
# # Point plots comparision ----
# poits_threshold_comparision <- lapply(seq_along(this.tss), FUN = function(i) {
#   # variables ----
#   var.name = "hex.standardised.leastcost.eca" # small hex area scales for
#   col.lim.var.name = "hex.leastcost.eca" # hex.leastcost.eca solves issue where small portion-hexes with high cover where skewing colour scale. If they are high they get high colour, but dont mess with the scale
#   
#   main.title = paste("Bigger, better, more joined up: functional connectivity of native woodland")
#   sub.title = paste0("Change in functional connectivity in ", nice.names[this.tss == change.hexgrids[[i]]$name], " ", 
#                      min(change.hexgrids[[i]]$year), " to ", max(change.hexgrids[[i]]$year), 
#                      "\nTotal change in ECA(PC) = ",landscape.metrics.all.change$leastcost.ECA[i]  %>%  round( digits = 0) %>% format( big.mark = ","), " ha"  )
#   
#   
#   #####
#   var = as.data.frame(grid)[,names(grid) == var.name]
#   map(all.hexgrids, "hexgrid")$hex.standardised.leastcost.eca
#   
#   
#   col.lim.weights.same.landscape = as.data.frame(bind_rows(
#     map(all.hexgrids, "hexgrid")[var.name]# all the hexgrids
#     map(all.hexgrids, "name") == all.hexgrids[[i]]$name ]))[, "hex.ha"] # that have same name as this one (n = n.years), extract the hex area and concat into single vector
#   
#   
#   ######
#   
#   grid = change.hexgrids[[i]]$hexgrid %>% st_simplify(dTolerance = 100)
#   
#   
#   map(all.hexgrids, "name") == this.tss[[i]] ]))[,col.lim.var.name]
#   
#   
#   # colour limits & plot ----
#   # set variable and weights to make colour scale
#   col.lim.var.same.landscape = as.data.frame(grid)[,names(grid) == col.lim.var.name] %>% 
#     .[. != 0]    # remove all 0
#   col.lim.weights.same.landscape = as.data.frame(grid)[,names(grid) == "hex.ha"][col.lim.var.same.landscape!=0]
#   col.lim.var.same.landscape = col.lim.var.same.landscape[col.lim.var.same.landscape!=0]
#   
#   colour.limits = c(-min(mean(abs(col.lim.var.same.landscape)) + 4* sd(abs(col.lim.var.same.landscape)),
#                          wtd.quantile (abs(col.lim.var.same.landscape), q = 0.99, na.rm = FALSE, weight= col.lim.weights.same.landscape)),
#                     min(mean(abs(col.lim.var.same.landscape)) + 4* sd(abs(col.lim.var.same.landscape)),
#                         wtd.quantile (abs(col.lim.var.same.landscape), q = 0.99, na.rm = FALSE
    
  
       
# CHANGE PLOTS - ggplot list - change in eca  ----

change.eca.hexmap <- lapply(seq_along(this.tss), FUN = function(i) {
  # variables ----

  var.name = "hex.standardised.leastcost.eca" # small hex area scales for
  col.lim.var.name = "hex.leastcost.eca" # this solves issue where small portion-hexes with high cover where skewing colour scale. If they are high they get high colour, but dont mess with the scale
  hist.var.name = "hex.leastcost.eca" # variable to use for histogram
  
  
  main.title = paste("Bigger, better, more joined up: functional connectivity of native woodland")
  sub.title = paste0("Change in functional connectivity in ", nice.names[this.tss == change.hexgrids[[i]]$name], " ", 
                     min(change.hexgrids[[i]]$year), " to ", max(change.hexgrids[[i]]$year), 
                     "\nTotal change in ECA(PC) = ",landscape.metrics.all.change$leastcost.ECA[i]  %>%  round( digits = 0) %>% format( big.mark = ","), " ha"  )
  fill.scale.title = expression(paste(Delta, "ECA (ha)")) 
  
  grid = change.hexgrids[[i]]$hexgrid %>% st_simplify(dTolerance = 100)
  var = as.data.frame(grid)[,names(grid) == var.name]
  hist_var = as.data.frame(grid)[,names(grid) == hist.var.name]
  
  low.col = E.cols$connectiv.low
  high.col = E.cols$connectiv.high
  bigmark.round.pos <- function(x, digits = 0) {
    y <- round(x, digits)
    z <- format(y, big.mark = ",")
    z[!is.na(y) & y > 0] <- paste0("+", z[!is.na(y) & y > 0])
    z
  }
  
  
  # colour limits ----
  # set variable and weights to make colour scale
  col.lim.var.same.landscape = as.data.frame(grid)[,names(grid) == col.lim.var.name] %>% 
    .[. != 0]    # remove all 0
  col.lim.weights.same.landscape = as.data.frame(grid)[,names(grid) == "hex.ha"][col.lim.var.same.landscape!=0]
  col.lim.var.same.landscape = col.lim.var.same.landscape[col.lim.var.same.landscape!=0]
  
  
  abs.col.lim = min(mean(abs(col.lim.var.same.landscape)) + 4* sd(abs(col.lim.var.same.landscape)),
                    wtd.quantile (abs(col.lim.var.same.landscape), q = 0.99, na.rm = FALSE, weight= col.lim.weights.same.landscape))

  colour.limits = c(-abs.col.lim, abs.col.lim) 
  dividor = 1
  
  # plot ----
  
  plot = ggplot(data=grid) + 
          geom_sf(data = world %>%  st_transform(27700), size = 0.1) +
    geom_sf(mapping = aes(fill = hex.standardised.leastcost.eca, 
                          text =  map((paste0("<b>Change in:</b> <br>",
                                              "ECA(PC) = ", bigmark.round.pos(hex.leastcost.eca), " ha <br>",
                                              "N patches = ", bigmark.round.pos(n.clumps), " <br>",
                                              "Total habitat = ", bigmark.round.pos(tot.patch.ha), " ha (",
                                              bigmark.round.pos(tot.aw.patch.ha), " ha ancient woodland ) <br>",
                                              "Core habitat = ", bigmark.round.pos(tot.patch.ha-tot.edge.patch.ha), " ha (",
                                              bigmark.round.pos(tot.aw.patch.ha-tot.awedge.patch.ha), " ha ancient woodland ) <br>",
                                              "Negative edge habitat = ", bigmark.round.pos(tot.edge.patch.ha), " ha (", 
                                              bigmark.round.pos(tot.awedge.patch.ha), " ha ancient woodland ) <br>",
                                              "Mean landscape permiability = ", bigmark.round.pos(mean.scaled.ecolog.cost.not.sea * constants$cost.scale.factor, 2)
                          )), HTML)
    ), colour = "black", size = hex.line.size) +
    scale_fill_gradient2(name = fill.scale.title,
                         low= low.col, high= high.col ,
                         breaks = colour.brks(lims = colour.limits, n = n_breaks_in_legend_change, round_to = round_to_for_legend_change),
                         labels = colour.lable(x = col.lim.var.same.landscape , n = n_breaks_in_legend_change,
                                               lims = colour.limits , round_to = round_to_for_legend_change),
                         limits = range(colour.limits, colour.brks(lims = colour.limits, n = n_breaks_in_legend_change, round_to = round_to_for_legend_change)),
                         oob = scales::squish,
                         guide = guide_colorbar(
                           direction = "horizontal", barheight = unit(2, units = "mm"),
                           barwidth = unit(50, units = "mm"), draw.ulim = F,
                           title.position = 'top', title.hjust = 0.5, label.hjust = 0.5))+
    labs(x = NULL, y = NULL , title = main.title, subtitle = sub.title#, caption = ""
    )+
    coord_sf(xlim = pad.lim(st_bbox(grid) [c(1,3)]), ylim = pad.lim(st_bbox(grid) [c(2,4)]), expand = FALSE) +
    theme_map() +
    theme(legend.position = "bottom"
    ) 
  
  # save plot ----
  
  # ggsave(plot, filename = paste0(func.conect.path, 
  #                                "\\analysis outputs\\", change.hexgrids[[i]]$name, "\\", paste(change.hexgrids[[i]]$year, collapse = "_"), "_hex.ECAchange.pdf"),
  #        height = 5  , width = stand.plot.width *  max(1,1/(change.hexgrids[[i]]$height.width.ratio)))

  # hotfix 12.07.22 - to make these more accessible.. hope it works, if no revert to above
  ggsave(plot, filename = paste0(func.conect.path,
                                 "\\analysis outputs\\.maps\\", all.hexgrids[[i]]$name, "\\change - ", change.hexgrids[[i]]$name,"_", paste(change.hexgrids[[i]]$year, collapse = "_"), "_hex.ECAchange.pdf"),
         height = stand.plot.height  , width = stand.plot.width *  max(1,1/(change.hexgrids[[i]]$height.width.ratio)))
  
  ggsave(plot, filename = paste0(func.conect.path, 
                                 "\\analysis outputs\\.maps\\", all.hexgrids[[i]]$name, "\\change - ", change.hexgrids[[i]]$name,"_", paste(change.hexgrids[[i]]$year, collapse = "_"), "_hex.ECAchange.png"),
         height = stand.plot.height  , width = stand.plot.width *  max(1,1/(change.hexgrids[[i]]$height.width.ratio)), bg = "white", dpi = 900)
  
  
  
  # change in eca histogram ----
  # Precompute histogram
  hist_data <- hist(hist_var, breaks = 50, plot = FALSE)
  # Create data frame
  hist_df <- data.frame(
    count = hist_data$counts,
    bin = hist_data$mids
  )
  
  # Plot using geom_col and gradient fill
  change_hex_hist <- ggplot(hist_df, aes(x = bin, y = count, fill = bin)) +
    geom_col(color = "black", width = diff(hist_data$breaks)[1]) +
    scale_y_continuous(name = "Count", 
                       limits = c(0, max(hist_df$count[abs(hist_df$bin)>min(abs(hist_df$bin))])) * 1.1) +
    scale_fill_gradient2(
      low = low.col,
      high = high.col,
      breaks = colour.brks(lims = colour.limits, n = n_breaks_in_legend_change, round_to = round_to_for_legend_change),
      labels = colour.lable(
        x = col.lim.var.same.landscape,
        n = n_breaks_in_legend_change,
        lims = colour.limits,
        round_to = round_to_for_legend_change
      ),
      limits = range(
        colour.limits,
        colour.brks(lims = colour.limits, n = n_breaks_in_legend_change, round_to = round_to_for_legend_change)
      ),
      oob = scales::squish,
      guide = "none"
    ) +
    theme_pubr()
  
  # save histogram ----
  
  
      # plot ----
  list(map = plot,
       change_hex_hist = change_hex_hist)
  })
  

## plotly change maps ----  
# change.eca.hexmap.plotly <- lapply(seq_along(this.tss), FUN = function(i) {
#   plotly_build(ggplotly(change.eca.hexmap[[i]], tooltip = "text", 
#                         dynamicTicks = T) %>%
#                  config(displayModeBar = FALSE) %>% layout(hoverlabel = list(align = "left")))
  # 
  # all.hexgrids[[i]] = append(all.hexgrids[[i]], print(eca.hexmap) ) %>%
  #   append(., print(eca.hexmap.plotly) )
# })

# SAVE objects ----
save(landscape.metrics.all,
     all.hexgrids,
     landscape.metrics.all.change,
     change.hexgrids,
     eca.hexmap,
     comparison.eca.hexmap,
     change.eca.hexmap,
     # eca.hexmap.plotly,
     # change.eca.hexmap.plotly, 
     file = paste0(func.conect.path, 
              "\\analysis outputs\\", all.hexgrids[[i]]$name, "\\r_plots_stats_change.RData")
)

