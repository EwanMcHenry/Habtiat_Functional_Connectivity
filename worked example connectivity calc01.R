## Ewan McHenry
##------ Wed Aug 31 14:12:23 2022 ------##

# to make a toy walk-through of functioanl connectivity calcs and produce illustrative figures 
# copied a lot of code from master

landscape.name = "Illustrative.disp.dist"

# working directories ----
maindrive = "D:\\Users\\Ewan McHenry\\OneDrive - the Woodland Trust"
#maindrive = "C:\\Users\\emc2\\OneDrive - The Woodland Trust"
ts.wd = paste0(maindrive , "\\Treescapes analysis")
gis.wd = paste0( maindrive, "\\GIS")
func.conect.path = paste0(gis.wd, "\\Connectivity\\Functional connectivity\\functional conectivity metric dev")
sub.code.path = paste0(func.conect.path, "\\code\\subparts of calculation")

    # read subscript 01 - load libraries and functions ----
    source(paste0(sub.code.path, "\\sub01- loading libraries and functions 01.R"))
    source("D:\\Users\\Ewan McHenry\\OneDrive - the Woodland Trust\\GIS\\Ewans functions.R")
    source("D:\\Users\\Ewan McHenry\\OneDrive - the Woodland Trust\\GIS\\Ewans gis specifications.R")
    
# DEFINE LANDSCAPE(S) ----
# sf of landscapes for whcih connectivitty is to be calcualted
Tscapes01 = st_read(paste0(gis.wd, "\\Data\\Treescape boundaries\\Ewan TS_priority_v1.01gbgrid01.shp")) %>% st_transform( 27700) %>% 
  arrange(name) 
Tscapes01 = Tscapes01[0,]# cut out all
countries = st_read(paste0(gis.wd,"\\Data\\administrative boundaries\\Countries\\R.5countries.simp100m.shp"))%>% st_transform( 27700) %>% arrange(name)


landscape.h = 5000
landscape.w = 5000

easting = 567705 
northing = 169258
# easting[2] = easting[1] + landscape.w
# northing[2] = northing[1] + landscape.h

hex.grid0 = st_make_grid(countries, c(landscape.h, landscape.w), what = "polygons", square = F)
hex.grid = st_sf(hex.grid0) %>%
  mutate(grid_id = 1:length(lengths(hex.grid0)))  # add grid ID

Tscapes01[1,"geometry"] = hex.grid %>% 
  filter(st_intersects(hex.grid, 
                       data.frame(easting, northing) %>%
                         st_as_sf(coords = c("easting", "northing"), crs = 27700),
                       sparse = F
  ))

# data.frame(easting, northing) %>%
# st_as_sf(coords = c("easting", "northing"), crs = 27700) %>% 
# st_bbox() %>% 
# st_as_sfc()

Tscapes01$fa_id = 1
# Tscapes01$name = ts.lcm.names = "Illustrative.restricted"
Tscapes01$name = ts.lcm.names = landscape.name


this.tss = Tscapes01$name[1] # vector of names of all landscapes to be calculated over
this.years = c( 2019, 1990)[1] # vector of years to be calcualted over -- must be LCM data availible and comparible for these years

# SET MODEL CONSTANTS  including buffer ----
source(paste0(sub.code.path, "\\sub01.2- setting model constants 01.R"))

buffer.roundLandscape = dispersal.dist.95 * 2 # buffer around landscape, to account for contribution of neighboring habitat -- illustrative set to 0
buffer.for.figure.landscape = landscape.h*0.1 #

landscape.buffer.simplification.tolerance = c(0,0) # simlification of buffer

hexdist.v = 1000000000#
hexdist.h = 1000000000#


this.ts.for.loop = this.tss
this.year = this.years
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
    
    # ## load curated data
    # load(paste0(func.conect.path, "\\analysis outputs\\", ts.lcm.names[this.ts.num], "\\", this.year, "\\r_curated data_.RData"))
    
## define landscape data ----

if(year == 2019){
  lcm.landscape = tsbuff.lcm19.rast25 #eg.lcm19.rast25
}
if(year == 1990){
  lcm.landscape = tsbuff.lcm90.rast25 #eg.lcm19.rast25
}

awi.landscape = tsbuff.awi#eg.awi
bl.lcm = lcm.landscape#eg.lcm19.rast25 
if(this.country == "Scotland"){
  nwss.landscape = tsbuff.nwss.raster    }else{ nwss.landscape = NA}

# awi.rast.landscape = tsbuff.awi.raster#eg.awi.rast
# hex.landscape = tsbuff.hexgrid#eg.hex

# bl.lcm = eg.lcm19.rast25
# awi.landscape = eg.awi
# lcm.landscape = eg.lcm19.rast25
# nwss.landscape = eg.nwss.raster
# hex.landscape = eg.hex
# # awi.rast.landscape = eg.awi.rast

    ## read subscript 04 - defining patch attributes ----
    ## define patches - Broadleaf (or conifer in >50% native NWSS) cells contiguous or within small buffer (buffer.for.patchid x 2) 
    ### polygonise and intersect with grid
    ### area and awi area of clumps, patches and patch fragments within each hex
    ## negative edge effects
    ### polygonise lcm, buffer by edge effect, dissolve, and find area of edge and awi edge in each hex-patch fragment
    ## patch centroids
    
    source(paste0(sub.code.path, "\\sub04- patchwork01.R"))
    
    ## read subscript 05 - dispersal cost between patches ----
    ## count ncells of within each hex of: broadleaf, land not coastal and all cells
    ## make cost layer using lcm and dispersal dataframe
    ## mean and median landscape cost per hex
    ## euclidian distance between patch centroids
    ## warning if row names dont equal ID, will mess up indexing
    ## least cost distance and rescale
    source(paste0(sub.code.path, "\\sub05- matric and patch isolation01.R"))
    
# ## Load cost distance data ----
load(paste0(func.conect.path, "\\analysis outputs\\", ts.lcm.names[this.ts.num], "\\", this.year, "\\r_funcconnect_MatrixCostDists.RData")
)

    ## read subscript 06 - effective area and eca calculation ----
    ## effective patch area, usign edge and awi area
    ## ECA calculation per patch, per hex and per landscape
    ### hex == contribution of all patches in hex to all patches in buffered landscape
    #### standardised to patch area of non-coastal land
    ### patch = contribution to all in buffered landscape
    ### landscape = to all in buffered
    ### euclidian, least cost and least cost scaled (so that mean cost == 1)
    source(paste0(sub.code.path, "\\sub06- efective patch area and eca calculation01.R"))
    
    ## load edge and awi polygons ----
    load(paste0(gis.wd, 
                "\\Connectivity\\Functional connectivity\\functional conectivity metric dev\\analysis outputs\\", 
                ts.lcm.names[this.ts.num], "\\", this.year, "\\edge_awi.polys.RData")
    )
    # awi.bl.patch.hexid
    # patch.edge
    # awi.edge
    # awi.landscape,
    # edge,
    ## load least cost path objects ----
    load(paste0(gis.wd, 
                "\\Connectivity\\Functional connectivity\\functional conectivity metric dev\\analysis outputs\\", 
                ts.lcm.names[this.ts.num], "\\", this.year, "\\cost.objects.RData"))
    
    
# curation for figures -----
    ## figure buffer ----
    # becasue showing entire buffered landscape takes up too much space
    ts.fig.buff = Tscapes01 %>% st_buffer(dist = buffer.for.figure.landscape)
    # and cut down all important spatial data to that for later use
    
    lcm.landscape = lcm.landscape  %>% crop(ts.fig.buff)%>% mask(ts.fig.buff)
    hab.cost.lcm = hab.cost.lcm %>% crop(ts.fig.buff)%>% mask(ts.fig.buff)
    patch.edge = patch.edge %>% st_intersection(ts.fig.buff)
    awi.edge = awi.edge %>% st_intersection(ts.fig.buff)
    bl.patch.hexid.centroids = bl.patch.hexid.centroids %>% st_intersection(ts.fig.buff)
    bl.patch.id.poly.hexid = bl.patch.id.poly.hexid[bl.patch.id.poly.hexid$hex_splitpatch_ID %in% bl.patch.hexid.centroids$hex_splitpatch_ID,]
    bl.patch.hexid.centroids.sp = bl.patch.hexid.centroids.sp[bl.patch.hexid.centroids.sp$hex_splitpatch_ID %in% bl.patch.hexid.centroids$hex_splitpatch_ID,]
    awi.bl.patch.hexid = awi.bl.patch.hexid %>% st_intersection(ts.fig.buff)
    effective.distance = effective.distance[bl.patch.hexid.centroids$hex_splitpatch_ID, bl.patch.hexid.centroids$hex_splitpatch_ID]
    
    ## lcm ----
    lcm.to.df = function(x){
      x %>% 
        as.data.frame( xy = TRUE) %>%
        mutate(cost = dispers.costs$scaled.ecolog.cost[layer],
               layer = factor(layer),
               'Land cover' = factor(ceh.full.habtype)[values(lcm.landscape)],
               ceh.colour = ceh.col.pallette[values(lcm.landscape)] ) %>%
        subset(!is.na(layer))
    }
    
    lcm.df = lcm.landscape %>%
      lcm.to.df()
    lcm.df.focal = lcm.landscape %>% mask(Tscapes01) %>% 
      lcm.to.df()
    
    ## separate quality polygons ----
    nonawi.edge = st_difference(patch.edge, st_union(awi.bl.patch.hexid)) %>% st_simplify(dTolerance = 10)
    nonawi.core = st_difference(bl.patch.id.poly.hexid, st_union(awi.bl.patch.hexid)) %>% st_difference(., st_union(patch.edge)) %>% st_simplify(dTolerance = 10)
    awi.core =   st_difference(awi.bl.patch.hexid, st_union(awi.edge))    
    nonawi.edge$qual.score = relative.edge.quality
    nonawi.core$qual.score = non.awi.qual.eff
    awi.core$qual.score = awi.qual.eff
    awi.edge$qual.score = relative.edge.quality*awi.qual.eff
    col.of.interest = c("geometry","qual.score")
    
    patch.qual.sf = rbind(nonawi.edge[, col.of.interest], 
                          nonawi.core[, col.of.interest], 
                          awi.core[, col.of.interest], 
                          awi.edge[, col.of.interest]
    ) %>% group_by(qual.score) %>%  summarize()
    
    ## mean patch.quality ----
    bl.patch.id.poly.hexid$effective.ha = bl.patch.hexid.centroids$effective.ha
    bl.patch.id.poly.hexid$mean.quality = bl.patch.hexid.centroids$effective.ha/bl.patch.hexid.centroids$clump.ha
    
    ## wider landscape raster of fading value ----
    # to give plots a fade out
    
    extra.fade.distance = 100 #m to block out harsh outline of some polys that falls outside origianl raster
    alpha.inter = 0.5 # alpha intercept - the initial fade level once out of focal landscape
    
    dist.to.center.df =   rasterize(Tscapes01, rasterFromXYZ(lcm.df)) %>% buffer(width = extra.fade.distance) %>% 
      distance() %>% mask(ts.fig.buff %>% st_buffer(dist = extra.fade.distance)) %>% 
      as.data.frame(., xy = TRUE)  %>%
      mutate(fade = layer/max(layer,na.rm = T)) %>%
      mutate(fade = (fade+alpha.inter)/max(fade+alpha.inter, na.rm = T)) %>%
      subset(!is.na(fade))
    
    dist.to.center.df$fade[dist.to.center.df$layer == 0] = 0
    
    bl.patch.hexid.centroids$fade = extract(rasterFromXYZ(dist.to.center.df[,c("x", "y", "fade")]), bl.patch.hexid.centroids )
    
    ## least cost path calculations  ----

values(hab.cost.lcm) = dispers.costs$scaled.ecolog.cost[values(lcm.landscape)]*cost.scale.factor
#replace gaps with high-cost landscape - these are normally sea, but can be beyond edge of landscape, this shouldnt be a problem, becasue landscape is buffered
values(hab.cost.lcm)[is.na(values(hab.cost.lcm))] = max(values(hab.cost.lcm),na.rm=T)

gridbuff_Tcost <- transition(hab.cost.lcm,function(x) 1/mean(x),8)
gridbuff_Tcost.C <- geoCorrection(gridbuff_Tcost, type="c", multpl=FALSE, scl=F)

bigger.patches = bl.patch.hexid.centroids$effective.ha> median (bl.patch.hexid.centroids$effective.ha)
eg.patch = which.max(as.numeric(1/st_distance(bl.patch.hexid.centroids
                                              , st_centroid(Tscapes01), by_element = T))*bigger.patches)
#which(bl.patch.hexid.centroids$effective.ha == max(bl.patch.hexid.centroids$effective.ha ))

eg.shortest.paths = shortestPath(gridbuff_Tcost.C, bl.patch.hexid.centroids.sp[eg.patch,1], bl.patch.hexid.centroids.sp, 
                                 output="SpatialLines")
eg.shortest.paths.focal = shortestPath(gridbuff_Tcost.C, bl.patch.hexid.centroids.sp[eg.patch,1], bl.patch.hexid.centroids.sp[bl.patch.hexid.centroids.sp$in.landscape,], 
                                       output="SpatialLines")

plot(hab.cost.lcm)
lines(eg.shortest.paths)

eg.accum.cost <- accCost(gridbuff_Tcost.C, bl.patch.hexid.centroids.sp[eg.patch,1])
eg.accum.cost <- accCost(gridbuff_Tcost.C, bl.patch.hexid.centroids.sp[eg.patch,1])
plot(eg.accum.cost)

eg.accum.cost.paths = as(eg.shortest.paths, "SpatialPoints") %>%
  rasterize(., eg.accum.cost, mask=TRUE) 
eg.accum.cost.paths.focal = as(eg.shortest.paths.focal, "SpatialPoints") %>%
  rasterize(., eg.accum.cost, mask=TRUE) %>% crop(Tscapes01) %>% mask(Tscapes01)

plot(eg.accum.cost.paths.focal)

bl.patch.hexid.centroids$eg.cost.toget = effective.distance[eg.patch,]

# could do fainter and fainter as p (connected) decreases

    ## create sf of ECA circle and effective area of patches ----
# pi*r^2 = ECA, sqrt(ECA/pi) = r

effective.area.circle.patch = st_buffer(bl.patch.hexid.centroids, dist = sqrt(bl.patch.hexid.centroids$effective.ha*10000/pi))# *10000 becasue ha to m
effective.area.circle.patch$qual.score = 1

eca.circle.patch = st_buffer(st_centroid(Tscapes01), dist = sqrt(landscape.metrics$leastcost.ECA*10000/pi))# *10000 becasue ha to m
eca.circle.patch$qual.score = 1

# plots ----
main.title = NULL
sub.title = NULL
quality.scale.title = "Relative habtiat quality"

patch.line.size = 0.25
patch.point.size = 3

patch.dat = bl.patch.id.poly.hexid
focal.centroid = bl.patch.hexid.centroids[eg.patch,]
eg.toget.patches = bl.patch.hexid.centroids[-eg.patch,]
eg.toget.patches.focal = eg.toget.patches[eg.toget.patches$in.landscape,]


  ## geoms and scales ----

    ### land area ----
    
    wider.fade.geom = ts.fig.buff %>% 
      geom_sf(data = ., mapping = aes(), fill = "white", alpha = 0.2, colour = NA, size = 0 )
    area.fade.geom = Tscapes01 %>% 
      geom_sf(data = ., mapping = aes(), fill = "white", alpha = 0.2, colour = NA, size = 0 )
    
    area.boarder.geom = Tscapes01 %>% 
      geom_sf(data = ., mapping = aes(), fill = "transparent", colour = "black", size = 1 )
    
    outter.boarder.geom = ts.fig.buff %>% 
      geom_sf(data = ., mapping = aes(), fill = "transparent", colour = "white", size = patch.line.size )
    
    
    ### land cover ----
    
    lcm.geom = geom_raster(data = lcm.df, aes(x = x, y = y, fill = layer ))
    lcm.geom.focal = geom_raster(data = lcm.df.focal, aes(x = x, y = y, fill = layer ))
    
    lcm.colscale = scale_fill_manual(values = ceh.col.pallette[lcm.df$layer  %>% levels() %>% as.numeric()], 
                                     labels = ceh.full.habtype[lcm.df$layer  %>% levels() %>% as.numeric()],
                                     name = "Land cover type", na.value="transparent",
                                     guide = guide_legend(title.position = "top",  title.hjust = 0))
    lcm.focal.colscale = scale_fill_manual(values = ceh.col.pallette[lcm.df.focal$layer  %>% levels() %>% as.numeric()], 
                                           labels = ceh.full.habtype[lcm.df.focal$layer  %>% levels() %>% as.numeric()],
                                           name = "Land cover type", na.value="transparent",
                                           guide = guide_legend(title.position = "top",  title.hjust = 0))
    lcm.colscale.noguide = scale_fill_manual(values = ceh.col.pallette[lcm.df$layer  %>% levels() %>% as.numeric()], 
                                             labels = ceh.full.habtype[lcm.df$layer  %>% levels() %>% as.numeric()],
                                             name = "Land cover type", na.value="transparent",
                                             guide = "none")
    lcm.focal.colscale.noguide = scale_fill_manual(values = ceh.col.pallette[lcm.df.focal$layer  %>% levels() %>% as.numeric()], 
                                                   labels = ceh.full.habtype[lcm.df.focal$layer  %>% levels() %>% as.numeric()],
                                                   name = "Land cover type", na.value="transparent",
                                                   guide = "none")
    
    ### Distance to landscape surround fade ----
    
    dist.to.center.fade.geom = geom_raster(data = dist.to.center.df, aes(x = x, y = y, alpha = fade ), fill = "white")
    dist.to.center.fade.scale = scale_alpha_continuous(range = c(0,1),
                                                       guide = "none")
    
    ### land cover cost ----
    cost.geom = function(scape = "focal") {
      if(scape == "focal"){lcm.df = lcm.df.focal } 
      if(scape == "wider"){lcm.df = lcm.df}
      geom_raster(data = lcm.df, aes(x = x, y = y, fill = cost*cost.scale.factor ))
    }
    
    cost.colscale = scale_fill_gradient(low = "white", high = "grey20", name = "Landscape\ndispersal cost",
                                        limits = c(0, NA), breaks = pretty_breaks(3),
                                        guide = guide_colorbar(title.position = "top",  title.hjust = 0))
    
    ### paths ----
    # line paths
    eg.shortest.paths.geom = eg.shortest.paths %>% st_as_sf() %>% 
      geom_sf(data = ., mapping = aes(), colour = "blue")
    eg.shortest.paths.geom.focal = eg.shortest.paths.focal %>% st_as_sf() %>% 
      geom_sf(data = ., mapping = aes(), colour = "blue")
    
    # path accum cost to get
    eg.accum_cost_path.geom = function(scape = "focal") {
      if(scape == "focal"){accum.cost.path = eg.accum.cost.paths.focal %>% as.data.frame(., xy = TRUE) } 
      if(scape == "wider"){accum.cost.path = eg.accum.cost.paths %>% as.data.frame(., xy = TRUE)}
      geom_raster(data = accum.cost.path, aes(x = x, y = y, fill = layer/1000 ))
    }
    
    eg.cost_path.colscale = scale_fill_gradient(low = "yellow" ,high = "#FF0000",
                                                name = "Accumulated\ndispersal cost (km)",
                                                limits = c(0,NA),
                                                na.value = NA,
                                                breaks = pretty_breaks(3),
                                                guide = "none")
    
    ### accum cost/p connect/contribution ----
    #### point cost/p to get/contribution ----
    
    eg.toget.geom = function(scape = "focal", toget = "cost") {
      if(scape == "focal"){cost.points = eg.toget.patches.focal
      cost.points$fade = NA} 
      if(scape == "wider"){cost.points = eg.toget.patches}
      if(toget == "cost"){fill.col = cost.points$eg.cost.toget/1000 } 
      if(toget == "prob"){fill.col = exp(-alpha * cost.points$eg.cost.toget)}
      if(toget == "contrib"){fill.col = sqrt(exp(-alpha * cost.points$eg.cost.toget) * 
                                               cost.points$effective.ha * 
                                               cost.points$effective.ha)}
      geom_sf(data = cost.points,
              mapping = aes(fill = fill.col, alpha = 1-fade),
              size = patch.point.size, shape = 21)
      # +
      #   scale_alpha(limits = c(0,1),guide = "none")
      #  # if(scape == "focal"){aa } 
      # if(scape == "wider"){geom_sf(data = cost.points, 
      #                              mapping = aes(fill = fill.col, alpha = fade),
      #                              size = patch.point.size, shape = 21) + aa}
    }
    
    eg.cost.toget.scale = scale_fill_gradient(low = "yellow" ,high = "#FF0000",
                                              name = "Accumulated\ndispersal cost (km)",
                                              limits = c(0,NA),
                                              na.value = NA,
                                              breaks = pretty_breaks(3),
                                              guide = guide_colorbar(title.position = "top",  title.hjust = 0)) 
    
    # point prob to get
    eg.connectivity_prob.scale = scale_fill_gradient(low = "blue", high = "grey94", limits = c(0,max(exp(-alpha * eg.toget.patches$eg.cost.toget))),
                                                     name = "Probability of\nconnectivity",
                                                     breaks = pretty_breaks(3),
                                                     guide = guide_colorbar(title.position = "top",  title.hjust = 0))
    
    eg.contribution.geom = geom_sf(data = eg.toget.patches, 
                                   mapping = aes(# eg.cost.toget, 
                                     fill = sqrt(exp(-alpha * eg.toget.patches$eg.cost.toget) * 
                                                   eg.toget.patches$effective.ha * 
                                                   focal.centroid$effective.ha)),
                                   size = patch.point.size, shape = 21
    )
    eg.contribution.scale = scale_fill_gradient(low = "grey90", high = "darkgreen", #limits = c(0.00001,max(sqrt(exp(-alpha * eg.toget.patches$eg.cost.toget) * eg.toget.patches$effective.ha * focal.centroid$effect ))),
                                                name = "Contribution to example patch",
                                                # trans = "log2",
                                                breaks = pretty_breaks(5),
                                                guide = guide_colorbar(title.position = "top",  title.hjust = 0.5))
    
    
    
    
    #### graph edge weight - p connect ----
    m.connect.prob = exp(-alpha * effective.distance)
    rownames(m.connect.prob) = 1:dim(m.connect.prob)[1]
    colnames(m.connect.prob) = 1:dim(m.connect.prob)[2]
    
    m.eg.connect.prob = m.connect.prob
    m.eg.connect.prob [-eg.patch, -eg.patch] = 0
    m.eg.connect.prob.focal = m.eg.connect.prob[bl.patch.hexid.centroids$in.landscape, bl.patch.hexid.centroids$in.landscape]
    
    m.connect.prob.focal = m.connect.prob[bl.patch.hexid.centroids$in.landscape, bl.patch.hexid.centroids$in.landscape]
    m.eg.connect.prob.focal = m.eg.connect.prob[bl.patch.hexid.centroids$in.landscape, bl.patch.hexid.centroids$in.landscape]
    
    df.connect.prob = m.connect.prob %>% as.data.frame %>% tibble::rownames_to_column() %>% 
      mutate(from.x = st_coordinates(bl.patch.hexid.centroids)[,1],
             from.y = st_coordinates(bl.patch.hexid.centroids)[,2],
             from.in.landscape = bl.patch.hexid.centroids$in.landscape,
             from.effective.ha = bl.patch.hexid.centroids$effective.ha) %>% 
      tidyr::pivot_longer(c(-rowname, -from.x, -from.y, -from.in.landscape, -from.effective.ha)) %>% 
      mutate(to.x = rep(st_coordinates(bl.patch.hexid.centroids)[,1], times = dim(bl.patch.hexid.centroids)[1]),
             to.y = rep(st_coordinates(bl.patch.hexid.centroids)[,2], times = dim(bl.patch.hexid.centroids)[1])) %>% 
      rename(from.patch = rowname,
             to.patch = name,
             p.connect = value)
    df.connect.prob$from.eg.patch = F
    df.connect.prob$from.eg.patch[df.connect.prob$from.patch == eg.patch] = T
    df.connect.prob$to.landscape = F
    df.connect.prob$to.landscape[df.connect.prob$to.patch %in% unique(df.connect.prob$from.patch[df.connect.prob$from.in.landscape])] = T
    
    p.connect.graphedge.geom = function(scape = "focal", connectors = "single", edge.col = "gray40", edge.alpha = 0.5) {
      if(scape == "focal"){patch.connections = df.connect.prob[df.connect.prob$from.in.landscape & df.connect.prob$to.landscape,]  } 
      if(scape == "wider"){patch.connections = df.connect.prob}
      if(connectors == "single"){patch.connections = patch.connections[patch.connections$from.eg.patch,]  } 
      if(connectors == "all"){patch.connections = patch.connections}
      
      patch.connections %>% 
        geom_edges(color=edge.col, alpha = edge.alpha,mapping = aes(x = from.x, y = from.y, xend = to.x, yend = to.y, size = p.connect),
                   curvature = 0)
    }
    p.connect.graphedge.scale = scale_size_continuous(range = c(0,3),
                                                      name = "Probability of\nconnectivity", limits = c(0,max(df.connect.prob$p.connect, na.rm = T)),
                                                      breaks = pretty_breaks(3),
                                                      guide = guide_legend(title.position = "top",  title.hjust = 0))
    
    
    ### edge ----
    edge.geom = function(scape = "focal", pattern_spacing = 0.01, size = patch.line.size) {
      if(scape == "focal"){patches = patch.edge %>% st_intersection(Tscapes01) } 
      if(scape == "wider"){patches = patch.edge %>% st_intersection(ts.fig.buff)}
      geom_sf_pattern(data = patches,  mapping = aes(), fill = NA, 
                      colour = "black", pattern = 'stripe',
                      pattern_spacing = pattern_spacing,
                      pattern_fill = "black",
                      size = size )
    }
    
    
    ### awi ----
    awi.geom <- function(scape = "focal") {
      if(scape == "focal"){patches = awi.landscape %>% st_intersection(Tscapes01) } 
      if(scape == "wider"){patches = awi.landscape %>% st_intersection(ts.fig.buff)}
      geom_sf(data = patches, mapping = aes(), fill = "yellow", alpha = 0.5,
              colour = "black", 
              size = patch.line.size )
    }
    
    ### patches ----
    patch_point.geom = geom_sf(data = bl.patch.hexid.centroids, mapping = aes(),  
                               colour = "black", size = patch.line.size )
    
    all.patch.geom <- function(scape = "focal", fill.col = ceh.col.pallette[ceh.full.habtype == "Broadleaf" ]) {
      if(scape == "focal"){patches = bl.patch.id.poly.hexid %>% st_intersection(Tscapes01) }
      if(scape == "wider"){patches = bl.patch.id.poly.hexid %>% st_intersection(ts.fig.buff) } 
      geom_sf(data = patches, mapping = aes(), fill = fill.col
              ,  colour = "black", size = patch.line.size )
    }
    
    # focal patch star
    focal.patch.geom =  focal.centroid %>% st_coordinates() %>% as.data.frame() %>% 
      geom_star(data = ., 
                aes(x = X, y = Y), fill = "yellow", colour = "black", size = patch.point.size*2
      )
    
    #### patch area ----
    patch_area.geom <- function(scape = "focal", fill.col = "subpatch.hex.ha") {
      if(scape == "focal"){patches = bl.patch.id.poly.hexid %>% st_intersection(Tscapes01) } 
      if(scape == "wider"){patches = bl.patch.id.poly.hexid %>% st_intersection(ts.fig.buff)}
      geom_sf(data = patches,  mapping = aes(fill = as.data.frame(patches)[,fill.col]), colour = "black", size = patch.line.size )
    }
    
    patch_area.colscale = scale_fill_gradient(low = "white" ,high = ceh.col.pallette[ceh.full.habtype == "Broadleaf" ],#option = "plasma", 
                                              limits = c(0, max(bl.patch.id.poly.hexid$subpatch.hex.ha)),
                                              name = "Patch area (ha)",
                                              breaks = pretty_breaks(3),
                                              guide = guide_colorbar(title.position = "top",  title.hjust = 0.5))
    #### patch quality ----
    # habtiat
    hab.quality.geom <- function(scape = "focal", fill.col = "qual.score") {
      if(scape == "focal"){patches = patch.qual.sf %>% st_intersection(Tscapes01) } 
      if(scape == "wider"){patches = patch.qual.sf %>% st_intersection(ts.fig.buff)}
      geom_sf(data = patches, mapping = aes(fill = as.data.frame(patches)[,fill.col]),
              colour = "black", size = patch.line.size )
    }
    
    hab.quality.colscale = scale_fill_gradient(low = "white" ,high = brewer.pal(n=9,"Blues")[9],#option = "plasma", 
                                               limits = c(0, max(patch.qual.sf$qual.score)),
                                               name = "Habitat quality",
                                               breaks = pretty_breaks(3),
                                               guide = guide_colorbar(title.position = "top",  title.hjust = 0)) 
    # mean for patch
    av_patch_quality.geom = function(scape = "focal", fill.col = "mean.quality") {
      if(scape == "focal"){patches = bl.patch.id.poly.hexid %>% st_intersection(Tscapes01) } 
      if(scape == "wider"){patches = bl.patch.id.poly.hexid %>% st_intersection(ts.fig.buff)}
      geom_sf(data = patches, mapping = aes(fill = as.data.frame(patches)[,fill.col]),
              colour = "black", size = patch.line.size )
    }
    
    av_patch_quality.colscale = scale_fill_gradient(low = "white" ,high = brewer.pal(n=9,"Blues")[9],#option = "plasma", 
                                                    limits = c(0, NA),
                                                    name = "Mean habitat quality",
                                                    breaks = pretty_breaks(3),
                                                    guide = guide_colorbar(title.position = "top",  title.hjust = 0))
    
    #### patch quality weighted/effective area ----
    effective_area.colscale = scale_fill_gradient(low = "white" ,high = brewer.pal(n=9,"Purples")[9],#option = "plasma", 
                                                  limits = c(0, NA),
                                                  name = "Quality-weighted area (ha)",
                                                  breaks = pretty_breaks(3),
                                                  guide = guide_colorbar(title.position = "top",  title.hjust = 0.5))
    
    effective_area_circlepatch.geom.focal = effective.area.circle.patch[effective.area.circle.patch$in.landscape,] %>% 
      st_intersection(Tscapes01) %>% 
      geom_sf(data = ., mapping = aes(fill = qual.score),
              colour = "black", size = patch.line.size )
    
    effective_area_circlepatch.geom.wider = effective.area.circle.patch %>% 
      st_intersection(ts.fig.buff)  %>% 
      geom_sf(data = ., mapping = aes(fill = qual.score, alpha = 1-fade),
              colour = NA, size = patch.line.size )
    
    
    ### patch contribution to landscape ECA ----
    # decided not to do this, feel that it misses a point
    # the whole idea is that this isnt on a per-patch level
    # looked at contribution per unit area, but small patches right next to big patches have massive, becasue they are probabilistically almost the same patch, but corrected for their tinty area
    
    ### ECA ----
    eca.geom =  geom_sf(data = eca.circle.patch, mapping = aes(fill = qual.score),
                        colour = "black", 
                        size = patch.line.size )
    eca.colscale = scale_fill_viridis_c(#option = "plasma", 
      name = "Habitat quality",
      breaks = pretty_breaks(3),
      guide = guide_colorbar(title.position = "top",  title.hjust = 0))
    
  # plots ------
p.default =   ggplot() +
  theme_map(leg.tit.size = 11,
            legend.position = "right")

      ## p01-02 lcm and patch id ----
      p01.01.lcm.focal = p.default +
        lcm.geom.focal + lcm.focal.colscale  + new_scale_fill() +
        area.boarder.geom +
        labs(x = NULL, y = NULL , title = paste("CEH Land Cover Map", this.year), subtitle = sub.title#, caption = ""
        )+
        coord_sf(expand = FALSE)
      
      p01.02.lcm.wider =  p.default +
        lcm.geom + lcm.colscale  + new_scale_fill() +
        dist.to.center.fade.geom + dist.to.center.fade.scale +
        area.boarder.geom +
        labs(x = NULL, y = NULL , title = paste("CEH Land Cover Map", this.year), subtitle = sub.title#, caption = ""
        )+
        coord_sf(expand = FALSE)
      
      p02.01.lcm_patch.focal = p.default +
        lcm.geom.focal + lcm.focal.colscale  + new_scale_fill() +
        area.fade.geom +
        all.patch.geom() +
        area.boarder.geom +
        labs(x = NULL, y = NULL , title = paste("Native woodland patches and surrounding matrix\nfrom CEH Land Cover Map"), subtitle = sub.title#, caption = ""
        )+
        coord_sf(expand = FALSE)
      
      p02.02.lcm_patch.wider = p.default +
        lcm.geom + lcm.colscale  + new_scale_fill() +
        wider.fade.geom +
        all.patch.geom(scape = "wider") +
        dist.to.center.fade.geom + dist.to.center.fade.scale +
        outter.boarder.geom + area.boarder.geom +
        labs(x = NULL, y = NULL , title = paste("Native woodland patches"), subtitle = sub.title#, caption = ""
        )+
        coord_sf(expand = FALSE)
      
      ## p03 patch area ----
      p03.01.lcm_patcharea.focal = p.default +
        # lcm.geom.focal + lcm.focal.colscale.noguide  + new_scale_fill() +
        patch_area.geom(scape = "focal") + patch_area.colscale + new_scale_fill() +
        area.boarder.geom +
        labs(x = NULL, y = NULL , title = paste("Patch area"), subtitle = sub.title#, caption = ""
        )+
        coord_sf(expand = FALSE)
      
      p03.02.lcm_patcharea.wider = p.default +
        # lcm.geom.focal + lcm.focal.colscale.noguide  + new_scale_fill() +
        patch_area.geom(scape = "wider") + patch_area.colscale + new_scale_fill() +
        dist.to.center.fade.geom + dist.to.center.fade.scale + new_scale_fill() +
        outter.boarder.geom + area.boarder.geom +
        labs(x = NULL, y = NULL , title = paste("Patch area"), subtitle = sub.title#, caption = ""
        )+
        coord_sf(expand = FALSE)
      
      ## p04-05 edge ancient woodland ----
      p04.01.01.lcm_patch_edge.focal = p.default +
        lcm.geom.focal + lcm.focal.colscale  + new_scale_fill() +
        area.fade.geom + 
        all.patch.geom(scape = "focal") +
        edge.geom(scape = "focal", pattern_spacing = 0.02, size = patch.line.size/2) + 
        area.boarder.geom +
        labs(x = NULL, y = NULL , title = paste("Unfavourable patch edge\nimpacted by neighbouring land use"), subtitle = sub.title    )+#, caption = ""
        coord_sf(expand = FALSE) 
      
      p04.02.01.lcm_patch_aw.focal = p.default +
        lcm.geom.focal + lcm.focal.colscale  + new_scale_fill() +
        area.fade.geom + 
        all.patch.geom(scape = "focal") +
        awi.geom(scape = "focal") +
        area.boarder.geom +
        labs(x = NULL, y = NULL , title = paste("Ancient woodland"), subtitle = sub.title#, caption = ""
        ) +
        coord_sf(expand = FALSE)
      
      p05.01.lcm_patch_edge_awi.focal = p04.02.01.lcm_patch_aw.focal +
        edge.geom(scape = "focal") + 
        area.boarder.geom +
        labs(x = NULL, y = NULL , title = paste("Unfavourable edge + ancient woodland"), subtitle = sub.title#, caption = ""
        )+
        coord_sf(expand = FALSE)
      
      ## p06-07 habitat quality ----
      p06.01.lcm_patchHab.focal = p.default +
        # lcm.geom + lcm.colscale.noguide + new_scale_fill() +
        hab.quality.geom(scape = "focal",fill.col = "qual.score" ) + hab.quality.colscale + new_scale_fill() +  
        area.boarder.geom +
        labs(x = NULL, y = NULL , title = paste("Quality of patch components"), subtitle = sub.title#, caption = ""
        )+
        coord_sf(expand = FALSE)
      
      p07.01.lcm_patchMeanHab.focal = p.default +
        # lcm.geom + lcm.colscale.noguide + new_scale_fill() +
        av_patch_quality.geom(scape = "focal" )   + hab.quality.colscale + new_scale_fill() +    
        area.boarder.geom +
        labs(x = NULL, y = NULL , title = paste("Mean patch quality"), subtitle = sub.title#, caption = ""
        )+
        coord_sf(expand = FALSE)
      
      ## p08 quality-weighted equivelent patch area ----
      p08.01.01.lcm_patchQualAreaFill.focal = p.default +
        # lcm.geom + lcm.colscale.noguide + new_scale_fill() +
        av_patch_quality.geom(scape = "focal", fill.col = "effective.ha" )  + effective_area.colscale + new_scale_fill() +   
        area.boarder.geom +
        labs(x = NULL, y = NULL , title = paste("Quality-weighted patch area"), subtitle = sub.title#, caption = ""
        )+
        coord_sf(expand = FALSE) 
      
      p08.02.01.patchQualAreaSize.focal = p.default +
        # lcm.geom + lcm.colscale.noguide + new_scale_fill() +
        effective_area_circlepatch.geom.focal + hab.quality.colscale  + new_scale_fill() +   
        area.boarder.geom +
        labs(x = NULL, y = NULL , title = paste("Quality-weighted patch area"), subtitle = sub.title#, caption = ""
        )+
        coord_sf(expand = FALSE) 
      
      p08.02.02.patchQualAreaSize.wider = p.default +
        # lcm.geom + lcm.colscale.noguide + new_scale_fill() +
        effective_area_circlepatch.geom.wider + hab.quality.colscale + scale_alpha(guide = "none")  + new_scale_fill() +   
        effective_area_circlepatch.geom.focal + hab.quality.colscale  + new_scale_fill() +   
        # dist.to.center.fade.geom + dist.to.center.fade.scale + new_scale_fill() +
        outter.boarder.geom + area.boarder.geom +
        labs(x = NULL, y = NULL , title = paste("Quality-weighted patch area"), subtitle = sub.title#, caption = ""
        )+
        coord_sf(expand = FALSE) 
      
      ## p09-10 dispersal cost ----
      p09.01.cost.focal = p.default +
        cost.geom(scape = "focal") + cost.colscale + new_scale_fill() +
        all.patch.geom(scape = "focal", fill.col = "transparent")  + 
        area.boarder.geom +
        labs(x = NULL, y = NULL , title = paste("Dispersal cost layer"), subtitle = sub.title#, caption = ""
        )+
        coord_sf(expand = FALSE) 
      
      p09.02.cost.wider = p.default +
        cost.geom(scape = "wider") + cost.colscale + new_scale_fill() +
        all.patch.geom(scape = "wider", fill.col = "transparent")  + 
        dist.to.center.fade.geom + dist.to.center.fade.scale + new_scale_fill() +
        outter.boarder.geom + area.boarder.geom +
        labs(x = NULL, y = NULL , title = paste("Dispersal cost layer"), subtitle = sub.title#, caption = ""
        )+
        coord_sf(expand = FALSE) 
      
      # p10.01.cost.patch.focal = p.default +
      #   cost.geom(scape = "focal") + cost.colscale + new_scale_fill() +
      #   effective_area.geom.focal + effective_area.colscale + new_scale_fill() +
      #   area.boarder.geom +
      #   labs(x = NULL, y = NULL , title = paste("Dispersal cost layer"), subtitle = sub.title#, caption = ""
      #   )+
      #   coord_sf(expand = FALSE)  
      # p10.02.cost.patch.wider = p.default +
      #   cost.geom(scape = "wider") + cost.colscale + new_scale_fill() +
      #   effective_area.geom + effective_area.colscale + new_scale_fill() +
      #   dist.to.center.fade.geom + dist.to.center.fade.scale + new_scale_fill() +
      #   area.boarder.geom +
      #   labs(x = NULL, y = NULL , title = paste("Dispersal cost layer"), subtitle = sub.title#, caption = ""
      #   )+
      #   coord_sf(expand = FALSE)  
      ## p11 eg patch accum cost -----
      p11.01.01.cost.patch.accumCost.focal = p.default +
        cost.geom(scape = "focal") + cost.colscale + new_scale_fill() +
        all.patch.geom(scape = "focal", fill.col = "transparent")  +
        # # effective_area.geom.focal + effective_area.colscale + new_scale_fill() +
        eg.accum_cost_path.geom(scape = "focal") + eg.cost_path.colscale + guides(fill = "none") + new_scale_fill() +
        eg.toget.geom(scape = "focal", toget = "cost") + guides(alpha = "none") + eg.cost.toget.scale + new_scale_fill() +
        # dist.to.center.fade.geom + dist.to.center.fade.scale + new_scale_fill() +
        area.boarder.geom +
        focal.patch.geom +
        labs(x = NULL, y = NULL , title = paste("Dispersal through the matrix from an example patch"), subtitle = sub.title#, caption = ""
        )+
        coord_sf(expand = FALSE) 
      
      p11.01.02.cost.patch.accumCost.wider = p.default +
        cost.geom(scape = "wider") + cost.colscale + new_scale_fill() +
        all.patch.geom(scape = "wider", fill.col = "transparent")  +
        # # effective_area.geom.focal + effective_area.colscale + new_scale_fill() +
        dist.to.center.fade.geom + dist.to.center.fade.scale + new_scale_fill() +
        eg.accum_cost_path.geom(scape = "wider") + eg.cost_path.colscale + guides(fill = "none") + new_scale_fill() +
        eg.toget.geom(scape = "wider", toget = "cost") + guides(alpha = "none") + eg.cost.toget.scale + new_scale_fill() +
        outter.boarder.geom + area.boarder.geom +
        focal.patch.geom +
        labs(x = NULL, y = NULL , title = paste("Dispersal through the matrix from an example patch"), subtitle = sub.title#, caption = ""
        )+
        coord_sf(expand = FALSE) 
      
      p11.02.01.patch.accumCost.focal = p.default +
        # cost.geom(scape = "focal") + cost.colscale + new_scale_fill() +
        # all.patch.geom(scape = "focal", fill.col = "transparent")  +
        # # effective_area.geom.focal + effective_area.colscale + new_scale_fill() +
        eg.accum_cost_path.geom(scape = "focal") + eg.cost_path.colscale + guides(fill = "none") + new_scale_fill() +
        eg.toget.geom(scape = "focal", toget = "cost") + guides(alpha = "none") + eg.cost.toget.scale + new_scale_fill() +
        # dist.to.center.fade.geom + dist.to.center.fade.scale + new_scale_fill() +
        area.boarder.geom +
        focal.patch.geom +
        labs(x = NULL, y = NULL , title = paste("Dispersal cost from an example patch"), subtitle = sub.title#, caption = ""
        )+
        coord_sf(expand = FALSE) 
      
      p11.03.01.cost.accumCost.pconnect.focal = p.default +
        # cost.geom(scape = "focal") + cost.colscale + new_scale_fill() +
        # all.patch.geom(scape = "focal", fill.col = "transparent")  +
        # # effective_area.geom.focal + effective_area.colscale + new_scale_fill() +
        p.connect.graphedge.geom(scape = "focal", connectors = "single") + p.connect.graphedge.scale +
        eg.toget.geom(scape = "focal", toget = "cost") + guides(alpha = "none") + eg.cost.toget.scale + new_scale_fill() +
        # dist.to.center.fade.geom + dist.to.center.fade.scale + new_scale_fill() +
        area.boarder.geom +
        focal.patch.geom + 
        labs(x = NULL, y = NULL , title = paste("Probability of connectivity between example patch\nand other patches in the landscape"), subtitle = sub.title#, caption = ""
        )+
        coord_sf(expand = FALSE) 
      
      
      ## p12-13 p connect ----
      p12.01.01.cost.patch.pconnect.focal.eg = p.default +
        # lcm.geom + lcm.colscale.noguide + new_scale_fill() +
        p.connect.graphedge.geom(scape = "focal", connectors = "single") + p.connect.graphedge.scale +
        effective_area_circlepatch.geom.focal + hab.quality.colscale  + new_scale_fill() +   
        focal.patch.geom +
        area.boarder.geom +
        labs(x = NULL, y = NULL , title = paste("Probability of connectivity between example patch\n and others in the landscape"), subtitle = sub.title#, caption = ""
        )+
        coord_sf(expand = FALSE)
      
      p12.01.02.cost.patch.pconnect.focal.all = p.default +
        # lcm.geom + lcm.colscale.noguide + new_scale_fill() +
        p.connect.graphedge.geom(scape = "focal", connectors = "all") + p.connect.graphedge.scale +
        effective_area_circlepatch.geom.focal + hab.quality.colscale  + new_scale_fill() +   
        # focal.patch.geom +
        area.boarder.geom +
        labs(x = NULL, y = NULL , title = paste("Probability of connectivity\nbetween all patches in the landscape"), subtitle = sub.title#, caption = ""
        )+
        coord_sf(expand = FALSE)
      
      p12.02.01.cost.patch.pconnect.wider.eg =     p.default +
        # lcm.geom + lcm.colscale.noguide + new_scale_fill() +
        p.connect.graphedge.geom(scape = "wider", connectors = "single", edge.col = "gray80") + p.connect.graphedge.scale +
        p.connect.graphedge.geom(scape = "focal", connectors = "single") + p.connect.graphedge.scale + scale_alpha(guide = "none") +
        
        effective_area_circlepatch.geom.wider + hab.quality.colscale + scale_alpha(guide = "none")  + new_scale_fill() +   
        effective_area_circlepatch.geom.focal + hab.quality.colscale  + new_scale_fill() +   
        # dist.to.center.fade.geom + dist.to.center.fade.scale + new_scale_fill() +
        focal.patch.geom +
        outter.boarder.geom + area.boarder.geom +
        
        labs(x = NULL, y = NULL , title = paste("Probability of connectivity between example patch\n and others in the wider landscape with different quality-weighted area"), subtitle = sub.title#, caption = ""
        )+
        coord_sf(expand = FALSE)
      
      
      p12.02.02.cost.patch.pconnect.wider.all =     p.default +
        # lcm.geom + lcm.colscale.noguide + new_scale_fill() +
        p.connect.graphedge.geom(scape = "wider", connectors = "all", edge.col = "gray80") + p.connect.graphedge.scale +
        p.connect.graphedge.geom(scape = "focal", connectors = "all") + p.connect.graphedge.scale + scale_alpha(guide = "none") +
        
        effective_area_circlepatch.geom.wider + hab.quality.colscale + scale_alpha(guide = "none")  + new_scale_fill() +   
        effective_area_circlepatch.geom.focal + hab.quality.colscale  + new_scale_fill() +   
        # dist.to.center.fade.geom + dist.to.center.fade.scale + new_scale_fill() +
        outter.boarder.geom + area.boarder.geom +
        
        labs(x = NULL, y = NULL , title = paste("Probability of connectivity of all patches in the landscape to each other\nand to patches in the surrounding landscape"), subtitle = sub.title#, caption = ""
        )+
        coord_sf(expand = FALSE)
      
      
      ## p14 contribution ----
      
      # 
      # 
      # p14.contribution = p.default +
      #     # cost.geom +
      #     # cost.colscale +
      #     # new_scale_fill() +
      #     
      #     effective_area.geom +
      #     effective_area.colscale +      
      #     new_scale_fill() +
      #   area.boarder.geom +
      #   
      #     eg.shortest.paths.geom +
      #     new_scale_fill() +
      #     
      #     eg.contribution.geom +
      #     eg.contribution.scale +
      #     # new_scale_fill() +
      #     
      #     focal.patch.geom +
      #     
      # 
      #     labs(x = NULL, y = NULL , title = paste("Equivelent connected area of landscape"), subtitle = sub.title#, caption = ""
      #     )+
      #     coord_sf(expand = FALSE)  
      ## p15 eca ----
      p15.eca.focal = p12.02.02.cost.patch.pconnect.wider.all +
        # cost.geom +
        # cost.colscale +
        # new_scale_fill() +
        geom_sf(ts.fig.buff, mapping = aes(), colour = "white", fill = "white", alpha = 0.8) +
        eca.geom + hab.quality.colscale +
        area.boarder.geom +
        labs(x = NULL, y = NULL , title = paste("Equivelent Connected Area (Probability of Connectivity)\nof native woodland in the landscape"), subtitle = sub.title#, caption = ""
        )+
        coord_sf(expand = FALSE) 
      
      
#save plots ----

h.a4.mm = 80
w.a4.mm = 50

ggsave(p04.01.01.lcm_patch_edge.focal+theme(legend.position="none",
                                            plot.title = element_blank()), 
       file = paste0("analysis outputs\\", ts.lcm.names, "\\" , this.year, "\\for a4\\just figs\\", 
                     "p04.01.01.lcm_patch_edge.focal.png"), 
       width = w.a4.mm, height = h.a4.mm, units = "mm")


plots.a4.list = list(p02.01.lcm_patch.focal ,
                     p04.01.01.lcm_patch_edge.focal  ,
                     p04.02.01.lcm_patch_aw.focal ,
                     p05.01.lcm_patch_edge_awi.focal,
                     p06.01.lcm_patchHab.focal ,
                     p07.01.lcm_patchMeanHab.focal ,
                     p08.02.01.patchQualAreaSize.focal ,
                     p09.01.cost.focal ,
                     p11.01.01.cost.patch.accumCost.focal ,
                     p11.02.01.patch.accumCost.focal ,
                     p11.03.01.cost.accumCost.pconnect.focal ,
                     p12.01.01.cost.patch.pconnect.focal.eg ,
                     p12.01.02.cost.patch.pconnect.focal.all ,
                     p12.02.02.cost.patch.pconnect.wider.all ,
                     p15.eca.focal + annotation_scale(style = "ticks"))

plots.for.a4fig = align_plots(plotlist = lapply(plots.a4.list, function(x){
  x +theme(legend.position= "right",
           plot.title = element_blank()
           )}),
                              align="hv", axis="tblr")

legends.for.a4fig = lapply(plots.a4.list, function(x){x %>% get_legend() %>% ggdraw()})

plots.for.a4fig.names = c("p02.01.lcm_patch.focal",
                          "p04.01.01.lcm_patch_edge.focal",
                          "p04.02.01.lcm_patch_aw.focal",
                          "p05.01.lcm_patch_edge_awi.focal",
                          "p06.01.lcm_patchHab.focal",
                          "p07.01.lcm_patchMeanHab.focal",
                          "p08.02.01.patchQualAreaSize.focal",
                          "p09.01.cost.focal",
                          "p11.01.01.cost.patch.accumCost.focal",
                          "p11.02.01.patch.accumCost.focal",
                          "p11.03.01.cost.accumCost.pconnect.focal",
                          "p12.01.01.cost.patch.pconnect.focal.eg",
                          "p12.01.02.cost.patch.pconnect.focal.all",
                          "p12.02.02.cost.patch.pconnect.wider.all",
                          "p15.eca.focal")    

dir.create(paste0("analysis outputs\\", ts.lcm.names, "\\" , this.year, "\\for a4"))
dir.create(paste0("analysis outputs\\", ts.lcm.names, "\\" , this.year, "\\for a4\\just figs"))
dir.create(paste0("analysis outputs\\", ts.lcm.names, "\\" , this.year, "\\for a4\\legends"))

for(i in 1:length(plots.for.a4fig.names)){
  ggsave(ggdraw(plots.for.a4fig[[i]]), file = paste0("analysis outputs\\", ts.lcm.names, "\\" , this.year, 
                                                     "\\for a4\\",  plots.for.a4fig.names[i], ".png"), width = 10, height = 6)
  
  # just fig
  ggsave(ggdraw(plots.for.a4fig[[i]]) , 
         file = paste0("analysis outputs\\", ts.lcm.names, "\\" , this.year, "\\for a4\\just figs\\",
                       plots.for.a4fig.names[i], ".svg"), 
         width = w.a4.mm, height = h.a4.mm, units = "mm")
  # just legend
  ggsave(legends.for.a4fig[[i]] , 
         file = paste0("analysis outputs\\", ts.lcm.names, "\\" , this.year, "\\for a4\\legends\\",
                       plots.for.a4fig.names[i], ".svg"), 
         width = w.a4.mm, height = h.a4.mm, units = "mm")
}

# to align space for legend /plot in all
plots.for.slides = align_plots(p01.lcm,
                               p02.lcm_patch,
                               p03.lcm_patcharea,
                               p04.01.lcm_patch_edge,
                               p04.02.lcm_patch_aw,
                               p05.lcm_patch_edge_awi,
                               p06.lcm_patchHab,
                               p07.lcm_patchMeanHab,
                               p08.lcm_patchQualArea,
                               p09.cost,
                               p10.cost.patch,
                               p11.cost.patch.accumCost,
                               p12.cost.patch.pconnect,
                               p13.patchQualArea.pconnect,
                               p14.contribution,
                               p15.eca,
                               align="hv", axis="tblr")
plot.names = c("p01.lcm"
               ,"p02.lcm_patch"
               ,"p03.lcm_patcharea"
               ,"p04.01.lcm_patch_edge"
               ,"p04.02.lcm_patch_aw"
               ,"p05.lcm_patch_edge_awi"
               ,"p06.lcm_patchHab"
               ,"p07.lcm_patchMeanHab"
               ,"p08.lcm_patchQualArea"
               ,"p09.cost"
               ,"p10.cost.patch"
               ,"p11.cost.patch.accumCost"
               ,"p12.cost.patch.pconnect"
               ,"p13.patchQualArea.pconnect"
               ,"p14.contribution"
               ,"p15.eca"
)

for(i in 1:length(plot.names)){
  ggsave(ggdraw(all.plots[[i]]), file = paste0("analysis outputs\\", ts.lcm.names, "\\" , 
                                               this.year, "\\",  plot.names[i], ".png"), 
         width = 14, height = 9)
}
