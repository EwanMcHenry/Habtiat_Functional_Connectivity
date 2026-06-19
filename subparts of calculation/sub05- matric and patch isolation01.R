  # Ewan McHenry
  ##------ Fri Feb 25 09:19:50 2022 ------##
  # functional connectivity metric dev
  
  # script 05 - matrix and patch isolation
  
  time.now = Sys.time()
  print("cost distance calculation")
  
  # ## load curated data ----
  load(paste0(func.conect.path, "\\analysis outputs\\", this.tss[this.ts.num], "\\r_curated data_.RData"))
  ## LOAD PATCH DATA ----
  load(paste0(func.conect.path, "\\analysis outputs\\", this.tss[this.ts.num], "\\", this.year, "\\r_funcconnect_patchwork.RData"))
  
  r_stack <- terra::rast(
    paste0(func.conect.path, "\\analysis outputs\\", 
           this.tss[this.ts.num], "\\", this.year, "\\r_stack.tif"))
  
  
  
  # BUILD COST LAYER ----
  
  hab.cost.lcm = r_stack$lcm  
  values(hab.cost.lcm) = dispers.costs$scaled.ecolog.cost[values(r_stack$lcm)]
  #replace gaps with high-cost landscape - these are normally sea, but can be beyond edge of landscape, this shouldnt be a problem, becasue landscape is buffered
  values(hab.cost.lcm)[is.na(values(hab.cost.lcm))] = max(dispers.costs$scaled.ecolog.cost[values(r_stack$lcm)],na.rm=T)
  trouble_plot(hab.cost.lcm, "habitat cost layer")
  writeRaster(hab.cost.lcm, paste0(func.conect.path, "\\analysis outputs\\", this.tss[this.ts.num], "\\", this.year, "\\habitat_cost_layer.tif"), overwrite = T)
  
  # aggragate cost raster by mean to reduce resolution and computational power ----
  if(constants$agg_cost_surface == T){
    hab.cost.agg.rast = aggregate(hab.cost.lcm, constants$cost_agg_n, max, na.rm = T) # 4 x 4 mean aggregation
    trouble_plot(hab.cost.agg.rast, "habitat cost layer aggregated")
    writeRaster(hab.cost.agg.rast, paste0(func.conect.path, "\\analysis outputs\\", this.tss[this.ts.num], "\\", this.year, "\\habitat_cost_layer_aggregated.tif"), overwrite = T)
    cost_raster = hab.cost.agg.rast
    rm(hab.cost.agg.rast)
  } else{
    cost_raster = hab.cost.lcm}
  trouble_plot(cost_raster, "habitat cost layer for cost distance calculation")
  
  ## save LANDSCAPE cost INFO ----
  # Mean landscape cost
  # mean cost traveling through "not sea" in landscape
  landscape_stats[[this.ts.num]]$year_stats$landscape.mean.ecolog.cost.not.sea = mean(values(hab.cost.lcm)[!(values(r_stack$lcm  ) %in% coastal_water_lcm) & !is.na(values(r_stack$lcm  ))])*constants$cost.scale.factor
  # landscape.median.scaled.ecolog.cost.not.sea = median(values(hab.cost.lcm)[!(values(r_stack$lcm) %in%coastal_water_lcm) & !is.na(values(r_stack$lcm))])*constants$cost.scale.factor
  
  # mean hex cost
  # mean costs of hexes "not sea" in land scape
  not.sea.cost = hab.cost.lcm
  values(not.sea.cost)[(values(r_stack$lcm) %in% coastal_water_lcm) | is.na(values(r_stack$lcm))] = NA
  landscape_stats[[this.ts.num]]$year_stats$mean.ecolog.cost.not.sea = exact_extract(not.sea.cost, ts.hexgrid, "mean" )*constants$cost.scale.factor
  # ts.hexgrid$median.ecolog.cost.not.sea = exact_extract(not.sea.cost, ts.hexgrid, "median" )*constants$cost.scale.factor
  
  # CANDIDATE CONNECTIONS ----
  
  patches <- patch_centroid_info %>% select(row_id, uid, grid_id, focal_landscape )
  coords <- sf::st_coordinates(patches)
  
  ## build candidate pairs - global neighbourhood graph ----
  focal_idx <- which(patches$focal_landscape == 1)
  
  nn <- sf::st_is_within_distance(
    patches[focal_idx, ], # 
    patches,
    dist = constants$max.dispersal.considered
  )
  
  candidate_pairs <- data.frame(
    from = rep(focal_idx, lengths(nn)),
    to   = unlist(nn)
  ) %>%
    # dplyr::filter(from != to) %>%  # remove self connections
    dplyr::mutate(
      uid_from = patch_centroid_info$uid[from],
      uid_to   = patch_centroid_info$uid[to],
      
      grid_from = patch_centroid_info$grid_id[from],
      grid_to   = patch_centroid_info$grid_id[to]
    )
  
  ## add euclidean distance between patches to candidate pairs
  candidate_pairs$euclid_dist <-
    sqrt(
      (coords[candidate_pairs$from,1] -
         coords[candidate_pairs$to,1])^2 +
        (coords[candidate_pairs$from,2] -
           coords[candidate_pairs$to,2])^2
    )
  
  ## build edge list for cost distance calculation ----
  edge_list <- split(
    candidate_pairs$uid_to,
    candidate_pairs$uid_from
  )
  
  # hex_id <- 151888
  # i = 1
  # hex_layer = ts.hexgrid
  # cost_raster = cost_raster
  # hex_geom <- hex_layer[hex_layer$grid_id == hex_id, ] # hexid should only be those within focal landscape
  #   A <- "1011_151888_1"
  # B <- "999_151888_1"
  # e.g.patches <- patch_centroid_info[patch_centroid_info$uid %in% c(A, B), ]
  # trouble_plot(e.g.patches, "test_patches")

  
  
  compute_lcd_from_hex <- function(hex_id, # shoudl only be those in focal landscape
                                   hex_layer,
                                   patches,
                                   cost_raster,
                                   edge_list,
                                   constants) {
    
    
    cost_raster <- terra::unwrap(cost_raster)
    
  #  hex_id, # should only be those in focal landscape
    
    # ---------------------------
    # 1. hex geometry
    # ---------------------------
    hex_geom <- hex_layer[hex_layer$grid_id == hex_id, ] # hexid should only be those within focal landscape
    
    if (nrow(hex_geom) == 0) return(NULL)
    
    # ---------------------------
    # 2. buffer hex (movement window)
    # ---------------------------
    hex_buffer <- st_buffer(
      hex_geom,
      dist = constants$max.dispersal.considered,
      nQuadSegs = 3
    )
    trouble_plot(hex_buffer, paste0("hex ", hex_id, " buffered"))
    
    # ---------------------------
    # 3. crop cost raster once
    # ---------------------------
    cost_crop <- terra::crop(cost_raster, hex_buffer) # note this is not masked, so the "corners" might have an impact, allowing corridors etc, but they are deemed too far to have an impact 
    trouble_plot(cost_crop, paste0("cost raster cropped for hex ", hex_id))
    
    if(constants$cost.dist.method == "gdistance"){
      cost_raster_r <- raster::raster(cost_crop)
      # transition matrix and correct for diagonal/ortho -- note scale = F, to get real cost distance measure, rather than scaled
      gridbuff_Tcost <- transition(cost_raster_r,function(x) 1/mean(x),8)
      gridbuff_Tcost.C <- geoCorrection(gridbuff_Tcost, type="c", multpl=FALSE, scl=F)
      
    }
    
    # ---------------------------
    # 4. patches inside buffered area
    # ---------------------------
    patches_in_buffer <- patches[hex_buffer, ]  # spatial subset
    trouble_plot(patches_in_buffer, paste0("patches in buffer for hex ", hex_id))
    
    if (nrow(patches_in_buffer) == 0) return(NULL)
    
    # ---------------------------
    # 5. focal patches (sources)
    # ---------------------------
    source_patches <- patches_in_buffer[patches_in_buffer$grid_id == hex_id &
                                          patches_in_buffer$focal_landscape == 1
                                          , ]
    trouble_plot(source_patches, paste0("source patches for hex ", hex_id))
    source_ids <- source_patches$uid
    
    if (nrow(source_patches) == 0) return(NULL)
    
    # ---------------------------
    # 6. loop - run cost accumulation over source patches
    # ---------------------------
    
    if(constants$cost.dist.method == "terra"){
      lcd_terra_list <- lapply(seq_len(nrow(source_patches)), function(i) {
      
      from_patch <- source_patches[i, ]
      from_uid <- from_patch$uid
      
      # candidate targets from edge list
      targets_uids <- edge_list[[from_uid]]
      
      if (is.null(targets_uids) || length(targets_uids) == 0) {
        return(NULL)
      }
  
      # subset target geometries
      to_patch <- patches_in_buffer[patches_in_buffer$uid %in% targets_uids, ]
      trouble_plot(to_patch, paste0("target patches for hex ", hex_id, " from patch ", from_uid))
      # convert to terra vectors
      from_vect <- terra::vect(from_patch)
      to_vect   <- terra::vect(to_patch)
      
      # make cost_matrix with target_id for focal patch cell (target)
      target_id <- 0
      if(min(values(cost_crop))==target_id){print("WARNING: multiple targets in cost matrix, consider changing target_id")}
      cost_matrix <- cost_crop
      from_cell_id <- terra::cellFromXY(cost_matrix, terra::crds(from_vect))
      cost_matrix[from_cell_id] <- target_id
  
      # compute accumulated cost distance matrix
      lcd_mat <- terra::costDist(
        cost_matrix, target = target_id
      )
      trouble_plot(lcd_mat, paste0("cost distance matrix for hex ", hex_id, " from patch ", from_uid))

      
      # extract lcd values for target patches
      lcd <- terra::extract(lcd_mat, to_vect)[, 2]
      trouble_plot(to_vect, paste0("target patches for hex ", hex_id, " from patch ", from_uid))
      # # add cost to pairs info and output
      data.frame(
        from = from_uid,
        to   = to_patch$uid,
        lcd  = as.numeric(lcd)
      )
    })
      lcd_list <- lcd_terra_list}
    
    if(constants$cost.dist.method == "gdistance"){
      lcd_gdist_list <- lapply(seq_len(nrow(source_patches)), function(i) {
      
      from_patch <- source_patches[i, ]
      from_uid <- from_patch$uid
      
      # candidate targets from edge list
      targets_uids <- edge_list[[from_uid]]
      
      if (is.null(targets_uids) || length(targets_uids) == 0) {
        return(NULL)
      }
      
      # subset target geometries
      to_patch <- patches_in_buffer[patches_in_buffer$uid %in% targets_uids, ]
      
      # convert to terra vectors
      from_vect <- terra::vect(from_patch)
      to_vect   <- terra::vect(to_patch)
  
      # compute accumulated cost distance matrix using gdistance
      lcd_mat_gdist <- gdistance::costDistance(
        gridbuff_Tcost.C,
        fromCoords = terra::crds(from_vect),
        toCoords   = terra::crds(to_vect)
      )
      
      # format output
      data.frame(
        from = from_uid,
        to   = to_patch$uid,
        lcd  = as.numeric(lcd_mat_gdist)
      )
    })
      lcd_list <- lcd_gdist_list}
    
    lcd_list
  }
  
  # run lcd engine ----
  hex_ids <- unique(patches$grid_id[patches$focal_landscape == 1])#[20]
  cost_raster_safe <- terra::wrap(cost_raster)
  
  handlers("txtprogressbar")
  
  with_progress({
    
    p <- progressor(along = hex_ids)
    
    lcd_results <- if (constants$run_parallel) {
      
      plan(multisession, workers = constants$cores_for_parallel)
      
      future_lapply(
        hex_ids,
        function(h) {
          
          p()
          
          compute_lcd_from_hex(
            h,
            hex_layer   = tsbuff.hexgrid,
            patches     = patches,
            cost_raster = cost_raster_safe,
            edge_list   = edge_list,
            constants   = constants
          )
        },
        future.seed = TRUE,
        future.scheduling = 1
      )
      
    } else {
      
      lapply(
        hex_ids,
        function(h) {
          
          p()
          
          compute_lcd_from_hex(
            h,
            hex_layer   = tsbuff.hexgrid,
            patches     = patches,
            cost_raster = cost_raster_safe,
            edge_list   = edge_list,
            constants   = constants
          )
        }
      )
    }
  })  
  
  # terra runs x6-7 times faster than gdist
  # attach lcds to candidate pairs ----
  
  lcd_edges <- lcd_results %>%
    unlist(recursive = FALSE) %>%  # flatten nested lists
    dplyr::bind_rows() %>% 
    # SCALE lcd to effective distance
    dplyr::mutate(
      lcd = lcd * constants$cost.scale.factor
    )
  
  candidate_pairs <- candidate_pairs %>%
    dplyr::left_join(
      lcd_edges,
      by = c(
        "uid_from" = "from",
        "uid_to"   = "to"
      )
    )
  

  print(Sys.time() -time.now)
  # rescale cost distance to actual effective distance (before costs scaled to be close to 1 for computational efficiency)
  print("cost distance done")
  
  
  # SAVE PATCH DATA ----
  save(candidate_pairs,
     #  big.cost.dist,
     landscape_stats,
       patch_centroid_info, 
       file = 
       paste0(func.conect.path, 
              "\\analysis outputs\\", this.ts.for.loop[this.ts.num], "\\", this.year, "\\r_funcconnect_MatrixCostDists.RData")
  )
  
  rm(list = c("hab.cost.lcm", "cost_raster", "patches", "patch_centroid_info", "candidate_pairs", "lcd_edges"))
  
  
  
  if(grepl("Illustrative", this.ts.for.loop[this.ts.num]) ){
    save(gridbuff_Tcost.C,
         file = 
           paste0(func.conect.path, 
                  "\\analysis outputs\\", 
                       this.ts.for.loop[this.ts.num], "\\", this.year, "\\cost.objects.RData")
    )
    
  }
  
  rm()
  gc()
       
  print("Landscape matrix and patch linkage (script05) done")
  
  
  
  
  
  #   EXPLORATION ----
  # ###exploring different cost-distance calculation methods ----
  # # plot an example points coloured by cost distance from a focal patch
  # e.g.focal_patch_uid <- "803_149085_1"
  # e.g.lcd_terra <- lcd_terra_list[[1]]
  # e.g.cost_terra <- e.g.lcd_terra[e.g.lcd_terra$from == e.g.focal_patch_uid, ]
  # e.g.lcd_gdist <- lcd_gdist_list[[1]]
  # e.g.cost_gdist <- e.g.lcd_gdist[e.g.lcd_gdist$from == e.g.focal_patch_uid, ]
  # patches_costs <- patches_in_buffer %>%
  #   dplyr::mutate(
  #     type = ifelse(uid == source_patches$uid[1], "source", "target"),
  #     lcd_terra = e.g.cost_terra$lcd[match(uid, e.g.cost_terra$to)],
  #     lcd_gdist = e.g.cost_gdist$lcd[match(uid, e.g.cost_gdist$to)]
  #   ) %>% 
  #   filter(!is.na(lcd_terra)| !is.na(lcd_gdist))
  # patches_costs$cost_ratio <- patches_costs$lcd_terra/ patches_costs$lcd_gdist
  # 
  # # # plot accum cost surface (cost_matrix) and points (patches_costs[,7])
  # ggplot() +
  #   geom_spatraster(data = lcd_mat) +
  #   scale_fill_viridis_c(direction = -1) +
  #   geom_sf(data = patches_costs, aes(color = lcd_gdist), size = 1) +
  #   #plasma viridis
  #   scale_color_viridis_c(option ="magma" ) +
  #   labs(color = "gDist") +
  #   theme_minimal()
  # 
  # ggplot() +
  #   geom_spatraster(data = lcd_mat) +
  #   scale_fill_viridis_c(direction = -1) +
  #   geom_sf(data = patches_costs, aes(color = lcd_terra), size = 1) +
  #   #plasma viridis
  #   scale_color_viridis_c(option ="magma" ) +
  #   labs(color = "terra") +
  #   theme_minimal()
  # 
  # 
  # plot(patches_costs$lcd_terra, patches_costs$cost_ratio)
  # 
  # patches_costs[patches_costs$lcd_terra<1000,] %>% arrange(lcd_terra)
  # patches_costs$is.low_ratio <- patches_costs$cost_ratio < 0.9 | patches_costs$cost_ratio > 1.1
  # 
  # ggplot() +
  #   geom_spatraster(data = lcd_mat) +
  #   scale_color_viridis_c(option ="magma" ) +
  #   geom_sf(data = patches_costs, aes(color = is.low_ratio), size = 1) +
  #   #plasma viridis
  #   scale_color_viridis_c(option ="magma" ) +
  #   theme_minimal()
  # 
  # 
  
#   LESSON - WHEN COST LAYER NOT AGGREGATED THEY ARE THE SAME
  # median looks pretty close and max looks even closer, which kinda makes sense! when agging 25 over 100 i.e. 4, and youve gotta cross over at least one of them
  # A <- "905_149435_1"
  # B <- "946_150486_1"
  # A <- "525_148035_1"
  # B <- "629_149787_1"
  # A <- "678_152240_1"
  # B <- "999_151888_1"
  # 
  # candidate_pairs %>%
  #   filter(
  #     (uid_from == A & uid_to == B) |
  #       (uid_from == B & uid_to == A)
  #   ) %>%
  #   select(uid_from, uid_to, grid_from, grid_to, lcd)
  # #