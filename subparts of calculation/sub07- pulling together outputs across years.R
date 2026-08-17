##------ Wed Jul 20 12:32:04 2022 ------##
# Ewan McHenry

multi_year_landscape_summary <- list()
multi_year_ts.hexgrid <- list()
multi_year_candidate_pairs <- list()
multi_year_source_patch_centroid_info <- list()



y = 1
for (y in seq_along(years.considered)) {
  this.year <- years.considered[y]
  
  load(paste0(func.conect.path, 
              "\\analysis outputs\\", this.tss[this.ts.num], "\\", this.year, 
              "\\lcmres", constants$lcm.res, 
              "\\06_r_funcconnect_EffectiveAreas_ECAobs_.RData"))
  
  
  year_name <- as.character(this.year)
  
  multi_year_landscape_summary[[year_name]] <- 
    landscape_stats[[this.tss[this.ts.num]]][[this.year]]$landscape.summary
  
  multi_year_ts.hexgrid[[year_name]] <- 
    landscape_stats[[this.tss[this.ts.num]]][[this.year]]$ts.hexgrid
  
  multi_year_source_patch_centroid_info[[year_name]] <- 
    landscape_stats[[this.tss[this.ts.num]]][[this.year]]$source_patch_centroid_info
  
  multi_year_candidate_pairs[[year_name]] <- 
    landscape_stats[[this.tss[this.ts.num]]][[this.year]]$candidate_pairs
  
  model_params <- 
    landscape_stats[[this.tss[this.ts.num]]][[this.year]]$model_params
}

multi_year_landscape_summary <- do.call(rbind, multi_year_landscape_summary)

multi_year_results <- list(
  landscape_summary = multi_year_landscape_summary,
  ts.hexgrid = multi_year_ts.hexgrid,
  candidate_pairs = multi_year_candidate_pairs,
  source_patch_centroid_info = multi_year_source_patch_centroid_info,
  model_params = model_params
)

save(
  multi_year_results,
  file = paste0(
    func.conect.path,
    "\\analysis outputs\\",
    this.tss[this.ts.num],
    "\\07_multi_year_landscape_output_lcmres",
    constants$lcm.res,
    "_years_",
    paste(years.considered, collapse = "-"),
    ".RData"
  )
)


print(paste(this.tss[this.ts.num], "Pulling across years (script07) done"))
