# pulling over landscapes

all_landscape_summary <- data.frame()

all_ts.hexgrid <- list()
all_candidate_pairs <- list()
all_source_patch_centroid_info <- list()
all_model_params <- list()


this.ts.num <- 1

for (this.ts.num in 1:length(pulling_landscapes)) {
  ts <- pulling_landscapes[this.ts.num]
  load(
    paste0(
      func.conect.path,
      "\\analysis outputs\\",
      ts,
      "\\07_multi_year_landscape_output_lcmres",
      constants$lcm.res,
      "_years_",
      paste(years.considered, collapse = "-"),
      ".RData"
    )
  )

  all_landscape_summary <-
    rbind(
      all_landscape_summary,
      multi_year_results$landscape_summary
    )
  all_ts.hexgrid[[ts]] <-
    multi_year_results$ts.hexgrid

  all_candidate_pairs[[ts]] <-
    multi_year_results$candidate_pairs

  all_source_patch_centroid_info[[ts]] <-
    multi_year_results$source_patch_centroid_info

  all_model_params[[ts]] <-
    multi_year_results$model_params
}

out_dir <- paste0(
  func.conect.path,
  "\\analysis outputs\\",
  "all_treescapes_all_years_lcmres",
  constants$lcm.res
)

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

save(
  all_landscape_summary,
  all_ts.hexgrid,
  all_model_params,
  file = paste0(
    func.conect.path,
    "\\analysis outputs\\",
    "all_treescapes_all_years_lcmres", constants$lcm.res,
    "\\", paste(years.considered, collapse = "-"), "_scapes_",
    paste(substr(pulling_landscapes, 1, 2), collapse = "-"),
    ".RData"
  )
)

save(
  all_candidate_pairs,
  all_source_patch_centroid_info,
  file = paste0(
    func.conect.path,
    "\\analysis outputs\\",
    "all_treescapes_all_years_lcmres", constants$lcm.res,
    "\\",
    paste(years.considered, collapse = "-"), "_scapes_",
    paste(substr(pulling_landscapes, 1, 2), collapse = "-"),
    "_patch_pairs_centroids.RData"
  )
)


rm(multi_year_results,
   all_source_patch_centroid_info,
   all_candidate_pairs,
   all_ts.hexgrid
   )