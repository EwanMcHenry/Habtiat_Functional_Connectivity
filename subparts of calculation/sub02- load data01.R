# Ewan McHenry
# functional connectivity metric dev
# script 02 - loading data

# data ----

irish_isles <- gisco_get_countries(
  country = c("UK", "IE", "IM"),
  resolution = "20"
) %>%
  st_transform(27700) %>%
  summarise()
#plot(st_geometry(irish_isles))

countries <-gisco_get_nuts(
  country = "UK",
  nuts_level = 1,
  resolution = "01"
) %>%
  st_transform(27700) %>%
  mutate(
    nation = case_when(
      grepl("SCOTLAND", NAME_LATN) ~ "Scotland",
      grepl("WALES", NAME_LATN) ~ "Wales",
      grepl("NORTHERN IRELAND", NAME_LATN) ~ "Northern Ireland",
      grepl("ENGLAND", NAME_LATN) ~ "England",
      grepl("LONDON", NAME_LATN) ~ "England",
      grepl("YORKSHIRE AND THE HUMBER", NAME_LATN) ~ "England",
      TRUE ~ NAME_LATN
    ))%>%
  group_by(nation) %>%
  summarise(geometry = st_union(geometry), .groups = "drop")
#plot(st_geometry(countries))


# full lcms ----

# woodland info -----
nwss = st_read(nwss.dir) %>% st_transform( 27700)

awi = st_read(awi.dir) %>% st_transform( 27700)

# big roads data ----

# roads_uk_major <- oe_get(
#   place = "United Kingdom",
#   layer = "lines",
#   query = "
#   SELECT highway, geometry
#   FROM lines
#   WHERE highway IN (
#     'motorway',
#     'trunk',
#     'primary'
#   )
#   "
# )
# st_write(
#   roads_uk_major,
#   roads.dir,
#   delete_dsn = TRUE
# )

#source(paste0(func.conect.path, "\\code\\subparts of calculation\\whats up with roads.R"))

roads_uk_major <- st_read(roads.dir) %>% st_transform(27700) %>% 
  mutate(buffer_size = case_when(
    highway == "motorway" ~ 9.08,
    highway == "trunk" ~ 5.70,
    highway == "primary" ~ 6.82,
    TRUE ~ NA_real_
  )) 



### save ----
save(roads_uk_major,
     irish_isles,
     countries,
     lcm,
     nwss,
     awi, 
     file = paste0(func.conect.path, 
                   "\\analysis outputs\\r_global_data_.RData")
)

print("script02 done")
