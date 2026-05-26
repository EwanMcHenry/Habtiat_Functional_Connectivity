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

### save ----
save(irish_isles,
     countries,
     lcm,
     nwss,
     awi, 
     file = paste0(func.conect.path, 
                   "\\analysis outputs\\r_global_data_.RData")
)

print("script02 done")
