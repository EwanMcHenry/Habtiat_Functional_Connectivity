# Ewan McHenry
# functional connectivity metric dev
# script 02 - loading data

# data ----
# hab cost and edge etc
countries = st_read(paste0(gis.wd, "\\Data\\administrative boundaries\\Countries\\R.5countries.simp100m.shp"))%>% st_transform( 27700) %>% arrange(name)

# lcm19.rast = raster(paste0(gis.wd, "CEH LCM\\LCM 2019\\lcm2019_v1.01.tif"))
# full lcms
lcm19.rast25.gb = raster(paste0(gis.wd, "\\Data\\LCM\\LCM2019\\25m land parcel\\gb2019lcm25m.tif"))
lcm90.rast25.gb = raster(paste0(gis.wd, "\\Data\\LCM\\LCM1990\\25m land parcel\\gb1990lcm25m.tif"))

lcm19.rast25.ni = raster(paste0(gis.wd, "\\Data\\LCM\\LCM2019\\25m land parcel\\ni2019lcm25m.tif"))
lcm90.rast25.ni = raster(paste0(gis.wd, "\\Data\\LCM\\LCM1990\\25m land parcel\\ni1990lcm25m.tif"))
lcm90.rast25.ni[lcm90.rast25.ni == 0] = 13 # # strange 0 habitat type in NI 1990s LCM, this makes it sea

nwss = st_read(paste0(gis.wd, "\\Data\\NWSS\\Native_Woodland_Survey_of_Scotland.shp")) %>% st_transform( 27700)

awi = st_read(paste0(gis.wd, "\\Data\\ancient woodland\\original\\AWI_joined_v4.02.shp")) %>% st_transform( 27700)

### save ----
save(countries,
     lcm19.rast25.gb,
     lcm90.rast25.gb,
     lcm19.rast25.ni,
     lcm90.rast25.ni,
     lcm90.rast25.ni,
     nwss,
     awi, 
     file = paste0(func.conect.path, 
                   "\\analysis outputs\\r_global_data_.RData")
)

print("script02 done")
