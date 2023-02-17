{
  library(dataRetrieval); library(nhdplusTools);
  library(raster); library(sf); library(fasterize);
  library(resample); library(wrfhydroSubsetter); library(AOI);
  library(tidyverse); library(dplyr);
}


#library(devtools)
#install_github("mikejohnson51/HydroData")
library(HydroData)
#detach("package:HydroData", unload=TRUE); remove.packages("HydroData");


locs = data.frame(comids = c(191739, 1631587, 5894384, 5781369, 19389766, 23762661),
                  siteID = c('06709000', '08173000', '01616500', '08159000' ,'03118500' , '14190500'))

#test_loc = data.frame(comids = c(5894384), siteID = c('01616500')) # In case using a single watershed

namelist = c()
for (i in 1:nrow(locs)){
  state = st_transform(AOI::aoi_get(state = "all", county = "all"), 4269)[st_transform(findNLDI(comid = locs$comids[i]), crs = 4269),]
  name = gsub(" ", "-", tolower(paste(state$name, state$state_name, locs$comids[i], locs$siteID[i], sep = "_")))
  namelist = c(namelist, name)
  rm(name)
}

a = findNWIS(site_no = locs$siteID[3])

a = findNLDI(comid = locs$comids[3])

for (i in 1:nrow(locs)){
  basin = findNLDI(comid = locs$comids[i], find = c("basin"))$basin  #%>%  as('Spatial')
  x <- nrow(basin)
  d <- data.frame(row.names = 1:x, id = 1:x)
  spd <- SpatialPolygonsDataFrame(basin %>% as('Spatial'), data = d)
  
  huc12 = get_huc12(AOI = basin) %>% 
    st_intersection(basin) %>% as('Spatial')
  
  flowline = get_nhdplus (AOI = findNLDI(comid = locs$comids[i], find = c("basin"))$basin) %>% 
    filter(streamorde >= 3) %>% 
    as('Spatial')
  
  site = readNWISsite(locs$siteID[i])
  site = sf::st_as_sf(site, coords = c("dec_long_va", "dec_lat_va"), crs=site$coord_datum_cd) %>% 
    as('Spatial')
  
  name = strsplit(namelist[i], '_')[[1]][1]
  dsn = paste0("./ShpFiles/", name)
  dir.create("./ShpFiles/")
  
  rgdal::writeOGR(spd, dsn = paste0(dsn, '_basin.shp'), layer = 'basin', driver = "ESRI Shapefile")
  rgdal::writeOGR(huc12, dsn = paste0(dsn, '_huc12.shp'), layer = 'huc12', driver = "ESRI Shapefile")
  rgdal::writeOGR(flowline, dsn = paste0(dsn, '_flowline.shp'), layer = 'flowline', driver = "ESRI Shapefile")
  rgdal::writeOGR(site, dsn = paste0(dsn, '_site.shp'), layer = 'site', driver = "ESRI Shapefile")
}


basin = findNLDI(comid = locs$comids[3], find = c("basin"))$basin %>% as('Spatial') #%>% fortify();
huc12 = get_huc12(AOI = basin) %>% st_intersection(basin) %>% as('Spatial')
#flowline = findNLDI(comid = locs$comids[3], find = c("flowline"))#$flowlines
flowline = get_nhdplus (AOI = findNLDI(comid = locs$comids[3], find = c("basin"))$basin) %>% 
  filter(streamorde >= 3) %>% as('Spatial') #%>% fortify();
#site = whatNWISsites(sites =locs$siteID[3]) %>% as('Spatial') %>% fortify();
site = readNWISsite(locs$siteID[3])
site = sf::st_as_sf(site, coords = c("dec_long_va", "dec_lat_va"), crs=site$coord_datum_cd) %>% as('Spatial') #%>% fortify();
#%>% sf::as_Spatial()# %>% fortify()
  #as('Spatial') %>% fortify();

dsn

x <- nrow(basin)
d <- data.frame(row.names = 1:x, id = 1:x)
spd <- SpatialPolygonsDataFrame(basin,data = d)

rgdal::writeOGR(spd, dsn = 'basin.shp', layer = 'basin', driver = "ESRI Shapefile")
rgdal::writeOGR(huc12, dsn = 'huc12.shp', layer = 'huc12', driver = "ESRI Shapefile")
rgdal::writeOGR(flowline, dsn = 'flowline.shp', layer = 'flowline', driver = "ESRI Shapefile")
rgdal::writeOGR(site, dsn = 'site.shp', layer = 'site', driver = "ESRI Shapefile")

#sp = sf::st_as_sf(x = df,  coords = c("lon", "lat"), crs = as.character(AOI::aoiProj)) %>% sf::as_Spatial()


g <- ggplot()+
  geom_path(data=basin, aes(long, lat, group=group), color = 'black')+
  geom_path(data=flowline, aes(long, lat, group=group), color = 'blue')+
  geom_point(data=site, aes(dec_long_va, dec_lat_va))+
  theme(
    plot.background = element_blank())
g
