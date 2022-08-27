{
  library(dataRetrieval)
  library(raster)
  library(sf)
  library(resample)
  library(wrfhydroSubsetter)
  library(tidyverse)
  library(AOI)
}
# This script is prone to crash during debugging.
.rs.restartR()


locs = data.frame(comids = c(23762661, 1631587, 5894384, 19389766, 191739, 5781369),
                  siteID = c('14190500', '14190500', '01616500', '03118500' ,'6709000' , '08159000'))


# trouble shooting subset_wrf_hydro_domain first
test_loc = data.frame(comids = c(5894384),
                      siteID = c('01616500'))

outDir = "/mnt/d/subsetDOMAINS/"
FULLDOMAINDIR = "/mnt/d/FULLDOMAIN/nwmCONUS-v216/"


## In process of package development.
#library(devtools)
#detach("package:wrfhydroSubsetter", unload=TRUE)
#remove.packages("wrfhydroSubsetter")
#devtools::install("/mnt/d/GitHub/wrfhydroSubsetter", repos = NULL, type="source")
#library(wrfhydroSubsetter)


# Sourcing modified version of subsetter. Because it is still in development phase.
source("/mnt/d/GitHub/wrfhydroSubsetter/steering_file/subsetter_mod.R")
source("/mnt/d/GitHub/wrfhydroSubsetter/steering_file/resampleDataMod.R")
source("/mnt/d/GitHub/wrfhydroSubsetter/steering_file/subset_sequence.R")


# Current version of NWM is likely to use NLCD 2016 LC, according to (https://www.weather.gov/media/notification/pdf2/scn20-119nwm_v2_1aad.pdf)
nlcdDir <- "/mnt/d/nlcd_2016_land_cover_l48_20210604/nlcd_2016_land_cover_l48_20210604.img"
nlcdObj <- raster(nlcdDir)
methodList <- c("nn", "maj","rawarea","area")


# Running the following wrapper do tasks as follow:
# 1. Subset NWM domain files based on  given COMIDs & siteIDs
# 2. If nlcdDir and methodList given, it resamples nlcd into nwm LULC using selected methods
for(i in 1:nrow(locs)){
  subset_sequence(locs$comids[i], locs$siteID[i], FULLDOMAINDIR, outDir, nlcdDir, methodList)
}


# Watershed Shuffle wrapper
for(i in 1:nrow(locs)){
  state = st_transform(AOI::aoi_get(state = "all", county = "all"), 4269)[st_transform(findNLDI(comid = locs$comids[i]), crs = 4269),]
  name = gsub(" ", "-", tolower(paste(state$name, state$state_name, locs$comids[i], locs$siteID[i], sep = "_")))
  basin = findNLDI(comid = locs$comids[i], find = c("basin"))$basin
  
  area_resample_Dir = list.files(paste0('/mnt/d/subsetDOMAINS/', name), "geogrid_area_resample", full = TRUE)
  area_resample = wrfhydroSubsetter::make_empty_geogrid_raster(area_resample_Dir, "LU_INDEX")
  
  basin_coord = st_transform(basin, st_crs(area_resample))
  
  # Regular raster::mask function only selects raster cells those have their centroid falls into polygon boundary
  # Using rasterize to work around
  # https://stackoverflow.com/questions/68781519/errors-when-clipping-raster-stacks-with-sfst-crop-and-rastercrop
  # https://gis.stackexchange.com/questions/255025/r-raster-masking-a-raster-by-polygon-also-remove-cells-partially-covered
  basin_coord_ras <- rasterize(basin_coord, area_resample, getCover=TRUE)
  basin_coord_ras[basin_coord_ras==0] = NA
  area_resample_mask_df <- stack(raster::mask(area_resample, basin_coord_ras)) %>% as.data.frame(xy = TRUE)

  # now shuffle LULC inside WS
  shuffled = transform(na.omit(area_resample_mask_df), layer = sample(layer))
  area_resample_shuffle = area_resample
  values(area_resample_shuffle)[as.numeric(rownames(shuffled))] <- shuffled$layer
  
  out_file = paste0('/mnt/d/subsetDOMAINS/', name, "/geogrid_shuffle.nc")
  file.copy(area_resample_Dir, out_file, overwrite = TRUE)
  
  nc = ncdf4::nc_open(out_file, write = TRUE)
  ncdf4::ncvar_put(nc, "LU_INDEX",
                   vals = as.vector(t(apply(as.matrix(area_resample_shuffle), 2, rev))))
  ncdf4::nc_close(nc)
  print(paste0(name, ": Done"))
}



# Watershed Shuffle
state = st_transform(AOI::aoi_get(state = "all", county = "all"), 4269)[st_transform(findNLDI(comid = locs$comids[4]), crs = 4269),]
name = gsub(" ", "-", tolower(paste(state$name, state$state_name, locs$comids[4], locs$siteID[4], sep = "_")))
basin = findNLDI(comid = locs$comids[4], find = c("basin"))$basin

area_resample_Dir = list.files('/mnt/d/subsetDOMAINS/polk_oregon_23762661_14190500', "geogrid_area_resample", full = TRUE)
area_resample = wrfhydroSubsetter::make_empty_geogrid_raster(area_resample_Dir, "LU_INDEX")

basin = findNLDI(comid = locs$comids[1], find = c("basin"))$basin
crs.here <- paste0(st_crs(area_resample)[1])
basin_coord = st_transform(basin, crs.here) %>% st_as_sf()
basin_coord2 = st_transform(basin, st_crs(area_resample))%>% st_as_sf()

# Regular raster::mask function only selects raster cells those have their centroid falls into polygon boundary
# Using rasterize to work around
# https://stackoverflow.com/questions/68781519/errors-when-clipping-raster-stacks-with-sfst-crop-and-rastercrop
# https://gis.stackexchange.com/questions/255025/r-raster-masking-a-raster-by-polygon-also-remove-cells-partially-covered
basin_coord_ras <- rasterize(basin_coord, area_resample, getCover=TRUE)
basin_coord_ras[basin_coord_ras==0] = NA

area_resample_mask_df <- stack(raster::mask(area_resample, basin_coord_ras)) %>% as.data.frame(xy = TRUE)
colnames(area_resample_mask_df) =c('x', 'y', 'layer')

# Simply checking out.
# https://datacarpentry.org/r-raster-vector-geospatial/02-raster-plot/
ggplot() +
  geom_raster(data = area_resample_mask_df, aes(x = x, y = y, fill = as.factor(layer)))+
  geom_sf(data = basin_coord, colour = "black", fill = "NA")+
  scale_fill_manual(name = "grp",values = myColors) +
  coord_sf()

# now shuffle LULC inside WS
shuffled = transform(na.omit(area_resample_mask_df), layer = sample(layer))
#area_resample_mask_df[as.numeric(rownames(non_na_shuffled)),] <- non_na_shuffled$layer

area_resample_shuffle = area_resample
values(area_resample_shuffle)[as.numeric(rownames(shuffled))] <- shuffled$layer

# Checking out if it is shuffled nicely.
{
  s = stack()
  s = addLayer(area_resample, area_resample_shuffle)
  names(s) = c("OG", "Shuffle")
  
  test = as(basin_coord, 'Spatial') %>% fortify()
  
  # https://stackoverflow.com/questions/68865682/overlay-polygon-layer-on-a-raster-stack-qplot
  rasterVis::gplot(s) +
    geom_tile(aes(fill = as.factor(value))) +
    geom_path(data=test, aes(long, lat, group=group), color = 'black')+
    facet_wrap(~ variable, ncol = 2) +
    scale_fill_manual(name = "grp",values = myColors) +
    coord_equal()
  #theme(legend.position = "none")
}


# subset sequence that includes LU resampling.
#subset_sequence(test_loc$comids[1], test_loc$siteID[1], FULLDOMAINDIR, outDir, nlcdDir, methodList) 

{### test========================================================================
# Defining the area
state = st_transform(AOI::aoi_get(state = "all", county = "all"), 4269)[st_transform(findNLDI(comid = test_loc$comids[1]), crs = 4269),]
name = gsub(" ", "-", tolower(paste(state$name, state$state_name, test_loc$comids[1], test_loc$siteID[1], sep = "_")))
basin = findNLDI(comid = test_loc$comids[1], find = c("basin"))$basin

# Subsetting NetCDF DOMAIN files
subset_files = paste0(outDir, name , "/")

# This data look-up table is most problematic at current point 08/10/2022
data = readxl::read_excel('/mnt/d/GitHub/wrfhydroSubsetter/nwm_land_use_mappings.xlsx') %>%
  dplyr::select(nlcd = Class, 
                description = Description,
                nwm = NWM)
geo  = list.files(subset_files, "geo_em.nc", full.names = TRUE)

nlcdDir <- "/mnt/d/nlcd_2016_land_cover_l48_20210604/nlcd_2016_land_cover_l48_20210604.img"
nlcdObj <- raster(nlcdDir)
#nlcdObj2 <- projectRaster(nlcdObj, crs=crs.here)

output_geo = spex::qm_rasterToPolygons(wrfhydroSubsetter::make_empty_geogrid_raster(geo))
output_geo_raster = wrfhydroSubsetter::make_empty_geogrid_raster(geo) # 43x28
#o1 =  output_geo %>%
#  AOI::bbox_get()
o2 = output_geo %>%
  AOI::bbox_get() %>%
  sf::st_transform(5070) %>%
  sf::st_buffer(40)

#crs.here <- paste0(st_crs(output_geo_raster))
#plot(o1)
#plot(o2)
#
#o3 = o2%>% 
#  sf::st_transform(st_crs(output_geo))
#
## So, 


#o0_to_grid <- output_grid(output_geo, 1000) # 46x31
#o1_to_grid <- output_grid(o1, 1000) # 46x31
crs.here <- paste0(st_crs(output_geo)[1])
#crs.here2 <- paste0(st_crs(output_geo)[1])
#
#o2_to_grid <- output_grid(o2, 1000) # 46x31
#o3_to_grid <- output_grid(o3, 1000) # 43x28
#
#test = o2 %>% 
#  sf::st_transform(sf::st_crs(output_geo))

nlcdCrop = raster::crop(nlcdObj, o2)
#nlcdCrop2 = projectRaster(nlcdCrop, crs.here)

nlcd_nwm = raster::reclassify(nlcdCrop, rcl = dplyr::select(data, from = nlcd, to = nwm))

#https://stackoverflow.com/questions/47415451/change-raster-extent
#
#nlcd_nwm2 = nlcd_nwm
#bb <- extent(output_geo_raster)
#extent(nlcd_nwm2) <- bb
#crs(nlcd_nwm2) <- paste0(st_crs(output_geo_raster))[1]
#

#
#
#nlcd_nwm2 = projectRaster(nlcd_nwm, crs = crs.here)
#nlcd_nwm3 = raster::crop(nlcd_nwm2, o1)

#crs.here <- paste0(st_crs(output_geo_raster)[1])
#pts = st_sfc(st_point(c(lng_max,lat_max)), crs = 4326) %>%
#  st_transform(crs.here) %>%
#  st_coordinates()


plot(nlcd_nwm)
plot(nlcd_nwm2)


nlcd_nwm2 = nlcd_nwm
bb <- extent(output_geo)
extent(nlcd_nwm2) <- bb
crs(nlcd_nwm2) <- paste0(st_crs(output_geo))[1]

nlcd_nwm3 <- projectRaster(nlcd_nwm, crs = 4326)
nlcd_nwm3 <- nlcd_nwm
crs(nlcd_nwm3) <- 4326
nlcd_nwm3 <- projectRaster(nlcd_nwm3, crs = 4326)

testout = call_gdal_mod(nlcd_nwm2, "mode", cellsize=1000)
plot(testout)

#x1 = as_tibble(nlcd_nwm)
#x2 = as_tibble(nlcd_nwm2)

methodList <- c("area", "rawarea", "maj", "nn")
methodList <- c("area", "maj")
methodList <- c("area")
source("/mnt/d/GitHub/wrfhydroSubsetter/steering_file/resampleDataMod.R")

for(j in 1:length(methodList)){
  #new_lu = resample::resampleData(input = nlcd_nwm, output = output_geo, method = methodList[j])
  new_lu = resampleDataMod(input = nlcd_nwm, output_geo = output_geo, method = methodList[j])
  out_file = paste0(subset_files, "geogrid_", methodList[j], "_resample.nc")
  file.copy(geo, out_file, overwrite = TRUE)
  nc = ncdf4::nc_open(out_file, write = TRUE)
  ncdf4::ncvar_put(nc, "LU_INDEX",
                   vals = as.vector(t(apply(as.matrix(new_lu), 2, rev))))
  ncdf4::nc_close(nc)
}

### testing croping new_lu with basin polygon
library(sf)
library(raster)
#library(terrainr)

#coordProj = geo_grid_proj(geo)
basin_coord = st_transform(basin, crs.here) %>% st_as_sf()
new_lu_df <- as.data.frame(new_lu, xy = TRUE) 

library(ggplot2)

# https://datacarpentry.org/r-raster-vector-geospatial/02-raster-plot/
# https://stackoverflow.com/questions/68781519/errors-when-clipping-raster-stacks-with-sfst-crop-and-rastercrop

ggplot() +
  geom_raster(data = new_lu_df, aes(x = x, y = y, fill = layer))+
  geom_sf(data = basin_coord, colour = "red", fill = "NA")+
  coord_sf()


my_stack <- stack(raster::crop(new_lu, basin_coord))
my_stack <- stack(raster::mask(my_stack, basin_coord))
my_stack_df <- my_stack %>% as.data.frame(xy=TRUE)

ggplot() +
  geom_raster(data = my_stack_df, aes(x = x, y = y, fill = layer))+
  geom_sf(data = basin_coord, colour = "red", fill = "NA")+
  coord_sf()

cropping_test <- raster::mask(new_lu2, basin_coord)

#new_lu2 = st_crop(new_lu, o0_to_grid)
#test = resample::raster.from.vector(new_lu2)
#
#
#new_lu2 = raster::crop(new_lu, output_geo_raster)
#new_lu3 = raster::crop(new_lu2, output_geo_raster)
#new_lu3 = raster(ext=extent(output_geo), res=c(1000,1000))
#plot(new_lu)
#plot(new_lu2)

}### test end========================================================================

### test plot========================================================================
data = readxl::read_excel('/mnt/d/GitHub/wrfhydroSubsetter/nwm_land_use_mappings.xlsx') %>%
  dplyr::select(nlcd = Class, 
                description = Description,
                nwm = NWM)

col_lu <- data.frame(
  nlcd.code = c(-.1, 0, 11, 12, 21, 22, 23, 24, 31, 41, 42, 43, 51, 52, 71, 72, 73, 74, 81, 82, 90, 95),
  nwm.code  = c(-.1, 0, 16, 23, NA, 1,  NA, NA, 19, 11, 14, 15, 22, 8,  7,  20, NA, NA,  5,  3,  18, 17),
  
  color = c("#000000",
            "#476BA0", "#D1DDF9",
            "#DDC9C9", "#D89382", "#ED0000", "#AA0000",
            "#B2ADA3",
            "#68AA63", "#1C6330", "#B5C98E",
            "#A58C30", "#CCBA7C",
            "#E2E2C1", "#C9C977", "#99C147", "#77AD93",
            "#DBD83D", "#AA7028",
            "#BAD8EA", "#70A3BA", NA) ,
  
  name = c(NA, "EMPTY", "Open Water", "Ice/Snow", "Developed (Open)", "Developed (Low)", 'Developed (Medium)', 'Developed (High)', "Barren",
           "Deciduous Forest", "Evergreen Forest", "Mixed Forest", "Dwarf Scrub", "Shurb", "Grassland", "Sedge", 'Lichens', "Moss",
           "Pasture", "Culitivated Crops", "Woody Wetlands", "Herbaceous Wetlands"),
  
  stringsAsFactors = FALSE)

tt = dplyr::left_join(data, col_lu, by = c("nlcd" = "nlcd.code"))
myColors <- tt$color
names(myColors) <- levels(tt$nwm)
colScale <- scale_colour_manual(name = "grp",values = myColors)



files = list.files('/mnt/d/subsetDOMAINS/polk_oregon_23762661_14190500', "geo", full = TRUE)
s = stack()
for(i in files){
  print(i)
  s=addLayer(s, wrfhydroSubsetter::make_empty_geogrid_raster(i, "LU_INDEX"))
}
names(s) = basename(files)
#plot(s)

tmp = values(s) %>% data.frame() %>% mutate(cell = 1:n()) %>%  tidyr::pivot_longer(-cell)

x = layerStats(s, stat = 'pearson')$pearson
x



#library(rasterVis)
#library(ggplot2)

rasterVis::gplot(s) +
  geom_tile(aes(fill = as.factor(value))) +
  facet_wrap(~ variable, ncol = 2) +
  scale_fill_manual(name = "grp",values = myColors) +
  coord_equal()
  #theme(legend.position = "none")

ggplot(data = tmp) +
  geom_histogram(aes(x = as.factor(value), fill= name), stat = "count", position='dodge') +
  #ggpubr::fill_palette('aaas') +
  theme_linedraw() +
  theme(legend.position = "bottom")

hist(s, breaks = 20)
### test plot end========================================================================




#### Below is wrapper function=================================================
# subsest_sequence is a wrapper that automates the procedure in build_subset.R
subset_sequence = function(comid, siteID, FULLDOMAINDIR, outDir, nlcdDir = NA, methodList = NA){
  
  # Defining the area
  state = st_transform(AOI::aoi_get(state = "all", county = "all"), 4269)[st_transform(findNLDI(comid = comid), crs = 4269),]
  name = gsub(" ", "-", tolower(paste(state$name, state$state_name, comid, siteID, sep = "_")))
  basin = findNLDI(comid = comid, find = c("basin"))$basin

  # Subsetting NetCDF DOMAIN files
  subset_files = paste0(outDir, "/", name , "/")
  subset_wrf_hydro_domain_mod(AOI = basin,  domain_files = FULLDOMAINDIR,  outDir = subset_files, config = 'LongRange')
  
  
  # Resampling sequence
  if (exists("methodList")==T & exists("nlcdObj")==T) {
    data = readxl::read_excel('/mnt/d/GitHub/wrfhydroSubsetter/nwm_land_use_mappings.xlsx') %>%
      dplyr::select(nlcd = Class, 
                    description = Description,
                    nwm = NWM)
    geo  = list.files(subset_files, "geo_em.nc", full.names = TRUE)
    
    output_geo = spex::qm_rasterToPolygons(wrfhydroSubsetter::make_empty_geogrid_raster(geo))
    o2 = output_geo %>%
      AOI::bbox_get() %>%
      sf::st_transform(5070) %>%
      sf::st_buffer(40)
    
    nlcdObj = raster(nlcdDir)
    methodList = methodList
    nlcdCrop = raster::crop(nlcdObj, o2)
    nlcd_nwm = raster::reclassify(nlcdCrop, rcl = dplyr::select(data,from = nlcd, to = nwm ))
    
    for(j in 1:length(methodList)){
      new_lu = resample::resampleData(input = nlcd_nwm, output = output, method = methodList[j])
      out_file = paste0(subset_files, "geogrid_", method[j], "_resample.nc")
      
      file.copy(geo, out_file, overwrite = TRUE)
      nc = ncdf4::nc_open(out_file, write = TRUE)
      ncdf4::ncvar_put(nc, "LU_INDEX",
                       vals = as.vector(t(apply(as.matrix(new_lu), 2, rev))))
      ncdf4::nc_close(nc)
    }
    
  }

  
}
