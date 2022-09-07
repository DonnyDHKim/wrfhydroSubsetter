{
  library(dataRetrieval); library(nhdplusTools);
  library(raster); library(sf); library(fasterize);
  library(resample); library(wrfhydroSubsetter); library(AOI);
  library(tidyverse); library(dplyr);
}

# This script is prone to crash during debugging. In such cases, just restart R session. It solves most (if not all) of the problems.
.rs.restartR()


### STEP 1======================================================================
# Setting up locations
locs = data.frame(comids = c(23762661, 1631587, 5894384, 19389766, 191739, 5781369),
                  siteID = c('14190500', '14190500', '01616500', '03118500' ,'6709000' , '08159000'))

#test_loc = data.frame(comids = c(5894384), siteID = c('01616500')) # In case using a single watershed

namelist = c()
for (i in 1:nrow(locs)){
  state = st_transform(AOI::aoi_get(state = "all", county = "all"), 4269)[st_transform(findNLDI(comid = locs$comids[i]), crs = 4269),]
  name = gsub(" ", "-", tolower(paste(state$name, state$state_name, locs$comids[i], locs$siteID[i], sep = "_")))
  namelist = c(namelist, name)
  rm(name)
  }

# Setting directories
outDir = "/mnt/d/subsetDOMAINS/"
FULLDOMAINDIR = "/mnt/d/FULLDOMAIN/nwmCONUS-v216/"
### STEP 1 END======================================================================




### STEP 2======================================================================
## Package development purpose
#library(devtools)
#detach("package:wrfhydroSubsetter", unload=TRUE); remove.packages("wrfhydroSubsetter");
#devtools::install("/mnt/d/GitHub/wrfhydroSubsetter", repos = NULL, type="source")
#library(wrfhydroSubsetter)

# Sourcing modified version of subsetter. Because it is still in development phase.
source("/mnt/d/GitHub/wrfhydroSubsetter/steering_file/subsetter_mod.R")
source("/mnt/d/GitHub/wrfhydroSubsetter/steering_file/resampleDataMod.R")
source("/mnt/d/GitHub/wrfhydroSubsetter/steering_file/subset_sequence.R")
### STEP 2 END======================================================================




### STEP 3======================================================================
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
### STEP 3 END======================================================================




### From step 4, it is project-specific code: not for typical use of wrfhydrosubsetter.
### STEP 4======================================================================
# Watershed Shuffle Loop
for(i in 1:nrow(locs)){
  # Specifying WS area
  state = st_transform(AOI::aoi_get(state = "all", county = "all"), 4269)[st_transform(findNLDI(comid = locs$comids[i]), crs = 4269),]
  name = gsub(" ", "-", tolower(paste(state$name, state$state_name, locs$comids[i], locs$siteID[i], sep = "_")))
  basin = findNLDI(comid = locs$comids[i], find = c("basin"))$basin
  
  # Using Area Preserve sampled GEOGRID to shuffle LC within watershed boundary
  area_resample_Dir = list.files(paste0('/mnt/d/subsetDOMAINS/', name), "geogrid_area_resample", full = TRUE)
  area_resample = wrfhydroSubsetter::make_empty_geogrid_raster(area_resample_Dir, "LU_INDEX")
  area_resample_shuffle = area_resample # RasterLayer
  area_resample_mask_df = mask_geogrid_byWS(area_resample, basin) #mask_geogrid_byWS is in subset_sequence.R
  
  # Shuffle LULC inside WS. Cells that are not in WS boundary has value of NA, and omitting them.
  # Then, perturb/shuffle cell values using "sample" function.
  shuffled = transform(na.omit(area_resample_mask_df), layer = sample(layer)) 
  values(area_resample_shuffle)[as.numeric(rownames(shuffled))] <- shuffled$layer # Overwriting shuffled value over RasterLayer.
  # This overwriting relies on "rownames" which is technically a row index.
  
  # Setting output file directory
  out_file = paste0('/mnt/d/subsetDOMAINS/', name, "/geogrid_shuffle.nc")
  file.copy(area_resample_Dir, out_file, overwrite = TRUE)
  
  #using ncdf4 package to write separate GEOGRID with shuffled LU_INDEX
  nc = ncdf4::nc_open(out_file, write = TRUE)
  ncdf4::ncvar_put(nc, "LU_INDEX", vals = as.vector(t(apply(as.matrix(area_resample_shuffle), 2, rev))))
  ncdf4::nc_close(nc)
  print(paste0(name, ": Done"))
}
### STEP 4 END======================================================================




### STEP 4-1 ======================================================================
# HUC 12 Single LU Loop
for(i in 1:nrow(locs)){
  name = namelist[i]; basin = findNLDI(comid = locs$comids[i], find = c("basin"))$basin;
  
  # Load every geogrid netcdf files
  files = list.files(paste0('/mnt/d/subsetDOMAINS/', name), "geo", full = TRUE);
  ss = stack()
  for(j in files){
    ss=addLayer(ss, wrfhydroSubsetter::make_empty_geogrid_raster(j, "LU_INDEX"))
  }
  names(ss) = basename(files)
  
  # crop them to geogrid box
  geogrid_box = spex::qm_rasterToPolygons(ss[[1]]) %>% AOI::bbox_get() %>% sf::st_transform(5070)# %>% sf::st_buffer(40)
  
  ### Resampling scheme that assigns single LULC values to  each HUC12 ws within study basins
  basin_geo_coord = st_transform(basin, st_crs(ss[[2]])) #ss[[2]] = area preserving geogrid.
  huc12s = get_huc12(AOI = geogrid_box) %>% st_transform(st_crs(basin_geo_coord))
  huc12s = huc12s %>% st_intersection(basin_geo_coord) #only when you want to "clip" huc12 polygons by basin-polygon
  
  out_geo = ss[[2]]; out_geo_df = as.data.frame(ss[[2]], xy=TRUE);
  
  # Application of step 4 & 5.
  for (n in 1:nrow(huc12s)){
    #huc12s_df =  mask_geogrid_byWS(out_geo, huc12s[n,]) %>% na.omit() #n stands for each HUC12 ws
    huc12s_df =  raster::mask(out_geo, huc12s[n,]) %>% as.data.frame(xy = TRUE) %>% na.omit()
    huc_dom_LU = huc12s_df %>%  group_by_at(3) %>%  tally() %>% filter (n==max(n)) %>% select_at(1)
    
    if (nrow(huc_dom_LU)>1) {
      huc12s_df =  mask_geogrid_byWS(out_geo, huc12s[n,]) %>% na.omit() #n stands for each HUC12 ws
      huc_dom_LU = huc12s_df %>%  group_by_at(3) %>%  tally() %>% filter (n==max(n)) %>% select_at(1) %>% as.integer()
    } else {
      huc_dom_LU = huc_dom_LU %>% as.integer() #finding the most dominant LULC
    }
    
    raster::values(out_geo)[as.numeric(rownames(huc12s_df))] = huc_dom_LU # Updating the raster cell values
    
  }
  
  
  # Setting output file directory
  out_file = paste0('/mnt/d/subsetDOMAINS/', name, "/geogrid_huc12uniform.nc")
  
  file.copy(paste0('/mnt/d/subsetDOMAINS/', name, "/geogrid_area_resample.nc"), 
            out_file, 
            overwrite = TRUE)
  
  #using ncdf4 package to write separate GEOGRID with HUC12 ws having single LU representation
  nc = ncdf4::nc_open(out_file, write = TRUE)
  ncdf4::ncvar_put(nc, "LU_INDEX", vals = as.vector(t(apply(as.matrix(out_geo), 2, rev))))
  ncdf4::nc_close(nc)
  print(paste0(name, ": HUC12 uniform LC sampling done"))
}


###### Plot
out_geo_mapping = as.data.frame(out_geo, xy=TRUE) %>% rename(LU = 3)
ggplot()+
  #geom_raster(as.data.frame(ss[[2]], xy=T),mapping = aes(x=x, y=y, fill=as.factor(geogrid_area_resample.nc)))+
  geom_raster(out_geo_mapping,mapping = aes(x=x, y=y, fill=as.factor(LU)))+
  #geom_sf(basin_geo_coord, mapping = aes(fill=NULL), color = 'red', alpha = 0)+
  geom_sf(huc12s, mapping = aes(fill=NULL), alpha = 0)#
###### Plot end
### STEP 4-1 END======================================================================




### STEP 5======================================================================
# WS LULC statistics
lookup_table = readxl::read_excel('/mnt/d/GitHub/wrfhydroSubsetter/nwm_land_use_mappings.xlsx') %>%
  dplyr::select(nlcd = Class, description = Description, nwm = NWM) %>% select(nwm, nlcd, description) %>% arrange(nwm, nlcd)

for(i in 1:nrow(locs)){
  name = namelist[i]; basin = findNLDI(comid = locs$comids[i], find = c("basin"))$basin;
  
  # Load every geogrid netcdf files
  files = list.files(paste0('/mnt/d/subsetDOMAINS/', name), "geo", full = TRUE)
  ss = stack()
  for(j in files){
    ss=addLayer(ss, wrfhydroSubsetter::make_empty_geogrid_raster(j, "LU_INDEX"))
  }
  names(ss) = basename(files)

  # crop them to geogrid box
  geogrid_box = spex::qm_rasterToPolygons(ss[[1]]) %>%
    AOI::bbox_get() %>%
    sf::st_transform(5070) %>%
    sf::st_buffer(40)
  nlcd_crop = raster::crop(nlcdObj, geogrid_box) 
  
  # Masking (clipping) cropped nlcd into watershed boundary
  # I am not using mask_geogrid_byWS function, as it leads to error.
  basin_coord = st_transform(basin, st_crs(nlcd_crop))
  nlcd_mask_df = raster::mask(nlcd_crop, fasterize::fasterize(basin_coord, nlcd_crop)) %>% getValues() %>% as.data.frame(xy = TRUE)
  nlcd_basin_stat = nlcd_mask_df %>% na.omit() %>% group_by_at(1) %>% tally() %>% rename(nlcd = 1) %>% rename(nlcd_2016 = 2)
  nlcd_basin_stat[,2] = nlcd_basin_stat[,2]/(sum(nlcd_basin_stat[,2])) * 100 # percent
  
  # I didn't want to spend time to write code that uses map or apply 
  mask_LU_stat1 = mask_geogrid_byWS(ss[[1]], basin) %>% na.omit() %>%  group_by_at(3) %>% tally() %>% rename(nwm = 1) %>%  rename(!!basename(files)[1] := 2)
  mask_LU_stat2 = mask_geogrid_byWS(ss[[2]], basin) %>% na.omit() %>%  group_by_at(3) %>% tally() %>% rename(nwm = 1) %>%  rename(!!basename(files)[2] := 2)
  mask_LU_stat3 = mask_geogrid_byWS(ss[[3]], basin) %>% na.omit() %>%  group_by_at(3) %>% tally() %>% rename(nwm = 1) %>%  rename(!!basename(files)[3] := 2)
  mask_LU_stat4 = mask_geogrid_byWS(ss[[4]], basin) %>% na.omit() %>%  group_by_at(3) %>% tally() %>% rename(nwm = 1) %>%  rename(!!basename(files)[4] := 2)
  mask_LU_stat5 = mask_geogrid_byWS(ss[[5]], basin) %>% na.omit() %>%  group_by_at(3) %>% tally() %>% rename(nwm = 1) %>%  rename(!!basename(files)[5] := 2)
  mask_LU_stat6 = mask_geogrid_byWS(ss[[6]], basin) %>% na.omit() %>%  group_by_at(3) %>% tally() %>% rename(nwm = 1) %>%  rename(!!basename(files)[6] := 2)
  mask_LU_stat7 = mask_geogrid_byWS(ss[[7]], basin) %>% na.omit() %>%  group_by_at(3) %>% tally() %>% rename(nwm = 1) %>%  rename(!!basename(files)[7] := 2)
  
  #merge all data frames in list
  df_list <- list(mask_LU_stat1, mask_LU_stat2,mask_LU_stat3, mask_LU_stat4, mask_LU_stat5, mask_LU_stat6, mask_LU_stat7) %>% reduce(full_join, by='nwm') %>% as.data.frame()
  rm(mask_LU_stat1, mask_LU_stat2,mask_LU_stat3, mask_LU_stat4, mask_LU_stat5, mask_LU_stat6, mask_LU_stat7)
  
  # make it percent
  for (k in 2:ncol(df_list)){
    df_list[,k] <- df_list[,k]/sum(df_list[,k], na.rm=TRUE) * 100
  }
  
  # join nlcd and different geogrid LULC stats
  LU_stat_merged = lookup_table %>% right_join(nlcd_basin_stat) %>% right_join(df_list)
  LU_stat_merged[,4:ncol(LU_stat_merged)] = round(LU_stat_merged[,4:ncol(LU_stat_merged)], digits = 1)

  # Let's also look into pearson spatial correlation between different resampling methods
  x = layerStats(ss, stat = 'pearson')$pearson %>% as.data.frame()
  
  # Save it as csv files
  write_csv(LU_stat_merged, file = paste0(outDir, "/", name, "/LU_WS_stats.csv"))
  write.csv(x, file = paste0(outDir, "/", name, "/Pearson.csv"), row.names= TRUE)
  print(paste0(name, ": Basin stats and Pearson spatial correlation saved as CSV"))
}
### STEP 5 END======================================================================




# Plotting
## Simply checking out.
## https://datacarpentry.org/r-raster-vector-geospatial/02-raster-plot/
#ggplot() +
#  geom_raster(data = area_resample_mask_df, aes(x = x, y = y, fill = as.factor(layer)))+
#  geom_sf(data = basin_coord, colour = "black", fill = "NA")+
#  scale_fill_manual(name = "grp",values = myColors) +
#  coord_sf()

# Plotting wrapper: messy but works
{
  
}


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



files = list.files(paste0('/mnt/d/subsetDOMAINS/',namelist[3]), "geo", full = TRUE)
s = stack()
for(i in files){
  print(i)
  s=addLayer(s, wrfhydroSubsetter::make_empty_geogrid_raster(i, "LU_INDEX"))
}
names(s) = basename(files)
#plot(s)

tmp = values(s) %>% data.frame() %>% mutate(cell = 1:n()) %>%  tidyr::pivot_longer(-cell)


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

