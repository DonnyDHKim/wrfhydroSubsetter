{
  library(sf)
  library(lubridate)
  #library(RNetCDF)
  library(dplyr)
  #library(ggplot2)
  library(raster)
  #library(stars)
  #library(gdalUtils)
  library(tidyverse)
  library(dataRetrieval)
  library(rgdal) # extremely important to use rgdal
  library(terra) # extremely important to use rgdal
  library(tictoc)
}

#library(devtools)
#devtools::install("/mnt/d/GitHub/wrfhydroSubsetter", repos = NULL, type="source")
#devtools::install("D:/GitHub/wrfhydroSubsetter", repos = NULL, type="source")
#library(wrfhydroSubsetter)

###
# Import comid list csv
getwd()
#setwd('/mnt/d/')
setwd('D:/')
site_tbl = read.csv('./GitHub/wrfhydroSubsetter/steering_file/site_selection_stage1_COMIDs.csv')
stage2 =  read.csv('./GitHub/wrfhydroSubsetter/steering_file/site_selection_stage2_COMIDs.csv')

site_tbl = left_join(stage2, site_tbl)


stage_add =  read.csv('./GitHub/wrfhydroSubsetter/steering_file/comids_dates_for_donny.csv')
stage_add$startDT = stage_add$startDT %>% as.Date("%m/%d/%Y")
stage_add$endDT = stage_add$endDT %>% as.Date("%m/%d/%Y")
site_tbl = stage_add


site_tbl = site_tbl %>% 
  dplyr::mutate(startWY = ifelse(lubridate::month(startDT) >= 10, lubridate::year(startDT) + 1, lubridate::year(startDT)),
                endWY = ifelse(lubridate::month(endDT) >= 10, lubridate::year(endDT) + 1, lubridate::year(endDT))
  )

NLDAS_dir = './NLDAS'
out_dir = './subsetNLDAS_add'; dir.create(out_dir);
#comid = 9453855
#startWY = 1999
#endWY = 2000

library(parallel); library(doParallel); library(foreach);

comid_list = site_tbl$COMID
startWY_list = site_tbl$startWY
endWY_list = site_tbl$endWY

in_list = list(comid_list=comid_list, startWY_list=startWY_list, endWY_list=endWY_list)

NLDAS_to_ts = function(comid, startWY, endWY){
  tic()
  require(tidyverse); require(foreach); require(raster); require(sf); require(rgdal);
  #require(parallel); require(doParallel); 

  # build a df that contains NLDAS file directory and its date/time
  if (startWY < 2000){
    fileList_pt1 = data.frame(files = list.files(paste0(NLDAS_dir, "/WY1980_WY1999"), pattern = ".grb", full.name = TRUE)) %>% 
      dplyr::mutate(date  = lubridate::ymd_hm(gsub("\\.","", str_split(str_split(files, "_FORA0125_H.A", simplify = TRUE)[,2], ".002.grb", simplify = TRUE)[,1])),
                    wy = ifelse(lubridate::month(date) >= 10, lubridate::year(date) + 1, lubridate::year(date)))
    fileList_pt1 = fileList_pt1 %>% filter(wy >= startWY)
    
    fileList_pt2 = data.frame(files = list.files(paste0(NLDAS_dir, "/WY2000_"), pattern = ".grb", full.name = TRUE)) %>% 
      dplyr::mutate(date  = lubridate::ymd_hm(gsub("\\.","", str_split(str_split(files, "_FORA0125_H.A", simplify = TRUE)[,2], ".002.grb", simplify = TRUE)[,1])),
                    wy = ifelse(lubridate::month(date) >= 10, lubridate::year(date) + 1, lubridate::year(date)))
    fileList_pt2 = fileList_pt2 %>% filter(wy <= endWY)
    fileList_final = rbind(fileList_pt1, fileList_pt2)
  } else {
    fileList_final = data.frame(files = list.files(paste0(NLDAS_dir, "/WY2000_"), pattern = ".grb", full.name = TRUE)) %>% 
      dplyr::mutate(date  = lubridate::ymd_hm(gsub("\\.","", str_split(str_split(files, "_FORA0125_H.A", simplify = TRUE)[,2], ".002.grb", simplify = TRUE)[,1])),
                    wy = ifelse(lubridate::month(date) >= 10, lubridate::year(date) + 1, lubridate::year(date)))
    fileList_final %>% filter(wy >= startWY) %>% filter(wy <= endWY)
  }
  rm(fileList_pt1, fileList_pt2)
  print(paste0(comid, ": NLDAS directory identified"))
  
  basin = dataRetrieval::findNLDI(comid = comid, find = c("basin"))$basin
  x = raster::stack(fileList_final$files[1])
  basin_coord = st_transform(basin, st_crs(x))
  maskBas = raster::rasterize(basin_coord, x, getCover=TRUE)
  rm(x)
  maskMatrix = raster::as.matrix(maskBas) 
  maskMatrix[maskMatrix==0] = NA
  print(paste0(comid, ": Now NLDAS subset and matrix calculation..."))
  
  ## snow progress bar
  #pb <- txtProgressBar(max = 100, style = 3)
  #progress <- function(n) setTxtProgressBar(pb, n)
  #opts <- list(progress = progress)
  
  out = foreach (
    file = fileList_final$files, #.options.snow = opts,
    .combine='rbind', 
    .packages = c("dplyr", "terra"))  %dopar% {
    
    r = terra::rast(file)
    data.frame(sum(matrix(r[[1]], nrow=224, ncol=464, byrow=T) * maskMatrix, na.rm = TRUE)/sum(maskMatrix[which(!is.na(maskMatrix))]), 
               sum(matrix(r[[9]], nrow=224, ncol=464, byrow=T) * maskMatrix, na.rm = TRUE)/sum(maskMatrix[which(!is.na(maskMatrix))]), 
               sum(matrix(r[[10]], nrow=224, ncol=464, byrow=T) * maskMatrix, na.rm = TRUE)/sum(maskMatrix[which(!is.na(maskMatrix))])
               ) %>% 
      setNames(c("T", "PET","P"))
    }
  print(paste0(comid, ": Done. Writing CSV."))

  o2    = cbind(fileList_final, dplyr::bind_rows(out))[,2:6]
  write.csv(o2, paste0(out_dir, "/", comid,"_",startWY,"_",endWY, ".csv"))
  #close(pb)
  toc()
}

#### test1: pmap
#library(furrr); library(purrr);
#plan(multisession, workers = 30)
#tic(); future_pmap(in_list, ~NLDAS_to_ts(comid=..1, startWY=..2, endWY=..3)); toc();
#future:::ClusterRegistry("stop")


# Use this one: for loop with inner %dopar%
cl <- parallel::makeCluster(4) #Somehow 3 is the fastest in my computer setting
registerDoParallel(cl)
{
  tic()
  for (i in 1:5){
    NLDAS_to_ts(comid_list[i], startWY_list[i], endWY_list[i])
  }
  toc()
}
parallel::stopCluster(cl) # Important to do this.


# dosnow variation test
library(doSNOW)
library(itertools)
cl <- makePSOCKcluster(3)
registerDoSNOW(cl)
{
  tic()
  NLDAS_to_ts(comid_list[2], startWY_list[2], endWY_list[2])
  toc()
}
stopCluster(cl) 



# test3: nested foreach loop with inner %dopar%
in_df = data.frame(comid_list=comid_list, startWY_list=startWY_list, endWY_list=endWY_list)
{
  cl <- parallel::makeCluster(30) # don't set makeCluster(n) too high on laptops.
  registerDoParallel(cl)
  tic()
  foreach(i = 1:nrow(in_df)) %dopar% {
    NLDAS_to_ts(in_df[i, 1], in_df[i, 2], in_df[i, 3])
  } 
  toc()
  parallel::stopCluster(cl) # Important to do this.
}
#### trash


# single-case test
cl <- parallel::makeCluster(3) # don't set makeCluster(n) too high on laptops.
registerDoParallel(cl)
{
  tic()
  NLDAS_to_ts(comid_list[2], startWY_list[2], endWY_list[2])
  toc()
}
parallel::stopCluster(cl) # Important to do this.




################# Junk below

test_matrix = matrix(raster::values(r[[1]]), nrow=224, ncol=464, byrow=T)
test_matrix = matrix(r[[1]], nrow=224, ncol=464, byrow=T)
test_matrix = matrix(r[[1]], nrow=224, ncol=464, byrow=T)
test_matrix2 = t(test_matrix)

maskMatrix[which(!is.na(maskMatrix))]
test_matrix[which(!is.na(maskMatrix))]
sum(test_matrix * maskMatrix, na.rm = TRUE) / sum(maskMatrix[which(!is.na(maskMatrix))])
mean(test_matrix2 * maskMatrix, na.rm=TRUE)

file = fileList_final$files[1]


{
  tic()
  a = mean(values(brick(file))[[1]]* maskMatrix, na.rm = TRUE)* sum(!is.na(maskMatrix))
  a = mean(values(brick(file))[[1]]* maskMatrix, na.rm = TRUE)* sum(!is.na(maskMatrix))
  a = mean(values(brick(file))[[1]]* maskMatrix, na.rm = TRUE)* sum(!is.na(maskMatrix))
  a = mean(values(brick(file))[[1]]* maskMatrix, na.rm = TRUE)* sum(!is.na(maskMatrix))
  a = mean(values(brick(file))[[1]]* maskMatrix, na.rm = TRUE)* sum(!is.na(maskMatrix))

  r = terra::rast(file)
  test = matrix(values(r[[1]]), dim(r)[1], dim(r)[2], byrow=T)
  test = matrix(values(r[[1]]), dim(r)[1], dim(r)[2], byrow=T)
  #matrix(values(r), d[1], d[2], byrow=TRUE)
  
  b = mean(test* maskMatrix, na.rm = TRUE)* sum(!is.na(maskMatrix))
  b = sum(test* maskMatrix, na.rm = TRUE)
  #b = mean(values(terra::rast(file))[[1]]* maskMatrix, na.rm = TRUE)* sum(!is.na(maskMatrix))
  #b = mean(values(terra::rast(file))[[1]]* maskMatrix, na.rm = TRUE)* sum(!is.na(maskMatrix))
  #b = mean(values(terra::rast(file))[[1]]* maskMatrix, na.rm = TRUE)* sum(!is.na(maskMatrix))
  #b = mean(values(terra::rast(file))[[1]]* maskMatrix, na.rm = TRUE)* sum(!is.na(maskMatrix))
  
  toc()
}
test[which(!is.na(maskMatrix))];  maskMatrix[which(!is.na(maskMatrix))];



#https://gis.stackexchange.com/questions/360547/era5-grib-file-how-to-know-what-each-band-means

x = brick(file)
y = as.matrix(x[[1]])

y[which(!is.na(maskMatrix))];  maskMatrix[which(!is.na(maskMatrix))];
y[91049] = 0.39

mean(y[which(!is.na(maskMatrix))], na.rm = TRUE)
mean(y * maskMatrix, na.rm = TRUE)* sum(!is.na(maskMatrix)) # this one works the best
sum(y * maskMatrix, na.rm = TRUE)

sum(as.numeric(values(x[[1]])) * maskMatrix, na.rm = TRUE)*sum(!is.na(maskMatrix))
mean(values(x[[9]])[c(90601,90602,90825,90826,91049,91050)], na.rm = TRUE)
mean(values(x[[1]])[which(!is.na(maskMatrix))], na.rm = TRUE)
mean(values(x[[1]])[which(!is.na(maskMatrix))] * maskMatrix[which(!is.na(maskMatrix))], na.rm = TRUE)

mean((y[y=="NaN"] == NA) * maskMatrix, na.rm = TRUE)
y[91049]

nldas_test_dir = "/mnt/d/NLDAS/1999_2010/NLDAS_FORA0125_H.A19991001.0000.002.grb"
meta = gdalinfo(nldas_test_dir)

file= fileList_final$files[1]

bands <- sf::gdal_utils(source = file,
                        options = c("-json"), 
                        quiet = T) %>%
  jsonlite::fromJSON()

bands$bands$metadata # So T:1, PET: 9, P:10

#s= stars::read_stars(nldas_test_dir)
#plot(s)

x<-stack(nldas_test_dir) 

extract_T = x[[1]]
extract_PET = x[[9]]
extract_P = x[[10]]

basin_coord = st_transform(basin, st_crs(x))
maskBas = raster::rasterize(basin_coord, extract_T, getCover=TRUE)
maskMatrix = raster::as.matrix(maskBas) 
maskMatrix[maskMatrix==0] = NA


  # Read in variable, transpose, multiply by basin mask, and average
  # Result is average *unit* of *variable* in the basin at time of file
T_avg =   mean(values(extract_T) * maskMatrix, na.rm = TRUE)


#nlcd_mask_df = raster::mask(nlcd_crop, fasterize::fasterize(basin_coord, nlcd_crop)) %>% getValues() %>% as.data.frame(xy = TRUE)
#mask_nldas_T = raster::mask(extract_T, fasterize::fasterize(basin_coord, extract_T)) %>% getValues() %>% as.data.frame(xy = TRUE)

#mask_nldas_T = raster::mask(extract_T,basin_coord) %>% na.omit()



{
  test = as(basin_coord, 'Spatial') %>% fortify()
  
  # https://stackoverflow.com/questions/68865682/overlay-polygon-layer-on-a-raster-stack-qplot
  rasterVis::gplot(extract_T) +
    geom_tile(aes(fill = as.factor(value))) +
    geom_path(data=test, aes(long, lat, group=group), color = 'black')+
    facet_wrap(~ variable, ncol = 2) +
    scale_fill_brewer() +
    coord_equal()
  #theme(legend.position = "none")
}

