### script for calculating spatial stats from LULC geogrid

{
  library(raster); library(sf); library(fasterize);
  library(resample); library(wrfhydroSubsetter); library(AOI);
  library(tidyverse); library(dplyr);
  library(spatialEco);
}

source("/mnt/d/GitHub/wrfhydroSubsetter/steering_file/subset_sequence.R")


locs = data.frame(comids = c(23762661, 1631587, 5894384, 19389766, 191739, 5781369),
                  siteID = c('14190500', '08159000', '01616500', '03118500' ,'6709000' , '08159000'))

#test_loc = data.frame(comids = c(5894384), siteID = c('01616500')) # In case using a single watershed

namelist = c()
for (i in 1:nrow(locs)){
  state = st_transform(AOI::aoi_get(state = "all", county = "all"), 4269)[st_transform(findNLDI(comid = locs$comids[i]), crs = 4269),]
  name = gsub(" ", "-", tolower(paste(state$name, state$state_name, locs$comids[i], locs$siteID[i], sep = "_")))
  namelist = c(namelist, name)
  rm(name)
}

scheme_list = c('NWM', 'AP',  'HUC12L', 'PERT', 'SHUF', 'maj', 'nn')

outDir = "/mnt/d/ANALYSIS/"


lookup_table = readxl::read_excel('/mnt/d/GitHub/wrfhydroSubsetter/nwm_land_use_mappings.xlsx') %>%
  dplyr::select(nlcd = Class, description = Description, nwm = NWM) %>% dplyr::select(nwm, nlcd, description) %>% arrange(nwm, nlcd)


for(i in 1:nrow(locs)){ #nrow(locs)
  name = namelist[i]; basin = findNLDI(comid = locs$comids[i], find = c("basin"))$basin;
  
  print(paste0(name, ": Starting..."))
  # Load every geogrid netcdf files
  files = list.files(paste0('/mnt/d/subsetDOMAINS/', name), "geo", full = TRUE)
  ss = stack()
  for(j in files){
    ss=addLayer(ss, wrfhydroSubsetter::make_empty_geogrid_raster(j, "LU_INDEX"))
  }
  names(ss) = basename(files)
  
  basin_coord = st_transform(basin, st_crs(ss[[1]]))
  
  raw_NWM    = ss[[2]] %>% as('SpatialGridDataFrame')
  raw_AP     = ss[[1]] %>% as('SpatialGridDataFrame')
  raw_HUC12L = ss[[3]] %>% as('SpatialGridDataFrame')
  raw_PERT   = ss[[4]] %>% as('SpatialGridDataFrame')
  raw_SHUF   = ss[[5]] %>% as('SpatialGridDataFrame')
  raw_maj    = ss[[6]] %>% as('SpatialGridDataFrame')
  raw_nn     = ss[[7]] %>% as('SpatialGridDataFrame')
  
  raw_list = list(raw_NWM, raw_AP, raw_HUC12L, raw_PERT, raw_SHUF, raw_maj, raw_nn)
  corr_df = data.frame(matrix(nrow = length(raw_list), ncol = length(raw_list)));   rownames(corr_df) = scheme_list;   colnames(corr_df) = scheme_list;
  cell_diff_df = data.frame(matrix(nrow = length(raw_list), ncol = length(raw_list)));
  rownames(cell_diff_df) = scheme_list;   colnames(cell_diff_df) = scheme_list;
  
  row_cnt = 0; col_cnt =0;
  tic()
  for (row in raw_list){
    row_cnt = row_cnt + 1
    col_cnt = 0
    for (col in raw_list){
      col_cnt = col_cnt +1
      
      # Difference between cell-counts
      cell_diff_temp = row@data != col@data; cell_diff_temp = sum(cell_diff_temp==TRUE)/length(cell_diff_temp);
      cell_diff_df[row_cnt, col_cnt] = cell_diff_temp
      
      ## Hoping d=9000 is enough. Except Douglas (n=5). 15000
      #corr_temp = spatialEco::raster.modified.ttest(row, col, d=9000) %>% 
      #  raster() %>% mask_geogrid_byWS(basin) %>% 
      #  na.omit() %>% dplyr::select(corr)
      #corr_temp = mean(corr_temp$corr)
      #corr_df[row_cnt, col_cnt] = corr_temp
    }
  }
  print(paste0(name, ": Done."))
  toc()
  
  # Save it as csv files
  #write.csv(corr_df, file = paste0(outDir, "/", name, "/", locs$comids[i],"_Raster_Correlation.csv"), row.names= TRUE)
  write.csv(cell_diff_df, file = paste0(outDir, "/", name, "/", locs$comids[i],"_Cell_Diff.csv"), row.names= TRUE)
}


### Stupid but yeah, let's summarize now
corr_1 = read.csv(paste0(outDir, "/", namelist[1],"/", locs$comids[1], "_Raster_Correlation.csv"), header = TRUE, row.names = 1) %>% matrix()
corr_2 = read.csv(paste0(outDir, "/", namelist[2],"/", locs$comids[2], "_Raster_Correlation.csv"), header = TRUE, row.names = 1) %>% matrix()
corr_3 = read.csv(paste0(outDir, "/", namelist[3],"/", locs$comids[3], "_Raster_Correlation.csv"), header = TRUE, row.names = 1) %>% matrix()
corr_4 = read.csv(paste0(outDir, "/", namelist[4],"/", locs$comids[4], "_Raster_Correlation.csv"), header = TRUE, row.names = 1) %>% matrix()
corr_5 = read.csv(paste0(outDir, "/", namelist[5],"/", locs$comids[5], "_Raster_Correlation.csv"), header = TRUE, row.names = 1) %>% matrix()
corr_6 = read.csv(paste0(outDir, "/", namelist[6],"/", locs$comids[6], "_Raster_Correlation.csv"), header = TRUE, row.names = 1) %>% matrix()

corr_list <- list(corr_1, corr_2, corr_3, corr_4, corr_5, corr_6)

corr_mean= Reduce("+", corr_list) / length(corr_list)
corr_max = Reduce(pmax, corr_list)
corr_min = Reduce(pmin, corr_list)
rm(corr_1, corr_2, corr_3, corr_4, corr_5, corr_6)

cell_diff_1 = read.csv(paste0(outDir, "/", namelist[1],"/", locs$comids[1], "_Cell_Diff.csv"), header = TRUE, row.names = 1)
cell_diff_2 = read.csv(paste0(outDir, "/", namelist[2],"/", locs$comids[2], "_Cell_Diff.csv"), header = TRUE, row.names = 1)
cell_diff_3 = read.csv(paste0(outDir, "/", namelist[3],"/", locs$comids[3], "_Cell_Diff.csv"), header = TRUE, row.names = 1)
cell_diff_4 = read.csv(paste0(outDir, "/", namelist[4],"/", locs$comids[4], "_Cell_Diff.csv"), header = TRUE, row.names = 1)
cell_diff_5 = read.csv(paste0(outDir, "/", namelist[5],"/", locs$comids[5], "_Cell_Diff.csv"), header = TRUE, row.names = 1)
cell_diff_6 = read.csv(paste0(outDir, "/", namelist[6],"/", locs$comids[6], "_Cell_Diff.csv"), header = TRUE, row.names = 1)

cell_diff_list <- list(cell_diff_1, cell_diff_2, cell_diff_3, cell_diff_4, cell_diff_5, cell_diff_6)
cell_diff_mean= 1 -(Reduce("+", cell_diff_list) / length(cell_diff_list))
cell_diff_min = 1- Reduce(pmax, cell_diff_list)
cell_diff_max = 1- Reduce(pmin, cell_diff_list)
rm(cell_diff_1, cell_diff_2, cell_diff_3, cell_diff_4, cell_diff_5, cell_diff_6)

write.csv(corr_mean, file=paste0(outDir, "/Raster_Correlation_mean.csv"))
write.csv(cell_diff_mean, file=paste0(outDir, "/Cell_Diff_mean.csv"))



### RMSE between NLCD vs AP,NWM, HUC12L
rmse_1 = read.csv(paste0("/mnt/d/subsetDOMAINS", "/", namelist[1],"/LU_WS_stats.csv"), header = TRUE) %>% mutate(NWM = modelr::rmse(geo_em.nc, nlcd_2016))
rmse_2 = read.csv(paste0("/mnt/d/subsetDOMAINS", "/", namelist[2],"/LU_WS_stats.csv"), header = TRUE, row.names = 1) %>% matrix()
rmse_3 = read.csv(paste0("/mnt/d/subsetDOMAINS", "/", namelist[3],"/LU_WS_stats.csv"), header = TRUE, row.names = 1) %>% matrix()
rmse_4 = read.csv(paste0("/mnt/d/subsetDOMAINS", "/", namelist[4],"/LU_WS_stats.csv"), header = TRUE, row.names = 1) %>% matrix()
rmse_5 = read.csv(paste0("/mnt/d/subsetDOMAINS", "/", namelist[5],"/LU_WS_stats.csv"), header = TRUE, row.names = 1) %>% matrix()
rmse_6 = read.csv(paste0("/mnt/d/subsetDOMAINS", "/", namelist[6],"/LU_WS_stats.csv"), header = TRUE, row.names = 1) %>% matrix()
