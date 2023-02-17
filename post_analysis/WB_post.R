data = readxl::read_excel('/mnt/d/GitHub/wrfhydroSubsetter/nwm_land_use_mappings.xlsx') %>%
  dplyr::select(nlcd = Class, 
                description = Description,
                nwm = NWM)

locs = data.frame(comids = c(191739, 1631587, 5894384, 5781369, 19389766, 23762661),
                  siteID = c('6709000', '08173000', '01616500', '08159000' ,'03118500' , '14190500'))

namelist = c()
for (i in 1:nrow(locs)){
  state = st_transform(AOI::aoi_get(state = "all", county = "all"), 4269)[st_transform(findNLDI(comid = locs$comids[i]), crs = 4269),]
  name = gsub(" ", "-", tolower(paste(state$name, state$state_name, locs$comids[i], locs$siteID[i], sep = "_")))
  namelist = c(namelist, name)
  rm(name)
}


col_lu <- data.frame(
  nlcd = c(0, 11, 12, 21, 22, 23, 24, 31, 41, 42, 43, 51, 52, 71, 72, 73, 74, 81, 82, 90, 95),
  nwm.code  = c(0, 16, 23, NA, NA,  1, NA, 19, 11, 14, 15, 22, 8,  7,  20, NA, NA,  5,  3,  18, 17),
  #nwm.code  = c(-.1, 0, 16, 23,  7, 1,  1,  1,  19, 11, 14, 15, 22, 8,  7,  20, NA, NA,  5,  3,  18, 17),
  
  color = c("#000000",
            "#476BA0", "#D1DDF9",
            "#DDC9C9", "#D89382", "#ff0000", "#AA0000",
            "#B2ADA3",
            "#68AA63", "#1C6330", "#B5C98E",
            "#A58C30", "#CCBA7C",
            "#E2E2C1", "#C9C977", "#99C147", "#77AD93",
            "#DBD83D", "#AA7028",
            "#BAD8EA", "#70A3BA") ,
  
  description = c(NA, 
           "Open Water", "Ice/Snow", 
           "Developed (Open)", "Developed (Low)", 'Developed (Medium)', 'Developed (High)', 
           "Barren",
           "Deciduous Forest", "Evergreen Forest", "Mixed Forest", 
           "Dwarf Scrub", "Shurb", 
           "Grassland", "Sedge",'Lichens', "Moss",
           "Pasture", "Culitivated Crops", 
           "Woody Wetlands", "Herbaceous Wetlands"),
  
  stringsAsFactors = FALSE)

lookup_table = na.omit(col_lu) %>% arrange(nwm.code) %>% dplyr::select(nwm.code, description)




for (name in namelist){
  print(name)
  test = read.csv(paste0("/mnt/d/subsetDOMAINS/", name, "/LU_WS_stats.csv")) %>% 
    rename(nwm.code = nwm)
  
  if  (nrow(test %>% filter(nwm.code==1)) > 1){
    x = test %>% filter(nwm.code==1) %>% dplyr::select(nlcd_2016) %>% sum()
    test = test %>% 
      mutate(nlcd_2016 = ifelse (nwm.code == 1, x, nlcd_2016))
  }
  
  if  (nrow(test %>% filter(nwm.code==7)) > 1){
    y = test %>% filter(nwm.code==7) %>% dplyr::select(nlcd_2016) %>% sum()
    test = test %>% 
      mutate(nlcd_2016 = ifelse (nwm.code == 7, y, nlcd_2016)) %>% 
      filter(duplicated(nwm.code) == FALSE)
  }
  
  test = test %>% 
    dplyr::select(-c('nlcd', 'description', 'geogrid_maj_resample.nc', 'geogrid_nn_resample.nc', "geogrid_rawarea_resample.nc", "geogrid_shuffle.nc"))

  
  county_name = str_split(name, "_")[[1]][1]
  colnames(test) = c('nwm.code', paste0(county_name,"_NLCD"), paste0(county_name,"_NWM"), paste0(county_name,"_AP"), paste0(county_name,"_HUC12L"))
  #colnames(test) = c('nwm.code', paste0("NLCD"), paste0("NWM"), paste0("AP"), paste0("HUC12L"))
  
  lookup_table = left_join(lookup_table, test, by =c('nwm.code'))
}


write.csv(lookup_table, file = '/mnt/d/ANALYSIS/LU_appendix_tbl2.csv')

