{
  library(sf)
  library(lubridate)
  library(RNetCDF)
  library(dplyr)
  library(ggplot2)
  library(raster)
}



## Base Material
#=======================================================
name = 'douglas_colorado_191739_6709000'
run_timestamp = '20220923_1945'
scheme_list = c('NWM', 'AP', 'SHUF', 'HUC12L', 'PERT')


#These are the variable names, descriptions and units we need
df = data.frame(variableNames = c("ACCECAN", "ACCEDIR", "ACCETRAN", "ACCPRCP",
                                  "CANICE", "CANLIQ", "SFCRUNOFF", "UDRUNOFF", "SNEQV"),
                
                description = c("Accumulated canopy evaporation",
                                "Accumulated direct soil evaporation",
                                "Accumulated transpiration",
                                "Accumulated precipitation",
                                "Canopy ice water content",
                                "Canopy liquid water content",
                                "Accumulated surface runoff",
                                "Accumulated underground runoff",
                                "Snow water equivalent"),
                
                units = c("mm", "mm", "mm", "mm", "mm", "mm", "mm", "mm", "kgm2"))

# USER PROVIDES SOIL DEPTH
soil_depths_mm = c(100, 300, 600, 1000)



# USER PROVIDES DIRECTORY, TO OUTPUT
{
  subPath = paste0('/mnt/d/subsetDOMAINS/', name, '/geo_em.nc')
  
  geogrid.raster = wrfhydroSubsetter::make_empty_geogrid_raster(subPath)
  basin = dataRetrieval::findNLDI(comid = str_split(name, "_")[[1]][3], find = c("basin"))$basin
  
  #maskBas = fasterize::fasterize(st_transform(basin, st_crs(geogrid.raster)), geogrid.raster) # may be this needs to be changed?
  maskBas = raster::rasterize(st_transform(basin, st_crs(geogrid.raster)), geogrid.raster, getCover=TRUE)
  ## Mask of 1, NA depending on if cell is in basin
  maskMatrix = raster::as.matrix(maskBas) 
  
  dir_list = c(paste0('/mnt/d/OUTPUTS/', name, '/RUN_', run_timestamp, '/'),
               paste0('/mnt/d/OUTPUTS/', name, '/RUN_AP_', run_timestamp, '/'),
               paste0('/mnt/d/OUTPUTS/', name, '/RUN_SHUF_', run_timestamp, '/'),
               paste0('/mnt/d/OUTPUTS/', name, '/RUN_HUC12L_', run_timestamp, '/'),
               paste0('/mnt/d/OUTPUTS/', name, '/RUN_PERT_', run_timestamp, '/')
  )
}



#=======================================================
## FUNCTION WILL START HERE:
## INPUTS ARE FILELIST (dir), soil depths (vector, mm) and maskBas (raster):

total_wb_insepction = function(dir, maskMatrix, sfcrt=FALSE){
  
  ########
  # Read all files in directory, and pull out date and water year
  fileList = data.frame(
    files = list.files(dir, pattern = "RESTART", full.names = TRUE)) %>%
    dplyr::mutate(date  = lubridate::ymd_h(gsub("_DOMAIN1", "", gsub("^.*\\.","", files)), tz = "UTC"),
                  wy = ifelse(lubridate::month(date) >= 10,
                              lubridate::year(date) + 1,
                              lubridate::year(date)))
  
  # The current water budget formulation is the first minus the last timestep.
  # So, lets only process those...
  fileList = fileList[c(1,nrow(fileList)),]
  
  #Empty list to populate
  out = list()
  
  #=======================================================
    for(i in 1:nrow(fileList)){
    # open the current netcdf here
    tmp.nc = RNetCDF::open.nc(fileList$files[i])
    
    # Begin surface variable extraction
    surface_extract = lapply(seq_along(df$variableNames), function(i){
      # Read in variable, transpose, multiply by basin mask, and average
      # Result is average *unit* of *variable* in the basin at time of file
      mean(apply(t(RNetCDF::var.get.nc(tmp.nc, df$variableNames[i])), 2, rev) *
             maskMatrix, na.rm = TRUE)
    })
    
    # Open 3D SMC variable
    soil = RNetCDF::var.get.nc(tmp.nc, "SMC")
    
    # Begin surface variable extraction
    soil_extract = lapply(1:length(soil_depths_mm), function(i){
      # Read in soil layer, transpose, multiply by basin mask, average,
      # multiply by soil depth of layer
      # Result is average mm in each layer in the basin at time of file
      mean(apply(t(soil[,i,]), 2, rev) *
             maskMatrix, na.rm = TRUE) *
        soil_depths_mm[i]
    })
    
    # 1 row data.frame, 1 date, all variables
    # add these to the 'out' list
    out[[i]] = data.frame(t(c(unlist(surface_extract),
                              unlist(soil_extract)))) %>%
      setNames(c(df$variableNames,
                 paste0('SOIL_', 1:length(soil_depths_mm))))
    RNetCDF::close.nc(tmp.nc)
  }
  
  # bind all rows in out list, column bind to fileList
  o2          = cbind(fileList, dplyr::bind_rows(out))
  
  # Not used in water balance calc
  # # Build delta columns with a lag of 1
  # for(i in seq_along(df$variableNames)){
  #   var = df$variableNames[i]
  #   o2[[paste0("DEL_", var)]] = dplyr::lag(o2[[var]], n = 1)
  # }
  
  # Here we are summing all soil columns
  TOT_SOIL = rowSums(dplyr::select(o2, dplyr::starts_with("SOIL")))
  
  # For all states subtract the inital from the final
  # For canopy, ice and liquid water content are combined
  if (!sfcrt){
    WB_SFCRNOFF  = o2$SFCRUNOFF[nrow(o2)] -  o2$SFCRUNOFF[1]
  } else {
    hyd_fileList = data.frame(files = list.files(dir, pattern = "HYDRO_RST", full.names = TRUE)) %>%
      dplyr::mutate(date  = lubridate::ymd_hm(gsub("_DOMAIN1", "", gsub("^.*\\.","", files)), tz = "UTC"),
                    wy = ifelse(lubridate::month(date) >= 10, lubridate::year(date) + 1, lubridate::year(date)))

    hyd_fileList = hyd_fileList[c(1,nrow(hyd_fileList)),]; hyd_out = list();
    
    for(i in 1:nrow(hyd_fileList)){
      tmp.nc = RNetCDF::open.nc(hyd_fileList$files[i])
      HYD_QSTRMVOL = mean(apply(t(RNetCDF::var.get.nc(tmp.nc, "qstrmvolrt")), 2, rev) *
                            maskMatrix, na.rm = TRUE)
      
      hyd_out[[i]] = data.frame(t(c(unlist(HYD_QSTRMVOL)))) %>% setNames("qstrmvolrt")
      
      RNetCDF::close.nc(tmp.nc)
      }
    
    hyd_o2          = cbind(hyd_fileList, dplyr::bind_rows(hyd_out))
    WB_SFCRNOFF = hyd_o2$qstrmvolrt[nrow(hyd_o2)] - hyd_o2$qstrmvolrt[1]
    }
  
  
  wbDf = data.frame(
    # RAINFALL
    LSM_PRCP      = o2$ACCPRCP[nrow(o2)]  - o2$ACCPRCP[1],
    # Canopy evaporation
    LSM_ECAN      = o2$ACCECAN[nrow(o2)]  - o2$ACCECAN[1],
    # Canopy transporations
    LSM_ETRAN     = o2$ACCETRAN[nrow(o2)] - o2$ACCETRAN[1],
    # Soil evaporation
    LSM_EDIR      = o2$ACCEDIR[nrow(o2)]  - o2$ACCEDIR[1],
    # Snow water
    LSM_DELSWE    = o2$SNEQV[nrow(o2)]    - o2$SNEQV[1],
    # Canopy water/ice
    LSM_DELCANWAT = (o2$CANICE[nrow(o2)]  + o2$CANLIQ[nrow(o2)]) - (o2$CANICE[1] + o2$CANLIQ[1]),
    # Surface runoff
    #LSM_SFCRNOFF  = o2$SFCRUNOFF[nrow(o2)] -  o2$SFCRUNOFF[1],
    WB_SFCRNOFF,
    # Underground runoff
    LSM_UGDRNOFF  = o2$UDRUNOFF[nrow(o2)]  -  o2$UDRUNOFF[1],
    # Soil moisture
    LSM_DELSOILM = TOT_SOIL[length(TOT_SOIL)] - TOT_SOIL[1]
  ) %>% dplyr::mutate(
    # Since we are only using the LSM we set the WB to the LSM
    #WB_SFCRNOFF =  LSM_SFCRNOFF,
    WB_GWOUT    =  LSM_UGDRNOFF,
    # I dont get this but its not used so who cares
    #WB_DELGWSTOR = (o2$UDRUNOFF[nrow(o2)] - o2$UDRUNOFF[1]) - WB_GWOUT,
    ERROR = (LSM_PRCP -
      LSM_ECAN - LSM_ETRAN - LSM_EDIR -
      WB_SFCRNOFF - WB_GWOUT - LSM_DELSOILM - LSM_DELSWE - LSM_DELCANWAT)/LSM_PRCP,
    RUN_FRAC = (WB_SFCRNOFF + WB_GWOUT)/LSM_PRCP,
    EVAP_FRAC = (LSM_ECAN + LSM_ETRAN + LSM_EDIR)/LSM_PRCP,
    STOR_FRAC = (LSM_DELSOILM + LSM_DELSWE + LSM_DELCANWAT)/LSM_PRCP)
  
  
  #=======================================================
  ## To REPLICATE WB BARPLOT:
  #=======================================================
  
  
  
  to_plot  = data.frame(
    class = c(
      #"Canopy Evap",
      #"Transpiration",
      #"Surface Evap",
      "Total ET",
      "Surface Runoff",
      "Groundwater Outflow",
      "Change in Storage",
      "Error"),
    pcts = with(
      wbDf,
      c(#LSM_ECAN / LSM_PRCP * 100,
        #LSM_ETRAN / LSM_PRCP * 100,
        #LSM_EDIR / LSM_PRCP * 100,
        (LSM_ECAN + LSM_ETRAN + LSM_EDIR) / LSM_PRCP * 100,
        WB_SFCRNOFF / LSM_PRCP * 100,
        WB_GWOUT / LSM_PRCP * 100,
        (LSM_DELSOILM + LSM_DELSWE + LSM_DELCANWAT) / LSM_PRCP * 100,
        ERROR * 100
        )
      )
  ) %>%
    dplyr::mutate(labels = paste0(class, "\n", round(pcts, 1), "%"), basin = "1") #%>% dplyr::arrange(pcts)
  
  return(to_plot)
}
###### FUNCTION end
#=======================================================


nwm_to_plot = total_wb_insepction(dir_list[1], maskMatrix, sfcrt = TRUE)
ap_to_plot = total_wb_insepction(dir_list[2], maskMatrix, sfcrt = TRUE)
shuf_to_plot = total_wb_insepction(dir_list[3], maskMatrix, sfcrt = TRUE)
huc12l_to_plot = total_wb_insepction(dir_list[4], maskMatrix, sfcrt = TRUE)
pert_to_plot = total_wb_insepction(dir_list[5], maskMatrix, sfcrt = TRUE)

combined = data.frame(
  scheme = c("NWM",
             "AP",
             "SHUF",
             "HUC12L",
             "PERT"),
  ET = c(nwm_to_plot[1,2], 
         ap_to_plot[1,2],
         shuf_to_plot[1,2],
         huc12l_to_plot[1,2],
         pert_to_plot[1,2]),
  SFC_Q = c(nwm_to_plot[2,2], 
         ap_to_plot[2,2],
         shuf_to_plot[2,2],
         huc12l_to_plot[2,2],
         pert_to_plot[2,2]),
  GW_Q = c(nwm_to_plot[3,2], 
            ap_to_plot[3,2],
            shuf_to_plot[3,2],
            huc12l_to_plot[3,2],
            pert_to_plot[3,2]),
  dS = c(nwm_to_plot[4,2], 
            ap_to_plot[4,2],
            shuf_to_plot[4,2],
            huc12l_to_plot[4,2],
            pert_to_plot[4,2]),
  Error = c(nwm_to_plot[5,2], 
         ap_to_plot[5,2],
         shuf_to_plot[5,2],
         huc12l_to_plot[5,2],
         pert_to_plot[5,2])
)
combined_plot <- reshape2::melt(combined, id.vars = "scheme")


ggplot(combined_plot, aes(scheme, value, fill = variable)) +
  geom_col(position = "dodge")





ggplot(data = huc12l_to_plot, aes(y = pcts, x = basin, fill = labels)) +
  geom_col(position="fill") +
  theme_bw() +
  scale_fill_viridis_d() +
  labs(fill = "", x = "", y = "% of Water Budget")



####====================================================

#### Loop for saving files

locs = data.frame(comids = c(23762661, 1631587, 5894384, 19389766, 191739, 5781369),
                  siteID = c('14190500', '14190500', '01616500', '03118500' ,'6709000' , '08159000'))

namelist = c()
for (i in 1:nrow(locs)){
  state = st_transform(AOI::aoi_get(state = "all", county = "all"), 4269)[st_transform(findNLDI(comid = locs$comids[i]), crs = 4269),]
  name = gsub(" ", "-", tolower(paste(state$name, state$state_name, locs$comids[i], locs$siteID[i], sep = "_")))
  namelist = c(namelist, name)
  rm(name)
}

scheme_list = c('NWM', 'AP', 'SHUF', 'HUC12L', 'PERT')

#######################
tic()
for (i in 1:nrow(locs)) {
  comid = locs$comids[i]
  name = namelist[i]
  print(paste0("Processing: ", name))
  subPath = paste0('/mnt/d/subsetDOMAINS/', name, '/geo_em.nc')
  
  geogrid.raster = wrfhydroSubsetter::make_empty_geogrid_raster(subPath)
  basin = dataRetrieval::findNLDI(comid = comid, find = c("basin"))$basin
  maskBas = raster::rasterize(st_transform(basin, st_crs(geogrid.raster)), geogrid.raster, getCover=TRUE)
  
  ## Mask of 1, NA depending on if cell is in basin
  maskMatrix = raster::as.matrix(maskBas) 
  
  NWM_OUT_path     = list.files(paste0("/mnt/d/OUTPUTS/", name), "RUN_", full.names = TRUE)[1]
  AP_OUT_path      = list.files(paste0("/mnt/d/OUTPUTS/", name), "RUN_", full.names = TRUE)[2]
  HUC12L_OUT_path  = list.files(paste0("/mnt/d/OUTPUTS/", name), "RUN_", full.names = TRUE)[3]
  SHUF_OUT_path    = list.files(paste0("/mnt/d/OUTPUTS/", name), "RUN_", full.names = TRUE)[5]
  PERT_OUT_path    = list.files(paste0("/mnt/d/OUTPUTS/", name), "RUN_", full.names = TRUE)[4]
  
  
  dir_list = c(NWM_OUT_path,   
               AP_OUT_path,    
               HUC12L_OUT_path,
               SHUF_OUT_path,
               PERT_OUT_path
               )
  
  nwm_to_plot = total_wb_insepction(dir_list[1], maskMatrix, sfcrt = TRUE)
  ap_to_plot = total_wb_insepction(dir_list[2], maskMatrix, sfcrt = TRUE)
  huc12l_to_plot = total_wb_insepction(dir_list[3], maskMatrix, sfcrt = TRUE)
  shuf_to_plot = total_wb_insepction(dir_list[4], maskMatrix, sfcrt = TRUE)
  pert_to_plot = total_wb_insepction(dir_list[5], maskMatrix, sfcrt = TRUE)
  
  combined = data.frame(
    scheme = c("NWM", 
               "AP",
               "HUC12L",
               "SHUF",
               "PERT"
               ),
    ET = c(nwm_to_plot[1,2], 
           ap_to_plot[1,2],
           huc12l_to_plot[1,2],
           shuf_to_plot[1,2],
           pert_to_plot[1,2]
           ),
    SFC_Q = c(nwm_to_plot[2,2], 
              ap_to_plot[2,2],
              huc12l_to_plot[2,2],
              shuf_to_plot[2,2],
              pert_to_plot[2,2]
              ),
    GW_Q = c(nwm_to_plot[3,2], 
             ap_to_plot[3,2],
             huc12l_to_plot[3,2],
             shuf_to_plot[3,2],
             pert_to_plot[3,2]
             ),
    dS = c(nwm_to_plot[4,2], 
           ap_to_plot[4,2],
           huc12l_to_plot[4,2],
           shuf_to_plot[4,2],
           pert_to_plot[4,2]
           ),
    Error = c(nwm_to_plot[5,2], 
              ap_to_plot[5,2],
              huc12l_to_plot[5,2],
              shuf_to_plot[5,2],
              pert_to_plot[5,2]
              )
  )
  
  combined_plot <- reshape2::melt(combined, id.vars = "scheme")
  
  write.csv(combined   , paste0("/mnt/d/ANALYSIS/", name, "/WB_combined.csv"), row.names=T)
  
}
toc()
