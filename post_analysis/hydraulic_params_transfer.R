#' @title Hydraulic soilproperties transfer
#' @description Copy hydraulic soilproperties from NWM soilproperties.nc
#' @param comid comid of outlet reach
#' @param siteID siteID for NWIS
#' @param FULLDOMAINDIR the value representing "NO DATA", default is NA
#' @param outDir output directory (e.g. /home/subsetDOMAINS/)
#' @param nlcdDir directory for nlcd dataset
#' @param methodList resampling algorithm to use for resampling
#' @importFrom sf st_transform st_buffer
#' @importFrom AOI aoi_get bbox_get
#' @importFrom dataRetrieval findNLDI
#' @importFrom spex qm_rasterToPolygons
#' @importFrom raster crop
#' @importFrom ncdf4 nc_open ncvar_put nc_close
#' @importFrom resampleDataMod.R resampleDataMod

{
  library(sf)
  library(RNetCDF)
  library(dplyr)
  library(wrfhydroSubsetter)
}

locs = data.frame(comids = c(23762661, 1631587, 5894384, 19389766, 191739, 5781369),
                  siteID = c('14190500', '14190500', '01616500', '03118500' ,'6709000' , '08159000'))

namelist = c()
for (i in 1:nrow(locs)){
  state = st_transform(AOI::aoi_get(state = "all", county = "all"), 4269)[st_transform(findNLDI(comid = locs$comids[i]), crs = 4269),]
  name = gsub(" ", "-", tolower(paste(state$name, state$state_name, locs$comids[i], locs$siteID[i], sep = "_")))
  namelist = c(namelist, name)
  rm(name)
}

varlist = c('LR', 'AP', 'SHUF', 'HUC12L', 'PERT')

#name = 'berkeley_west-virginia_5894384_01616500'
#var = 'AP'
range = 'LongRange'

for (i in namelist){
  for (j in varlist){
    hyd_params_transfer(i, j, range)
  }
}


hyd_params_transfer = function(name, var, range){
  
  nwm_soil = paste0('/mnt/d/subsetDOMAINS/', name, "/soilproperties_", range, ".nc")
  print(paste0("Hydraulic paramters will be copied from NWM soilproperties in: ", name, ", /", var))
  
  # Open NC
  nc = RNetCDF::open.nc(nwm_soil)
  slopeVal = RNetCDF::var.get.nc(nc, 'slope') %>% as.data.frame(xy=TRUE)
  refkdtVal = RNetCDF::var.get.nc(nc, 'refkdt') %>% as.data.frame(xy=TRUE)
  #refdkVal = RNetCDF::var.get.nc(nc, 'refdk') %>% as.data.frame(xy=TRUE)

  dir.create(paste0('/mnt/d/subsetDOMAINS/', name, '/archive/'), showWarnings = FALSE)
  
  file.copy(paste0('/mnt/d/subsetDOMAINS/', name, "/soil_properties_",var,".nc"), 
            paste0('/mnt/d/subsetDOMAINS/', name, "/archive/soil_properties_",var,".nc"), overwrite = TRUE)
  
  #out_file = paste0('/mnt/d/subsetDOMAINS/', name, "/soil_properties_", var,"_test.nc")
  #file.copy(paste0('/mnt/d/subsetDOMAINS/', name, "/soil_properties_",var,".nc"), out_file, overwrite = TRUE)
  out_file = paste0('/mnt/d/subsetDOMAINS/', name, "/soil_properties_", var,".nc")
  
  nc = ncdf4::nc_open(out_file, write = TRUE)
  ncdf4::ncvar_put(nc, "slope", vals = as.vector(as.matrix(slopeVal)))
  ncdf4::ncvar_put(nc, "refkdt", vals = as.vector(as.matrix(refkdtVal)))
  #ncdf4::ncvar_put(nc, "refdk", vals = as.vector(as.matrix(refdkVal)))
  ncdf4::nc_close(nc)
  print(paste0(name, "_", var, ": Hydraulic params in soilproperties.nc transferred"))
}
