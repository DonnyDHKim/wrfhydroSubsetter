subset_sequence = function(comid, siteID, FULLDOMAINDIR, outDir, nlcdDir = NA, methodList = NA){
  
  # Defining the area * This could go out and just stay in the loop.
  state = st_transform(AOI::aoi_get(state = "all", county = "all"), 4269)[st_transform(findNLDI(comid = comid), crs = 4269),]
  name = gsub(" ", "-", tolower(paste(state$name, state$state_name, comid, siteID, sep = "_")))
  basin = findNLDI(comid = comid, find = c("basin"))$basin
  
  # Subsetting NetCDF DOMAIN files
  subset_files = paste0(outDir, "/", name , "/")
  subset_wrf_hydro_domain_mod(AOI = basin,  domain_files = FULLDOMAINDIR,  outDir = subset_files, config = 'LongRange')
  
  
  # Resampling sequence
  if (exists("methodList")==T & exists("nlcdDir")==T) {
    
    # Reading in look-up table. This is somewhat problematic at this point.
    data = readxl::read_excel('/mnt/d/GitHub/wrfhydroSubsetter/nwm_land_use_mappings.xlsx') %>%
      dplyr::select(nlcd = Class, 
                    description = Description,
                    nwm = NWM)
    
    # Cropping nlcd to geogrid resolution
    geo  = list.files(subset_files, "geo_em.nc", full.names = TRUE)
    output_geo = spex::qm_rasterToPolygons(wrfhydroSubsetter::make_empty_geogrid_raster(geo))
    o2 = output_geo %>%
      AOI::bbox_get() %>%
      sf::st_transform(5070) %>%
      sf::st_buffer(40)
    nlcdObj = raster(nlcdDir)
    methodList = methodList
    nlcdCrop = raster::crop(nlcdObj, o2)
    
    # Reclassifying nlcd LULC class into nwm class *DK: is it necessary to do this at this phase tho?
    nlcd_nwm = raster::reclassify(nlcdCrop, rcl = dplyr::select(data,from = nlcd, to = nwm ))
    
    # Resampling sequence.
    for(j in 1:length(methodList)){
      new_lu = resampleDataMod(input = nlcd_nwm, output = output_geo, method = methodList[j])
      out_file = paste0(subset_files, "geogrid_", methodList[j], "_resample.nc")
      
      file.copy(geo, out_file, overwrite = TRUE)
      nc = ncdf4::nc_open(out_file, write = TRUE)
      ncdf4::ncvar_put(nc, "LU_INDEX",
                       vals = as.vector(t(apply(as.matrix(new_lu), 2, rev))))
      ncdf4::nc_close(nc)
    }
  }
  
  
}