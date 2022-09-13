{
  library(sf)
  library(lubridate)
  library(RNetCDF)
  library(dplyr)
  library(ggplot2)
  library(rwrfhydro) # DO I need it?
}

#devtools::install_github("NCAR/rwrfhydro")


rtlinkFile <- ReadRouteLink("/mnt/d/subsetDOMAINS/berkeley_west-virginia_5894384_01616500/RouteLink.nc")
NWM_OUT_path <- "/mnt/d/OUTPUTS/berkeley_west-virginia_5894384_01616500/RUN_20220912_1401"
AP_OUT_path <- "/mnt/d/OUTPUTS/berkeley_west-virginia_5894384_01616500/RUN_AP_20220912_1401"
SHUF_OUT_path <- "/mnt/d/OUTPUTS/berkeley_west-virginia_5894384_01616500/RUN_SHUF_20220912_1401"
HUC12L_OUT_path <- "/mnt/d/OUTPUTS/berkeley_west-virginia_5894384_01616500/RUN_HUC12L_20220912_1401"

head(rtlinkFile)

#list gages that are in the route_link file
listgages <- rtlinkFile$site_no
listgages <- na.omit(as.numeric(listgages))
listgages <- paste0("0",listgages)
cat("The USGS gages found in the Route_Link.nc file are:", listgages)

#pull out the rows in route_link with gages and print metadata
gage_comid = arrange(rtlinkFile, desc(site_no))
gage_comid = gage_comid[1:8,]

gage_USGS = "01616500"


# find the NHDPlus COMID that corresponds to the above USGS gage ID within the Route_Link.nc file
correspondCOMID <- rtlinkFile$link[which(rtlinkFile$site_no == gage_USGS)]
correspondCOMID <- 5894384

# check to make sure that USGS gage is in fact contained in the current Route_Link.nc file
if (length(GAGEindex)==0) {
  message("!!! ERROR: USGS Gage not found in Route_Link.nc file !!!")
} else {
  message("GAGE FOUND! COMID '",correspondCOMID,"' corresponds to USGS gage '", gage_USGS,"'")
}

#channelRead_List <- future_pmap(OUT_path_list, ReadChrtout(path=~.x, idList = 5894384, parallel = TRUE))

OUT_path_list = c(NWM_OUT_path, AP_OUT_path, SHUF_OUT_path,HUC12L_OUT_path)

#if (!require("furrr")) install.packages("furrr"); library("furrr") # furrr package is essential.
library(parallel); library(doParallel); library(foreach); library(tictoc)
# Make cluster and register it.
cl <- parallel::makeCluster(detectCores()-4) # don't set makeCluster(n) too high on laptops.

registerDoParallel(cl)
#registerDoParallel(cores=detectCores()-4)



#future::plan(multisession, workers=4)
#future_map(OUT_path_list, print((~.x)))



{
  tic()
  channelRead_NWM = ReadChrtout(NWM_OUT_path, correspondCOMID, parallel=TRUE) #%>% select(q_cms) %>% as.numeric()
  channelRead_AP = ReadChrtout(AP_OUT_path, correspondCOMID, parallel=TRUE) #%>% select(q_cms) %>% as.numeric()
  channelRead_SHUF = ReadChrtout(SHUF_OUT_path, correspondCOMID, parallel=TRUE) #%>% select(q_cms) %>% as.numeric()
  channelRead_HUC12L = ReadChrtout(HUC12L_OUT_path, correspondCOMID, parallel=TRUE) #%>% select(q_cms) %>% as.numeric()
  toc()
}


#future:::ClusterRegistry("stop") # very important to kill the nodes. If you don't kill it, it will just eat up computer's memory.
parallel::stopCluster(cl) # Important to do this.


channelRead_NWM$q_cms <- as.numeric(channelRead_NWM$q_cms)
message("Reading discharge output from Channel Routing files complete")
print(names(channelRead.NWM))



fdc_df_builder = function(dat, nam) {
  dat <- sort(dat, decreasing = T)
  df  <- data.frame(x = 100/length(dat) * 1:length(dat), y = dat)
  colnames(df) = c(paste0(nam, ".x"), paste0(nam, ".y"))
  return(df)
}

fdc_nwm = fdc_df_builder(channelRead_NWM$q_cms, "NWM")
fdc_AP = fdc_df_builder(channelRead_AP$q_cms, "AP")
fdc_SHUF = fdc_df_builder(channelRead_SHUF$q_cms, "SHUF")
fdc_HUC12L = fdc_df_builder(channelRead_HUC12L$q_cms, "HUC12L")

fdc_merged = cbind(fdc_nwm, fdc_AP, fdc_SHUF, fdc_HUC12L)

fdc_merged_plot = ggplot(fdc_merged)+
  geom_line(aes(x=NWM.x, y=NWM.y, color = "NWM"), cex=0.5) +
  geom_line(aes(x=AP.x, y=AP.y, color = "AP"), cex=0.5) +
  geom_line(aes(x=SHUF.x, y=SHUF.y, color = "SHUF"), cex=0.5) +
  geom_line(aes(x=HUC12L.x, y=HUC12L.y, color = "HUC12L"), cex=0.5) +
  scale_y_continuous(limits = c(min= 1, max=200),  
                     breaks = c(5, 10, 50, 100, 200), trans="log2")+ 
  scale_x_continuous(breaks = c(0, 20, 40, 60, 80, 100), limits = c(min=0, max=100),expand=c(0.02,0.02))+
  theme(#axis.line=element_blank(),
    axis.line = element_line(size = 0.5, colour = "black", linetype=1),
    plot.margin=unit(c(0,0.2,0,0),"cm"),
    #axis.text.x=element_blank(),
    #axis.text.y=element_blank(),
    #axis.ticks=element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    #legend.position="none",
    panel.grid.minor=element_blank(),
    #panel.grid=element_blank(),
    axis.ticks = element_line(size = 1, color="black"),
    axis.ticks.length=unit(-0.1, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.2,0.2,0.2,0.2), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.2,0.2,0.2,0.2), "cm")),
    #legend.text=element_text(size=10),
    #legend.justification=c(1.1,-0.1),
    #legend.position="none",
    legend.title=element_blank(),
    panel.background = element_blank())
  
fdc_merged_plot

{
  fdc12plot =
    ggplot(fdc12)+
    geom_ribbon(aes(x=opt1.x, ymin=lo1.y, ymax=hi1.y), fill = fdc12.colors[3], alpha=0.15)+
    geom_ribbon(aes(x=opt2.x, ymin=lo2.y, ymax=hi2.y), fill = fdc12.colors[4], alpha=0.15)+
    geom_line(aes(x=opt1.x, y=opt1.y, color = "opt1"), cex=0.5) +
    geom_line(aes(x=opt2.x, y=opt2.y, color = "opt2"), cex=0.5) +
    geom_line(aes(x=NWM.x, y=NWM.y, color = "NWM"),linetype=6, cex=0.5) +
    geom_line(aes(x=OBS.x, y=OBS.y, color = "OBS"), linetype=1,cex=1) +
    scale_y_continuous(limits = c(min= 0.01, max=40),  
                       breaks = c(0.01, 0.1, 1, 5, 10, 20, 40), trans="log2")+
    scale_x_continuous(breaks = c(0, 20, 40, 60, 80, 100), limits = c(min=0, max=100),expand=c(0.02,0.02))+
    theme(#axis.line=element_blank(),
      axis.line = element_line(size = 0.5, colour = "black", linetype=1),
      plot.margin=unit(c(0,0.2,0,0),"cm"),
      #axis.text.x=element_blank(),
      #axis.text.y=element_blank(),
      #axis.ticks=element_blank(),
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      #legend.position="none",
      panel.grid.minor=element_blank(),
      #panel.grid=element_blank(),
      axis.ticks = element_line(size = 1, color="black"),
      axis.ticks.length=unit(-0.1, "cm"), 
      axis.text.x = element_text(margin=unit(c(0.2,0.2,0.2,0.2), "cm")), 
      axis.text.y = element_text(margin=unit(c(0.2,0.2,0.2,0.2), "cm")),
      #legend.text=element_text(size=10),
      #legend.justification=c(1.1,-0.1),
      legend.position="none",
      legend.title=element_blank(),
      panel.background = element_blank())+
    scale_color_manual(values=fdc12.colors, breaks=c("OBS", "NWM","opt1", "opt2"))
}
fdc12plot









