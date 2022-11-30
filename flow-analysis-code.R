{
  library(sf)
  library(lubridate)
  library(RNetCDF)
  library(dplyr)
  library(ggplot2)
  library(rwrfhydro) # DO I need it?
}

#devtools::install_github("NCAR/rwrfhydro")

locs = data.frame(comids = c(191739, 1631587, 5894384, 5781369, 19389766, 23762661),
                  siteID = c('6709000', '08173000', '01616500', '08159000' ,'03118500' , '14190500'))

namelist = c()
for (i in 1:nrow(locs)){
  state = st_transform(AOI::aoi_get(state = "all", county = "all"), 4269)[st_transform(findNLDI(comid = locs$comids[i]), crs = 4269),]
  name = gsub(" ", "-", tolower(paste(state$name, state$state_name, locs$comids[i], locs$siteID[i], sep = "_")))
  namelist = c(namelist, name)
  rm(name)
}

scheme_list = c('NWM', 'AP', 'SHUF', 'HUC12L', 'PERT')


list.files("/mnt/d/OUTPUTS/douglas_colorado_191739_6709000/", "RUN_", full.names = TRUE)
##### Loop to save csv files

{
  library(parallel); library(doParallel); library(foreach); library(tictoc)
  cl <- parallel::makeCluster(detectCores()-4) # don't set makeCluster(n) too high on laptops.
  registerDoParallel(cl)
  
  tic()
  for (i in 1:nrow(locs)) {
    name = namelist[i]
    rtlinkFile <- ReadRouteLink(paste0("/mnt/d/subsetDOMAINS/", name, "/RouteLink.nc"))
    
    NWM_OUT_path     = list.files(paste0("/mnt/d/OUTPUTS/", name), "RUN_", full.names = TRUE)[1]
    AP_OUT_path      = list.files(paste0("/mnt/d/OUTPUTS/", name), "RUN_", full.names = TRUE)[2]
    SHUF_OUT_path    = list.files(paste0("/mnt/d/OUTPUTS/", name), "RUN_", full.names = TRUE)[5]
    HUC12L_OUT_path  = list.files(paste0("/mnt/d/OUTPUTS/", name), "RUN_", full.names = TRUE)[3]
    PERT_OUT_path    = list.files(paste0("/mnt/d/OUTPUTS/", name), "RUN_", full.names = TRUE)[4]
    
    correspondCOMID = locs$comids[i]
    
    channelRead_NWM    = ReadChrtout(NWM_OUT_path, correspondCOMID, parallel=TRUE) #%>% select(q_cms) %>% as.numeric()
    channelRead_AP     = ReadChrtout(AP_OUT_path, correspondCOMID, parallel=TRUE) #%>% select(q_cms) %>% as.numeric()
    channelRead_SHUF   = ReadChrtout(SHUF_OUT_path, correspondCOMID, parallel=TRUE) #%>% select(q_cms) %>% as.numeric()
    channelRead_HUC12L = ReadChrtout(HUC12L_OUT_path, correspondCOMID, parallel=TRUE) #%>% select(q_cms) %>% as.numeric()
    channelRead_PERT   = ReadChrtout(PERT_OUT_path, correspondCOMID, parallel=TRUE) #%>% select(q_cms) %>% as.numeric()
    
    dir.create(paste0("/mnt/d/ANALYSIS/", name))
    write.csv(channelRead_NWM   , paste0("/mnt/d/ANALYSIS/", name, "/RawFlow_NWM.csv"), row.names=T)
    write.csv(channelRead_AP    , paste0("/mnt/d/ANALYSIS/", name, "/RawFlow_AP.csv"), row.names=T)
    write.csv(channelRead_SHUF  , paste0("/mnt/d/ANALYSIS/", name, "/RawFlow_SHUF.csv"), row.names=T)
    write.csv(channelRead_HUC12L, paste0("/mnt/d/ANALYSIS/", name, "/RawFlow_HUC12L.csv"), row.names=T)
    write.csv(channelRead_PERT  , paste0("/mnt/d/ANALYSIS/", name, "/RawFlow_PERT.csv"), row.names=T)
  }
  toc()
  
  parallel::stopCluster(cl) # Important to do this.
}




###### below is original
name = 'douglas_colorado_191739_6709000'
run_timestamp = '20220923_1945'
scheme_list = c('NWM', 'AP', 'SHUF', 'HUC12L', 'PERT')


rtlinkFile <- ReadRouteLink(paste0("/mnt/d/subsetDOMAINS/", name, "/RouteLink.nc"))
NWM_OUT_path   <-  paste0("/mnt/d/OUTPUTS/", name, "/RUN_", run_timestamp)
AP_OUT_path   <-   paste0("/mnt/d/OUTPUTS/", name, "/RUN_AP_", run_timestamp)
SHUF_OUT_path   <- paste0("/mnt/d/OUTPUTS/", name, "/RUN_SHUF_", run_timestamp)
HUC12L_OUT_path <- paste0("/mnt/d/OUTPUTS/", name, "/RUN_HUC12L_", run_timestamp)
PERT_OUT_path <- paste0("/mnt/d/OUTPUTS/", name, "/RUN_PERT_", run_timestamp)


head(rtlinkFile)

#list gages that are in the route_link file
listgages <- rtlinkFile$site_no
listgages <- na.omit(as.numeric(listgages))
listgages <- paste0("0",listgages)
cat("The USGS gages found in the Route_Link.nc file are:", listgages)

#pull out the rows in route_link with gages and print metadata
gage_comid = arrange(rtlinkFile, desc(site_no))
gage_comid = gage_comid[1:8,]

gage_USGS = "06709000"


# find the NHDPlus COMID that corresponds to the above USGS gage ID within the Route_Link.nc file
correspondCOMID <- rtlinkFile$link[which(rtlinkFile$site_no == gage_USGS)]
correspondCOMID = str_split(name, "_")[[1]][3]
#correspondCOMID <- 191739


#if (!require("furrr")) install.packages("furrr"); library("furrr") # furrr package is essential.
library(parallel); library(doParallel); library(foreach); library(tictoc)
# Make cluster and register it.
cl <- parallel::makeCluster(detectCores()-4) # don't set makeCluster(n) too high on laptops.

registerDoParallel(cl)

#future::plan(multisession, workers=4)
#future_map(OUT_path_list, print((~.x)))

{
  tic()
  channelRead_NWM = ReadChrtout(NWM_OUT_path, correspondCOMID, parallel=TRUE) #%>% select(q_cms) %>% as.numeric()
  channelRead_AP = ReadChrtout(AP_OUT_path, correspondCOMID, parallel=TRUE) #%>% select(q_cms) %>% as.numeric()
  channelRead_SHUF = ReadChrtout(SHUF_OUT_path, correspondCOMID, parallel=TRUE) #%>% select(q_cms) %>% as.numeric()
  channelRead_HUC12L = ReadChrtout(HUC12L_OUT_path, correspondCOMID, parallel=TRUE) #%>% select(q_cms) %>% as.numeric()
  channelRead_PERT = ReadChrtout(PERT_OUT_path, correspondCOMID, parallel=TRUE) #%>% select(q_cms) %>% as.numeric()
  
  toc()
}


#future:::ClusterRegistry("stop") # very important to kill the nodes. If you don't kill it, it will just eat up computer's memory.
parallel::stopCluster(cl) # Important to do this.


#channelRead_NWM$q_cms <- as.numeric(channelRead_NWM$q_cms)
message("Reading discharge output from Channel Routing files complete")
print(names(channelRead.NWM))

# Well why don't I improt the USGS gage streamflow data
# loop

locs$siteID = zeroPad(locs$siteID, 8)

library(xts)

{
  tic()
  for (i in 2:2) {
    name = namelist[i]
    gageid = locs$siteID[i]
    
    USGS_flow = readNWISuv(siteNumbers = gageid, parameterCd = "00060", "2000-09-30", "2020-09-30") %>% 
      renameNWISColumns() #%>% 
    
    USGS_flow2 = USGS_flow %>% mutate(Flow_Inst = .[[4]] * 0.028316846592)

    #USGS_flow2 = USGS_flow2 %>% 
    #  mutate(WY = ifelse(lubridate::month(dateTime) >= 10, lubridate::year(dateTime)+1, lubridate::year(dateTime)))
    #
    #USGS_flow2 = USGS_flow2 %>% 
    #  filter(WY > lubridate::year(USGS_flow2$dateTime[1]))
    
    
    if (lubridate::month(USGS_flow2$dateTime[1]) < 10){
      USGS_flow2 = USGS_flow2 %>% 
        mutate(WY = ifelse(lubridate::month(dateTime) >= 10, lubridate::year(dateTime)+1, lubridate::year(dateTime))) %>% 
        filter(WY > lubridate::year(USGS_flow2$dateTime[1]))
      USGS_flow_trim = USGS_flow2 %>% 
        filter(dateTime >= as.POSIXct("2000-09-30 23:59:00", tz = "utc") & dateTime <= as.POSIXct("2020-10-01 00:01:00", tz = "utc"))
    } else {    
      USGS_flow_trim = USGS_flow2 %>% 
      filter(dateTime >= as.POSIXct("2000-09-30 23:59:00", tz = "utc") & dateTime <= as.POSIXct("2020-10-01 00:01:00", tz = "utc"))
      }
    
    

    
    USGS_flow_hrly <- aggregate(USGS_flow_trim["Flow_Inst"], 
                     list(hour=cut(as.POSIXct(USGS_flow_trim$dateTime), "6 hour")),
                     mean)
    #USGS_flow_hrly = USGS_flow_hrly[18:173090,]
    
    #flow.xts = xts(USGS_flow_hrly$Flow_Inst, as.POSIXct(USGS_flow_hrly$hour))
    #ends = as.integer(endpoints(flow.xts, on='hours', k=6)+1)[1:length(ends)-1]
    #flow.xts2 = period.apply(flow.xts, ends, mean)
    #
    #USGS_flow <- readNWISdv(siteNumbers = gageid, parameterCd = "00060", "2000-10-01", "2020-09-30") %>% 
    #  renameNWISColumns() %>% 
    #  mutate(Flow = Flow * 0.028316846592)

    
    print(paste0(name, ": Compete, without problem."))
    dir.create(paste0("/mnt/d/ANALYSIS/", name))
    write.csv(USGS_flow_hrly   , paste0("/mnt/d/ANALYSIS/", name, "/USGS_6hrly.csv"), row.names=T)
  }
  toc()
}

USGS_flow <- readNWISdv(siteNumbers = gage_USGS, parameterCd = "00060", "2000-10-01", "2020-09-30") %>% 
  renameNWISColumns() %>% 
  mutate(Flow = Flow * 0.028316846592)







make_daily = function(df){
  x = df %>% 
    mutate(date = floor_date(as_date(POSIXct))) %>%
    group_by(date) %>% 
    summarize (q_cms = mean(q_cms)) %>% 
    filter(date >= as.POSIXct("2000-09-30") & date <= as.POSIXct("2020-09-30")) %>% 
    dplyr::select(q_cms) %>% as.matrix() %>% as.numeric()
  return(x)
}

test = channelRead_AP %>% make_daily()
test2 = test %>% dplyr::select(q_cms)

fdc_df_builder = function(dat, nam) {
  dat <- sort(dat, decreasing = T)
  df  <- data.frame(x = 100/length(dat) * 1:length(dat), y = dat)
  colnames(df) = c(paste0(nam, ".x"), paste0(nam, ".y"))
  return(df)
}


fdc_USGS = fdc_df_builder(USGS_flow$Flow, "USGS")
fdc_nwm = fdc_df_builder(channelRead_NWM %>% make_daily(), "NWM")
fdc_AP = fdc_df_builder(channelRead_AP %>% make_daily(), "AP")
fdc_SHUF = fdc_df_builder(channelRead_SHUF %>% make_daily(), "SHUF")
fdc_HUC12L = fdc_df_builder(channelRead_HUC12L %>% make_daily(), "HUC12L")
fdc_PERT = fdc_df_builder(channelRead_PERT %>% make_daily(), "PERT")

fdc_merged = cbind(fdc_USGS, fdc_nwm, fdc_AP, fdc_SHUF, fdc_HUC12L, fdc_PERT)
max_val = plyr::round_any(max(fdc_merged[1,]), 10, f=ceiling)

fdc_merged_plot = ggplot(fdc_merged)+
  geom_line(aes(x=USGS.x, y=USGS.y, color = "USGS"), cex=0.5) +
  geom_line(aes(x=NWM.x, y=NWM.y, color = "NWM"), cex=0.5) +
  geom_line(aes(x=AP.x, y=AP.y, color = "AP"), cex=0.5) +
  geom_line(aes(x=SHUF.x, y=SHUF.y, color = "SHUF"), cex=0.5) +
  geom_line(aes(x=HUC12L.x, y=HUC12L.y, color = "HUC12L"), cex=0.5) +
  geom_line(aes(x=PERT.x, y=PERT.y, color = "PERT"), cex=0.5) +
  scale_y_continuous(limits = c(min= 0.1, max=max_val), trans='log10',
                     breaks = c(1, max_val*0.1, max_val*0.2, max_val*0.5, max_val))+ 
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









#####
#####
##### FDC slope
# https://github.com/IvanHeriver/hsa-R-package/blob/master/R/flow_duration_curve.R
fdc_slope <- function(Q, probs = c(0.33, 0.66)) {
  Qp <- quantile(Q, probs = probs, na.rm = TRUE, names = FALSE)
  - (log10(Qp[1L]) - log10(Qp[2L])) / diff(probs)
}

fdc_midslop_df = cbind(fdc_slope(USGS_flow$Flow),
                        fdc_slope(channelRead_NWM$q_cms),
                        fdc_slope(channelRead_AP$q_cms),
                        fdc_slope(channelRead_SHUF$q_cms),
                        fdc_slope(channelRead_HUC12L$q_cms),
                        fdc_slope(channelRead_PERT$q_cms)
                        )
colnames(fdc_midslop_df) = c("USGS", scheme_list)




#####
#####
##### 7Q10
# https://vt-hydroinformatics.github.io/lfas.html
library(dplyr)
library(lubridate)
library(zoo)
library(moments)

Q7_calc = function(flow, Xday = 7, YrecInt = 10){
  require(moments); library(zoo); library(dplyr); library(lubridate);
  
  Xday = Xday;   YrecInt = YrecInt;
  
  x = flow %>% 
    mutate(date = floor_date(as_date(POSIXct))) %>%
    group_by(date) %>% 
    summarize (daily_q_cms = mean(q_cms)) %>% 
    mutate(WY = ifelse(lubridate::month(date) >= 10, lubridate::year(date) + 1, lubridate::year(date)))
  
  #X day rolling mean, don't fill the ends of the timeseries,
  #don't ignore NAs, use a backward-looking window (right align)
  x2 = x %>% mutate(xdaymean = rollmean(daily_q_cms, Xday, fill = NA, na.rm = F, align = "right"))
  
  QyearlyMins = x2 %>% 
    group_by(WY) %>%
    summarize(minQ = min(xdaymean, na.rm = T), lenDat = length(daily_q_cms), lenNAs = sum(is.na(xdaymean))) %>%
    filter(lenDat > 328 & lenNAs / lenDat < 0.1) %>% 
    mutate(rank = rank(minQ, ties.method = "first")) %>%
    mutate(ReturnInterval = (length(rank) + 1)/rank) %>%
    mutate(ExceedProb = 1 / ReturnInterval)
  
  #Measures of the distribution
  Xbar <- mean(log10(QyearlyMins$minQ))
  S    <- sd(log10(QyearlyMins$minQ))
  g    <- skewness(log10(QyearlyMins$minQ))
  
  QyearlyMins = QyearlyMins %>% 
    mutate(z = 4.91 * ((1 / ReturnInterval) ^ 0.14 - (1 - 1 / ReturnInterval) ^ 0.14)) %>%
    mutate(K = (2 / g) * (((1 + (g * z) / 6 - (g ^ 2) / 36) ^ 3) - 1) ) %>%
    mutate(Qfit = 10^(Xbar + (K * S)))
  
  #xQy ei: 7Q10
  y = YrecInt
  
  #Find these values based on established relationships
  z    <- 4.91 * ((1 / y) ^ 0.14 - (1 - 1 / y) ^ 0.14)
  K    <- (2 / g) * (((1 + (g * z) / 6 - (g ^ 2) / 36) ^ 3) - 1) 
  PearsonxQy <- 10^(Xbar + K * S)
  
  
  return(PearsonxQy)
}

Q7_10_df = cbind(Q7_calc(channelRead_NWM),
                         Q7_calc(channelRead_AP),
                         Q7_calc(channelRead_SHUF),
                         Q7_calc(channelRead_HUC12L),
                         Q7_calc(channelRead_PERT))
colnames(Q7_10_df) = scheme_list



#####
#####
##### Annual Peak
Annual_peak_calc = function(flow){
  
  x = flow %>% 
    mutate(date = floor_date(as_date(POSIXct))) %>%
    group_by(date) %>% 
    mutate(WY = ifelse(lubridate::month(date) >= 10, lubridate::year(date) + 1, lubridate::year(date))) %>% 
    group_by(WY) %>%
    summarize(maxQ = max(q_cms, na.rm = T), lenDat = length(q_cms), lenNAs = sum(is.na(q_cms))) %>%
    filter(lenDat > 328*4 & lenNAs / lenDat < 0.1) %>% 
    mutate(rank = rank(maxQ, ties.method = "last")) %>%
    mutate(ReturnInterval = (length(rank) + 1)/(length(rank) + 1 - rank)) %>%
    mutate(ExceedProb = 1 / ReturnInterval)
    
  return(x[,1:2])
}

Annual_peak_df = cbind(Annual_peak_calc(channelRead_NWM),
                  Annual_peak_calc(channelRead_AP)[,2],
                  Annual_peak_calc(channelRead_SHUF)[,2],
                  Annual_peak_calc(channelRead_HUC12L)[,2],
                  Annual_peak_calc(channelRead_PERT)[,2])
names(Annual_peak_df) = c("WY", scheme_list)





#####
#####
##### Rising limb density
#https://github.com/TOSSHtoolbox/TOSSH/blob/e76cf3217f37cac079ad31173a89652f2034a0a5/TOSSH_code/signature_functions/sig_RisingLimbDensity.m

RLD_varList = c(1, 0, 0)
varargin = RLD_varList

sig_RisingLimbDensity = function(Q, varargin, ts_hr = 6){
  
  Q = Q
  
  rising_limb_length = varargin[1];   eps = varargin[2];   minimum_peak = varargin[3];
  
  len_increase = rising_limb_length/(ts_hr/24) #days(t[2]-t[1]) = 0.25
  
  increasing_flow = Q[2:length(Q)]>(Q[1:length(Q)-1]-eps) #%>% as.data.frame()#it must be picking timesteps where flow is increasing
  start_point = min(which(increasing_flow == FALSE)) 
  increasing_flow = increasing_flow[start_point:length(increasing_flow)] 
  
  flow_change = which(increasing_flow[1:length(increasing_flow)-1] != increasing_flow[2:length(increasing_flow)])#indices for start and end of increasing section
  flow_change2 = as.data.frame(cbind(flow_change[seq(1, length(flow_change), 2)], flow_change[seq(2, length(flow_change), 2)]))
  
  flow_section = flow_change2 %>% 
    dplyr::filter((flow_change2[,2] - flow_change2[,1]) >= len_increase)
  flow_section = flow_section + start_point
  flow_section = flow_section %>% 
    filter(!row_number() %in% which((Q[flow_section[,2]]- Q[flow_section[,1]]) < minimum_peak))
  
  RLD = 1/mean(flow_section[,2] - flow_section[,1])
  
  return(RLD)
}

RLD_df = cbind(sig_RisingLimbDensity(USGS_flow$Flow, RLD_varList, ts_hr = 24),
               sig_RisingLimbDensity(   channelRead_NWM$q_cms,     RLD_varList),
               sig_RisingLimbDensity(    channelRead_AP$q_cms,     RLD_varList),
               sig_RisingLimbDensity(  channelRead_SHUF$q_cms,     RLD_varList),
               sig_RisingLimbDensity(channelRead_HUC12L$q_cms,     RLD_varList),
               sig_RisingLimbDensity(  channelRead_PERT$q_cms,     RLD_varList))
colnames(RLD_df) = c("USGS", scheme_list)






#####
#####
##### Master Recession Curve
#https://github.com/TOSSHtoolbox/TOSSH/blob/e76cf3217f37cac079ad31173a89652f2034a0a5/TOSSH_code/signature_functions/sig_RisingLimbDensity.m

#recession_length, n_start, eps, filter_par, seg_test
MRC_varlist =c(15, 0 , 0, 0.925, 0.75)

sig_MRC_SlopeChanges = function(Q, varargin, ts_hr = 6){}

library(devtools)
install_github("cran/lfstat")
library(lfstat)

ts_build = ts(channelRead_NWM %>% 
                filter(POSIXct >= as.POSIXct("2000-09-30") & POSIXct <= as.POSIXct("2020-09-30")) %>% 
                dplyr::select(c("q_cms"))              )

ts_build = ts(channelRead_AP %>% make_daily())
ts_build = ts(USGS_flow$Flow)

my_lfobj = createlfobj(ts_build, startdate = "01/10/2000",hyearstart = 10, dateformat = "%d/%m/%Y")
seglenplot(my_lfobj, threslevel = 70)
my_recession = recession(my_lfobj, method = "MRC",seglen = 7,threshold = 70)

