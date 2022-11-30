{
  library(sf)
  library(lubridate)
  library(tidyverse)
  library(ggplot2)
  library(cowplot)
  library(wrfhydroSubsetter)
  library(reshape2)
  library(stringr)
}


getwd()

outDir = paste0("./ANALYSIS/")

#library(AOI)
#library(nhdplusTools)
#x= get_nwis(AOI = st_transform(findNLDI(comid = 1631587), crs = 4269))

locs = data.frame(comids = c(23762661, 1631587, 5894384, 19389766, 191739, 5781369),
                  siteID = c('14190500', '08173000', '01616500', '03118500' ,'6709000' , '08159000'))

locs = data.frame(comids = c(191739, 1631587, 5894384, 5781369, 19389766, 23762661),
                  siteID = c('6709000', '08173000', '01616500', '08159000' ,'03118500' , '14190500'))
namelist = c()
for (i in 1:nrow(locs)){
  state = st_transform(AOI::aoi_get(state = "all", county = "all"), 4269)[st_transform(findNLDI(comid = locs$comids[i]), crs = 4269),]
  name = gsub(" ", "-", tolower(paste(state$name, state$state_name, locs$comids[i], locs$siteID[i], sep = "_")))
  namelist = c(namelist, name)
  rm(name)
}


scheme_list = c('NWM', 'AP',  'HUC12L', 'SHUF(NWM)', 'SHUF(AP)')

fdc_df_builder = function(dat, nam) {
  dat <- sort(dat, decreasing = T)
  df  <- data.frame(x = 100/length(dat) * 1:length(dat), y = dat)
  colnames(df) = c(paste0(nam, ".x"), paste0(nam, ".y"))
  return(df)
}

make_daily = function(df){
  x = df %>% 
    mutate(date = floor_date(as_date(POSIXct))) %>%
    group_by(date) %>% 
    summarize (q_cms = mean(q_cms)) %>% 
    filter(date >= as.POSIXct("2000-10-01 00:00:00") & date <= as.POSIXct("2020-10-01 00:00:00")) %>% 
    dplyr::select(q_cms) %>% as.matrix() %>% as.numeric()
  return(x)
}

#library(scales)


myplots <- vector('list', length(namelist))
for (i in 1:length(namelist)){
  
  myplots[[i]] <- local({
    i = i
    
    name = namelist[i]
    print(name)
    
    USGS_flow  = read.csv(paste0(outDir, name,"/USGS_6hrly.csv"), header = TRUE, row.names = 1)
    colnames(USGS_flow) = c("POSIXct", "q_cms2")
    nwm_flow   = read.csv(paste0(outDir, name,"/RawFlow_NWM.csv"), header = TRUE, row.names = 1) %>% dplyr::filter(POSIXct >= as.POSIXct("2000-10-01 00:00:00") & POSIXct <= as.POSIXct("2020-10-01 00:00:00"))
    na_vals_USGS = which(is.na(left_join(nwm_flow, USGS_flow)$q_cms2))
    
    nwm_flow = nwm_flow %>% filter(!row_number() %in% as.vector(na_vals_USGS))
    
    colnames(USGS_flow) = c("POSIXct", "q_cms")
    AP_flow    = read.csv(paste0(outDir, name,"/RawFlow_AP.csv"), header = TRUE, row.names = 1) %>% 
      dplyr::filter(POSIXct >= as.POSIXct("2000-10-01 00:00:00") & POSIXct <= as.POSIXct("2020-10-01 00:00:00")) %>%
      filter(!row_number() %in% as.vector(na_vals_USGS))
    
    SHUF_flow  = read.csv(paste0(outDir, name,"/RawFlow_SHUF.csv"), header = TRUE, row.names = 1) %>% 
      dplyr::filter(POSIXct >= as.POSIXct("2000-10-01 00:00:00") & POSIXct <= as.POSIXct("2020-10-01 00:00:00")) %>% 
      filter(!row_number() %in% as.vector(na_vals_USGS))
    
    HUC12_flow = read.csv(paste0(outDir, name,"/RawFlow_HUC12L.csv"), header = TRUE, row.names = 1) %>% 
      dplyr::filter(POSIXct >= as.POSIXct("2000-10-01 00:00:00") & POSIXct <= as.POSIXct("2020-10-01 00:00:00")) %>% 
      filter(!row_number() %in% as.vector(na_vals_USGS))
    
    PERT_flow  = read.csv(paste0(outDir, name,"/RawFlow_PERT.csv"), header = TRUE, row.names = 1) %>% 
      dplyr::filter(POSIXct >= as.POSIXct("2000-10-01 00:00:00") & POSIXct <= as.POSIXct("2020-10-01 00:00:00")) %>% 
      filter(!row_number() %in% as.vector(na_vals_USGS))
    
    
    fdc_USGS = fdc_df_builder(USGS_flow %>% make_daily(), "USGS")
    fdc_nwm = fdc_df_builder(nwm_flow %>% make_daily(), "NWM")
    fdc_AP = fdc_df_builder(AP_flow %>% make_daily(), "AP")
    fdc_SHUF = fdc_df_builder(SHUF_flow %>% make_daily(), "SHUF(AP)")
    fdc_HUC12L = fdc_df_builder(HUC12_flow %>% make_daily(), "HUC12L")
    fdc_PERT = fdc_df_builder(PERT_flow %>% make_daily(), "SHUF(NWM)")
    
    fdc_merged = cbind(fdc_USGS, fdc_nwm, fdc_AP, fdc_SHUF, fdc_HUC12L, fdc_PERT) %>% mutate(x = USGS.x) %>% 
      select(-c('USGS.x', 'NWM.x', 'AP.x', 'HUC12L.x', 'SHUF(NWM).x', 'SHUF(AP).x'))
    colnames(fdc_merged) = c('USGS', 'NWM', 'AP', 'SHUF(AP)', 'HUC12L', 'SHUF(NWM)', 'x')

    #min_val = plyr::round_any(min(fdc_merged[1,]), 0.01, f=ceiling)
    min_val = plyr::round_any(pmax(min(fdc_merged[nrow(fdc_merged),]), 0.1), 0.1, f=ceiling)
    if (max(fdc_merged[1,]) < 50){
      max_val = plyr::round_any(max(fdc_merged[1,]), 10, f=ceiling)
    } else if (max(fdc_merged[1,]) < 500) {
      max_val = plyr::round_any(max(fdc_merged[1,]), 50, f=ceiling)
    } else {
      max_val = plyr::round_any(max(fdc_merged[1,]), 100, f=ceiling)
    }
    #max_val = plyr::round_any(max(fdc_merged[1,]), 10, f=ceiling)
    #coeff = 1/max_val
    
    fdc_merged = fdc_merged %>% gather(key = "Scheme", value = "value", -x)
    #fdc_merged$variable <- factor(fdc_merged$variable, levels = c('NWM', 'AP', 'HUC12L', 'SHUF(AP)', 'SHUF(NWM)', 'USGS'))
    fdc_merged$Scheme <- factor(fdc_merged$Scheme, levels = c('AP', 'SHUF(AP)', 'NWM', 'SHUF(NWM)', 'HUC12L',  'USGS'))
    
    #https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
    fdc_merged_plot = ggplot(fdc_merged, aes(x=x, y=value))+
      geom_line(aes(color = Scheme, linetype = Scheme), cex=0.5) + 
      #geom_line(aes(x=USGS.x, y=USGS.y, color = "USGS"), cex=0.5) +
      #geom_line(aes(x=NWM.x, y=NWM.y, color = "NWM"), cex=0.5) +
      #geom_line(aes(x=AP.x, y=AP.y, color = "AP"), cex=0.5) +
      #geom_line(aes(x=SHUF.x, y=SHUF.y, color = "SHUF"), cex=0.5) +
      #geom_line(aes(x=HUC12L.x, y=HUC12L.y, color = "HUC12L"), cex=0.5) +
      #geom_line(aes(x=SHUF(NWM).x, y=SHUF(NWM).y, color = "SHUF(NWM)"), cex=0.5) +
      #scale_y_continuous(limits = c(min=min_val, max=max_val), trans='log2',
      #                   breaks = c(min_val, plyr::round_any((min_val+max_val)*0.025, 5, f=floor), max_val*0.1, max_val*0.2, max_val*0.5, max_val))+
      scale_y_log10(limits = c(min=min_val, max=max_val),labels= ~ format(.x, digits=0, scientific = FALSE))+
      scale_x_continuous(breaks = c(0, 20, 40, 60, 80, 100), limits = c(min=0, max=100),expand=c(0.02,0.02))+
      scale_linetype_manual(values=c("solid", "twodash","solid","twodash", "solid","dotted"))+
      scale_color_manual(values=c('#b60a1c', "#e03531", '#117733', "#009E73", "#E69F00", "black"))+
      ggtitle(paste0(str_to_title(str_split(name, "_")[[1]][1]), " (", str_to_title(str_split(name, "_")[[1]][2]), ")"))+
      theme(#axis.line=element_blank(),
        axis.line = element_line(size = 0.5, colour = "black", linetype=1),
        #plot.margin=unit(c(0,0.2,0,0),"cm"),
        #axis.text.x=element_blank(),
        #axis.text.y=element_blank(),
        #axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        panel.grid.minor=element_blank(),
        #panel.grid=element_blank(),
        axis.ticks = element_line(size = 1, color="black"),
        axis.ticks.length=unit(-0.1, "cm"), 
        axis.text.x = element_text(margin=unit(c(0.2,0.2,0.2,0.2), "cm")), 
        axis.text.y = element_text(margin=unit(c(0.2,0.2,0.2,0.2), "cm")),
        #legend.text=element_text(size=10),
        #legend.justification=c(1.1,-0.1),
        #legend.position="none",
        legend.title=element_text(),
        panel.background = element_blank(),
        plot.title = element_text(),
        plot.margin = margin(6, 6, 6, 6)
        )
    
    #print(fdc_merged_plot)
    #rm(fdc_merged, min_val, max_val)
  })
}

p_grid = plot_grid(myplots[[1]], myplots[[2]], myplots[[3]], myplots[[4]], myplots[[5]], myplots[[6]], labels = c('A.', 'B.', 'C.', 'D.' ,'E.', 'F.'), ncol=3, align="v")

legend_b <- get_legend(myplots[[1]] + 
                         guides(color = guide_legend(nrow = 1)) +
                         theme(legend.position = "bottom"))

p_grid2 = plot_grid(p_grid, legend_b, rel_heights = c(1, .1), scale = c(1, 3), ncol=1)

p_grid2
