detach("package:wrfhydroSubsetter", unload=TRUE)
remove.packages("wrfhydroSubsetter")

devtools::install("/mnt/d/GitHub/wrfhydroSubsetter", repos = NULL, type="source")
library(wrfhydroSubsetter)
version_ck <- latest_nwm_version()
#download_conus_nwm("/mnt/d/FULLDOMAIN")

library(devtools)
detach("package:resample", unload=TRUE)
remove.packages("resample")
devtools::install_github("mikejohnson51/resample")
package_list <-utils::installed.packages()


