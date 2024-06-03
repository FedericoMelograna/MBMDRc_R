
# install dependencies
webr::install("data.table")
webr::install("tidyverse")

# download the R script for MBMDRc
download.file("https://raw.githubusercontent.com/FedericoMelograna/MBMDRc_R/main/Code/MBMDRc_script.R", "MBMDRc_script.R")
source("MBMDRc_script.R")


Wrapper_MBMDRc = function(p_threshold, is_0_considered, model,data,  chisq, d_level, seed = 7){
  d_level = as.numeric(d_level)
  if (d_level == 1){
    individual_risk = Calc_MBMDR_Risk_1d(p_threshold, is_0_considered, model, data, chisq, seed)
    print("MBMDRc-1d successfully run")
  } else if (d_level == 2) {
    individual_risk = Calc_MBMDR_Risk_2d(p_threshold, is_0_considered, model, data, chisq, seed)
    print("MBMDRc-2d successfully run")
  }
  filename = paste0("Individual_risk_", d_level, "d.txt")
  write(x = individual_risk, file = paste0("Individual_risk_", d_level, "d.txt"), ncolumns = 1)
  
  return(filename)
}

