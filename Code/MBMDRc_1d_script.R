
library(data.table)
library(tidyverse)

rm(list = ls())

create_MDR_model_1d = function( raw_model, chisq_df, p_threshold = 0.05, row_1d = 7){
  
  mdr_model <- list()
  nrow_sign = sum(chisq_df$p_val < p_threshold)
  
  for (i in seq(nrow_sign)){
    idx <- row_1d * (i-1) + 1
    ma_names <- as.character(raw_model[idx, 1])
    mdr_model[[ma_names]] <- list(
      affected = as_tibble(raw_model[idx + 2, ]),
      avg_cont_trait = as_tibble(raw_model[idx + 4, ]),
      HLO = as_tibble(raw_model[idx + 6,  ])%>%
        map_df(~ dplyr::recode(.x, L = -1, H = 1, O = 0, N = 0)) %>%
        as.matrix() )
  }
  
  return(mdr_model)
}

Risk_calculation_1d = function(snp_dat, chisq_df, mdr_model, is_0_considered = "0", p_threshold = 0.05 ){
  risk_1d <- vector(mode = 'numeric', length = nrow(snp_dat))
  for (subj_idx in seq(nrow(snp_dat))){
    sums = 0
    for (j in seq(nrow(chisq_df))){
      if (chisq_df$p_val[j] < p_threshold){
        ma <- chisq_df$ma1[j] # e.g. X4
        chi_sq <- chisq_df$chi_sq[j] # chi-squared values of that SNP comb
        snp1_val <- as.integer(snp_dat[subj_idx, ma])
        hlo <- mdr_model[[ma]]$HLO
        if ((snp1_val != -9)){ # if SNP1 is not missing
          if (!is.na(hlo[snp1_val + 1])){
            avg_trait_i = as.numeric(as.character( mdr_model[[ma]]$avg_cont_trait[[snp1_val + 1]] ) ) 
            hlo_i = hlo[snp1_val + 1] # hlo of our individual
            # sum <- sum + hlo[snp1_val + 1]*chi_sq # check order here   
            sums = sums + ifelse(is_0_considered == "0",  abs(hlo_i) * avg_trait_i ,  avg_trait_i)  # 0 if there is 0 in HLO?
            
          }
        }
      }
    }
    risk_1d[subj_idx] <- sums
  }
  
  return(risk_1d)
}





data_path = "C:/Users/fmelo/Documents/Github/MBMDRc_R/Data/1d"


result_path = "C:/Users/fmelo/Documents/Github/MBMDRc_R/Result/1d"

args = commandArgs(trailingOnly=TRUE)
p_threshold <- as.numeric(args[1]) #0.95
is_0_considered <- (args[2]) # 0 = considered ; 1 Not_considered

p_threshold <- 0.05
is_0_considered = "0"
set.seed(7)
row_1d <- 7

setwd(data_path)
checkStrict <- function(f, silent=FALSE) {
  vars <- codetools::findGlobals(f)
  found <- !vapply(vars, exists, logical(1), envir=as.environment(2))
  if (!silent && any(found)) {
    warning("global variables used: ", paste(names(found)[found], collapse=', '))
    return(invisible(FALSE))
  }
  
  !any(found)
}

raw_model = read.table(file = paste0("mbmdr_model_1d", ".txt"),fill = T, header = F)
chisq_df = fread(file = paste0("mbmdr_output_1d", ".txt"), header = F, skip = 0,col.names = c('ma1', 'chi_sq', 'p_val'))
# chisq_df = chisq_df %>% mutate(ma_names = paste(ma1,ma2, sep = "_"))

setwd(data_path)

data <- as.tibble(fread("MockData_1d.txt"))
snp_dat = select(data, - c(D))




# Create MDR model 1d -----------------------------------------------------

mdr_model = create_MDR_model_1d( raw_model, chisq_df, p_threshold = p_threshold)


setwd(result_path)
saveRDS(mdr_model, file = paste0("MDR_model_1d_",p_threshold,".rds"))


# Calculate risk ----------------------------------------------------------

risk_1d = Risk_calculation_1d(snp_dat, chisq_df, mdr_model, is_0_considered,p_threshold)


saveRDS(risk_1d, file = paste0("MBMDRc_risk_1d_",p_threshold,"with",is_0_considered,".rds"))


result_df <- data.frame(risk_1d, data)
result_df %>% 
  rownames_to_column('Subj') %>%
  dplyr::select(Subj, D, risk_1d) %>%
  fwrite(paste0("MBMDRc_risk_file_MBMDR_1d_",p_threshold,"with",is_0_considered ,".txt"))
