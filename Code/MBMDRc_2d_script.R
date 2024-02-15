

library(data.table)
library(tidyverse)


rm(list = ls())



# Functions ---------------------------------------------------------------

create_MDR_model_2d = function( raw_model, chisq_df, p_threshold = 0.05, row_2d = 14){
  # """
  # Input: p_threshold : a double indicating the p-value threshold you are considering
  #        raw_model : the model output from MBMDR, usually a file named "MBMDRinput_models"
  #        chisq_df : the output of MBMDR, a dataframe with the SNP pair, the chi square,  the p value and the snp1_snp2 names
  # Output: a list for each SNP significant (below the thresold)  pair, 
  #  with the # individuals in each cell, the HLO matrix and the average outcome per cell
  # """
  mdr_model <- list()
  nrow_sign = sum(chisq_df$p_val < p_threshold)
  
  for (i in seq(nrow_sign)){
    idx <- row_2d * (i-1) + 1
    ma_names <- paste(raw_model[idx, 1:2], collapse = '_')
    
    mdr_model[[ma_names]] <- list(
      affected = as_tibble(sapply(raw_model[idx + (2:4), ], as.numeric)),
      avg_cont_trait = as_tibble(raw_model[idx + (6:8), ]),
      HLO = as_tibble(raw_model[idx + (10:12), ]) %>%
        map_df(~ recode(.x, L = -1, H = 1, O = 0)) %>%
        as.matrix()
    )
  }
  return(mdr_model)
}


Risk_calculation_2d = function(snp_dat, chisq_df, mdr_model, is_0_considered = "0", p_threshold = 0.05 ){
  # """
  # Input: snp_dat: the base data (without the phenotype) where the rows are individuals and the columns are SNPs
  #        chisq_df : the output of MBMDR, a dataframe with the SNP pair, the chi square,  the p value and the snp1_snp2 names
  #        mdr_model :a list for each SNP significant (below the thresold)  pair, 
  #                   with the # individuals in each cell, the HLO matrix and the average outcome per cell
  #        is_0_considered : "0" or "1" indicating if we take into account the 0 of the HL0 matrix ("0") or not ("1")
  #        p_threshold : a double indicating the p-value threshold you are considering
  # Output: a vector with the estimated individual risk score per each indivudal
  # """
  risk_2d <- vector(mode = 'numeric', length = nrow(snp_dat)) # For each individual
  
  for (subj_idx in seq(nrow(snp_dat))){ # For every subject
    sums = 0
    for (j in seq(nrow(chisq_df))){ # For every interaction -->
      if (chisq_df$p_val[j] < p_threshold){ # If the interaction is significant 
        
        ma_pair <- chisq_df$ma_names[j] # e.g. X1_X2
        chi_sq <- chisq_df$chi_sq[j] # chi-squared values of that SNP comb
        mas <- strsplit(ma_pair, '_') %>% unlist # SNPs names, e.g. c('X1', 'X2')
        snp1_val <- snp_dat[[subj_idx, mas[1]]]
        snp2_val <- snp_dat[[subj_idx, mas[2]]]
        
        snp2_val = 9
        s1 = snp1_val == 9 | snp1_val == -9 ; s2 = snp2_val == 9 | snp2_val == -9 
        
        if (! (s1 || s2)){
          
          hlo <- mdr_model[[ma_pair]]$HLO
          hlo_i = hlo[snp1_val + 1, snp2_val + 1]
          avg_trait_i = as.numeric( mdr_model[[ma_pair]]$avg_cont_trait[snp1_val + 1, snp2_val + 1] )
          
          sums = sums + ifelse(is_0_considered == "0", abs(hlo_i) * avg_trait_i ,  avg_trait_i)  # 0 if there is 0 in HLO?
          
        } else if (s1 && !s2){ # s1 is NA
          
          conditional_prop = mdr_model[[ma_pair]]$affected[, snp2_val+1] / sum(mdr_model[[ma_pair]]$affected[, snp2_val+1] )
          avg_trait_i = mdr_model[[ma_pair]]$avg_cont_trait[, snp2_val + 1,drop=FALSE] %>%
            mutate_all(as.numeric)
          
          sums = sums + sum(avg_trait_i * conditional_prop ) 
          
        } else if (s2 && !s1){ ###s2 is NA
          
          conditional_prop = (mdr_model[[ma_pair]]$affected[snp1_val+1,] / sum(mdr_model[[ma_pair]]$affected[snp1_val+1,] ) ) %>% 
            t() %>% 
            as.tibble()
          
          avg_trait_i = (mdr_model[[ma_pair]]$avg_cont_trait[snp1_val+1, ,drop=FALSE] ) %>% 
            t() %>% 
            as.tibble() %>% 
            mutate_all(as.numeric)
          
          sums = sums + sum(avg_trait_i * conditional_prop ) 
          
        } else { ## BOTH MISSING
          
          unconditional_trait = mdr_model[[ma_pair]]$avg_cont_trait %>%
            mutate_all(as.numeric)
          
          unconditional_prop = mdr_model[[ma_pair]]$affected / sum(mdr_model[[ma_pair]]$affected )
          
          sums = sums + sum(unconditional_trait * unconditional_prop )
        }
        
        
        
      }
    }
    risk_2d[subj_idx] <- sums
  }
  return(risk_2d)
}


# Parameters --------------------------------------------------------------

data_path = "C:/Users/fmelo/Documents/Github/MBMDRc_R/Data/2d"
result_path = "C:/Users/fmelo/Documents/Github/MBMDRc_R/Result/2d/"

args = commandArgs(trailingOnly=TRUE)
p_threshold <- as.numeric(args[1]) #i.e., 0.95
is_0_considered <- args[2]  # 0 = considered ; 1 Not_considered

set.seed(7)
row_2d <- 14 # This is not a parameter to change: number of rows in each HLO matrix

setwd(data_path)


# Data input --------------------------------------------------------------



raw_model = fread(file = paste0("MBMDRinput_models", ".txt"),fill = T, header = F)
chisq_df = fread(file = paste0("mbmdr_output", ".txt"), header = F, skip = 0, col.names = c('ma1', 'ma2', 'chi_sq', 'p_val'))
chisq_df = chisq_df %>% mutate(ma_names = paste(ma1,ma2, sep = "_"))


data <- as.tibble(fread("MockData.txt"))
snp_dat = select(data, - c(D))



### MDR model IMPLEMENTATION --> it exploit the facts that the entries are ORDERED --> took only the top rows, that have significant p-value
setwd(result_path)



mdr_model = create_MDR_model_2d( raw_model, chisq_df, p_threshold = p_threshold)
saveRDS(mdr_model, file = paste0("MDR_model_",p_threshold,".rds"))

print(mdr_model[[1]])

# Calculate risk ----------------------------------------------------------


risk_2d = Risk_calculation_2d(snp_dat, chisq_df, mdr_model, is_0_considered,p_threshold)
saveRDS(risk_2d, file = paste0("MBMDRc_risk_",p_threshold,"with",is_0_considered,".rds"))



result_df <- data.frame(risk_2d, data) 

result_df %>%
  rownames_to_column('Subj') %>%
  dplyr::select(Subj, D, risk_2d) %>%
  fwrite(paste0("MBMDRc_risk_file_MBMDR_",p_threshold,"with",is_0_considered, ".txt"))


