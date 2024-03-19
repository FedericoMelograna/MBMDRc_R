
Wrapper_MBMDRc = function(p_threshold, is_0_considered, data_path, result_path, model, output, data, d_level){
  Calc_MBMDR_1d = function(p_threshold, is_0_considered, data_path, result_path, model, output, data){
    
    library(data.table)
    library(tidyverse)
    
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
    
    
    p_threshold <- as.numeric(p_threshold) #0.95
    
    
    set.seed(7)
    row_1d <- 7
    
    setwd(data_path)
    
    
    raw_model = read.table(file = model,fill = T, header = F)
    chisq_df = fread(file = output, header = F, skip = 0,col.names = c('ma1', 'chi_sq', 'p_val'))
    data <- as.tibble(fread(data))
    
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
    
  }
  Calc_MBMDR_2d = function(p_threshold, is_0_considered, data_path, result_path, model, output, data){
    
    library(data.table)
    library(tidyverse)
    
    
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
    
    
    p_threshold <- as.numeric(p_threshold) #i.e., 0.95
    set.seed(7)
    row_2d <- 14 # This is not a parameter to change: number of rows in each HLO matrix
    
    setwd(data_path)
    
    
    # Data input --------------------------------------------------------------
    
    raw_model = fread(file = model,fill = T, header = F)
    chisq_df = fread(file = output, header = F, skip = 0, col.names = c('ma1', 'ma2', 'chi_sq', 'p_val'))
    chisq_df = chisq_df %>% mutate(ma_names = paste(ma1,ma2, sep = "_"))
    data <- as.tibble(fread(data))
    snp_dat = select(data, - c(D))
    
    
    
    ### MDR model IMPLEMENTATION --> it exploit the facts that the entries are ORDERED --> took only the top rows, that have significant p-value
    setwd(result_path)
    mdr_model = create_MDR_model_2d( raw_model, chisq_df, p_threshold = p_threshold)
    saveRDS(mdr_model, file = paste0("MDR_model_",p_threshold,".rds"))
    
    
    # Calculate risk ----------------------------------------------------------
    
    risk_2d = Risk_calculation_2d(snp_dat, chisq_df, mdr_model, is_0_considered,p_threshold)
    saveRDS(risk_2d, file = paste0("MBMDRc_risk_",p_threshold,"with",is_0_considered,".rds"))
    result_df <- data.frame(risk_2d, data) 
    result_df %>%
      rownames_to_column('Subj') %>%
      dplyr::select(Subj, D, risk_2d) %>%
      fwrite(paste0("MBMDRc_risk_file_MBMDR_",p_threshold,"with",is_0_considered, ".txt"))
    
  }
  
  d_level = as.numeric(d_level)
  if (d_level == 1){
    Calc_MBMDR_1d(p_threshold, is_0_considered, data_path, result_path, model, output, data)
    print("MBMDRc-1d successfully run")
  } else if (d_level == 2) {
    Calc_MBMDR_2d(p_threshold, is_0_considered, data_path, result_path, model, output, data)
    print("MBMDRc-2d successfully run")
  }
  
}


