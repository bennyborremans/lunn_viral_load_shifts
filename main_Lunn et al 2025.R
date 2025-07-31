## Title: Periodic shifts in viral load increase risk of spillover from bats
## Code authors: Tamika Lunn, Griffith University & Benny Borremans, Wildlife Health Ecology Research Organization
## Version: V4, created 10th February 2025

rm(list=ls())

##------------------------------------------------------------------------------------------##
##---------------------------------------Load packages--------------------------------------##
##------------------------------------------------------------------------------------------##

## R version 4.2.2 ##
library(tidyverse) ## v2.0.0
library(broom) ## v1.0.3
library(ggplot2) ## v3.4.1
library(ggpubr) ## v0.6.0
library(cowplot) ## v1.1.1
library(RColorBrewer) ## v1.1-3
library(patchwork) ## v1.1.2
library(kableExtra) ## v1.3.4
library(mgcv)
library(splines)

##------------------------------------------------------------------------------------------##
##-------------------------------------Custom functions-------------------------------------##
##------------------------------------------------------------------------------------------##

##----------------------------Function(s): Interquartile range -----------------------------##
Liqr<-function(x) { 
  return(round(quantile(x,0.25,na.rm=TRUE),4)) 
}
Uiqr<-function(x) { 
  return(round(quantile(x,0.75,na.rm=TRUE),4)) 
}

##------------------------------- Function to prepare data ---------------------------------##
prepdata <- function(urv_data, threshold){ 
  ## Start data summaries:  
  results_w_meta <- dplyr::mutate(urv_data,
                                  buffer_pos_hev = dplyr::if_else(ct_hev > 0, preservative, ""))
  results_w_meta <- dplyr::filter(results_w_meta, site %in% continuous_ur_s$site) 
  results_w_meta <- results_w_meta %>% 
    mutate(nnonbff_over_sheet = ifelse(is.na(nnonbff_over_sheet)==TRUE, 0, nnonbff_over_sheet))## NA values should be 0; if there's a BFF count then sheets were counted
  results_w_meta <- dplyr::filter(results_w_meta, bff_present_update == "bff_present") #look only at sheets for which it was "known" at least one BFF present
  
  results_w_meta <- results_w_meta %>% ## Define cutoffs for thresholds:
    mutate(ct_hev = ifelse(ct_hev<= threshold, ct_hev, 0)) ## Change values above threshold value to zero 
  
  results_by_sample <- dplyr::group_by(results_w_meta,
                                       accession_update, site_code, 
                                       date, sampling_cluster, sampling_cluster_median_date,sampling_group_month, sampling_group_year,
                                       sample_id)
  results_by_sample <- dplyr::summarize(results_by_sample, 
                                        n = sum(!is.na(ct_hev)),
                                        ct_hev = max(ct_hev),
                                        Nml_hev = max(Nml_hev),
                                        mean_BFF = mean(nbff_over_sheet),
                                        mean_nonBFF = mean(nnonbff_over_sheet),
                                        buffer_pos_hev = paste0(unique(buffer_pos_hev), collapse = ""))
  results_by_sample <- dplyr::ungroup(results_by_sample) %>%
    dplyr::select(-c(Nml_hev))%>% ## remove before merge
    dplyr::left_join(genomeHev, by=c("ct_hev"))
  
  results_by_sample <- dplyr::mutate_at(results_by_sample,
                                        c("buffer_pos_hev"),
                                        function(x) dplyr::if_else(x %in% c("", "NA"), NA_character_, x))
  
  hendra_urine_ungrouped <- dplyr::mutate(results_by_sample, 
                                          hen_bin = dplyr::case_when(ct_hev == 0 ~ 0,
                                                                     ct_hev > 0 ~ 1),
                                          hen_ct_of_pos = dplyr::if_else(ct_hev > 0, ct_hev, NA_real_),
                                          hen_Nml_of_pos = dplyr::if_else(Nml_hev > 0, Nml_hev, NA_real_)) 
  hendra_urine_ungrouped <- dplyr::mutate(hendra_urine_ungrouped,
                                          site_code = dplyr::if_else(is.na(site_code), 
                                                                     substr(accession_update, 3, 5), 
                                                                     site_code))
  hendra_urine <- dplyr::group_by(hendra_urine_ungrouped,
                                  sampling_cluster_median_date, site_code)
  
  hendra_urine <- dplyr::summarize(hendra_urine, 
                                   year = format(sampling_cluster_median_date, "%Y"),
                                   month = format(sampling_cluster_median_date, "%m"),
                                   hen_prevalence = sum(hen_bin)/dplyr::n(),
                                   hen_num_pos = sum(hen_bin),
                                   mean_hen_ct_pos = mean(hen_ct_of_pos, na.rm = T),
                                   min_hen_ct_pos = min(hen_ct_of_pos, na.rm = T),
                                   mean_hen_Nml_pos = mean(hen_Nml_of_pos, na.rm = T),
                                   max_hen_Nml_pos = max(hen_Nml_of_pos, na.rm = T),
                                   LIQR_hen_Nml_pos = Liqr(hen_Nml_of_pos[hen_Nml_of_pos!="NA"]),
                                   UIQR_hen_Nml_pos = Uiqr(hen_Nml_of_pos[hen_Nml_of_pos!="NA"]),
                                       
                                   mean_BFF = mean(mean_BFF), 
                                   mean_nonBFF = mean(mean_nonBFF),
                                   number_of_tests = dplyr::n()) %>%
    mutate(
      mean_hen_ct_pos = ifelse(hen_prevalence==0, NA, mean_hen_ct_pos), ## replace NaN
      min_hen_ct_pos = ifelse(hen_prevalence==0, NA, min_hen_ct_pos), ## replace Inf
      mean_hen_Nml_pos = ifelse(hen_prevalence==0, NA, mean_hen_Nml_pos), ## replace NaN
      max_hen_Nml_pos = ifelse(hen_prevalence==0, NA, max_hen_Nml_pos), ## replace -Inf
    )
  hendra_urine <- dplyr::ungroup(hendra_urine)
  hendra_urine <- dplyr::mutate(hendra_urine, 
                                number_of_tests = dplyr::if_else(is.na(hen_prevalence) & number_of_tests == 1, 
                                                                 NA_integer_, 
                                                                 number_of_tests))
  
  site_code_factor_levels <- dplyr::count(hendra_urine, site_code)
  site_code_factor_levels <- dplyr::count(site_code_factor_levels, site_code)
  site_code_factor_levels <- dplyr::arrange(site_code_factor_levels, dplyr::desc(n))
  site_code_factor_levels <- dplyr::pull(site_code_factor_levels, site_code)
  site_code_factor_levels <- unique(site_code_factor_levels)
  
  hendra_urine_with_availability <- dplyr::mutate(hendra_urine,
                                                  tested = dplyr::if_else(!is.na(hen_prevalence), 
                                                                          "available", 
                                                                          "not_available"),
                                                  site_code = factor(site_code, levels = site_code_factor_levels)) %>%
    distinct()
  
  hen_with_ci <- dplyr::mutate(hendra_urine_with_availability,
                               hen_a_param = hen_num_pos + 0.5,
                               hen_b_param = number_of_tests - hen_num_pos + 0.5,
                               hen_lower = qbeta(0.025, hen_a_param, hen_b_param), 
                               hen_upper = qbeta(0.975, hen_a_param, hen_b_param)) %>% 
    distinct()
  
  ## Add hen_prevalence=0 for sampling sites where bats were absent (number_of_tests=0)
  ##         - Burleigh March-April 2018, Burleigh April 2019, Burleigh June 2019, Sunnybank May-July 2018, Clunes March-April 2018
  nobats <- as.data.frame(matrix(nrow=9, ncol=ncol(hen_with_ci)))
  colnames(nobats) <- colnames(hen_with_ci)
  nobats[1,c(1:2)] <- c("2018-03-15", "BUR") ## Burleigh, March 2018 ## need to remove when doing restricted date range
  nobats[2,c(1:2)] <- c("2018-04-15", "BUR") ## Burleigh, April 2018
  nobats[3,c(1:2)] <- c("2019-04-15", "BUR") ## Burleigh, April 2019
  nobats[4,c(1:2)] <- c("2019-06-15", "BUR") ## Burleigh, June 2019
  nobats[5,c(1:2)] <- c("2018-05-15", "SUN") ## Sunnybank, May 2018
  nobats[6,c(1:2)] <- c("2018-06-15", "SUN") ## Sunnybank, June 2018
  nobats[7,c(1:2)] <- c("2018-07-15", "SUN") ## Sunnybank, July 2018
  nobats[8,c(1:2)] <- c("2018-03-15", "CLU") ## Clunes, March 2018 ## need to remove when doing restricted date range
  nobats[9,c(1:2)] <- c("2018-04-15", "CLU") ## Clunes, April 2018
  nobats[,1] = as.Date(nobats[,1])
  nobats[,"hen_prevalence"] <- 0
  nobats[,"number_of_tests"] <- 0
  hen_with_ci <- rbind(hen_with_ci,nobats)
  
  season_info <- dplyr::tibble(month_start = rep(c(1, 6, 9), times = 5), 
                               month_end = rep(c(6, 9, 1), times = 5),
                               year_start = rep(2016:2020, each = 3),
                               year_end = c(rep(2016, 2), rep(2017:2020, each = 3), 2021))
  season_info <- dplyr::mutate(season_info, 
                               season = dplyr::if_else(month_start == 6, "winter", "not_winter"),
                               start_date = as.Date(paste(year_start, month_start, "01", sep = "-")),
                               end_date = as.Date(paste(year_end, month_end, "01", sep = "-"))) 
  
  
  hendra_urine_by_month <- dplyr::mutate(results_by_sample,
                                         hen_bin = dplyr::case_when(ct_hev == 0 ~ 0,
                                                                    ct_hev > 0 ~ 1)) 
  hendra_urine_by_month <- dplyr::mutate(hendra_urine_by_month, 
                                         site_code = dplyr::if_else(is.na(site_code), 
                                                                    substr(accession_update, 3, 5), 
                                                                    site_code))
  hendra_urine_by_month <- dplyr::group_by(hendra_urine_by_month,
                                           sampling_cluster_median_date) 
  hendra_urine_by_month <- dplyr::summarize(hendra_urine_by_month,
                                            year = format(sampling_cluster_median_date,"%Y"),
                                            month = format(sampling_cluster_median_date,"%m"),
                                            hen_prevalence = sum(hen_bin)/dplyr::n(),
                                            hen_num_pos = sum(hen_bin),
                                            number_of_tests = dplyr::n(), 
                                            min_hen_ct = min(ct_hev[ct_hev>0.01]), #but minimum of positive samples, otherwise all will be 0
                                            mean_hen_ct = mean(ct_hev[ct_hev>0.01]),
                                            max_hen_Nml = max(Nml_hev[Nml_hev>0.01], na.rm = T),
                                            mean_hen_Nml = mean(Nml_hev[Nml_hev>0], na.rm = T),
                                            LIQR_hen_Nml_pos = Liqr(Nml_hev[Nml_hev>0]),
                                            UIQR_hen_Nml_pos = Uiqr(Nml_hev[Nml_hev>0])) %>%
    mutate(min_hen_ct=ifelse(hen_prevalence==0,NA,min_hen_ct))  %>%
    mutate(mean_hen_ct=ifelse(hen_prevalence==0,NA,mean_hen_ct))  %>%
    mutate(max_hen_Nml=ifelse(hen_prevalence==0,NA,max_hen_Nml))  %>%
    mutate(mean_hen_Nml=ifelse(hen_prevalence==0,NA,mean_hen_Nml))

  hendra_urine_by_month <- dplyr::ungroup(hendra_urine_by_month)
  hendra_urine_by_month <- dplyr::mutate(hendra_urine_by_month, 
                                         number_of_tests = dplyr::if_else(is.na(hen_prevalence) & number_of_tests == 1, 
                                                                          NA_integer_, 
                                                                          number_of_tests))
  
  
  hendra_urine_with_availability_by_month <- dplyr::mutate(hendra_urine_by_month,
                                                           tested = dplyr::if_else(!is.na(hen_prevalence), 
                                                                                   "available", 
                                                                                   "not_available")) %>%
    distinct()
  
  hen_with_ci_by_month <- dplyr::mutate(hendra_urine_with_availability_by_month,
                                        hen_a_param = hen_num_pos + 0.5,
                                        hen_b_param = number_of_tests - hen_num_pos + 0.5,
                                        hen_lower = qbeta(0.025, hen_a_param, hen_b_param), 
                                        hen_upper = qbeta(0.975, hen_a_param, hen_b_param))  %>% 
    distinct()
  
  output <- list(hen_with_ci, hen_with_ci_by_month, season_info, hendra_urine_ungrouped, results_w_meta) ## store information to be used later in a list for recall
  return(output)
}

##------------------------------------------------------------------------------------------##
##-----------------------------------------Load data----------------------------------------##
##------------------------------------------------------------------------------------------##

continuous_ur_s <- read.csv("Data/Processed/continuous_ur_s.csv") ## Names and locations of sites
urv_data <- read.csv("Data/Processed/urv_data.csv")   ## Urine sampling meta-data with virus data
genomeHev <- read.csv("Data/Processed/Ct-to-genome-convert_HeV.csv") ## Ct to genome per ml conversion
spillover <- read.csv("Data/Processed/spillover.csv") ## spillovers within study

##-------------------------------------Numbers over sheet-----------------------------------##
bats_over_sheet_summary <- urv_data %>%
  filter(site %in% continuous_ur_s$site) %>%
  filter(bff_present_update == "bff_present") %>% #look only at sheets for which it was "known" at least one BFF present
  drop_na(nbff_over_sheet) %>% ## 9 rows dropped
  #drop_na(nnonbff_over_sheet) %>% ## additional NA values should be 0; if there's a BFF count then sheets were counted
  mutate(nnonbff_over_sheet = ifelse(is.na(nnonbff_over_sheet)==TRUE, 0, nnonbff_over_sheet)) %>% ## additional NA values should be 0; if there's a BFF count then sheets were counted
  mutate(bats_over_sheet = nbff_over_sheet+nnonbff_over_sheet) %>%
  group_by(site, sampling_cluster_median_date) %>%
  summarise(mean_all = mean(bats_over_sheet),
            Liqr_total = Liqr(bats_over_sheet),
            Uiqr_total = Uiqr(bats_over_sheet),
            ##BFF
            mean_bff = mean(nbff_over_sheet),
            Liqr_bff = Liqr(nbff_over_sheet),
            Uiqr_bff = Uiqr(nbff_over_sheet),
            ##other
            mean_other = mean(nnonbff_over_sheet),
            Liqr_other = Liqr(nnonbff_over_sheet),
            Uiqr_other = Uiqr(nnonbff_over_sheet)) %>%
  mutate(mean_all = round(mean_all, digits=1),
         all = paste(c(round(mean_all, digits=1))," (", round(Liqr_total, digits=1),"-", round(Uiqr_total, digits=1), ")", sep=""),
         bff = paste(c(round(mean_bff, digits=1))," (", round(Liqr_bff, digits=1),"-", round(Uiqr_bff, digits=1), ")", sep=""),
         non_bff = paste(c(round(mean_other, digits=1))," (", round(Liqr_other, digits=1),"-", round(Uiqr_other, digits=1), ")", sep=""))
  
kable(bats_over_sheet_summary[,c(1,2,12:14)], caption = "Mean (lower IQR - upper IQR)",  escape = F) %>%
  kable_styling("striped", full_width = F) %>% 
  scroll_box(height = "500px") ## Make scroller for table if it is long

#write.csv(bats_over_sheet_summary[,c(1,2,3,12:14)], "bats-over-sheet.csv")


## Functions to calculate virus concentrations in pooled samples ------------##


# calculate mean Ct value of a positive individual in a pooled sample, given the pooled Ct
ct_pos_fun = function(ct_pooled, n_pos, n_total) {
  
  if(ct_pooled > 0) {
  # convert Ct to virus N/ml
  cur_v_conc = exp((46.606 - ct_pooled) / 1.619)
  
  # virus per positive individual
  mean_conc = cur_v_conc * n_total / n_pos
  
  # convert virus N/ml to Ct. = Ct per positive individual, assuming they have the same Ct value
  ct_pos = 46.606 - 1.619  * log(mean_conc)
  
  } else { ct_pos = 0 }
  # round up because that is how the PCR works
  return(ceiling(ct_pos))
}

# e.g., if the Ct value of the pooled sample is 26, if there were 3 positive bats out of 5 total contributing bats over the sheet, the mean Ct value each positive individual bat has would be:
# ct_pos_fun(26, 3, 5)


# calculate pooled Ct value given the mean Ct value of a positive individual, and a number of total bats
ct_pooled_fun = function(ct_pos, n_pos, n_total) {
  
  if(ct_pos > 0) {
  # convert Ct to virus N/ml
  cur_v_conc = exp((46.606 - ct_pos) / 1.619)
  
  # pooled virus concentration given number of positives, and divided by the total number of individuals
  cur_v_pooled = cur_v_conc * n_pos / n_total
  
  # convert virus N/ml to Ct    
  ct_pooled = 46.606 - 1.619  * log(cur_v_pooled)
  } else { ct_pooled = 0 }
  
  # round up because that is how the PCR works
  return(ceiling(ct_pooled))
  
}


##------------------------------------------------------------------------------------------##
##-----------------------------------------Fit models---------------------------------------##
##------------------------------------------------------------------------------------------##

##------------------------------------- Run permutation ------------------------------------##
## Run one threshold for data format:
threshold <- 35
Ct_data <- prepdata(urv_data, threshold)
hen_with_ci_by_month <- Ct_data[[2]] %>%
  mutate(hen_prevalence_normalized = hen_prevalence/max(hen_prevalence))
hen_prev_all = data.frame(sampling_cluster_median_date = hen_with_ci_by_month$sampling_cluster_median_date,
                          prev = hen_with_ci_by_month$hen_prevalence_normalized) 
season_info <- Ct_data[[3]]

## Permutation:
set.seed(111)
permutations = 500
ct_thresholds = 28:40

run.again = F ## T if running for the first time (to save output), F if not (to read prior output)
run.date = "20250210"

if(run.again) {      
  
  # prepare dataframe for permutation prevalence data
  perm.out = data.frame(sampling_cluster_median_date = rep(hen_prev_all$sampling_cluster_median_date, permutations * length(ct_thresholds)),
                        permutation = rep(1:permutations, each = nrow(hen_prev_all) * length(ct_thresholds)),
                        hen_prevalence_normalized = numeric(length(nrow(hen_prev_all) * permutations * length(ct_thresholds))),
                        hen_prevalence_normalized_smooth = numeric(length(nrow(hen_prev_all) * permutations * length(ct_thresholds))),
                        ct_threshold = rep(rep(ct_thresholds, each = nrow(hen_prev_all)) ,permutations))
  
  # prepare dataframe for observed prevalence data
  obs.out = data.frame(sampling_cluster_median_date = rep(hen_prev_all$sampling_cluster_median_date,length(ct_thresholds)),
                       hen_prevalence_normalized = numeric(length(nrow(hen_prev_all) * length(ct_thresholds))),
                       hen_prevalence_normalized_smooth = numeric(length(nrow(hen_prev_all) * length(ct_thresholds))),
                       ct_threshold = rep(ct_thresholds, each = nrow(hen_prev_all)),
                       prev_quantile = numeric(length(nrow(hen_prev_all) * length(ct_thresholds))))
  
  # prepare dataframe for permutation statistics    
  perm.out.stats.gam.poisson = data.frame(permutation = rep(1:permutations,length(ct_thresholds)),
                                          ct_threshold = rep(ct_thresholds, each = permutations),
                                          AIC = NA)
  perm.out.stats.gam.binomial = data.frame(permutation = rep(1:permutations,length(ct_thresholds)),
                                           ct_threshold = rep(ct_thresholds, each = permutations),
                                           AIC = NA)
  perm.out.stats.glm.poisson = data.frame(permutation = rep(1:permutations,length(ct_thresholds)),
                                          ct_threshold = rep(ct_thresholds, each = permutations),
                                          AIC = NA)
  perm.out.stats.glm.binomial = data.frame(permutation = rep(1:permutations,length(ct_thresholds)),
                                           ct_threshold = rep(ct_thresholds, each = permutations),
                                           AIC = NA)
  
  # prepare dataframe for observed statistics    
  obs.out.stats.gam.poisson = data.frame(ct_threshold = ct_thresholds,
                                         AIC = NA,
                                         AIC.quantile = NA)
  
  obs.out.stats.gam.binomial = data.frame(ct_threshold = ct_thresholds,
                                          AIC = NA,
                                          AIC.quantile = NA)
  obs.out.stats.glm.poisson = data.frame(ct_threshold = ct_thresholds,
                                         AIC = NA,
                                         AIC.quantile = NA)
  
  obs.out.stats.glm.binomial = data.frame(ct_threshold = ct_thresholds,
                                          AIC = NA,
                                          AIC.quantile = NA)
  
  for(j in 1:length(ct_thresholds)) {
    
    cur.ct = ct_thresholds[j]  
    
    for(i in 1:permutations){
      
      # load observed data
      urv_data_perm = urv_data
      
      # add column total bats over sheet
      urv_data_perm$n_total = apply(urv_data[,c("nbff_over_sheet","nnonbff_over_sheet")], 1, sum, na.rm=T)
      # if n_total is zero (i.e., there has been at least one bat but not at the time of counting), set n_total to 1
      urv_data_perm$n_total[which(urv_data_perm$n_total==0)] = 1
      
      # calculate individual Ct values given a randomly selected number of positive bats out of the total number of bats over the sheet
      urv_data_perm$cur_n_pos = sapply(urv_data_perm$n_total, FUN = function(x) sample(1:x, size = 1))
      
      # calculate average Ct value of each positive individual
      urv_data_perm$cur_ct_pos = apply(urv_data_perm[,c("ct_hev","cur_n_pos","n_total")], 1, function(x) ct_pos_fun(ct_pooled = x[1], n_pos = x[2], n_total = x[3]))
      
      # permute Ct values, within buffer type
      urv_data_perm$ct_hev[which(urv_data_perm$ct_hev>0 & urv_data_perm$preservative=="AVL")] = sample(urv_data_perm$cur_ct_pos[which(urv_data_perm$ct_hev>0 & urv_data_perm$preservative=="AVL")], replace = F)
      urv_data_perm$ct_hev[which(urv_data_perm$ct_hev>0 & urv_data_perm$preservative=="VTM")] = sample(urv_data_perm$cur_ct_pos[which(urv_data_perm$ct_hev>0 & urv_data_perm$preservative=="VTM")], replace = F)
      urv_data_perm$ct_hev[which(urv_data_perm$ct_hev>0 & urv_data_perm$preservative=="NB")] = sample(urv_data_perm$cur_ct_pos[which(urv_data_perm$ct_hev>0 & urv_data_perm$preservative=="NB")], replace = F)
      
      # calculate pooled Ct value, given the current Ct value for a positive individual, and the current number of positive individuals
      urv_data_perm$ct_hev = apply(urv_data_perm[,c("ct_hev","cur_n_pos","n_total")], 1, function(x) ct_pooled_fun(ct_pos = x[1], n_pos = x[2], n_total = x[3]))
      
      
      
      # prepare data using permuted Ct values
      Ct_data_perm = prepdata(urv_data_perm, cur.ct)
      hen_with_ci_by_month_perm <- Ct_data_perm[[2]] %>%
        mutate(hen_prevalence_normalized = hen_prevalence/max(hen_prevalence))
      # add spillover data
      hen_with_ci_by_month_perm$spillover_count = 0
      hen_with_ci_by_month_perm$spillover_count[which(hen_with_ci_by_month_perm$sampling_cluster_median_date==as.Date("2017-07-17"))] = 1
      hen_with_ci_by_month_perm$spillover_count[which(hen_with_ci_by_month_perm$sampling_cluster_median_date==as.Date("2017-08-08"))] = 2
      hen_with_ci_by_month_perm$spillover_count[which(hen_with_ci_by_month_perm$sampling_cluster_median_date==as.Date("2018-08-30"))] = 1
      hen_with_ci_by_month_perm$spillover_count[which(hen_with_ci_by_month_perm$sampling_cluster_median_date==as.Date("2020-05-26"))] = 1
      
      hen_with_ci_by_month_perm$spillover_binary = as.numeric(hen_with_ci_by_month_perm$spillover_count>0)
      
      # smooth prevalence
      perm.out$hen_prevalence_normalized[which(perm.out$permutation==i & perm.out$ct_threshold==cur.ct)] = hen_with_ci_by_month_perm$hen_prevalence_normalized
      prev.smooth = gam(hen_prevalence_normalized ~ ns(sampling_cluster_median_date, df = 8), family = "binomial", data = hen_with_ci_by_month_perm)
      perm.out$hen_prevalence_normalized_smooth[which(perm.out$permutation==i & perm.out$ct_threshold==cur.ct)] = predict(prev.smooth,type = "response")
      
      hen_with_ci_by_month_perm$hen_prevalence_normalized_smooth = predict(prev.smooth,type = "response")
      
      # fit gam, using smoothed prevalence
      gam1 = gam(spillover_count ~ ns(hen_prevalence_normalized_smooth), data = hen_with_ci_by_month_perm, family = poisson)
      gam2 = gam(spillover_binary ~ ns(hen_prevalence_normalized_smooth), data = hen_with_ci_by_month_perm, family = binomial)
      
      # fit glm, using smoothed prevalence
      glm1 = glm(spillover_count ~ hen_prevalence_normalized_smooth, data = hen_with_ci_by_month_perm, family = poisson)
      glm2 = glm(spillover_binary ~ hen_prevalence_normalized_smooth, data = hen_with_ci_by_month_perm, family = binomial)
      
      # # fit gam, using non-smoothed prevalence
      # gam1 = gam(spillover_count ~ ns(hen_prevalence_normalized), data = hen_with_ci_by_month_perm, family = poisson)
      # gam2 = gam(spillover_binary ~ ns(hen_prevalence_normalized), data = hen_with_ci_by_month_perm, family = binomial)
      # # fit glm, using non-smoothed prevalence
      # glm1 = glm(spillover_count ~ hen_prevalence_normalized, data = hen_with_ci_by_month_perm, family = poisson)
      # glm2 = glm(spillover_binary ~ hen_prevalence_normalized, data = hen_with_ci_by_month_perm, family = binomial)
      
      # store data and statistics
      
      perm.out.stats.gam.poisson$AIC[which(perm.out.stats.gam.poisson$permutation==i & perm.out.stats.gam.poisson$ct_threshold==cur.ct)] = round(AIC(gam1),3)
      
      perm.out.stats.gam.binomial$AIC[which(perm.out.stats.gam.binomial$permutation==i & perm.out.stats.gam.binomial$ct_threshold==cur.ct)] = round(AIC(gam2),3)
      
      perm.out.stats.glm.poisson$AIC[which(perm.out.stats.glm.poisson$permutation==i & perm.out.stats.glm.poisson$ct_threshold==cur.ct)] = round(AIC(glm1),3)
      
      perm.out.stats.glm.binomial$AIC[which(perm.out.stats.glm.binomial$permutation==i & perm.out.stats.glm.binomial$ct_threshold==cur.ct)] = round(AIC(glm2),3)
      
      message("--------------")
      message("--------------")
      message(paste0("Ct threshold ",cur.ct," -- Permutation ",i))
      message("--------------")
      message("--------------")
    }
    # calculate the same stats for the observed data
    
    Ct_data = prepdata(urv_data, cur.ct)
    hen_with_ci_by_month <- Ct_data[[2]] %>%
      mutate(hen_prevalence_normalized = hen_prevalence/max(hen_prevalence))
    # add spillover data
    hen_with_ci_by_month$spillover_count = 0
    hen_with_ci_by_month$spillover_count[which(hen_with_ci_by_month$sampling_cluster_median_date==as.Date("2017-07-17"))] = 1
    hen_with_ci_by_month$spillover_count[which(hen_with_ci_by_month$sampling_cluster_median_date==as.Date("2017-08-08"))] = 2
    hen_with_ci_by_month$spillover_count[which(hen_with_ci_by_month$sampling_cluster_median_date==as.Date("2018-08-30"))] = 1
    hen_with_ci_by_month$spillover_count[which(hen_with_ci_by_month$sampling_cluster_median_date==as.Date("2020-05-26"))] = 1
    
    hen_with_ci_by_month$spillover_binary = as.numeric(hen_with_ci_by_month$spillover_count>0)
    
    # smoothed prevalence
    obs.out$hen_prevalence_normalized[which(obs.out$ct_threshold==cur.ct)] = hen_with_ci_by_month$hen_prevalence_normalized
    prev.smooth = gam(hen_prevalence_normalized ~ ns(sampling_cluster_median_date, df = 8), family = "binomial", data = hen_with_ci_by_month)
    obs.out$hen_prevalence_normalized_smooth[which(obs.out$ct_threshold==cur.ct)] = predict(prev.smooth,type = "response")
    
    hen_with_ci_by_month$hen_prevalence_normalized_smooth = predict(prev.smooth,type = "response")
    
    # fit gam, using smoothed prevalence
    gam1 = gam(spillover_count ~ ns(hen_prevalence_normalized_smooth), data = hen_with_ci_by_month, family = poisson)
    gam2 = gam(spillover_binary ~ ns(hen_prevalence_normalized_smooth), data = hen_with_ci_by_month, family = binomial)
    # fit glm, using smoothed prevalence
    glm1 = glm(spillover_count ~ hen_prevalence_normalized_smooth, data = hen_with_ci_by_month, family = poisson)
    glm2 = glm(spillover_binary ~ hen_prevalence_normalized_smooth, data = hen_with_ci_by_month, family = binomial)
    
    # # fit gam, using non-smoothed prevalence
    # gam1 = gam(spillover_count ~ ns(hen_prevalence_normalized), data = hen_with_ci_by_month, family = poisson)
    # gam2 = gam(spillover_binary ~ ns(hen_prevalence_normalized), data = hen_with_ci_by_month, family = binomial)
    # # fit glm, using non-smoothed prevalence
    # glm1 = glm(spillover_count ~ hen_prevalence_normalized, data = hen_with_ci_by_month, family = poisson)
    # glm2 = glm(spillover_binary ~ hen_prevalence_normalized, data = hen_with_ci_by_month, family = binomial)
    
    # store data and statistics
    
    # gam
    obs.out.stats.gam.poisson$AIC[which(obs.out.stats.gam.poisson$ct_threshold==cur.ct)] = round(AIC(gam1),3)
    
    obs.out.stats.gam.poisson$AIC.quantile[which(obs.out.stats.gam.poisson$ct_threshold==cur.ct)] = sum(perm.out.stats.gam.poisson$AIC[which(perm.out.stats.gam.poisson$ct_threshold==cur.ct)] < obs.out.stats.gam.poisson$AIC[which(obs.out.stats.gam.poisson$ct_threshold==cur.ct)])/permutations
    
    obs.out.stats.gam.binomial$AIC[which(obs.out.stats.gam.binomial$ct_threshold==cur.ct)] = round(AIC(gam2),3)
    
    obs.out.stats.gam.binomial$AIC.quantile[which(obs.out.stats.gam.binomial$ct_threshold==cur.ct)] = sum(perm.out.stats.gam.binomial$AIC[which(perm.out.stats.gam.binomial$ct_threshold==cur.ct)] < obs.out.stats.gam.binomial$AIC[which(obs.out.stats.gam.binomial$ct_threshold==cur.ct)])/permutations
    
    # glm
    obs.out.stats.glm.poisson$AIC[which(obs.out.stats.glm.poisson$ct_threshold==cur.ct)] = round(AIC(glm1),3)
    
    obs.out.stats.glm.poisson$AIC.quantile[which(obs.out.stats.glm.poisson$ct_threshold==cur.ct)] = sum(perm.out.stats.glm.poisson$AIC[which(perm.out.stats.glm.poisson$ct_threshold==cur.ct)] < obs.out.stats.glm.poisson$AIC[which(obs.out.stats.glm.poisson$ct_threshold==cur.ct)])/permutations
    
    obs.out.stats.glm.binomial$AIC[which(obs.out.stats.glm.binomial$ct_threshold==cur.ct)] = round(AIC(glm2),3)
    
    obs.out.stats.glm.binomial$AIC.quantile[which(obs.out.stats.glm.binomial$ct_threshold==cur.ct)] = sum(perm.out.stats.glm.binomial$AIC[which(perm.out.stats.glm.binomial$ct_threshold==cur.ct)] < obs.out.stats.glm.binomial$AIC[which(obs.out.stats.glm.binomial$ct_threshold==cur.ct)])/permutations
    
    # calculate quantiles for the prevalence curves
    # if the distribution of Ct values shifts towards low values during periods of high prevalence, 
    # we would expect observed prevalence during high prevalence to be higher than permuted prevalence (==> upper quantile), for lower Ct thresholds,
    # and at the same time we would expect observed prevalence during low prevalence to be lower than permuted prevalence (==> lower quantile), for lower Ct thresholds.
    # if the Ct distribution does not change, then we would expect similar quantiles (around 50%) during any period.
    # done for each time period
    # this is done using smoothed prevalence
    
    cur.obs.out.idx = which(obs.out$ct_threshold==cur.ct)
    for(k in 1:length(cur.obs.out.idx)){
      cur.prev.obs = obs.out$hen_prevalence_normalized_smooth[which(obs.out$ct_threshold==cur.ct)[k]]
      cur.date.obs = obs.out$sampling_cluster_median_date[which(obs.out$ct_threshold==cur.ct)[k]]
      cur.prev.perm = perm.out$hen_prevalence_normalized_smooth[which(perm.out$ct_threshold==cur.ct & perm.out$sampling_cluster_median_date==cur.date.obs)]
      
      obs.out$prev_quantile[cur.obs.out.idx[k]] = sum(cur.prev.perm < cur.prev.obs)/permutations
    }
    
  }
  
  # save output so it doesn't have the be run again
  write.csv(perm.out, paste0("perm.out.",run.date,".csv"))
  write.csv(perm.out.stats.gam.poisson, paste0("perm.out.stats.gam.poisson.",run.date,".csv"))
  write.csv(perm.out.stats.gam.binomial, paste0("perm.out.stats.gam.binomial.",run.date,".csv"))
  write.csv(perm.out.stats.glm.poisson, paste0("perm.out.stats.glm.poisson.",run.date,".csv"))
  write.csv(perm.out.stats.glm.binomial, paste0("perm.out.stats.glm.binomial.",run.date,".csv"))
  
  write.csv(obs.out, paste0("obs.out.",run.date,".csv"))
  write.csv(obs.out.stats.gam.poisson, paste0("obs.out.stats.gam.poisson.",run.date,".csv"))
  write.csv(obs.out.stats.gam.binomial, paste0("obs.out.stats.gam.binomial.",run.date,".csv"))
  write.csv(obs.out.stats.glm.poisson, paste0("obs.out.stats.glm.poisson.",run.date,".csv"))
  write.csv(obs.out.stats.glm.binomial, paste0("obs.out.stats.glm.binomial.",run.date,".csv"))
} else {
  perm.out = read.csv(paste0("perm.out.",run.date,".csv"))
  perm.out.stats.gam.poisson = read.csv(paste0("perm.out.stats.gam.poisson.",run.date,".csv"))
  perm.out.stats.gam.binomial = read.csv(paste0("perm.out.stats.gam.binomial.",run.date,".csv"))
  perm.out.stats.glm.poisson = read.csv(paste0("perm.out.stats.glm.poisson.",run.date,".csv"))
  perm.out.stats.glm.binomial = read.csv(paste0("perm.out.stats.glm.binomial.",run.date,".csv"))
  
  obs.out = read.csv(paste0("obs.out.",run.date,".csv"))
  obs.out.stats.gam.poisson = read.csv(paste0("obs.out.stats.gam.poisson.",run.date,".csv"))
  obs.out.stats.gam.binomial = read.csv(paste0("obs.out.stats.gam.binomial.",run.date,".csv"))
  obs.out.stats.glm.poisson = read.csv(paste0("obs.out.stats.glm.poisson.",run.date,".csv"))
  obs.out.stats.glm.binomial = read.csv(paste0("obs.out.stats.glm.binomial.",run.date,".csv"))
}

# add xmin and xmax dates for quantile plotting
obs.out.temp = obs.out[which(obs.out$ct_threshold==29),]
obs.out.temp$xmin = obs.out.temp$sampling_cluster_median_date
obs.out.temp$xmax = obs.out.temp$sampling_cluster_median_date

for(i in 1:nrow(obs.out.temp)){
  if(i==1) {
    obs.out.temp$xmin[i] = as.Date("2017-07-01")
    obs.out.temp$xmax[i] = obs.out.temp$sampling_cluster_median_date[i] + floor((obs.out.temp$sampling_cluster_median_date[i+1]-obs.out.temp$sampling_cluster_median_date[i])/2)
  } else {
    obs.out.temp$xmin[i] = obs.out.temp$xmax[i-1]
    obs.out.temp$xmax[i] = obs.out.temp$sampling_cluster_median_date[i] + floor((obs.out.temp$sampling_cluster_median_date[i+1]-obs.out.temp$sampling_cluster_median_date[i])/2)
  }
}
obs.out.temp$xmax[45] = as.Date("2020-09-09")+5

obs.out$xmin = rep(obs.out.temp$xmin, length(ct_thresholds))
obs.out$xmax = rep(obs.out.temp$xmax, length(ct_thresholds))

# add quantile colors
obs.out$fill.col = NA
obs.out$fill.col[which(obs.out$prev_quantile < 0.1)] = rev(brewer.pal(n = 10,"RdYlBu"))[1]
obs.out$fill.col[which(obs.out$prev_quantile >= 0.1 & obs.out$prev_quantile < 0.2)] = rev(brewer.pal(n = 10,"RdYlBu"))[2]
obs.out$fill.col[which(obs.out$prev_quantile >= 0.2 & obs.out$prev_quantile < 0.3)] = rev(brewer.pal(n = 10,"RdYlBu"))[3]
obs.out$fill.col[which(obs.out$prev_quantile >= 0.3 & obs.out$prev_quantile < 0.4)] = rev(brewer.pal(n = 10,"RdYlBu"))[4]
obs.out$fill.col[which(obs.out$prev_quantile >= 0.4 & obs.out$prev_quantile < 0.5)] = rev(brewer.pal(n = 10,"RdYlBu"))[5]
obs.out$fill.col[which(obs.out$prev_quantile >= 0.5 & obs.out$prev_quantile < 0.6)] = rev(brewer.pal(n = 10,"RdYlBu"))[6]
obs.out$fill.col[which(obs.out$prev_quantile >= 0.6 & obs.out$prev_quantile < 0.7)] = rev(brewer.pal(n = 10,"RdYlBu"))[7]
obs.out$fill.col[which(obs.out$prev_quantile >= 0.7 & obs.out$prev_quantile < 0.8)] = rev(brewer.pal(n = 10,"RdYlBu"))[8]
obs.out$fill.col[which(obs.out$prev_quantile >= 0.8 & obs.out$prev_quantile < 0.9)] = rev(brewer.pal(n = 10,"RdYlBu"))[9]
obs.out$fill.col[which(obs.out$prev_quantile >= 0.9)] = rev(brewer.pal(n = 10,"RdYlBu"))[10]

# add quantile categories
obs.out$quantile.cat = NA
obs.out$quantile.cat[which(obs.out$prev_quantile < 0.1)] = 1
obs.out$quantile.cat[which(obs.out$prev_quantile >= 0.1 & obs.out$prev_quantile < 0.2)] = 2
obs.out$quantile.cat[which(obs.out$prev_quantile >= 0.2 & obs.out$prev_quantile < 0.3)] = 3
obs.out$quantile.cat[which(obs.out$prev_quantile >= 0.3 & obs.out$prev_quantile < 0.4)] = 4
obs.out$quantile.cat[which(obs.out$prev_quantile >= 0.4 & obs.out$prev_quantile < 0.5)] = 5
obs.out$quantile.cat[which(obs.out$prev_quantile >= 0.5 & obs.out$prev_quantile < 0.6)] = 6
obs.out$quantile.cat[which(obs.out$prev_quantile >= 0.6 & obs.out$prev_quantile < 0.7)] = 7
obs.out$quantile.cat[which(obs.out$prev_quantile >= 0.7 & obs.out$prev_quantile < 0.8)] = 8
obs.out$quantile.cat[which(obs.out$prev_quantile >= 0.8 & obs.out$prev_quantile < 0.9)] = 9
obs.out$quantile.cat[which(obs.out$prev_quantile >= 0.9)] = 10


##------------------------------------- Supporting GLMs, Figure 4 ------------------------------------##
obs.out.f4 <- obs.out %>% 
  mutate(month = as.numeric(format(as.Date(sampling_cluster_median_date),"%m")),
         year = as.numeric(format(as.Date(sampling_cluster_median_date),"%Y"))) %>%
  mutate(season_onoff = case_when(month <5 ~ "off-season",
                                  month > 9 ~ "off-season",
                                  TRUE ~ "peak-season")) %>% 
  mutate(spillover = case_when(year == 2017 & month == 7 ~ "spillover", ## x1 - all_058
                               year == 2017 & month == 8 ~ "spillover", ## x2 - all_059 & all_060
                               year == 2018 & month == 8 ~ "spillover", ## x1 - all_061
                               year == 2020 & month == 5 ~ "spillover", ## x1 - all_063
                               TRUE~ "no_spillover")) %>% 
  mutate(spillover_season = paste(spillover, season_onoff, sep = "_"))

## When the time series was subdivided based on season and spillover event observation, prevalence was highest for peak seasons in which a spillover event was observed, and this was consistent for all Ct threshold values (effect size: peak-season & spillover=13.7±1.43, peak-season & no spillover=2.33±1.34, p<0.0001, df =2, χ2=53.2)
glm1 = glm(hen_prevalence_normalized_smooth ~ spillover_season, data = obs.out.f4, family = binomial)
glm2 = glm(hen_prevalence_normalized_smooth ~ 1, data = obs.out.f4, family = binomial)

anova(glm1,glm2, test = "Chisq")

sum.glm1 = summary(glm1)
exp(sum.glm1$coefficients)

## For peak season with spillover events, prevalence trended higher for lower Ct threshold values, excluding the two lowest values, although this was not statistically supported (effect size=0.90, p=0.30, df=1, χ2=1.12)
glm1 = glm(hen_prevalence_normalized_smooth ~ ct_threshold, data = obs.out.f4 %>% filter(spillover_season == "spillover_peak-season", ct_threshold %in% 30:40), family = binomial)
glm2 = glm(hen_prevalence_normalized_smooth ~ 1, data = obs.out.f4 %>% filter(spillover_season == "spillover_peak-season", ct_threshold %in% 30:40), family = binomial)
anova(glm1,glm2, test = "Chisq")

sum.glm1 = summary(glm1)
exp(sum.glm1$coefficients)


##------------------------------------------------------------------------------------------##
##--------------------------------------Generate figures------------------------------------##
##------------------------------------------------------------------------------------------##

##--------------------------------------- Figure 2 -----------------------------------------##
## Permutation analysis at different Ct thresholds. For each Ct threshold, Ct values were permuted 500 times, and prevalence was re-estimated (thin grey lines). 
## For each time-point, the quantile of observed prevalence (bold black lines) was calculated as the proportion of permuted prevalence values that was smaller than the observed value. 
## Background colors indicate prevalence quantiles for each time-point. When Ct distributions are similar over time, prevalence quantiles are expected to be around 50%, whereas changing Ct distributions are expected to result in low and high quantiles. 
## Prevalence is pooled for all sites and normalized to allow comparison between thresholds. 
## Cases of spillover are shown with crosses. 

## Plot functions:
plot.left.not.bottom = function() {
  
  cur.colors = colors[sort(unique(obs.out$quantile.cat[which(obs.out$ct_threshold==cur.ct)]))]
  
  ggplot() + 
    #geom_rect(data = obs.out[which(obs.out$ct_threshold==cur.ct),], aes(xmin = xmin, xmax = xmax, ymin = -0.5, ymax = Inf), fill = obs.out[which(obs.out$ct_threshold==cur.ct),"fill.col"]) +
    geom_rect(data = obs.out[which(obs.out$ct_threshold==cur.ct),], aes(xmin = xmin, xmax = xmax, ymin = -0.5, ymax = Inf, fill = factor(quantile.cat)), alpha = 0.8) +
    scale_fill_manual(values = cur.colors, name = "Prevalence\nQuantile", labels = c("0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%")) +
    geom_line(data=perm.out[which(perm.out$ct_threshold==cur.ct),],
              stat = "smooth",
              aes(y = hen_prevalence_normalized, x = sampling_cluster_median_date, group = permutation),
              method = "glm",
              method.args = list(family = binomial),
              formula = y ~ splines::ns(x, df = 8),
              se = F, colour = color.perm.prev, size = 0.5, span = 0.5, alpha = 0.2) +
    geom_smooth(data=obs.out[which(obs.out$ct_threshold==cur.ct),],
                aes(y = hen_prevalence_normalized, x = sampling_cluster_median_date),
                method = "glm",
                method.args = list(family = binomial),
                formula = y ~ splines::ns(x, df = 8),
                se = F, colour = color.obs.prev, size = 3, span = 0.1) +
    scale_x_date(date_labels = "%b\n%Y",
                 #date_breaks = "3 months",
                 breaks = seq.Date(from = as.Date("2017-08-01"), to = as.Date("2020-08-01"), by = "6 months"),
                 date_minor_breaks = "month") + 
    # theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    labs(x = "Sampling date",
         y = "Normalized prevalence") +
    theme_bw() +
    geom_point(data=spillover, aes(y = -0.02, x = dates_for_plots), alpha=1, colour="black", shape=4, size=4, stroke = 1) + 
    background_grid("none") +
    ggtitle(paste0("Ct threshold = ",cur.ct)) + 
    coord_cartesian(ylim=c(-0.1,1.0), xlim = as.Date(c("2017-07-01", "2020-10-01"))) + 
    theme(axis.text.y = element_text(size=axis.text.size), 
          axis.text.x = element_blank(), 
          axis.title.x = element_blank(), 
          axis.title.y = element_text(size=axis.title.size),
          legend.title = element_text(size=16), 
          legend.text = element_text(size=12), 
          legend.position = "none",
          plot.title = element_text(size=plot.title.size)) + 
    labs(subtitle = NULL) + 
    scale_y_continuous(expand=c(0,0))
}

plot.left.bottom = function() {
  
  cur.colors = colors[sort(unique(obs.out$quantile.cat[which(obs.out$ct_threshold==cur.ct)]))]
  
  ggplot() + 
    #geom_rect(data = obs.out[which(obs.out$ct_threshold==cur.ct),], aes(xmin = xmin, xmax = xmax, ymin = -0.5, ymax = Inf), fill = obs.out[which(obs.out$ct_threshold==cur.ct),"fill.col"]) +
    geom_rect(data = obs.out[which(obs.out$ct_threshold==cur.ct),], aes(xmin = xmin, xmax = xmax, ymin = -0.5, ymax = Inf, fill = factor(quantile.cat)), alpha = 0.8) +
    scale_fill_manual(values = cur.colors, name = "Prevalence\nQuantile", labels = c("0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%")) +
    geom_line(data=perm.out[which(perm.out$ct_threshold==cur.ct),],
              stat = "smooth",
              aes(y = hen_prevalence_normalized, x = sampling_cluster_median_date, group = permutation),
              method = "glm",
              method.args = list(family = binomial),
              formula = y ~ splines::ns(x, df = 8),
              se = F, colour = color.perm.prev, size = 0.5, span = 0.5, alpha = 0.2) +
    geom_smooth(data=obs.out[which(obs.out$ct_threshold==cur.ct),],
                aes(y = hen_prevalence_normalized, x = sampling_cluster_median_date),
                method = "glm",
                method.args = list(family = binomial),
                formula = y ~ splines::ns(x, df = 8),
                se = F, colour = color.obs.prev, size = 3, span = 0.1) +
    scale_x_date(date_labels = "%b\n%Y",
                 #date_breaks = "3 months",
                 breaks = seq.Date(from = as.Date("2017-08-01"), to = as.Date("2020-08-01"), by = "6 months"),
                 date_minor_breaks = "month") + 
    # theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    labs(x = "Sampling date",
         y = "Normalized prevalence") +
    theme_bw() +
    geom_point(data=spillover, aes(y = -0.02, x = dates_for_plots), alpha=1, colour="black", shape=4, size=4, stroke = 1) + 
    background_grid("none") +
    ggtitle(paste0("Ct threshold = ",cur.ct)) + 
    coord_cartesian(ylim=c(-0.1,1.0), xlim = as.Date(c("2017-07-01", "2020-10-01"))) + 
    theme(axis.text.y = element_text(size=axis.text.size), 
          axis.text.x = element_text(size=axis.text.size), 
          axis.title.x = element_text(size=axis.title.size), 
          axis.title.y = element_text(size=axis.title.size),
          legend.title = element_text(size=16), 
          legend.text = element_text(size=12), 
          legend.position = "none",
          plot.title = element_text(size=plot.title.size)) + 
    labs(subtitle = NULL) + 
    scale_y_continuous(expand=c(0,0))
}

plot.middle.not.bottom = function() {
  
  cur.colors = colors[sort(unique(obs.out$quantile.cat[which(obs.out$ct_threshold==cur.ct)]))]
  
  ggplot() + 
    #geom_rect(data = obs.out[which(obs.out$ct_threshold==cur.ct),], aes(xmin = xmin, xmax = xmax, ymin = -0.5, ymax = Inf), fill = obs.out[which(obs.out$ct_threshold==cur.ct),"fill.col"]) +
    geom_rect(data = obs.out[which(obs.out$ct_threshold==cur.ct),], aes(xmin = xmin, xmax = xmax, ymin = -0.5, ymax = Inf, fill = factor(quantile.cat)), alpha = 0.8) +
    scale_fill_manual(values = cur.colors, name = "Prevalence\nQuantile", labels = c("0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%")) +
    geom_line(data=perm.out[which(perm.out$ct_threshold==cur.ct),],
              stat = "smooth",
              aes(y = hen_prevalence_normalized, x = sampling_cluster_median_date, group = permutation),
              method = "glm",
              method.args = list(family = binomial),
              formula = y ~ splines::ns(x, df = 8),
              se = F, colour = color.perm.prev, size = 0.5, span = 0.5, alpha = 0.2) +
    geom_smooth(data=obs.out[which(obs.out$ct_threshold==cur.ct),],
                aes(y = hen_prevalence_normalized, x = sampling_cluster_median_date),
                method = "glm",
                method.args = list(family = binomial),
                formula = y ~ splines::ns(x, df = 8),
                se = F, colour = color.obs.prev, size = 3, span = 0.1) +
    scale_x_date(date_labels = "%b\n%Y",
                 #date_breaks = "3 months",
                 breaks = seq.Date(from = as.Date("2017-08-01"), to = as.Date("2020-08-01"), by = "6 months"),
                 date_minor_breaks = "month") + 
    # theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    labs(x = "Sampling date",
         y = "Normalized prevalence") +
    theme_bw() +
    geom_point(data=spillover, aes(y = -0.02, x = dates_for_plots), alpha=1, colour="black", shape=4, size=4, stroke = 1) + 
    background_grid("none") +
    ggtitle(paste0("Ct threshold = ",cur.ct)) + 
    coord_cartesian(ylim=c(-0.1,1.0), xlim = as.Date(c("2017-07-01", "2020-10-01"))) + 
    theme(axis.text.y = element_blank(), 
          axis.text.x = element_blank(), 
          axis.title.x = element_blank(), 
          axis.title.y = element_blank(),
          legend.title = element_text(size=16), 
          legend.text = element_text(size=12), 
          legend.position = "none",
          plot.title = element_text(size=plot.title.size)) + 
    labs(subtitle = NULL) + 
    scale_y_continuous(expand=c(0,0))
}

plot.right.not.bottom = function() {
  
  cur.colors = colors[sort(unique(obs.out$quantile.cat[which(obs.out$ct_threshold==cur.ct)]))]
  
  ggplot() + 
    #geom_rect(data = obs.out[which(obs.out$ct_threshold==cur.ct),], aes(xmin = xmin, xmax = xmax, ymin = -0.5, ymax = Inf), fill = obs.out[which(obs.out$ct_threshold==cur.ct),"fill.col"]) +
    geom_rect(data = obs.out[which(obs.out$ct_threshold==cur.ct),], aes(xmin = xmin, xmax = xmax, ymin = -0.5, ymax = Inf, fill = factor(quantile.cat)), alpha = 0.8) +
    scale_fill_manual(values = cur.colors, name = "Prevalence\nQuantile", labels = c("0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%")) +
    geom_line(data=perm.out[which(perm.out$ct_threshold==cur.ct),],
              stat = "smooth",
              aes(y = hen_prevalence_normalized, x = sampling_cluster_median_date, group = permutation),
              method = "glm",
              method.args = list(family = binomial),
              formula = y ~ splines::ns(x, df = 8),
              se = F, colour = color.perm.prev, size = 0.5, span = 0.5, alpha = 0.2) +
    geom_smooth(data=obs.out[which(obs.out$ct_threshold==cur.ct),],
                aes(y = hen_prevalence_normalized, x = sampling_cluster_median_date),
                method = "glm",
                method.args = list(family = binomial),
                formula = y ~ splines::ns(x, df = 8),
                se = F, colour = color.obs.prev, size = 3, span = 0.1) +
    scale_x_date(date_labels = "%b\n%Y",
                 #date_breaks = "3 months",
                 breaks = seq.Date(from = as.Date("2017-08-01"), to = as.Date("2020-08-01"), by = "6 months"),
                 date_minor_breaks = "month") + 
    # theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    labs(x = "Sampling date",
         y = "Normalized prevalence") +
    theme_bw() +
    geom_point(data=spillover, aes(y = -0.02, x = dates_for_plots), alpha=1, colour="black", shape=4, size=4, stroke = 1) + 
    background_grid("none") +
    ggtitle(paste0("Ct threshold = ",cur.ct)) + 
    coord_cartesian(ylim=c(-0.1,1.0), xlim = as.Date(c("2017-07-01", "2020-10-01"))) + 
    theme(axis.text.y = element_blank(), 
          axis.text.x = element_blank(), 
          axis.title.x = element_blank(), 
          axis.title.y = element_blank(),
          legend.title = element_text(size=legend.title.size), 
          legend.text = element_text(size=legend.text.size), 
          legend.position = "none",
          plot.title = element_text(size=plot.title.size)) + 
    labs(subtitle = NULL) + 
    scale_y_continuous(expand=c(0,0))
}

plot.right.bottom = function() {
  
  cur.colors = colors[sort(unique(obs.out$quantile.cat[which(obs.out$ct_threshold==cur.ct)]))]
  
  ggplot() + 
    #geom_rect(data = obs.out[which(obs.out$ct_threshold==cur.ct),], aes(xmin = xmin, xmax = xmax, ymin = -0.5, ymax = Inf), fill = obs.out[which(obs.out$ct_threshold==cur.ct),"fill.col"]) +
    geom_rect(data = obs.out[which(obs.out$ct_threshold==cur.ct),], aes(xmin = xmin, xmax = xmax, ymin = -0.5, ymax = Inf, fill = factor(quantile.cat)), alpha = 0.8) +
    scale_fill_manual(values = cur.colors, name = "Prevalence\nQuantile", labels = c("0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%")) +
    geom_line(data=perm.out[which(perm.out$ct_threshold==cur.ct),],
              stat = "smooth",
              aes(y = hen_prevalence_normalized, x = sampling_cluster_median_date, group = permutation),
              method = "glm",
              method.args = list(family = binomial),
              formula = y ~ splines::ns(x, df = 8),
              se = F, colour = color.perm.prev, size = 0.5, span = 0.5, alpha = 0.2) +
    geom_smooth(data=obs.out[which(obs.out$ct_threshold==cur.ct),],
                aes(y = hen_prevalence_normalized, x = sampling_cluster_median_date),
                method = "glm",
                method.args = list(family = binomial),
                formula = y ~ splines::ns(x, df = 8),
                se = F, colour = color.obs.prev, size = 3, span = 0.1) +
    scale_x_date(date_labels = "%b\n%Y",
                 #date_breaks = "3 months",
                 breaks = seq.Date(from = as.Date("2017-08-01"), to = as.Date("2020-08-01"), by = "6 months"),
                 date_minor_breaks = "month") + 
    # theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    labs(x = "Sampling date",
         y = "Normalized prevalence") +
    theme_bw() +
    geom_point(data=spillover, aes(y = -0.02, x = dates_for_plots), alpha=1, colour="black", shape=4, size=4, stroke = 1) + 
    background_grid("none") +
    ggtitle(paste0("Ct threshold = ",cur.ct)) + 
    coord_cartesian(ylim=c(-0.1,1.0), xlim = as.Date(c("2017-07-01", "2020-10-01"))) + 
    theme(axis.text.y = element_blank(), 
          axis.text.x = element_text(size = axis.text.size), 
          axis.title.x = element_text(size = axis.title.size), 
          axis.title.y = element_blank(),
          legend.title = element_text(size=legend.title.size), 
          legend.text = element_text(size=legend.text.size), 
          legend.position = "none",
          plot.title = element_text(size=plot.title.size)) + 
    labs(subtitle = NULL) + 
    scale_y_continuous(expand=c(0,0))
}

plot.legend.only = function() {
  
  cur.plot = ggplot() + 
    #geom_rect(data = obs.out[which(obs.out$ct_threshold==cur.ct),], aes(xmin = xmin, xmax = xmax, ymin = -0.5, ymax = Inf), fill = obs.out[which(obs.out$ct_threshold==cur.ct),"fill.col"]) +
    geom_rect(data = obs.out[which(obs.out$ct_threshold==cur.ct),], aes(xmin = xmin, xmax = xmax, ymin = -0.5, ymax = Inf, fill = factor(quantile.cat)), alpha = 0.8) +
    scale_fill_manual(values = colors, name = "Prevalence\nQuantile", labels = c("0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%")) +
    geom_line(data=perm.out[which(perm.out$ct_threshold==cur.ct),],
              stat = "smooth",
              aes(y = hen_prevalence_normalized, x = sampling_cluster_median_date, group = permutation),
              method = "glm",
              method.args = list(family = binomial),
              formula = y ~ splines::ns(x, df = 8),
              se = F, colour = color.perm.prev, size = 0.5, span = 0.5, alpha = 0.2) +
    geom_smooth(data=obs.out[which(obs.out$ct_threshold==cur.ct),],
                aes(y = hen_prevalence_normalized, x = sampling_cluster_median_date),
                method = "glm",
                method.args = list(family = binomial),
                formula = y ~ splines::ns(x, df = 8),
                se = F, colour = color.obs.prev, size = 2.5, span = 0.1) +
    scale_x_date(date_labels = "%b\n%Y",
                 #date_breaks = "3 months",
                 breaks = seq.Date(from = as.Date("2017-08-01"), to = as.Date("2020-08-01"), by = "6 months"),
                 date_minor_breaks = "month") + 
    # theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    labs(x = "Sampling date",
         y = "Normalized prevalence") +
    theme_bw() +
    geom_point(data=spillover, aes(y = -0.02, x = dates_for_plots), alpha=1, colour="black", shape=4, size=4, stroke = 1) + 
    background_grid("none") +
    ggtitle(paste0("Ct threshold = ",cur.ct)) + 
    coord_cartesian(ylim=c(-0.1,1.0), xlim = as.Date(c("2017-07-01", "2020-10-01"))) + 
    theme(axis.text.y = element_blank(), 
          axis.text.x = element_text(size = axis.text.size), 
          axis.title.x = element_text(size = axis.title.size), 
          axis.title.y = element_blank(),
          legend.title = element_text(size=legend.title.size), 
          legend.text = element_text(size=legend.text.size), 
          plot.title = element_text(size=plot.title.size)) + 
    labs(subtitle = NULL) + 
    scale_y_continuous(expand=c(0,0))
  
  return(as_ggplot(get_legend(cur.plot)))
  
}

## Set aesthetics:
plot.title.size = 32
axis.text.size = 28
axis.title.size = 30
legend.title.size = 32
legend.text.size = 30
color.perm.prev = "grey30"
color.obs.prev = "black"
    
## Create figure:
colors = rev(brewer.pal(n = 10,"RdYlBu"))
  
cur.ct = 28;plot.28 = plot.left.not.bottom()
cur.ct = 29;plot.29 = plot.middle.not.bottom()
cur.ct = 30;plot.30 = plot.right.not.bottom()
cur.ct = 31;plot.31 = plot.left.not.bottom()
cur.ct = 32;plot.32 = plot.middle.not.bottom()
cur.ct = 33;plot.33 = plot.right.not.bottom()
cur.ct = 34;plot.34 = plot.left.not.bottom()
cur.ct = 35;plot.35 = plot.middle.not.bottom()
cur.ct = 36;plot.36 = plot.right.not.bottom()
cur.ct = 37;plot.37 = plot.left.bottom()
cur.ct = 38;plot.38 = plot.right.bottom()
cur.ct = 39;plot.39 = plot.right.bottom()
cur.ct = 35;plot.legend = plot.legend.only()
  
layout = "
AAABBBCCC#
AAABBBCCC#
DDDEEEFFFM
DDDEEEFFFM
GGGHHHIIIM
GGGHHHIIIM
JJJKKKLLL#
JJJKKKLLL#
"
plot.out = plot.28 + plot.29 + plot.30 + plot.31 + plot.32 + plot.33 + plot.34 + plot.35 + plot.36 + plot.37 + plot.38 + plot.39 + plot.legend +
plot_layout(design = layout)

plot.out
# ggsave(plot = plot.out,filename = "Figure 2 20250210.png",width = 30, height = 24, dpi=300)
  

##--------------------------------------- Figure 3 -----------------------------------------##
## Normalized pooled prevalence for different Ct threshold values across off-season months (September-May) with no cases of spillover, and peak season (June-August) without and with cases of spillover. 
## Boxplots show the 25%, 50%, 75% percentiles, lines indicate the smallest and largest values within 1.5 times the interquartile range, dots indicate values beyond that. 

## Set aesthetics:
plot.title.size = 18
axis.text.size = 22
axis.title.size = 22
legend.title.size = 20
legend.text.size = 20

## Create figure:
plot.out = obs.out.f4 %>% 
  ggplot() +
  geom_boxplot(aes(x = spillover_season, y = hen_prevalence_normalized_smooth, fill = factor(ct_threshold)), color = "grey10") +
  scale_fill_manual(values = c("#6f112b",brewer.pal(n = 11,"RdYlBu"),"#000061")) +
  scale_x_discrete(labels = c("No spillover\nOff season","No spillover\nPeak season","Spillover\nPeak season")) +
  labs(fill = "Ct\nThreshold",
       x = "",
       y = "Normalized pooled prevalence") +
  theme_bw() +
  background_grid("none") +
  theme(axis.text.y = element_text(size = axis.text.size), 
        axis.text.x = element_text(size = axis.text.size), 
        axis.title.x = element_text(size = axis.title.size), 
        axis.title.y = element_text(size = axis.title.size),
        legend.title = element_text(size=legend.title.size), 
        legend.text = element_text(size=legend.text.size), 
        plot.title = element_text(size=plot.title.size))

plot.out
# ggsave(plot = plot.out,filename = "Figure 3 20250210.png",width = 16, height = 8, dpi=300)


##--------------------------------------- Figure 4 -----------------------------------------##
## Akaike Information Criterion (AIC) values for models with normalized pooled prevalence as predictor variable and number (A) or occurrence (B) of spillover events as outcome variable, for a range of Ct thresholds. 
## AIC values are shown in blue for permuted data and red for observed data. 
## Insets show the proportion of permuted AIC values that is larger than the observed value. Low quantile values provide statistical support for a shift in the distribution of Ct values that results in a better correlation between prevalence and spillover.

plot1 = ggplot() +
  geom_jitter(data = perm.out.stats.glm.poisson, aes(x = ct_threshold, y = AIC), alpha = 0.5, col = brewer.pal(n = 11,"RdYlBu")[9], shape = 21, width = 0.1) +
  geom_point(data = obs.out.stats.glm.poisson, aes(x = ct_threshold, y = AIC), col = brewer.pal(n = 11,"RdYlBu")[1], size = 3) +
  scale_x_continuous(breaks = 28:40) +
  xlab("Ct Threshold") +
  ylab("AIC") +
  theme_bw() +
  background_grid("none") +
  ggtitle("Number of Spillover Events") + 
  theme(axis.text.y = element_text(size=20), 
        axis.text.x = element_text(size=20), 
        axis.title.x = element_text(size=22), 
        axis.title.y = element_text(size=22),
        legend.title = element_text(size=20), 
        legend.text = element_text(size=20), 
        plot.title = element_text(size=22), 
        strip.text = element_text(size = 22), 
        legend.position = "none",
        plot.tag = element_text(size=22)) +
  labs(tag = "A")

plot2 = ggplot() +
  geom_point(data = obs.out.stats.glm.poisson[-which(obs.out.stats.glm.poisson$ct_threshold==40),], aes(x = ct_threshold, y = AIC.quantile), col = brewer.pal(n = 11,"RdYlBu")[1], size = 3) +
  scale_x_continuous(breaks = seq(28,40,2), limits = c(28,40)) +
  scale_y_continuous(limits = c(0,1)) +
  xlab("Ct Threshold") +
  ylab("") +
  theme_bw() +
  background_grid("none") +
  ggtitle("AIC Permutation Quantile") + 
  theme(axis.text.y = element_text(size=16), 
        axis.text.x = element_text(size=16), 
        axis.title.x = element_text(size=15), 
        axis.title.y = element_text(size=15),
        plot.title = element_text(size=15), 
        strip.text = element_text(size = 18), 
        legend.position = "none",
        plot.tag = element_text(size=22),
        panel.background = element_rect(fill = "grey90"))

plot3 = ggplot() +
  geom_jitter(data = perm.out.stats.glm.binomial, aes(x = ct_threshold, y = AIC), alpha = 0.2, col = brewer.pal(n = 11,"RdYlBu")[9], shape = 21, width = 0.1) +
  geom_point(data = obs.out.stats.glm.binomial, aes(x = ct_threshold, y = AIC), col = brewer.pal(n = 11,"RdYlBu")[1], size = 3) +
  scale_x_continuous(breaks = 28:40) +
  xlab("Ct Threshold") +
  ylab("AIC") +
  theme_bw() +
  background_grid("none") +
  ggtitle("Occurrence of Spillover Events") + 
  theme(axis.text.y = element_text(size=20), 
        axis.text.x = element_text(size=20), 
        axis.title.x = element_text(size=22), 
        axis.title.y = element_text(size=22),
        legend.title = element_text(size=20), 
        legend.text = element_text(size=20), 
        plot.title = element_text(size=22), 
        strip.text = element_text(size = 22), 
        legend.position = "none",
        plot.tag = element_text(size=22)) +
  labs(tag = "B")

plot4 = ggplot() +
  geom_point(data = obs.out.stats.glm.binomial[-which(obs.out.stats.glm.binomial$ct_threshold==40),], aes(x = ct_threshold, y = AIC.quantile), col = brewer.pal(n = 11,"RdYlBu")[1], size = 3) +
  scale_x_continuous(breaks = seq(28,40,2), limits = c(28,40)) +
  scale_y_continuous(limits = c(0,1)) +
  xlab("Ct Threshold") +
  ylab("") +
  theme_bw() +
  background_grid("none") +
  ggtitle("AIC Permutation Quantile") + 
  theme(axis.text.y = element_text(size=16), 
        axis.text.x = element_text(size=16), 
        axis.title.x = element_text(size=15), 
        axis.title.y = element_text(size=15),
        plot.title = element_text(size=15), 
        strip.text = element_text(size = 18), 
        legend.position = "none",
        plot.tag = element_text(size=22),
        panel.background = element_rect(fill = "grey90"))

plotA = plot1 + inset_element(
  plot2, 
  left = 0.45, 
  bottom = 0.45, 
  right = unit(1, 'npc') - unit(0.7, 'cm'), 
  top = unit(1, 'npc') - unit(0.5, 'cm')
)

plotB = plot3 + inset_element(
  plot4, 
  left = 0.45, 
  bottom = 0.45, 
  right = unit(1, 'npc') - unit(0.7, 'cm'), 
  top = unit(1, 'npc') - unit(0.5, 'cm')
)

plot.out = plotA + plotB

plot.out
# ggsave(plot = plot.out,filename = "Figure 4 20250210.png",width = 15, height = 6, dpi=600)

