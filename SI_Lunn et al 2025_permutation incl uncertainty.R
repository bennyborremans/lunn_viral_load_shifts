## Title: Periodic shifts in viral load increase risk of spillover from bats
##        Supplement - permutation analysis accounting for uncertainty in prevalence estimates  
## Code authors: Tamika Lunn, Griffith University & Benny Borremans, Wildlife Health Ecology Research Organization
## Version: V5, created 16th July 2025

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
prepdata <- function(urv_data, threshold, incl_prev_uncertainty = F){ 
  ## Start data summaries:  
  results_w_meta <- dplyr::mutate(urv_data,
                                  buffer_pos_hev = dplyr::if_else(ct_hev > 0, preservative, ""))
  results_w_meta <- dplyr::filter(results_w_meta, site %in% continuous_ur_s$site) 
  
  ### updated July 30 2024 to correct numbers over sheets ###
  results_w_meta <- results_w_meta %>% 
    mutate(nnonbff_over_sheet = ifelse(is.na(nnonbff_over_sheet)==TRUE, 0, nnonbff_over_sheet))## NA values should be 0; if there's a BFF count then sheets were counted
  ### end correction ###
  
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
    ## Nml_HeV throwing an error here (creating NA values instead of retaining actual values, can check with accession update ARTOW008 and the sample with Ct=40). Re-merge with conversion sheets:
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
                                   hen_prevalence = ifelse(incl_prev_uncertainty, rbeta(1, 2 + sum(hen_bin), 1 + n() - sum(hen_bin)), sum(hen_bin)/dplyr::n()),  # if uncertainty is included, prevalence is sampled from a beta distribution with a weakly informative prior (lower values slightly more likely). should only used for permutation analysis
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
                                            hen_prevalence = ifelse(incl_prev_uncertainty, rbeta(1, 2 + sum(hen_bin), 1 + n() - sum(hen_bin)), sum(hen_bin)/dplyr::n()),   # if uncertainty is included, prevalence is sampled from a beta distribution with a weakly informative prior (lower values slightly more likely). should only used for permutation analysis
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

continuous_ur_s <- readRDS("Data/Processed/continuous_ur_s.rds") ## Names and locations of sites
urv_data <- readRDS("Data/Processed/urv_data.rds")   ## Urine sampling meta-data with virus data
genomeHev <- read.csv("Data/Processed/Ct-to-genome-convert_HeV.csv") ## Ct to genome per ml conversion
spillover <- readRDS("Data/Processed/spillover.rds") ## spillovers within study

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
run.date = "20250716"


if(run.again) {      
  
  options(dplyr.summarise.inform = FALSE)  # silence summarise function warnings  
  
  # prepare dataframe for permutation prevalence data
  perm.out = data.frame(sampling_cluster_median_date = rep(hen_prev_all$sampling_cluster_median_date, permutations * length(ct_thresholds)),
                        permutation = rep(1:permutations, each = nrow(hen_prev_all) * length(ct_thresholds)),
                        hen_prevalence_normalized = numeric(length(nrow(hen_prev_all) * permutations * length(ct_thresholds))),
                        hen_prevalence_normalized_smooth = numeric(length(nrow(hen_prev_all) * permutations * length(ct_thresholds))),
                        ct_threshold = rep(rep(ct_thresholds, each = nrow(hen_prev_all)) ,permutations))
  
  # prepare dataframe for observed prevalence data inside permutation
  obs.out.perm = perm.out
  
  
  # prepare dataframe for observed prevalence data
  obs.out = data.frame(sampling_cluster_median_date = rep(hen_prev_all$sampling_cluster_median_date,length(ct_thresholds)),
                       hen_prevalence_normalized = numeric(length(nrow(hen_prev_all) * length(ct_thresholds))),
                       hen_prevalence_normalized_smooth = numeric(length(nrow(hen_prev_all) * length(ct_thresholds))),
                       ct_threshold = rep(ct_thresholds, each = nrow(hen_prev_all)),
                       prev_quantile = numeric(length(nrow(hen_prev_all) * length(ct_thresholds))),
                       hen_prevalence_normalized_lower = numeric(length(nrow(hen_prev_all) * length(ct_thresholds))),
                       hen_prevalence_normalized_upper = numeric(length(nrow(hen_prev_all) * length(ct_thresholds))),
                       hen_prevalence_normalized_smooth_lower = numeric(length(nrow(hen_prev_all) * length(ct_thresholds))),
                       hen_prevalence_normalized_smooth_upper = numeric(length(nrow(hen_prev_all) * length(ct_thresholds))))
  
  
  
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
      Ct_data_perm = prepdata(urv_data_perm, cur.ct, incl_prev_uncertainty = T)
      hen_with_ci_by_month_perm <- Ct_data_perm[[2]] %>%
        mutate(hen_prevalence_normalized = hen_prevalence/max(hen_prevalence))
      
      # smooth prevalence
      perm.out$hen_prevalence_normalized[which(perm.out$permutation==i & perm.out$ct_threshold==cur.ct)] = hen_with_ci_by_month_perm$hen_prevalence_normalized
      prev.smooth = gam(hen_prevalence_normalized ~ ns(sampling_cluster_median_date, df = 8), family = "binomial", data = hen_with_ci_by_month_perm)
      perm.out$hen_prevalence_normalized_smooth[which(perm.out$permutation==i & perm.out$ct_threshold==cur.ct)] = predict(prev.smooth,type = "response")
      
      
      
      # calculate the same stats for the observed data
      
      Ct_data_obs = prepdata(urv_data, cur.ct, incl_prev_uncertainty = T)
      hen_with_ci_by_month_obs <- Ct_data_obs[[2]] %>%
        mutate(hen_prevalence_normalized = hen_prevalence/max(hen_prevalence))
      
      # smoothed prevalence
      obs.out.perm$hen_prevalence_normalized[which(obs.out.perm$permutation==i & obs.out.perm$ct_threshold==cur.ct)] = hen_with_ci_by_month_obs$hen_prevalence_normalized
      prev.smooth = gam(hen_prevalence_normalized ~ ns(sampling_cluster_median_date, df = 8), family = "binomial", data = hen_with_ci_by_month_obs)
      obs.out.perm$hen_prevalence_normalized_smooth[which(obs.out.perm$permutation==i & obs.out.perm$ct_threshold==cur.ct)] = predict(prev.smooth,type = "response")
      
      
      if(i %% 50 == 0) cat(paste0("\r--- Ct threshold ",cur.ct,"  -  Permutation ",i))
    }
    
    
    
    Ct_data = prepdata(urv_data, cur.ct, incl_prev_uncertainty = F)
    hen_with_ci_by_month <- Ct_data[[2]] %>%
      mutate(hen_prevalence_normalized = hen_prevalence/max(hen_prevalence),
             hen_prevalence_normalized_lower = hen_lower/max(hen_prevalence),
             hen_prevalence_normalized_upper = hen_upper/max(hen_prevalence))
    # smoothed prevalence
    obs.out$hen_prevalence_normalized[which(obs.out$ct_threshold==cur.ct)] = hen_with_ci_by_month$hen_prevalence_normalized
    prev.smooth = gam(hen_prevalence_normalized ~ ns(sampling_cluster_median_date, df = 8), family = "binomial", data = hen_with_ci_by_month)
    obs.out$hen_prevalence_normalized_smooth[which(obs.out$ct_threshold==cur.ct)] = predict(prev.smooth,type = "response")
    
    
    # calculate proportion of observed prevalences < permuted prevalence
    # this is done using smoothed prevalence
    
    cur.obs.out.idx = which(obs.out$ct_threshold==cur.ct)
    for(k in 1:length(cur.obs.out.idx)){
      cur.date.obs = obs.out$sampling_cluster_median_date[which(obs.out$ct_threshold==cur.ct)[k]]
      cur.prev.obs = obs.out.perm$hen_prevalence_normalized_smooth[which(obs.out.perm$ct_threshold==cur.ct & perm.out$sampling_cluster_median_date==cur.date.obs)]
      cur.prev.perm = perm.out$hen_prevalence_normalized_smooth[which(perm.out$ct_threshold==cur.ct & perm.out$sampling_cluster_median_date==cur.date.obs)]
      
      obs.out$prev_quantile[cur.obs.out.idx[k]] = sum(cur.prev.perm < cur.prev.obs)/permutations
    }
    
  }
  
  # save output so it doesn't have the be run again
  saveRDS(perm.out, paste0("perm.out.incl.uncertainty",run.date,".RDS"))
  
  saveRDS(obs.out.perm, paste0("obs.out.perm.incl.uncertainty",run.date,".RDS"))
  
  saveRDS(obs.out, paste0("obs.out.incl.uncertainty",run.date,".RDS"))
} else {
  perm.out = readRDS(paste0("perm.out.incl.uncertainty",run.date,".RDS"))
  obs.out.perm = readRDS(paste0("obs.out.perm.incl.uncertainty",run.date,".RDS"))
  obs.out = readRDS(paste0("obs.out.incl.uncertainty",run.date,".RDS"))
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

obs.out$ci.low = NA
obs.out$ci.high = NA
# add observed prevalence uncertainty min/max values for bands
for(i in 1:nrow(obs.out)){
  cur.dat <- urv_data %>% filter(sampling_cluster_median_date == obs.out$sampling_cluster_median_date[i])  
  cur.pos = sum(cur.dat$ct_hev>0)
  cur.n = nrow(cur.dat)
  obs.out[i, c("ci.low", "ci.high")] = qbeta(p = c(0.05,0.95),2 + cur.pos, 1 + cur.n - cur.pos)
}


##--------------------------------------- Figure S3 -----------------------------------------##
## Analogue to Figure 2 in the main text
## Permutation analysis at different Ct thresholds. For each Ct threshold, Ct values were permuted 500 times, and prevalence was re-estimated (thin grey lines). 
## For each time-point, the quantile of observed prevalence (bold black lines) was calculated as the proportion of permuted prevalence values that was smaller than the observed value. 
## Background colors indicate prevalence quantiles for each time-point. When Ct distributions are similar over time, prevalence quantiles are expected to be around 50%, whereas changing Ct distributions are expected to result in low and high quantiles. 
## Prevalence is pooled for all sites and normalized to allow comparison between thresholds. 
## Cases of spillover are shown with crosses. 

obs.out$quantile.cat = factor(obs.out$quantile.cat, levels = c(1,2,3,4,5,6,7,8,9,10))

## Plot functions:
plot.left.not.bottom = function() {
  
  cur.colors = colors[sort(unique(obs.out$quantile.cat[which(obs.out$ct_threshold==cur.ct)]))]
  
  ggplot() + 
    #geom_rect(data = obs.out[which(obs.out$ct_threshold==cur.ct),], aes(xmin = xmin, xmax = xmax, ymin = -0.5, ymax = Inf), fill = obs.out[which(obs.out$ct_threshold==cur.ct),"fill.col"]) +
    geom_rect(data = obs.out[which(obs.out$ct_threshold==cur.ct),], aes(xmin = xmin, xmax = xmax, ymin = -0.5, ymax = Inf, fill = factor(quantile.cat)), alpha = 0.8) +
    scale_fill_manual(values = cur.colors, name = "Prevalence\nQuantile", labels = c("0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%"), drop = F) +
    geom_line(data=perm.out[which(perm.out$ct_threshold==cur.ct),],
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
    geom_point(data=spillover %>% filter(all_spillovers_ID == "all_059"), aes(y = -0.06, x = dates_for_plots), alpha=1, colour="black", shape=4, size=4, stroke = 1) + 
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
    scale_fill_manual(values = cur.colors, name = "Prevalence\nQuantile", labels = c("0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%"), drop = F) +
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
    geom_point(data=spillover %>% filter(all_spillovers_ID == "all_059"), aes(y = -0.06, x = dates_for_plots), alpha=1, colour="black", shape=4, size=4, stroke = 1) + 
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
    scale_fill_manual(values = cur.colors, name = "Prevalence\nQuantile", labels = c("0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%"), drop = F) +
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
    geom_point(data=spillover %>% filter(all_spillovers_ID == "all_059"), aes(y = -0.06, x = dates_for_plots), alpha=1, colour="black", shape=4, size=4, stroke = 1) + 
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
    scale_fill_manual(values = cur.colors, name = "Prevalence\nQuantile", labels = c("0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%"), drop = F) +
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
    geom_point(data=spillover %>% filter(all_spillovers_ID == "all_059"), aes(y = -0.06, x = dates_for_plots), alpha=1, colour="black", shape=4, size=4, stroke = 1) + 
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
    scale_fill_manual(values = cur.colors, name = "Prevalence\nQuantile", labels = c("0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%"), drop = F) +
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
    geom_point(data=spillover %>% filter(all_spillovers_ID == "all_059"), aes(y = -0.06, x = dates_for_plots), alpha=1, colour="black", shape=4, size=4, stroke = 1) + 
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




# create dummy data with all quantile levels, for correct legend plotting
dummy.cat = obs.out
dummy.cat$quantile.cat = sample(1:10, nrow(obs.out), replace = T)
dummy.cat$quantile.cat = factor(dummy.cat$quantile.cat, levels = c(1,2,3,4,5,6,7,8,9,10))


plot.legend.only = function() {
  
  cur.plot = ggplot() + 
    #geom_rect(data = obs.out[which(obs.out$ct_threshold==cur.ct),], aes(xmin = xmin, xmax = xmax, ymin = -0.5, ymax = Inf), fill = obs.out[which(obs.out$ct_threshold==cur.ct),"fill.col"]) +
    geom_rect(data = obs.out, aes(xmin = xmin, xmax = xmax, ymin = -0.5, ymax = Inf, fill = quantile.cat), alpha = 0.8) +
    geom_rect(data = dummy.cat, aes(xmin = xmin, xmax = xmax, ymin = -0.5, ymax = Inf, fill = quantile.cat), alpha = 0, show.legend = T) +
    scale_fill_manual(values = colors, name = "Prevalence\nQuantile", labels = c("0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%"), drop = F) +
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
    geom_point(data=spillover %>% filter(all_spillovers_ID == "all_059"), aes(y = -0.06, x = dates_for_plots), alpha=1, colour="black", shape=4, size=4, stroke = 1) + 
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
    scale_y_continuous(expand=c(0,0)) +
    guides(fill = guide_legend(override.aes = list(alpha = 1)))  # Ensure the legend colors are fully visible
  
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
# ggsave(plot = plot.out,filename = "Figure S3.png",width = 30, height = 24, dpi=300)

