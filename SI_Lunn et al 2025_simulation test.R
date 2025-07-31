## Title: Periodic shifts in viral load increase risk of spillover from bats
##        Supplement - permutation analysis test through simulation
## Code author: Benny Borremans, Wildlife Health Ecology Research Organization
## Version: V2, created 16th July 2025

rm(list=ls())

##------------------------------------------------------------------------------------------##
##---------------------------------------Load packages--------------------------------------##
##------------------------------------------------------------------------------------------##

## R version 4.3.2 ##
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

##----------------- Function: simulate data to test permutation analysis -------------------##
# for two possible prevalences per session, low or high


sim_for_perm = function(n_sessions, prev_low, prev_high, n_per_session, n_per_sample_low, n_per_sample_high, prev_ct_skew_threshold) {
        sim_out = data.frame(session = rep(1:n_sessions, each = n_per_session),
                             session_prev = rep(rnorm(n = n_sessions, sample(c(prev_low, prev_high),size = n_sessions, replace=T, prob = c(0.5,0.5)), sd = 0.02), each = n_per_session),
                             uniform_ct = NA,
                             skewed_ct = NA,
                             leftskewed_ct_all_prev = NA,
                             rightskewed_ct_all_prev = NA
        )
        
        sim_out$n_bats = sample(n_per_sample_low:n_per_sample_high, size = nrow(sim_out), replace = T)
        sim_out$pos_bats = apply(sim_out[,c("session_prev","n_bats")], 1, function(x) sum(rbinom(x[2], size = 1, prob = x[1])))
        
        for(i in 1:nrow(sim_out)){
                
                if(sim_out$pos_bats[i] > 0) {
                        # uniform
                        # log-uniform virus concentrations, independent of prevalence
                        cur_virus_conc = exp(runif(sim_out$pos_bats[i], log(60), log(7000000)))   # taking 60 and 7000000 as min and max concentrations
                        # calculate sample concentration, pos/all
                        cur_pooled_conc = sum(cur_virus_conc)/sim_out$n_bats[i]
                        # convert to Ct
                        sim_out$uniform_ct[i] = ceiling(46.606 - 1.619  * log(cur_pooled_conc))   # ceiling() to simulate PCR cycles
                        
                        
                        # left skewed Ct any prev
                        # log-skewed virus concentrations, independent of prevalence
                        cur_virus_conc = exp(log(60) + (log(7000000) - log(60)) * rbeta(sim_out$pos_bats[i], 1, 2))
                        # calculate sample concentration, pos/all
                        cur_pooled_conc = sum(cur_virus_conc)/sim_out$n_bats[i]
                        # convert to Ct
                        sim_out$leftskewed_ct_all_prev[i] = ceiling(46.606 - 1.619  * log(cur_pooled_conc))   # ceiling() to simulate PCR cycles
                        
                        
                        # right skewed Ct any prev
                        # log-skewed virus concentrations, independent of prevalence
                        cur_virus_conc = exp(log(60) + (log(7000000) - log(60)) * rbeta(sim_out$pos_bats[i], 2, 1))
                        # calculate sample concentration, pos/all
                        cur_pooled_conc = sum(cur_virus_conc)/sim_out$n_bats[i]
                        # convert to Ct
                        sim_out$rightskewed_ct_all_prev[i] = ceiling(46.606 - 1.619  * log(cur_pooled_conc))   # ceiling() to simulate PCR cycles
                        
                        
                        # skewed with prev
                        # log-skewed virus concentrations, prevalence-dependent
                        if(sim_out$session_prev[i] < prev_ct_skew_threshold) {
                                # prev low = lower concentrations more likely
                                cur_virus_conc = exp(log(60) + (log(7000000) - log(60)) * rbeta(sim_out$pos_bats[i], 1, 2))   # sampling from declining beta distribution, then converting to range of possible virus concentrations
                        } else {
                                # prev high = higher concentrations more likely
                                cur_virus_conc = exp(log(60) + (log(7000000) - log(60)) * rbeta(sim_out$pos_bats[i], 2, 1))   # sampling from declining beta distribution, then converting to range of possible virus concentrations
                        }
                        # calculate sample concentration, pos/all
                        cur_pooled_conc = sum(cur_virus_conc)/sim_out$n_bats[i]
                        # convert to Ct
                        sim_out$skewed_ct[i] = ceiling(46.606 - 1.619  * log(cur_pooled_conc))   # ceiling() to simulate PCR cycles
                } else {
                        sim_out$uniform_ct[i] = 0
                        sim_out$leftskewed_ct_all_prev[i] = 0
                        sim_out$rightskewed_ct_all_prev[i] = 0
                        sim_out$skewed_ct[i] = 0
                }
        }
        return(sim_out)
}

##------------------------- Function(s): viral load calculations ---------------------------##

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
##------------------------------------permutation analysis----------------------------------##
##------------------------------------------------------------------------------------------##



set.seed(11)
n_sessions <- 10
permutations <- 500
ct_thresholds <- 28:40
n_per_session <- 100


##------------------------------- uniform Ct distribution ----------------------------------##


# simulate data
sim_dat = sim_for_perm(n_sessions = n_sessions, prev_low = 0.05, prev_high = 0.2, n_per_session = n_per_session, n_per_sample_low = 1, n_per_sample_high = 10, prev_ct_skew_threshold = 0.12)




# uniform Ct distribution    ------

# results dataframes
perm.out <- data.frame(
        session                   = rep(seq_len(n_sessions),
                                        permutations * length(ct_thresholds)),
        permutation               = rep(seq_len(permutations),
                                        each = n_sessions * length(ct_thresholds)),
        hen_prevalence = NA_real_,
        hen_prevalence_normalized = NA_real_,
        hen_prevalence_normalized_smooth = NA_real_,
        ct_threshold              = rep(rep(ct_thresholds, each = n_sessions),
                                        permutations),
        distribution = NA
)

obs.out <- data.frame(
        session                   = rep(seq_len(n_sessions),
                                        length(ct_thresholds)),
        hen_prevalence = NA_real_,
        hen_prevalence_normalized = NA_real_,
        hen_prevalence_normalized_smooth = NA_real_,
        ct_threshold              = rep(ct_thresholds, each = n_sessions),
        prev_quantile             = NA_real_,
        distribution = NA
)


perm_all_out <- data.frame(permutation = rep(1:permutations,each = n_sessions),
                           session = rep(1:n_sessions, permutations),
                           ct_threshold = rep(rep(ct_thresholds, each = n_sessions),
                                              permutations),
                           prev = NA,
                           distribution = NA)


for(threshold in ct_thresholds) {
        
        
        
        for(i in 1:permutations) {
                
                perm_dat = sim_dat
                # calculate individual Ct values given a randomly selected number of positive bats out of the total number of bats over the sheet
                perm_dat$cur_n_pos = sapply(sim_dat$n_bats, FUN = function(x) sample(1:x, size = 1))
                
                # calculate average Ct value of each positive individual
                perm_dat$cur_ct_pos_avg = apply(perm_dat[,c("uniform_ct","cur_n_pos","n_bats")], 1, function(x) ct_pos_fun(ct_pooled = x[1], n_pos = x[2], n_total = x[3]))
                
                # permute Ct values
                perm_dat$perm_ct_avg = 0
                perm_dat$perm_ct_avg[which(perm_dat$uniform_ct>0)] = sample(perm_dat$cur_ct_pos_avg[which(perm_dat$uniform_ct>0)], replace = F)
                
                # calculate pooled Ct value, given the current average Ct value for a positive individual, and the current number of positive individuals
                perm_dat$perm_ct_pooled = apply(perm_dat[,c("perm_ct_avg","cur_n_pos","n_bats")], 1, function(x) ct_pooled_fun(ct_pos = x[1], n_pos = x[2], n_total = x[3]))
                
                # calculate prevalence per session, applying Ct threshold
                perm_sess_prev = perm_dat %>% 
                        group_by(session) %>% 
                        summarize(prev = sum(perm_ct_pooled < threshold & perm_ct_pooled > 0)/n()) %>% 
                        mutate(prev_normalized = prev/max(prev))
                
                perm.prev.smooth = gam(prev_normalized ~ ns(session, df = 6), family = "binomial", data = perm_sess_prev)
                perm_sess_prev$prev_normalized_smooth = predict(perm.prev.smooth,type = "response")
                
                perm.out[which(perm.out$ct_threshold == threshold & perm.out$permutation==i), "hen_prevalence"] = perm_sess_prev$prev
                perm.out[which(perm.out$ct_threshold == threshold & perm.out$permutation==i), "hen_prevalence_normalized"] = perm_sess_prev$prev_normalized
                perm.out[which(perm.out$ct_threshold == threshold & perm.out$permutation==i), "hen_prevalence_normalized_smooth"] = perm_sess_prev$prev_normalized_smooth
                
                perm_all_out$prev[which(perm_all_out$ct_threshold == threshold & perm_all_out$permutation == i)] = perm_sess_prev$prev
                
        }
        
        
        sim_sess_prev = sim_dat %>% 
                group_by(session) %>% 
                summarize(prev = sum(uniform_ct < threshold & uniform_ct > 0)/n())
        
        obs.out[which(obs.out$ct_threshold == threshold), "hen_prevalence"] = sim_sess_prev$prev
        obs.out[which(obs.out$ct_threshold == threshold), "hen_prevalence_normalized"] = sim_sess_prev$prev / max(sim_sess_prev$prev)
        
        obs.prev.smooth = gam(hen_prevalence_normalized ~ ns(session, df = 6), family = "binomial", data = obs.out[which(obs.out$ct_threshold == threshold),])
        obs.out$hen_prevalence_normalized_smooth[which(obs.out$ct_threshold==threshold)] = predict(obs.prev.smooth,type = "response")
        # plot(obs.out$hen_prevalence_normalized_smooth[which(obs.out$ct_threshold==threshold)], type = "l", ylim = c(0,1))
        # lines(obs.out$hen_prevalence_normalized[which(obs.out$ct_threshold==threshold)], col = "red")
        
        
        cur.obs.out.idx = which(obs.out$ct_threshold == threshold)
        for(k in 1:length(cur.obs.out.idx)){
                cur.prev.obs = obs.out$hen_prevalence[which(obs.out$ct_threshold==threshold)[k]]
                cur.session.obs = obs.out$session[which(obs.out$ct_threshold==threshold)[k]]
                cur.prev.perm = perm.out$hen_prevalence[which(perm.out$ct_threshold==threshold & perm.out$session==cur.session.obs)]
                
                obs.out$prev_quantile[cur.obs.out.idx[k]] = sum(cur.prev.perm < cur.prev.obs)/permutations
        }
        
        
        
        
        cat(paste0("\rthreshold ", threshold))
}



perm.out.uniform = perm.out
perm.out.uniform$distribution = "uniform"
obs.out.uniform = obs.out
obs.out.uniform$distribution = "uniform"
perm_all_out_uniform = perm_all_out
perm_all_out_uniform$distribution = "uniform"






set.seed(11)
# skewed Ct distribution    ------

# results dataframes
perm.out <- data.frame(
        session                   = rep(seq_len(n_sessions),
                                        permutations * length(ct_thresholds)),
        permutation               = rep(seq_len(permutations),
                                        each = n_sessions * length(ct_thresholds)),
        hen_prevalence = NA_real_,
        hen_prevalence_normalized = NA_real_,
        hen_prevalence_normalized_smooth = NA_real_,
        ct_threshold              = rep(rep(ct_thresholds, each = n_sessions),
                                        permutations),
        distribution = NA
)

obs.out <- data.frame(
        session                   = rep(seq_len(n_sessions),
                                        length(ct_thresholds)),
        hen_prevalence = NA_real_,
        hen_prevalence_normalized = NA_real_,
        hen_prevalence_normalized_smooth = NA_real_,
        ct_threshold              = rep(ct_thresholds, each = n_sessions),
        prev_quantile             = NA_real_,
        distribution = NA
)


perm_all_out <- data.frame(permutation = rep(1:permutations,each = n_sessions),
                           session = rep(1:n_sessions, permutations),
                           ct_threshold = rep(rep(ct_thresholds, each = n_sessions),
                                              permutations),
                           prev = NA,
                           distribution = NA)


for(threshold in ct_thresholds) {
        
        
        
        for(i in 1:permutations) {
                
                perm_dat = sim_dat
                # calculate individual Ct values given a randomly selected number of positive bats out of the total number of bats over the sheet
                perm_dat$cur_n_pos = sapply(sim_dat$n_bats, FUN = function(x) sample(1:x, size = 1))
                
                # calculate average Ct value of each positive individual
                perm_dat$cur_ct_pos_avg = apply(perm_dat[,c("skewed_ct","cur_n_pos","n_bats")], 1, function(x) ct_pos_fun(ct_pooled = x[1], n_pos = x[2], n_total = x[3]))
                
                # permute Ct values
                perm_dat$perm_ct_avg = 0
                perm_dat$perm_ct_avg[which(perm_dat$uniform_ct>0)] = sample(perm_dat$cur_ct_pos_avg[which(perm_dat$uniform_ct>0)], replace = F)
                
                # calculate pooled Ct value, given the current average Ct value for a positive individual, and the current number of positive individuals
                perm_dat$perm_ct_pooled = apply(perm_dat[,c("perm_ct_avg","cur_n_pos","n_bats")], 1, function(x) ct_pooled_fun(ct_pos = x[1], n_pos = x[2], n_total = x[3]))
                
                # calculate prevalence per session, applying Ct threshold
                perm_sess_prev = perm_dat %>% 
                        group_by(session) %>% 
                        summarize(prev = sum(perm_ct_pooled < threshold & perm_ct_pooled > 0)/n()) %>% 
                        mutate(prev_normalized = prev/max(prev))
                
                perm.prev.smooth = gam(prev_normalized ~ ns(session, df = 6), family = "binomial", data = perm_sess_prev)
                perm_sess_prev$prev_normalized_smooth = predict(perm.prev.smooth,type = "response")
                
                perm.out[which(perm.out$ct_threshold == threshold & perm.out$permutation==i), "hen_prevalence"] = perm_sess_prev$prev
                perm.out[which(perm.out$ct_threshold == threshold & perm.out$permutation==i), "hen_prevalence_normalized"] = perm_sess_prev$prev_normalized
                perm.out[which(perm.out$ct_threshold == threshold & perm.out$permutation==i), "hen_prevalence_normalized_smooth"] = perm_sess_prev$prev_normalized_smooth
                
                perm_all_out$prev[which(perm_all_out$ct_threshold == threshold & perm_all_out$permutation == i)] = perm_sess_prev$prev
                
        }
        
        
        sim_sess_prev = sim_dat %>% 
                group_by(session) %>% 
                summarize(prev = sum(skewed_ct < threshold & skewed_ct > 0)/n())
        
        obs.out[which(obs.out$ct_threshold == threshold), "hen_prevalence"] = sim_sess_prev$prev
        obs.out[which(obs.out$ct_threshold == threshold), "hen_prevalence_normalized"] = sim_sess_prev$prev / max(sim_sess_prev$prev)
        
        obs.prev.smooth = gam(hen_prevalence_normalized ~ ns(session, df = 6), family = "binomial", data = obs.out[which(obs.out$ct_threshold == threshold),])
        obs.out$hen_prevalence_normalized_smooth[which(obs.out$ct_threshold==threshold)] = predict(obs.prev.smooth,type = "response")
        # plot(obs.out$hen_prevalence_normalized_smooth[which(obs.out$ct_threshold==threshold)], type = "l", ylim = c(0,1))
        # lines(obs.out$hen_prevalence_normalized[which(obs.out$ct_threshold==threshold)], col = "red")
        
        
        cur.obs.out.idx = which(obs.out$ct_threshold == threshold)
        for(k in 1:length(cur.obs.out.idx)){
                cur.prev.obs = obs.out$hen_prevalence[which(obs.out$ct_threshold==threshold)[k]]
                cur.session.obs = obs.out$session[which(obs.out$ct_threshold==threshold)[k]]
                cur.prev.perm = perm.out$hen_prevalence[which(perm.out$ct_threshold==threshold & perm.out$session==cur.session.obs)]
                
                obs.out$prev_quantile[cur.obs.out.idx[k]] = sum(cur.prev.perm < cur.prev.obs)/permutations
        }
        
        
        
        
        cat(paste0("\rthreshold ", threshold))
}



perm.out.skewed = perm.out
perm.out.skewed$distribution = "skewed"
obs.out.skewed = obs.out
obs.out.skewed$distribution = "skewed"
perm_all_out_skewed = perm_all_out
perm_all_out_skewed$distribution = "skewed"






# plot output -----

plot1 = ggplot() +
        geom_jitter(data = perm_all_out_uniform, aes(x = session, y = prev), col = "grey30", alpha = 0.5, width = 0.2, height = 0) +
        geom_point(data = obs.out.uniform, aes(x = session, y = hen_prevalence), col = "orange", shape = 17, size = 4) +
        ylab("Prevalence") +
        xlab("Session") +
        scale_x_continuous(breaks = 1:10) +
        ggtitle("Uniform Distribution\n(Prevalence-Independent)") +
        theme_bw() +
        theme(axis.text.y = element_text(size=16), 
              axis.text.x = element_text(size=16), 
              axis.title.x = element_text(size=15), 
              axis.title.y = element_text(size=15),
              plot.title = element_text(size=15), 
              strip.text = element_text(size = 18), 
              legend.position = "none") +
        facet_wrap(~ct_threshold)




plot1 = ggplot() +
        geom_jitter(data = perm_all_out_skewed, aes(x = session, y = prev), col = "grey30", alpha = 0.5, width = 0.2, height = 0) +
        geom_point(data = obs.out.skewed, aes(x = session, y = hen_prevalence), col = "orange", shape = 17, size = 4) +
        ylab("Prevalence") +
        xlab("Session") +
        scale_x_continuous(breaks = 1:10) +
        ggtitle("Skewed Distribution\n(Prevalence-Dependent)") +
        theme_bw() +
        theme(axis.text.y = element_text(size=16), 
              axis.text.x = element_text(size=16), 
              axis.title.x = element_text(size=15), 
              axis.title.y = element_text(size=15),
              plot.title = element_text(size=15), 
              strip.text = element_text(size = 18), 
              legend.position = "none") +
        facet_wrap(~ct_threshold)






# quantile plots -------  

# load functions   


## Plot functions:
plot.left.not.bottom = function() {
        
        cur.colors = colors[sort(unique(obs.out$quantile.cat[which(obs.out$ct_threshold==cur.ct)]))]
        
        ggplot() + 
                #geom_rect(data = obs.out[which(obs.out$ct_threshold==cur.ct),], aes(xmin = xmin, xmax = xmax, ymin = -0.5, ymax = Inf), fill = obs.out[which(obs.out$ct_threshold==cur.ct),"fill.col"]) +
                geom_rect(data = obs.out[which(obs.out$ct_threshold==cur.ct),], aes(xmin = xmin, xmax = xmax, ymin = -0.5, ymax = Inf, fill = factor(quantile.cat)), alpha = 0.8) +
                scale_fill_manual(values = cur.colors, name = "Prevalence\nQuantile", labels = c("0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%")) +
                geom_line(data=perm.out[which(perm.out$ct_threshold==cur.ct),],
                          aes(y = hen_prevalence_normalized, x = session, group = permutation),
                          colour = color.perm.prev, size = 0.5, alpha = 0.2) +
                geom_line(data=obs.out[which(obs.out$ct_threshold==cur.ct),],
                          aes(y = hen_prevalence_normalized, x = session),
                          colour = color.obs.prev, size = 2.5) +
                # theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
                scale_x_continuous(limits = c(0.5,10.5)) +
                labs(x = "Session",
                     y = "Normalized prevalence") +
                theme_bw() +
                background_grid("none") +
                ggtitle(paste0("Ct threshold = ",cur.ct)) + 
                coord_cartesian(ylim=c(-0.1,1.0), xlim = c(1,10)) + 
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
                          aes(y = hen_prevalence_normalized, x = session, group = permutation),
                          colour = color.perm.prev, size = 0.5, alpha = 0.2) +
                geom_line(data=obs.out[which(obs.out$ct_threshold==cur.ct),],
                          aes(y = hen_prevalence_normalized, x = session),
                          colour = color.obs.prev, size = 2.5) +
                # theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
                scale_x_continuous(limits = c(0.5,10.5)) +
                labs(x = "Session",
                     y = "Normalized prevalence") +
                theme_bw() +
                background_grid("none") +
                ggtitle(paste0("Ct threshold = ",cur.ct)) + 
                coord_cartesian(ylim=c(-0.1,1.0), xlim = c(1,10)) + 
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
                          aes(y = hen_prevalence_normalized, x = session, group = permutation),
                          colour = color.perm.prev, size = 0.5, alpha = 0.2) +
                geom_line(data=obs.out[which(obs.out$ct_threshold==cur.ct),],
                          aes(y = hen_prevalence_normalized, x = session),
                          colour = color.obs.prev, size = 2.5) +
                # theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
                scale_x_continuous(limits = c(0.5,10.5)) +
                labs(x = "Session",
                     y = "Normalized prevalence") +
                theme_bw() +
                background_grid("none") +
                ggtitle(paste0("Ct threshold = ",cur.ct)) + 
                coord_cartesian(ylim=c(-0.1,1.0), xlim = c(1,10)) + 
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
                          aes(y = hen_prevalence_normalized, x = session, group = permutation),
                          colour = color.perm.prev, size = 0.5, alpha = 0.2) +
                geom_line(data=obs.out[which(obs.out$ct_threshold==cur.ct),],
                          aes(y = hen_prevalence_normalized, x = session),
                          colour = color.obs.prev, size = 2.5) +
                # theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
                scale_x_continuous(limits = c(0.5,10.5)) +
                labs(x = "Session",
                     y = "Normalized prevalence") +
                theme_bw() +
                background_grid("none") +
                ggtitle(paste0("Ct threshold = ",cur.ct)) + 
                coord_cartesian(ylim=c(-0.1,1.0), xlim = c(1,10)) + 
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
                          aes(y = hen_prevalence_normalized, x = session, group = permutation),
                          colour = color.perm.prev, size = 0.5, alpha = 0.2) +
                geom_line(data=obs.out[which(obs.out$ct_threshold==cur.ct),],
                          aes(y = hen_prevalence_normalized, x = session),
                          colour = color.obs.prev, size = 2.5) +
                # theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
                scale_x_continuous(limits = c(0.5,10.5)) +
                labs(x = "Session",
                     y = "Normalized prevalence") +
                theme_bw() +
                background_grid("none") +
                ggtitle(paste0("Ct threshold = ",cur.ct)) + 
                coord_cartesian(ylim=c(-0.1,1.0), xlim = c(1,10)) + 
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
        
        obs.out.temp = obs.out %>% filter(ct_threshold == 30)
        obs.out.temp$quantile.cat = 1:10
        
        
        cur.plot = ggplot() + 
                #geom_rect(data = obs.out[which(obs.out$ct_threshold==cur.ct),], aes(xmin = xmin, xmax = xmax, ymin = -0.5, ymax = Inf), fill = obs.out[which(obs.out$ct_threshold==cur.ct),"fill.col"]) +
                geom_rect(data = obs.out.temp, aes(xmin = xmin, xmax = xmax, ymin = -0.5, ymax = Inf, fill = factor(quantile.cat)), alpha = 0.8) +
                scale_fill_manual(values = colors, name = "Prevalence\nQuantile", labels = c("0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%")) +
                geom_line(data=perm.out[which(perm.out$ct_threshold==cur.ct),],
                          aes(y = hen_prevalence_normalized, x = session, group = permutation),
                          colour = color.perm.prev, size = 0.5, alpha = 0.2) +
                geom_line(data=obs.out[which(obs.out$ct_threshold==cur.ct),],
                          aes(y = hen_prevalence_normalized, x = session),
                          colour = color.obs.prev, size = 2.5) +
                # theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
                scale_x_continuous(limits = c(0.5,10.5)) +
                labs(x = "Session",
                     y = "Normalized prevalence") +
                theme_bw() +
                background_grid("none") +
                ggtitle(paste0("Ct threshold = ",cur.ct)) + 
                coord_cartesian(ylim=c(-0.1,1.0), xlim = c(1,10)) + 
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
colors = rev(brewer.pal(n = 10,"RdYlBu"))



# uniform distribution   ------

obs.out = obs.out.uniform

# add xmin and xmax session times for quantile plotting


obs.out$xmin = rep(0.5:9.5, length(ct_thresholds))
obs.out$xmax = rep(1.5:10.5, length(ct_thresholds))

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
plot.out.uniform = plot.28 + plot.29 + plot.30 + plot.31 + plot.32 + plot.33 + plot.34 + plot.35 + plot.36 + plot.37 + plot.38 + plot.39 + plot.legend +
        plot_layout(design = layout)

plot.out.uniform





# skewed distribution   ------

obs.out = obs.out.skewed

# add xmin and xmax session times for quantile plotting


obs.out$xmin = rep(0.5:9.5, length(ct_thresholds))
obs.out$xmax = rep(1.5:10.5, length(ct_thresholds))

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
plot.out.skewed = plot.28 + plot.29 + plot.30 + plot.31 + plot.32 + plot.33 + plot.34 + plot.35 + plot.36 + plot.37 + plot.38 + plot.39 + plot.legend +
        plot_layout(design = layout)

plot.out.skewed




##------------------------------------------------------------------------------------------##
##---------------------------------------Export figures-------------------------------------##
##------------------------------------------------------------------------------------------##

##-------------------------------------- Figure S2A ----------------------------------------##

ggsave(plot = plot.out.uniform,filename = "Permutation analysis simulation test uniform 20250615.png",width = 30, height = 24, dpi=300)


##-------------------------------------- Figure S2A ----------------------------------------##

ggsave(plot = plot.out.skewed,filename = "Permutation analysis simulation test skewed 20250615.png",width = 30, height = 24, dpi=300)
