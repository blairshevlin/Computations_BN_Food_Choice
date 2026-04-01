libs = c("rstan","tidyverse","fs","here","loo")
suppressPackageStartupMessages(lapply(libs, require, character.only = TRUE))

# Stan model options
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

resFolder = path(here()) / "results" / "stDDM" / "estimation"

# Check for convergence
for (m in c("M0","M1","M3")){
  for (gg in c("BN","HC")) {
    for (cc in c("Neutral","Negative")) {
      print(paste0("Checking R-Hat for Model ",m," for Group ", gg, " in ",cc," Condition"))
      
      load(file.path(resFolder, 
                     paste0("Feb2026_params_STAN_FIT_", m, 
                            "singleb_Dx-", gg, "_", "Cond-", cc, ".RData")))
      
      # extract R-hat column
      rhat_values <- summary(results)$summary[, "Rhat"]
      
      # Check if any R-hat > 1.1 (problematic)
      print(any(rhat_values > 1.1, na.rm = TRUE))
      
      # See which parameters have high R-hat
      print(rhat_values[rhat_values > 1.1])
    }
  }
}

fit.df = NULL
loo.df = list()
for (m in c("M0","M1","M3")){
  for (gg in c("BN","HC")) {
    for (cc in c("Neutral","Negative")) {
      load(file.path(resFolder, 
                     paste0("Feb2026_params_STAN_FIT_", m, 
                            "singleb_Dx-", gg, "_", "Cond-", cc, ".RData")))
      tmp.fit = data.frame(model = m,
                           Dx = gg,
                           cond = cc,
                           loo = looic$estimates[3],
                           waic = waic$estimates[3])
      fit.df = rbind(fit.df,tmp.fit)
      loo_name = paste0(m,"-",gg,"-",cc)
      loo.df[[loo_name]] = looic
      
    }
  }
}

# For loo, top model is better

# Compare BN Negative
loo_compare(loo.df[grepl(paste0("-","BN-Negative"), names(loo.df))])

# Compare BN Neutral
loo_compare(loo.df[grepl(paste0("-","BN-Neutral"), names(loo.df))])

# Compare HC Negative
loo_compare(loo.df[grepl(paste0("-","HC-Negative"), names(loo.df))])

# Compare HC Neutral
loo_compare(loo.df[grepl(paste0("-","HC-Neutral"), names(loo.df))])

# M1 wins for each

# Extract parameters
params.df = NULL
for (gg in c("BN","HC")) {
  for (cc in c("Neutral","Negative")) {
    load(file.path(resFolder, 
                   paste0("Feb2026_params_STAN_FIT_M3", 
                          "singleb_Dx-", gg, "_", "Cond-", cc, ".RData")))
    # Get summary for all subject-level parameters
    subject_summary <- summary(results, 
                               pars = c("theta_p", "alpha_p", "bias_p",  #"sv_p",
                                        "time_p", "b1_p", "b2_p"))$summary
    
    # Convert to data frame
    subject_summary_df <- as.data.frame(subject_summary)
    subject_summary_df$parameter <- rownames(subject_summary_df)
    # Idx
    subject_summary_df$idx <- as.numeric(str_extract(subject_summary_df$parameter, "(?<=\\[)\\d+"))
    subject_summary_df$param_type <- str_extract(subject_summary_df$parameter, "^[^\\[]+")
    subject_summary_df$foodType <- as.numeric(str_extract(subject_summary_df$parameter, "(?<=,)\\d+(?=\\])"))
    subject_summary_df  = subject_summary_df %>% 
      mutate(idxP = idx,
             cond = cc,
             Dx = gg,
             foodType = case_when(foodType == 1 ~ "Low-Fat", 
                                  foodType == 2 ~ "High-Fat", TRUE ~ NA)) %>%
      select(c(idxP,Dx,cond,param_type, foodType, mean, `2.5%`, `97.5%`,n_eff, Rhat))
    rownames(subject_summary_df) = NULL
    params.df = rbind(params.df,subject_summary_df)
  }
}

save(loo.df, fit.df, params.df, file = file.path(path(here()) / "results" / "stDDM" , 
                                                 paste0("Feb2026_STAN_params_Params_FIT_M3.RData")))




# Compare to JAGS
# Load subject-level parameters
df.fit.full <- NULL

for (cc in c("Neutral","Negative")) {
  for (gg in c("BN","HC")) {
    if (gg == "BN" & cc == "Negative"){
      # Load in file that required re-estimation to ensure convergence
      load(file.path(path(here()) / 'results',paste("/stDDM/estimation/params_HtSSM_FIT_M3_Dx-",gg,"_","Cond-",cc,"_rawRatings_converted_rest.RData",sep="")) )
    } else {
        load(file.path(path(here()) / 'results',paste("/stDDM/estimation/params_HtSSM_FIT_M3_Dx-",gg,"_","Cond-",cc,"_rawRatings_converted.RData",sep="")) )
    }
    idxP = unique(Data_partial$idxP)
    
    chain=rbind(results$mcmc[[1]], results$mcmc[[2]], results$mcmc[[3]])
    for (s in idxP){
      df.fit <- NULL
      df.fit$wt_lf = median(chain[,c( paste( c("b1.p[",toString(s),",1]"), collapse = ""))])
      df.fit$wh_lf = median(chain[,c( paste( c("b2.p[",toString(s),",1]"), collapse = ""))])
      df.fit$wt_hf = median(chain[,c( paste( c("b1.p[",toString(s),",2]"), collapse = ""))])
      df.fit$wh_hf = median(chain[,c( paste( c("b2.p[",toString(s),",2]"), collapse = ""))])
      df.fit$boundary = median(chain[,c( paste( c("alpha.p[",toString(s),"]"), collapse = ""))])
      df.fit$nDT= median(chain[,c( paste( c("theta.p[",toString(s),"]"), collapse = ""))])

      df.fit$tHin_lf = median(chain[,c( paste( c("time.p[",toString(s),",1]"), collapse = ""))])
      df.fit$tHin_hf = median(chain[,c( paste( c("time.p[",toString(s),",2]"), collapse = ""))])
      df.fit$bias = median(chain[,c( paste( c("bias[",toString(s),"]"), collapse = ""))])
      
      df.fit <-df.fit %>% as.data.frame() %>% pivot_longer(cols = c(wt_lf,wt_hf,wh_lf,wh_hf,boundary,
                                                                    nDT,
                                                                    tHin_lf,tHin_hf,bias),
                                                           names_to = "params",
                                                           values_to = "vals") %>%
        mutate(cond = cc,
               Dx = gg,
               idxP = s,
               idx = Data_partial$idx[Data_partial$idxP==s][1])
      
      df.fit.full <- rbind(df.fit.full,df.fit)
    }
  }
}

# Split params up by food-type
params_JAGS <- df.fit.full %>%
  mutate(foodType = ifelse(grepl("_lf",params),"Low-Fat",
                           ifelse(grepl("_hf",params),"High-Fat", "NA")),
         params = recode(params,
                         "wt_lf" = "wt",
                         "wt_hf" = "wt",
                         "wh_lf" = "wh",
                         "wh_hf" = "wh",
                         "tHin_lf" = "tHin",
                         "tHin_hf" = "tHin")
  ) %>% select(!c(idx))%>% rename(jags = vals)

params_STAN <- params.df %>%
  mutate(params = recode(param_type,
                         "b1_p" = "wt",
                         "b2_p" = "wh",
                         "alpha_p" = "boundary",
                         "theta_p" = "nDT",
                         "time_p" = "tHin",
                         "bias_p" = "bias")) %>%
  rename(stan = mean) %>% select(!c(`2.5%`,`97.5%`,n_eff,Rhat,param_type))

params_JAGS %>% head()
params_STAN %>% head()

params = merge(params_JAGS, params_STAN)

library(ggpubr)
ggplot(params, aes(x = jags, y = stan, group = interaction(cond,foodType), 
color = interaction(cond,foodType), fill = interaction(cond,foodType))) + 
facet_wrap(~params, scales = "free") +
theme_classic() +
geom_point() +
stat_cor(method = "pearson") 

params_STAN_v = params_STAN %>% rename(vals = stan) %>% mutate(method = "STAN")
params_JAGS_v = params_JAGS %>% rename(vals = jags) %>% mutate(method = "JAGS")

params_all = rbind(params_JAGS_v, params_STAN_v)

params_all %>% ggplot(aes(x = vals, fill = method, color = method)) + 
  facet_wrap(~params, scales = "free") +
  theme_classic() +
  geom_density(alpha = 0.5)

params_all %>% filter(!is.na(foodType), params %in% c("wt","wh","tHin")) %>%
  ggplot(aes(x = Dx, y = vals, fill = interaction(cond,foodType), color = interaction(cond,foodType))) + 
  facet_wrap(method ~ params, scales = "free") +
  theme_classic() +
  stat_summary(position = position_dodge2(width=0.5))
