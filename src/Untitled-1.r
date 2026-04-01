required_packages = c("DEoptim", "Rcpp", "parallel", "RcppParallel","stats4","pracma","tidyverse","ggpubr", "fs", "here", "lmerTest","ARTool")

install_if_missing <- function(p) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}

# Install missing packages
invisible(sapply(required_packages, install_if_missing))

# Load all packages
invisible(sapply(required_packages, library, character.only = TRUE))

RcppParallel::setThreadOptions(numThreads = 3) #this is critical for mclapply


scriptFolder = path(here()) / "src" / "stDDM"
dataFolder = path(here()) / "data"

# Code for simulating data
sourceCpp(scriptFolder / "2DDM_r_cpp_M1.cpp")

# Load data
Data<-read.csv(file=paste0(dataFolder, "/deidentified_ChoiceData.csv"))

# Simulate each condition separately
DDataSim <- NULL
for (gg in c("BN","HC")) {
  for (cc in c("Neutral","Negative")) {
  Data_partial <- Data %>%
    distinct() %>%
    filter(cond == cc,
           Dx == gg) %>%
    mutate(foodType = ifelse(foodType == "High Fat",2,1)
           )
  
  idxR = ordered(Data_partial$idx)
  idxP = as.numeric(ordered(Data_partial$idx)) #makes a sequentially numbered subj index
  
  idxDF = data.frame("idxR" = idxR,
                     "idxP" = idxP)
  
  Data_partial$idxP = idxDF$idxP[idxDF$idxR == Data_partial$idx]
  DDataSim_tmp<- data.frame(hd =Data_partial$hd, 
                        td = Data_partial$td, 
                        idx = Data_partial$idx, 
                        idxP = Data_partial$idxP,
                        cond = Data_partial$cond, 
                        foodType = Data_partial$foodType,
                        Dx = Data_partial$Dx)
  DDataSim = rbind(DDataSim,DDataSim_tmp)

  }
}

load(file = file.path(path(here()) / "results" / "stDDM" , 
                                                 paste0("Feb2026_STAN_params_Params_FIT_M3.RData")))
params_STAN <- params.df %>%
  mutate(params = recode(param_type,
                         "b1_p" = "wt",
                         "b2_p" = "wh",
                         "alpha_p" = "boundary",
                         "theta_p" = "nDT",
                         "time_p" = "tHin",
                         "bias_p" = "bias")) %>%
  rename(vals = mean) %>% select(!c(`2.5%`,`97.5%`,n_eff,Rhat,param_type))
# Model 3 (Time, b1, b2 depends on food type)
DataSim_m3<-NULL

for (gg in c("BN","HC")) {
  subj <- unique(params_STAN$idxP[params_STAN$Dx == gg])
  Dx_df = params_STAN %>% filter(Dx == gg)
  for (cc in c("Negative","Neutral")) {
    Cond_df = Dx_df %>% filter(cond == cc)
    #load(paste0(resFolder,"/Feb2026_params_HtSSM_FIT_M3bias_Dx-",gg,"_","Cond-",cc,".RData"))
    #chain=rbind(results$mcmc[[1]], results$mcmc[[2]], results$mcmc[[3]])
    for (s in subj){
      #wt_lf =  mean(chain[,c( paste( c("b1.p[",toString(s),",1]"), collapse = ""))])
      #wh_lf = #mean(chain[,c( paste( c("b2.p[",toString(s),",1]"), collapse = ""))])
      #wt_hf = #mean(chain[,c( paste( c("b1.p[",toString(s),",2]"), collapse = ""))])
      #wh_hf = #mean(chain[,c( paste( c("b2.p[",toString(s),",2]"), collapse = ""))])
      #nDT = #mean(chain[,c( paste( c("theta.p[",toString(s),"]"), collapse = ""))])

      #bias_lf = mean(chain[,c( paste( c("bias[",toString(s),",1]"), collapse = ""))])
      #bias_hf = mean(chain[,c( paste( c("bias[",toString(s),",2]"), collapse = ""))])

      #bound = mean(chain[,c( paste( c("alpha.p[",toString(s),"]"), collapse = ""))])

      #time_lf = mean(chain[,c( paste( c("time.p[",toString(s),",1]"), collapse = ""))])
      #time_hf = mean(chain[,c( paste( c("time.p[",toString(s),",2]"), collapse = ""))])

      wt_lf = Cond_df %>% filter(params == "wt", foodType == "Low-Fat", idxP == s) %>% pull(vals)
      wh_lf = Cond_df %>% filter(params == "wh", foodType == "Low-Fat", idxP == s) %>% pull(vals)
      wt_hf = Cond_df %>% filter(params == "wt", foodType == "High-Fat", idxP == s) %>% pull(vals)
      wh_hf = Cond_df %>% filter(params == "wh", foodType == "High-Fat", idxP == s) %>% pull(vals)
      nDT = Cond_df %>% filter(params == "nDT", idxP == s) %>% pull(vals)
      bias = Cond_df %>% filter(params == "bias", idxP == s) %>% pull(vals)
      bound = Cond_df %>% filter(params == "boundary", idxP == s) %>% pull(vals)
      time_lf = Cond_df %>% filter(params == "tHin", foodType == "Low-Fat", idxP == s) %>% pull(vals)
      time_hf = Cond_df %>% filter(params == "tHin", foodType == "High-Fat", idxP == s) %>% pull(vals)

      td_lf <- (DDataSim$td[DDataSim$idxP==s & DDataSim$foodType==1 & DDataSim$Dx == gg])
      td_hf <- (DDataSim$td[DDataSim$idxP==s & DDataSim$foodType==2 & DDataSim$Dx == gg])
      hd_lf <- (DDataSim$hd[DDataSim$idxP==s & DDataSim$foodType==1 & DDataSim$Dx == gg])
      hd_hf <- (DDataSim$hd[DDataSim$idxP==s & DDataSim$foodType==2 & DDataSim$Dx == gg])
      
      # pre-allocate matrix with combinations of vd and hd
      trial_mat_hf = expand.grid(unique(td_hf),unique(hd_hf))
      trial_mat_lf = expand.grid(unique(td_lf),unique(hd_lf))
      
      colnames(trial_mat_hf) = c("td","hd")
      colnames(trial_mat_lf) = c("td","hd")
      
      trial_mat_lf$foodType = 1
      trial_mat_hf$foodType = 2
      
      trial_mat_lf$count = 0
      trial_mat_hf$count = 0
      
      for (i in 1:length(td_lf)) {
        trial_mat_lf$count[trial_mat_lf$td == td_lf[i] & trial_mat_lf$hd == hd_lf[i]] = trial_mat_lf$count[trial_mat_lf$td == td_lf[i] & trial_mat_lf$hd == hd_lf[i]] + 1
      }
      if(sum(trial_mat_lf$count) != length(td_lf) ) {
        print(paste0("Subject ", s, "problem with lf trials"))
      }
      
      
      for (i in 1:length(td_hf)) {
        trial_mat_hf$count[trial_mat_hf$td == td_hf[i] & trial_mat_hf$hd == hd_hf[i]] = trial_mat_hf$count[trial_mat_hf$td == td_hf[i] & trial_mat_hf$hd == hd_hf[i]] + 1
      }
      if(sum(trial_mat_hf$count) != length(td_hf) ) {
        print(paste0("Subject ", s, "problem with negative trials"))
      }
      
      # Only keep trials WITH observations
      trial_mat_hf <- trial_mat_hf[trial_mat_hf$count > 0,]
      trial_mat_lf <- trial_mat_lf[trial_mat_lf$count > 0,]
      
      DataSim_Mat<-NULL
      for (i in 1:nrow(trial_mat_hf)) {
        rts = ddm2_parallel( wt_hf,wh_hf,bound,nDT,time_hf,bias,trial_mat_hf$td[i],trial_mat_hf$hd[i],1,trial_mat_hf$count[i] )
        for (j in 1:length(rts)) {
         
          DataSim_Mat_Temp = NULL
          DataSim_Mat_Temp$idx = unique(DDataSim$idx[DDataSim$idxP==s & DDataSim$Dx == gg])
          DataSim_Mat_Temp$idxP = s
          DataSim_Mat_Temp$Dx = gg
          DataSim_Mat_Temp$cond = cc
          DataSim_Mat_Temp$foodType = trial_mat_hf$foodType[i]
          DataSim_Mat_Temp$td = trial_mat_hf$td[i]
          DataSim_Mat_Temp$hd = trial_mat_hf$hd[i]
          DataSim_Mat_Temp$wt = wt_hf
          DataSim_Mat_Temp$wh = wh_hf
          DataSim_Mat_Temp$nDT = nDT
          DataSim_Mat_Temp$bias = bias
          DataSim_Mat_Temp$bound = bound
          DataSim_Mat_Temp$time = time_hf
          
          if (rts[j]<0){
            DataSim_Mat_Temp$choice <- 0
            DataSim_Mat_Temp$rt <- -rts[j]
            DataSim_Mat = rbind(DataSim_Mat, DataSim_Mat_Temp)
          } else {
            DataSim_Mat_Temp$choice <- 1
            DataSim_Mat_Temp$rt <- rts[j]
            DataSim_Mat = rbind(DataSim_Mat, DataSim_Mat_Temp)
          }
        }
      }
      
      for (i in 1:nrow(trial_mat_lf)) {
        rts = ddm2_parallel(wt_lf,wh_lf,bound,nDT,time_lf,bias,trial_mat_lf$td[i],trial_mat_lf$hd[i],1,trial_mat_lf$count[i])
        for (j in 1:length(rts)) {
          
          DataSim_Mat_Temp = NULL
          DataSim_Mat_Temp$idx = unique(DDataSim$idx[DDataSim$idxP==s & DDataSim$Dx == gg])
          DataSim_Mat_Temp$idxP = s
          DataSim_Mat_Temp$Dx = gg
          DataSim_Mat_Temp$cond = cc
          DataSim_Mat_Temp$foodType = trial_mat_lf$foodType[i]
          DataSim_Mat_Temp$td = trial_mat_lf$td[i]
          DataSim_Mat_Temp$hd = trial_mat_lf$hd[i]
          DataSim_Mat_Temp$wt = wt_lf
          DataSim_Mat_Temp$wh = wh_lf
          DataSim_Mat_Temp$nDT = nDT
          DataSim_Mat_Temp$bias = bias
          DataSim_Mat_Temp$bound = bound
          DataSim_Mat_Temp$time = time_lf
          
          if (rts[j]<0){
            DataSim_Mat_Temp$choice <- 0
            DataSim_Mat_Temp$rt <- -rts[j]
            DataSim_Mat = rbind(DataSim_Mat, DataSim_Mat_Temp)
          } else {
            DataSim_Mat_Temp$choice <- 1
            DataSim_Mat_Temp$rt <- rts[j]
            DataSim_Mat = rbind(DataSim_Mat, DataSim_Mat_Temp)
          }
        }
      }
      
      DataSim_T<-DataSim_Mat %>% as.data.frame()
      DataSim_m3<-rbind(DataSim_m3,DataSim_T)
      
    }
  }
}
# Convert from list
DataSim_m3 <- DataSim_m3 %>% 
  mutate(idx = unlist(idx),foodType = unlist(foodType),
         idxP = unlist(idxP), Dx = unlist(Dx), 
         td = unlist(td), hd = unlist(hd), wt = unlist(wt), wh = unlist(wh),
         nDT = unlist(nDT),bias = unlist(bias), 
         bound = unlist(bound), cond = unlist(cond),
         time = unlist(time), choice = unlist(choice),
         rt = unlist(rt),
         foodType = ifelse(foodType == 2, "High Fat","Low Fat"))
# Save
save(DataSim_m3, file=dataFolder / "DataSim_M3_STAN.RData")

DataSim_m3 = DataSim_m3 %>%
  mutate(cond = factor(cond, levels = c("Neutral","Negative")),
         foodType = factor(foodType, levels = c("Low Fat","High Fat")),
         Dx = factor(Dx, levels = c("HC","BN")),
         m = "M3ndt"
  )

emp_summary = Data %>%
  group_by(idx,Dx,cond,foodType) %>%
  summarise(mch = mean(choice),
            mrt = median(rt)) %>%
  mutate(type = "emp")
emp_summary_quantiles = Data %>%
  group_by(idx,Dx,cond,foodType,choice) %>%
  summarise(mQ1 = quantile(rt, probs = .1),
            mQ3 = quantile(rt, probs=.3),
            mQ5 = quantile(rt, probs=.5),
            mQ7 = quantile(rt, probs=.7),
            mQ9 = quantile(rt, probs = .9)) %>%
  mutate(type = "emp")
emp_summary = merge(emp_summary,emp_summary_quantiles) %>%
  mutate(mch = ifelse(choice == 0, 1 - mch, mch))

sim_summary_m3 = DataSim_m3 %>%
  group_by(idx,Dx,cond,foodType) %>%
  summarise(mch = mean(choice),
            mrt = median(rt)) %>%
  mutate(type = "sim")

sim_summary_m3_quantiles = DataSim_m3 %>%
  group_by(idx,Dx,cond,foodType,choice) %>%
  summarise(mQ1 = quantile(rt, probs = .1),
            mQ3 = quantile(rt, probs=.3),
            mQ5 = quantile(rt, probs=.5),
            mQ7 = quantile(rt, probs=.7),
            mQ9 = quantile(rt, probs = .9)) %>%
  mutate(type = "sim")
sim_summary_m3 = merge(sim_summary_m3,sim_summary_m3_quantiles) %>%
  mutate(mch = ifelse(choice == 0, 1 - mch, mch))
dat_summary_m3 = rbind(sim_summary_m3, emp_summary) %>%
  mutate(foodType = factor(foodType, levels = c("Low Fat","High Fat")),
         cond = factor(cond, levels = c("Neutral","Negative")),
         Dx = factor(Dx, levels = c("HC","BN","BED")),
         source = factor(type,levels=c("emp","sim"))
  )

dat_summary_m3 %>%
  mutate(source = ifelse(source == "emp", "Empirical","Simulation")) %>%
  group_by(Dx,cond, foodType, source, choice) %>%
  summarise(c = mean(mch),
            Q1 = mean(mQ1),
            Q3 = mean(mQ3),
            Q5 = mean(mQ5),
            Q7 = mean(mQ7),
            Q9 = mean(mQ9)
  ) %>% 
  ggplot()+
  theme_pubr(base_size=18)+
  facet_grid(cond~Dx)+
  geom_point(aes(x = c, y = Q1, group = interaction(cond,foodType), shape = source), color = "magenta")+
  geom_line(data = . %>% filter(source == "Simulation"), 
            aes(x = c, y = Q1), color = "magenta")+
  geom_point(aes(x = c, y = Q3, group = interaction(cond,foodType), shape = source), color = "blue")+
  geom_line(data = . %>% filter(source == "Simulation"), 
            aes(x = c, y = Q3), color = "blue")+
  geom_point(aes(x = c, y = Q5, group = interaction(cond,foodType), shape = source), color = "forestgreen")+
  geom_line(data = . %>% filter(source == "Simulation"), 
            aes(x = c, y = Q5), color = "forestgreen")+
  geom_point(aes(x = c, y = Q7, group = interaction(cond,foodType), shape = source), color = "gold")+
  geom_line(data = . %>% filter(source == "Simulation"), 
            aes(x = c, y = Q7), color = "gold")+
  geom_point(aes(x = c, y = Q9, group = interaction(cond,foodType), shape = source), color = "darkred")+
  geom_line(data = . %>% filter(source == "Simulation"), 
            aes(x = c, y = Q9), color = "darkred")+
  labs(x = "Choice Frequency", y = "Response Time [sec]", title = "Model 3 (STAN)", shape = "Data Type") +
  scale_shape_manual(values = c(0,4))
