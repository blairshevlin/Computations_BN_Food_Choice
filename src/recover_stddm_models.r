# recover_stddm_models.R - Code for recovering winning model from "Negative affect influences the computations underlying food choice in bulimia nervosa" 
#
# Copyright (C) 2024 Blair Shevlin, <blairshevlin@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Record of Revisions
#
# Date            Programmers                         Descriptions of Change
# ====         ================                       ======================
# 06/18/24      Blair Shevlin                         Wrote code for original manuscript

required_packages = c("DEoptim", "Rcpp", "plyr", "parallel", "RcppParallel","stats4","pracma","runjags","tidyverse","ggpubr", "fs", "here")

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

base_path = here()

dataFolder = path(base_path) / "data"
scriptFolder = path(base_path) / "src" / "stDDM"
resFolder = path(base_path) / "results" / "stDDM" / "estimation"
recFolder = path(base_path) / "results" / "stDDM"  / "recovery"
figPath <- path(base_path) / 'results' / 'figures'

# Simulation code
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
### Simulations ###

# parameters to recover with the fitting

DataSim<-NULL

for (gg in c("BN","HC")) {
  subj <- unique(DDataSim$idxP[DDataSim$Dx == gg])
  for (cc in c("Negative","Neutral")) {
    load(paste0(resFolder,"/params_HtSSM_FIT_M3_Dx-",gg,"_","Cond-",cc,"_rawRatings_converted.RData"))
    chain=rbind(results$mcmc[[1]], results$mcmc[[2]], results$mcmc[[3]])
    for (s in subj){
      wt_lf = mean(chain[,c( paste( c("b1.p[",toString(s),",1]"), collapse = ""))])
      wh_lf = mean(chain[,c( paste( c("b2.p[",toString(s),",1]"), collapse = ""))])
      
      wt_hf = mean(chain[,c( paste( c("b1.p[",toString(s),",2]"), collapse = ""))])
      wh_hf = mean(chain[,c( paste( c("b2.p[",toString(s),",2]"), collapse = ""))])
      
      nDT = mean(chain[,c( paste( c("theta.p[",toString(s),"]"), collapse = ""))])

      bias = mean(chain[,c( paste( c("bias[",toString(s),"]"), collapse = ""))])
      
      bound = mean(chain[,c( paste( c("alpha.p[",toString(s),"]"), collapse = ""))])
      time_lf = mean(chain[,c( paste( c("time.p[",toString(s),",1]"), collapse = ""))])
      time_hf = mean(chain[,c( paste( c("time.p[",toString(s),",2]"), collapse = ""))])
      
      td_lf <- (DDataSim$td[DDataSim$idxP==s & DDataSim$foodType==2 & DDataSim$Dx == gg])
      td_hf <- (DDataSim$td[DDataSim$idxP==s & DDataSim$foodType==1 & DDataSim$Dx == gg])
      hd_lf <- (DDataSim$hd[DDataSim$idxP==s & DDataSim$foodType==2 & DDataSim$Dx == gg])
      hd_hf <- (DDataSim$hd[DDataSim$idxP==s & DDataSim$foodType==1 & DDataSim$Dx == gg])
      
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
          good = 0
          # Make sure RT is within good boundaries, otherwise re-run
          while (good == 0) {
            if (abs(rts[j]) > 0.2 & abs(rts[j]) <= 7) {
              good = 1
            }
            rts[j] = ddm2_parallel(wt_lf,wh_lf,bound,nDT,time_hf,bias,trial_mat_hf$td[i],trial_mat_hf$hd[i],1,1)
          }
          
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
          good = 0
          # Make sure RT is within good boundaries, otherwise re-run
          while (good == 0) {
            if (abs(rts[j]) > 0.2 & abs(rts[j]) <= 7) {
              good = 1
            }
            rts[j] = ddm2_parallel(wt_lf,wh_lf,bound,nDT,time_lf,bias,trial_mat_lf$td[i],trial_mat_lf$hd[i],1,1)
          }
          
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
      DataSim<-rbind(DataSim,DataSim_T)
        
    }
  }
}
# Convert from list
DataSim <- DataSim %>% 
  mutate(idx = unlist(idx),foodType = unlist(foodType),
         idxP = unlist(idxP), Dx = unlist(Dx), 
         td = unlist(td), hd = unlist(hd), wt = unlist(wt), wh = unlist(wh),
         nDT = unlist(nDT),bias = unlist(bias), 
         bound = unlist(bound), cond = unlist(cond),
         time = unlist(time), choice = unlist(choice),
         rt = unlist(rt),
         foodType = ifelse(foodType == 2, "High Fat","Low Fat"))

# Extract parameters
params <- DataSim %>%
  select(idx,foodType,idxP,Dx,wt,wh,nDT,bias,bound,time,cond) %>%
  distinct()

params %>% group_by(Dx,foodType,cond) %>% summarise(m = mean(time), s = sd(time))

# Save
save(DataSim, file=dataFolder / "DataSim_M3.RData")

# Run recovery

for (gg in c("BN","HC")){
    for (cc in c("Neutral","Negative")) {
        Data_partial <- 
            DataSim %>%
            distinct() %>%
            filter(Dx == gg,
                    cond == cc) 
        
        idx = which(Data_partial$choice==0)
        Data_partial$RT <- Data_partial$rt
        Data_partial$RT[idx] = Data_partial$rt[idx] * -1
        
        idxP = as.numeric(ordered(Data_partial$idx)) #makes a sequentially numbered subj index
        
        idxDF = data.frame("idxR" = ordered(Data_partial$idx),
                            "idxP" = idxP)
        
        Data_partial$idxP = idxDF$idxP[idxDF$idxR == Data_partial$idx]
        
        td = Data_partial$td # taste (b1)
        
        hd = Data_partial$hd # health (b2)
        
        rtpos = Data_partial$rt
        
        y= Data_partial$RT
        
        N = length(y)
        
        ns = length(unique(idxP))
        
        fat = ifelse(Data_partial$foodType=="High Fat",2,1)
        
        nc = length(unique(fat))
        
        dat <- dump.format(list(N=N, y=y, idxP=idxP, hd=hd, td=td, rt=rtpos, ns=ns, cond=fat, nc = nc))
        
        
        inits3 <- dump.format(list( alpha.mu=2, time.mu=0.1, 
                                    alpha.pr=0.5, time.pr= 0.5, theta.mu=0.1,
                                    theta.pr=0.05,  b1.mu=0.3, b1.pr=0.05, b2.mu=0.01, b2.pr=0.05, 
                                    bias.mu=0.4,
                                    bias.kappa=1, y_pred=y,  .RNG.name="base::Super-Duper", .RNG.seed=99999))
        
        
        inits2 <- dump.format(list( alpha.mu=2.2, time.mu=-0.1, 
                                    alpha.pr=0.05, time.pr= 0.05, theta.mu=0.01,
                                    theta.pr=0.05, b1.mu=0.3, b1.pr=0.05, b2.mu=0.1, b2.pr=0.05, 
                                    bias.mu=0.4,
                                    bias.kappa=1, y_pred=y,  .RNG.name="base::Wichmann-Hill", .RNG.seed=1234))
        
        inits1 <- dump.format(list( alpha.mu=2.4, time.mu=0,  
                                    alpha.pr=0.05, time.pr= 0.05, theta.mu=0.15,
                                    theta.pr=0.05, b1.mu=0.1, b1.pr=0.05, b2.mu=0.05, b2.pr=0.05, 
                                    bias.mu=0.4,
                                    bias.kappa=1, y_pred=y, .RNG.name="base::Mersenne-Twister", .RNG.seed=6666 ))
        
        monitor = c(
        "time.mu","alpha.mu","theta.mu",
        "b1.mu","b2.mu","bias.mu",
        "time.p",
        "b1.p","b2.p", 
        "theta.p", 
        "bias",
        "alpha.p","log_lik",
        "deviance")
        
        model = file.path(scriptFolder / paste0("tssmHT_model_priors_M3.txt"))
        
        results <- run.jags(model=model, 
                            monitor=monitor, data=dat, n.chains=3, inits=c(inits1,inits2, inits3), 
                            plots = TRUE, method="parallel", module="wiener", burnin=85000, sample=15000, thin=10)
        save(results,
        Data_partial, file= paste0(recFolder,"/recovery_HtSSM_FIT_M3_Dx-",gg,"_","Cond-",cc,"_rawRatings_converted.RData"))   
    }
}

# RECOVERY ASSESSMENT

# Load estimated parameters
df.fit.full <- NULL
for (cc in c('Negative',"Neutral")) {
  for (gg in c("BN",'HC')) {
    load(paste0(resFolder,"/params_HtSSM_FIT_M3_Dx-",gg,"_","Cond-",cc,"_rawRatings_converted.RData"))
  chain=rbind(results$mcmc[[1]], results$mcmc[[2]], results$mcmc[[3]])

  subj <- unique(Data_partial$idxP)
  df.fit <- NULL
  for (s in subj){
    
    wt_lf = mean(chain[,c( paste( c("b1.p[",toString(s),",1]"), collapse = ""))])
    wh_lf = mean(chain[,c( paste( c("b2.p[",toString(s),",1]"), collapse = ""))])
    wt_hf = mean(chain[,c( paste( c("b1.p[",toString(s),",2]"), collapse = ""))])
    wh_hf = mean(chain[,c( paste( c("b2.p[",toString(s),",2]"), collapse = ""))])
    bound = mean(chain[,c( paste( c("alpha.p[",toString(s),"]"), collapse = ""))])
    nDT = mean(chain[,c( paste( c("theta.p[",toString(s),"]"), collapse = ""))])
    time_lf = mean(chain[,c( paste( c("time.p[",toString(s),",1]"), collapse = ""))])
    time_hf = mean(chain[,c( paste( c("time.p[",toString(s),",2]"), collapse = ""))])
    bias = mean(chain[,c( paste( c("bias[",toString(s),"]"), collapse = ""))])
    
    df.fit$bound = bound
    df.fit$wt_lf = wt_lf
    df.fit$wh_lf = wh_lf
    df.fit$wt_hf = wt_hf
    df.fit$wh_hf = wh_hf
    df.fit$nDT = nDT

    df.fit$time_lf = time_lf
    df.fit$time_hf = time_hf
    df.fit$bias = bias

    df.fit$idxP = s
    df.fit$Dx = gg
    df.fit$cond = cc
    df.fit$idx = Data_partial$idx[Data_partial$idxP==s][1]
  
    df.fit <- df.fit %>% as.data.frame() %>% distinct() 
    
    df.fit.full <- rbind(df.fit.full,df.fit)
  }
  rm(Data_partial,results)
  
}
}
df.fit.full$type = "Estimated"

# Load recovery params
df.rec.full <- NULL

for (cc in c('Negative',"Neutral")) {
  for (gg in c("BN",'HC')) {
    load(paste0(recFolder,"/recovery_HtSSM_FIT_M3_Dx-",gg,"_","Cond-",cc,"_rawRatings_converted.RData"))
    
    chain=rbind(results$mcmc[[1]], results$mcmc[[2]], results$mcmc[[3]])
    
    subj <- unique(Data_partial$idxP)
    df.rec <- NULL
    for (s in subj){
      
      wt_lf = mean(chain[,c( paste( c("b1.p[",toString(s),",1]"), collapse = ""))])
      wh_lf = mean(chain[,c( paste( c("b2.p[",toString(s),",1]"), collapse = ""))])
      wt_hf = mean(chain[,c( paste( c("b1.p[",toString(s),",2]"), collapse = ""))])
      wh_hf = mean(chain[,c( paste( c("b2.p[",toString(s),",2]"), collapse = ""))])
      bound = mean(chain[,c( paste( c("alpha.p[",toString(s),"]"), collapse = ""))])
      nDT = mean(chain[,c( paste( c("theta.p[",toString(s),"]"), collapse = ""))])
      time_lf = mean(chain[,c( paste( c("time.p[",toString(s),",1]"), collapse = ""))])
      time_hf = mean(chain[,c( paste( c("time.p[",toString(s),",2]"), collapse = ""))])
      bias = mean(chain[,c( paste( c("bias[",toString(s),"]"), collapse = ""))])

      df.rec$bound = bound
      df.rec$wt_lf = wt_lf
      df.rec$wh_lf = wh_lf
      df.rec$wt_hf = wt_hf
      df.rec$wh_hf = wh_hf
      df.rec$nDT = nDT

      df.rec$time_lf = time_lf
      df.rec$time_hf = time_hf
      df.rec$bias = bias

      df.rec$idxP = s
      df.rec$Dx = gg
      df.rec$cond = cc
      df.rec$idx = Data_partial$idx[Data_partial$idxP==s][1]
      
      df.rec <- df.rec %>% as.data.frame() %>% distinct() 
      
      df.rec.full <- rbind(df.rec.full,df.rec)
    }
    rm(Data_partial,results)
  }
}

df.rec.full$type = "Recovered"


df.full <- rbind(df.rec.full, df.fit.full)


df.long <- df.full %>% 
  pivot_longer(cols = c(wt_lf,wh_lf,wt_hf,wh_hf,bound,nDT,bias,time_hf,time_lf),
               names_to = "param") %>%
  pivot_wider(names_from = type,
              values_from = value)

# Combine
params <- df.long %>%
  mutate(foodType = ifelse(grepl("_lf",param),"Low-Fat",
                           ifelse(grepl("_hf",param),"High-Fat", "NA")),
         params = recode(param,
                         "wt_lf" = "omega[taste]",
                         "wt_hf" = "omega[taste]",
                         "wh_lf" = "omega[health]",
                         "wh_hf" = "omega[health]",
                         "nDT" = "tau[ND]",
                         "time_lf" = "tau[s]",
                         "time_hf" = "tau[s]",
                         "bound" = "alpha",
                         "bias" = "z"
         )
  )

figureS1 = params %>%
  ggplot(aes(x = Estimated, y = Recovered)) +
  theme_pubr(base_size = 18) +
  facet_wrap(~params, scales = "free",
             labeller = label_parsed) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "gray") +
  geom_point() + #aes( color = Dx)) +
  stat_cor(method = "pearson") +
  labs(x = "Estimated",
       y = "Recovered")

ggsave(file = figPath / "figureS1.png", plot = figureS1, width = 12, height = 8)
