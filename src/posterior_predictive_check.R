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
# 09/23/24      Blair Shevlin                         Wrote code to generate simulations and figures
# 09/24/24      Blair Shevlin                         Wrote code for regressions and ranked tests

required_packages = c("DEoptim", "Rcpp", "plyr", "parallel", "RcppParallel","stats4","pracma","runjags","tidyverse","ggpubr", "fs", "here", "lmerTest","ARTool")

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

base_path = path(here())
#here
dataFolder = path(base_path) / "data"
scriptFolder = path(base_path) / "src" / "stDDM"
resFolder = path(base_path) / "results" / "stDDM" / "estimation"
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

DataSim_m0<-NULL
DataSim_m1<-NULL
DataSim_m2<-NULL
DataSim_m3<-NULL

# Model 0 (Null Model)
for (gg in c("BN","HC")) {
  subj <- unique(DDataSim$idxP[DDataSim$Dx == gg])
  for (cc in c("Negative","Neutral")) {
    load(paste0(resFolder,"/params_HtSSM_FIT_M0_Dx-",gg,"_","Cond-",cc,"_rawRatings_converted.RData"))
    chain=rbind(results$mcmc[[1]], results$mcmc[[2]], results$mcmc[[3]])
    for (s in subj){
      wt = mean(chain[,c( paste( c("b1.p[",toString(s),"]"), collapse = ""))])

      wh = mean(chain[,c( paste( c("b2.p[",toString(s),"]"), collapse = ""))])

      nDT = mean(chain[,c( paste( c("theta.p[",toString(s),"]"), collapse = ""))])

      bias = mean(chain[,c( paste( c("bias[",toString(s),"]"), collapse = ""))])
      
      bound = mean(chain[,c( paste( c("alpha.p[",toString(s),"]"), collapse = ""))])
      time = mean(chain[,c( paste( c("time.p[",toString(s),"]"), collapse = ""))])

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
        rts = ddm2_parallel( wt,wh,bound,nDT,time,bias,trial_mat_hf$td[i],trial_mat_hf$hd[i],1,trial_mat_hf$count[i] )
        for (j in 1:length(rts)) {
          good = 0
          # Make sure RT is within good boundaries, otherwise re-run
          while (good == 0) {
            if (abs(rts[j]) > 0.2 & abs(rts[j]) <= 7) {
              good = 1
            }
            rts[j] = ddm2_parallel(wt,wh,bound,nDT,time,bias,trial_mat_hf$td[i],trial_mat_hf$hd[i],1,1)
          }
          
          DataSim_Mat_Temp = NULL
          DataSim_Mat_Temp$idx = unique(DDataSim$idx[DDataSim$idxP==s & DDataSim$Dx == gg])
          DataSim_Mat_Temp$idxP = s
          DataSim_Mat_Temp$Dx = gg
          DataSim_Mat_Temp$cond = cc
          DataSim_Mat_Temp$foodType = trial_mat_hf$foodType[i]
          DataSim_Mat_Temp$td = trial_mat_hf$td[i]
          DataSim_Mat_Temp$hd = trial_mat_hf$hd[i]
          DataSim_Mat_Temp$wt = wt
          DataSim_Mat_Temp$wh = wh
          DataSim_Mat_Temp$nDT = nDT
          DataSim_Mat_Temp$bias = bias
          DataSim_Mat_Temp$bound = bound
          DataSim_Mat_Temp$time = time
          
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
        rts = ddm2_parallel(wt,wh,bound,nDT,time,bias,trial_mat_lf$td[i],trial_mat_lf$hd[i],1,trial_mat_lf$count[i])
        for (j in 1:length(rts)) {
          good = 0
          # Make sure RT is within good boundaries, otherwise re-run
          while (good == 0) {
            if (abs(rts[j]) > 0.2 & abs(rts[j]) <= 7) {
              good = 1
            }
            rts[j] = ddm2_parallel(wt,wh,bound,nDT,time,bias,trial_mat_lf$td[i],trial_mat_lf$hd[i],1,1)
          }
          
          DataSim_Mat_Temp = NULL
          DataSim_Mat_Temp$idx = unique(DDataSim$idx[DDataSim$idxP==s & DDataSim$Dx == gg])
          DataSim_Mat_Temp$idxP = s
          DataSim_Mat_Temp$Dx = gg
          DataSim_Mat_Temp$cond = cc
          DataSim_Mat_Temp$foodType = trial_mat_lf$foodType[i]
          DataSim_Mat_Temp$td = trial_mat_lf$td[i]
          DataSim_Mat_Temp$hd = trial_mat_lf$hd[i]
          DataSim_Mat_Temp$wt = wt
          DataSim_Mat_Temp$wh = wh
          DataSim_Mat_Temp$nDT = nDT
          DataSim_Mat_Temp$bias = bias
          DataSim_Mat_Temp$bound = bound
          DataSim_Mat_Temp$time = time
          
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
      DataSim_m0<-rbind(DataSim_m0,DataSim_T)
        
    }
  }
}
# Convert from list
DataSim_m0 <- DataSim_m0 %>% 
  mutate(idx = unlist(idx),foodType = unlist(foodType),
         idxP = unlist(idxP), Dx = unlist(Dx), 
         td = unlist(td), hd = unlist(hd), wt = unlist(wt), wh = unlist(wh),
         nDT = unlist(nDT),bias = unlist(bias), 
         bound = unlist(bound), cond = unlist(cond),
         time = unlist(time), choice = unlist(choice),
         rt = unlist(rt),
         foodType = ifelse(foodType == 2, "High Fat","Low Fat"))
# Save
save(DataSim_m0, file=dataFolder / "DataSim_M0.RData")


# Model 1 (Attribute weights depend on food type)
for (gg in c("BN","HC")) {
  subj <- unique(DDataSim$idxP[DDataSim$Dx == gg])
  for (cc in c("Negative","Neutral")) {
    load(paste0(resFolder,"/params_HtSSM_FIT_M1_Dx-",gg,"_","Cond-",cc,"_rawRatings_converted.RData"))
    chain=rbind(results$mcmc[[1]], results$mcmc[[2]], results$mcmc[[3]])
    for (s in subj){
      wt_lf = mean(chain[,c( paste( c("b1.p[",toString(s),",1]"), collapse = ""))])
      wh_lf = mean(chain[,c( paste( c("b2.p[",toString(s),",1]"), collapse = ""))])
      
      wt_hf = mean(chain[,c( paste( c("b1.p[",toString(s),",2]"), collapse = ""))])
      wh_hf = mean(chain[,c( paste( c("b2.p[",toString(s),",2]"), collapse = ""))])
      
      nDT = mean(chain[,c( paste( c("theta.p[",toString(s),"]"), collapse = ""))])
      
      bias = mean(chain[,c( paste( c("bias[",toString(s),"]"), collapse = ""))])
      
      bound = mean(chain[,c( paste( c("alpha.p[",toString(s),"]"), collapse = ""))])
      time = mean(chain[,c( paste( c("time.p[",toString(s),"]"), collapse = ""))])

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
        rts = ddm2_parallel( wt_hf,wh_hf,bound,nDT,time,bias,trial_mat_hf$td[i],trial_mat_hf$hd[i],1,trial_mat_hf$count[i] )
        for (j in 1:length(rts)) {
          good = 0
          # Make sure RT is within good boundaries, otherwise re-run
          while (good == 0) {
            if (abs(rts[j]) > 0.2 & abs(rts[j]) <= 7) {
              good = 1
            }
            rts[j] = ddm2_parallel(wt_lf,wh_lf,bound,nDT,time,bias,trial_mat_hf$td[i],trial_mat_hf$hd[i],1,1)
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
          DataSim_Mat_Temp$time = time
          
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
        rts = ddm2_parallel(wt_lf,wh_lf,bound,nDT,time,bias,trial_mat_lf$td[i],trial_mat_lf$hd[i],1,trial_mat_lf$count[i])
        for (j in 1:length(rts)) {
          good = 0
          # Make sure RT is within good boundaries, otherwise re-run
          while (good == 0) {
            if (abs(rts[j]) > 0.2 & abs(rts[j]) <= 7) {
              good = 1
            }
            rts[j] = ddm2_parallel(wt_lf,wh_lf,bound,nDT,time,bias,trial_mat_lf$td[i],trial_mat_lf$hd[i],1,1)
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
          DataSim_Mat_Temp$time = time
          
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
      DataSim_m1<-rbind(DataSim_m1,DataSim_T)
      
    }
  }
}
# Convert from list
DataSim_m1 <- DataSim_m1 %>% 
  mutate(idx = unlist(idx),foodType = unlist(foodType),
         idxP = unlist(idxP), Dx = unlist(Dx), 
         td = unlist(td), hd = unlist(hd), wt = unlist(wt), wh = unlist(wh),
         nDT = unlist(nDT),bias = unlist(bias), 
         bound = unlist(bound), cond = unlist(cond),
         time = unlist(time), choice = unlist(choice),
         rt = unlist(rt),
         foodType = ifelse(foodType == 2, "High Fat","Low Fat"))
# Save
save(DataSim_m1, file=dataFolder / "DataSim_M1.RData")

# Model 2 (Time depends on food type)
for (gg in c("BN","HC")) {
  subj <- unique(DDataSim$idxP[DDataSim$Dx == gg])
  for (cc in c("Negative","Neutral")) {
    load(paste0(resFolder,"/params_HtSSM_FIT_M2_Dx-",gg,"_","Cond-",cc,"_rawRatings_converted.RData"))
    chain=rbind(results$mcmc[[1]], results$mcmc[[2]], results$mcmc[[3]])
    for (s in subj){
      wt = mean(chain[,c( paste( c("b1.p[",toString(s),"]"), collapse = ""))])
      wh = mean(chain[,c( paste( c("b2.p[",toString(s),"]"), collapse = ""))])

      nDT = mean(chain[,c( paste( c("theta.p[",toString(s),"]"), collapse = ""))])
      
      bias = mean(chain[,c( paste( c("bias[",toString(s),"]"), collapse = ""))])
      
      bound = mean(chain[,c( paste( c("alpha.p[",toString(s),"]"), collapse = ""))])
      time_lf = mean(chain[,c( paste( c("time.p[",toString(s),",1]"), collapse = ""))])
      time_hf = mean(chain[,c( paste( c("time.p[",toString(s),",2]"), collapse = ""))])
      
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
        rts = ddm2_parallel( wt,wh,bound,nDT,time_hf,bias,trial_mat_hf$td[i],trial_mat_hf$hd[i],1,trial_mat_hf$count[i] )
        for (j in 1:length(rts)) {
          good = 0
          # Make sure RT is within good boundaries, otherwise re-run
          while (good == 0) {
            if (abs(rts[j]) > 0.2 & abs(rts[j]) <= 7) {
              good = 1
            }
            rts[j] = ddm2_parallel(wt,wh,bound,nDT,time_hf,bias,trial_mat_hf$td[i],trial_mat_hf$hd[i],1,1)
          }
          
          DataSim_Mat_Temp = NULL
          DataSim_Mat_Temp$idx = unique(DDataSim$idx[DDataSim$idxP==s & DDataSim$Dx == gg])
          DataSim_Mat_Temp$idxP = s
          DataSim_Mat_Temp$Dx = gg
          DataSim_Mat_Temp$cond = cc
          DataSim_Mat_Temp$foodType = trial_mat_hf$foodType[i]
          DataSim_Mat_Temp$td = trial_mat_hf$td[i]
          DataSim_Mat_Temp$hd = trial_mat_hf$hd[i]
          DataSim_Mat_Temp$wt = wt
          DataSim_Mat_Temp$wh = wh
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
        rts = ddm2_parallel(wt,wh,bound,nDT,time_lf,bias,trial_mat_lf$td[i],trial_mat_lf$hd[i],1,trial_mat_lf$count[i])
        for (j in 1:length(rts)) {
          good = 0
          # Make sure RT is within good boundaries, otherwise re-run
          while (good == 0) {
            if (abs(rts[j]) > 0.2 & abs(rts[j]) <= 7) {
              good = 1
            }
            rts[j] = ddm2_parallel(wt,wh,bound,nDT,time_lf,bias,trial_mat_lf$td[i],trial_mat_lf$hd[i],1,1)
          }
          
          DataSim_Mat_Temp = NULL
          DataSim_Mat_Temp$idx = unique(DDataSim$idx[DDataSim$idxP==s & DDataSim$Dx == gg])
          DataSim_Mat_Temp$idxP = s
          DataSim_Mat_Temp$Dx = gg
          DataSim_Mat_Temp$cond = cc
          DataSim_Mat_Temp$foodType = trial_mat_lf$foodType[i]
          DataSim_Mat_Temp$td = trial_mat_lf$td[i]
          DataSim_Mat_Temp$hd = trial_mat_lf$hd[i]
          DataSim_Mat_Temp$wt = wt
          DataSim_Mat_Temp$wh = wh
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
      DataSim_m2<-rbind(DataSim_m2,DataSim_T)
      
    }
  }
}
# Convert from list
DataSim_m2 <- DataSim_m2 %>% 
  mutate(idx = unlist(idx),foodType = unlist(foodType),
         idxP = unlist(idxP), Dx = unlist(Dx), 
         td = unlist(td), hd = unlist(hd), wt = unlist(wt), wh = unlist(wh),
         nDT = unlist(nDT),bias = unlist(bias), 
         bound = unlist(bound), cond = unlist(cond),
         time = unlist(time), choice = unlist(choice),
         rt = unlist(rt),
         foodType = ifelse(foodType == 2, "High Fat","Low Fat"))
# Save
save(DataSim_m2, file=dataFolder / "DataSim_M2.RData")

# Model 3 (Time, b1, b2 depends on food type)
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
          good = 0
          # Make sure RT is within good boundaries, otherwise re-run
          while (good == 0) {
            if (abs(rts[j]) > 0.2 & abs(rts[j]) <= 7) {
              good = 1
            }
            rts[j] = ddm2_parallel(wt_hf,wh_hf,bound,nDT,time_hf,bias,trial_mat_hf$td[i],trial_mat_hf$hd[i],1,1)
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
save(DataSim_m3, file=dataFolder / "DataSim_M3.RData")

# Load Simulated Data
load(file=dataFolder / "DataSim_M0.RData")
load(file=dataFolder / "DataSim_M1.RData")
load(file=dataFolder / "DataSim_M2.RData")
load(file=dataFolder / "DataSim_M3.RData") 

DataSim_m0 = DataSim_m0 %>%
  mutate(cond = factor(cond, levels = c("Neutral","Negative")),
         foodType = factor(foodType, levels = c("Low Fat","High Fat")),
         Dx = factor(Dx, levels = c("HC","BN")),
         m = "M0"
  )
DataSim_m1 = DataSim_m1 %>%
  mutate(cond = factor(cond, levels = c("Neutral","Negative")),
         foodType = factor(foodType, levels = c("Low Fat","High Fat")),
         Dx = factor(Dx, levels = c("HC","BN")),
         m = "M1"
  )
DataSim_m2 = DataSim_m2 %>%
  mutate(cond = factor(cond, levels = c("Neutral","Negative")),
         foodType = factor(foodType, levels = c("Low Fat","High Fat")),
         Dx = factor(Dx, levels = c("HC","BN")),
         m = "M2"
  )
DataSim_m3 = DataSim_m3 %>%
  mutate(cond = factor(cond, levels = c("Neutral","Negative")),
         foodType = factor(foodType, levels = c("Low Fat","High Fat")),
         Dx = factor(Dx, levels = c("HC","BN")),
         m = "M3"
  )

# Recode Empirical Data
Data<-Data %>%
  mutate(cond = factor(cond, levels = c("Neutral","Negative")),
         foodType = factor(foodType, levels = c("Low Fat","High Fat")),
         Dx = factor(Dx, levels = c("HC","BN")),
         m = "Empirical"
         ) %>%
  select(!X)

# Set contrasts
contrasts(Data$cond) <- c(-1,1)
contrasts(Data$foodType) <- c(-1,1)
contrasts(Data$Dx) <- c(-1,1)

contrasts(DataSim_m0$cond) <- c(-1,1)
contrasts(DataSim_m0$foodType) <- c(-1,1)
contrasts(DataSim_m0$Dx) <- c(-1,1)

contrasts(DataSim_m1$cond) <- c(-1,1)
contrasts(DataSim_m1$foodType) <- c(-1,1)
contrasts(DataSim_m1$Dx) <- c(-1,1)

contrasts(DataSim_m2$cond) <- c(-1,1)
contrasts(DataSim_m2$foodType) <- c(-1,1)
contrasts(DataSim_m2$Dx) <- c(-1,1)

contrasts(DataSim_m3$cond) <- c(-1,1)
contrasts(DataSim_m3$foodType) <- c(-1,1)
contrasts(DataSim_m3$Dx) <- c(-1,1)

# Can we replicate GLM1?

glm.1.emp <- glmer(data = Data,
               formula = choice ~
                 Dx * cond * foodType  +
                 (1 + cond * foodType |idx),
               family=binomial,
               control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=20000)))
summary(glm.1.emp)
# Key features:
## Group (b = -0.48, p = 0.009)
## Food Type (b = -0.62, p < 0.001)
## Group x Food Type (b = -0.27, p = 0.047)
## MARGINAL: Group x Affect x Food Type (b = -0.10, p = 0.056)


# Model 0
glm.1.m0 <- glmer(data = DataSim_m0,
                   formula = choice ~
                     Dx * cond * foodType  +
                     (1 + cond * foodType |idx),
                   family=binomial,
                   control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=20000)))
summary(glm.1.m0)
## Group (b = -0.34, p = 0.006) - yes
## Food Type (b = -0.01, p = 0.848) - no!
## Group x Food Type (b = -0.05, p < 0.001) - yes!
## Group x Affect x Food Type (b = -0.06, p = 0.205) - no


# Model 1 (Attribute Weights)
glm.1.m1 <- glmer(data = DataSim_m1,
                  formula = choice ~
                    Dx * cond * foodType  +
                    (1 + cond * foodType |idx),
                  family=binomial,
                  control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=20000)))
summary(glm.1.m1)
## Group (b = -0.28, p = 0.053) - marginal!
## Food Type (b = 0.42, p < 0.001) - wrong direction
## Group x Food Type (b = -0.22, p = 0.001) - yes!
## Group x Affect x Food Type (b = 0.02, p = 0.552) - no
## NEW: Group x Affect (b = 0.17, p = 0.020) - not in original


# Model 2 (Onset Time)
glm.1.m2 <- glmer(data = DataSim_m2,
                  formula = choice ~
                    Dx * cond * foodType  +
                    (1 + cond * foodType |idx),
                  family=binomial,
                  control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=20000)))
summary(glm.1.m2)
## Group (b = -0.42, p = 0.002) - yes
## Food Type (b = -0.16, p = .128) - no!
## Group x Food Type (b = -0.36, p < 0.001) - yes!
## Group x Affect x Food Type (b = -0.07, p = 0.195) - yes


# Model 3 (FULL)
glm.1.m3 <- glmer(data = DataSim_m3,
                  formula = choice ~
                    Dx * cond * foodType  +
                    (1 + cond * foodType |idx),
                  family=binomial,
                  control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=20000)))
summary(glm.1.m3)
## Group (b = -0.36, p = 0.007) - yes
## Food Type (b = -0.36, p < 0.001) - yes!
## Group x Food Type (b = -0.27, p = 0.002) - yes!
## Group x Affect x Food Type (b = 0.03, p = 0.562) - no

# Summarise data
sim_summary_m0 = DataSim_m0 %>%
  group_by(idx,Dx,cond,foodType) %>%
  summarise(mch = mean(choice),
            mrt = median(rt),
            mQ1 = quantile(rt, probs=.1),
            mQ9 = quantile(rt, probs=.9)) %>%
  mutate(type = "sim")
sim_summary_m1 = DataSim_m1 %>%
  group_by(idx,Dx,cond,foodType) %>%
  summarise(mch = mean(choice),
            mrt = median(rt),
            mQ1 = quantile(rt, probs=.1),
            mQ9 = quantile(rt, probs=.9)) %>%
  mutate(type = "sim")
sim_summary_m2 = DataSim_m2 %>%
  group_by(idx,Dx,cond,foodType) %>%
  summarise(mch = mean(choice),
            mrt = median(rt),
            mQ1 = quantile(rt, probs=.1),
            mQ9 = quantile(rt, probs=.9)) %>%
  mutate(type = "sim")
sim_summary_m3 = DataSim_m3 %>%
  group_by(idx,Dx,cond,foodType) %>%
  summarise(mch = mean(choice),
            mrt = median(rt),
            mQ1 = quantile(rt, probs=.1),
            mQ9 = quantile(rt, probs=.9)) %>%
  mutate(type = "sim")

emp_summary = Data %>%
  group_by(idx,Dx,cond,foodType) %>%
  summarise(mch = mean(choice),
            mrt = median(rt),
            mQ1 = quantile(rt, probs = .1),
            mQ9 = quantile(rt, probs = .9)) %>%
  mutate(type = "emp")

dat_summary_m0 = rbind(sim_summary_m0, emp_summary) %>%
  mutate(cond = factor(cond, levels = c("Neutral","Negative")),
         foodType = factor(foodType, levels = c("Low Fat","High Fat")),
         Dx = factor(Dx, levels = c("HC","BN")),
         source = factor(type,levels=c("emp","sim"))
  )
dat_summary_m1 = rbind(sim_summary_m1, emp_summary) %>%
  mutate(cond = factor(cond, levels = c("Neutral","Negative")),
         foodType = factor(foodType, levels = c("Low Fat","High Fat")),
         Dx = factor(Dx, levels = c("HC","BN")),
         source = factor(type,levels=c("emp","sim"))
  )
dat_summary_m2 = rbind(sim_summary_m2, emp_summary) %>%
  mutate(cond = factor(cond, levels = c("Neutral","Negative")),
         foodType = factor(foodType, levels = c("Low Fat","High Fat")),
         Dx = factor(Dx, levels = c("HC","BN")),
         source = factor(type,levels=c("emp","sim"))
  )
dat_summary_m3 = rbind(sim_summary_m3, emp_summary) %>%
  mutate(cond = factor(cond, levels = c("Neutral","Negative")),
         foodType = factor(foodType, levels = c("Low Fat","High Fat")),
         Dx = factor(Dx, levels = c("HC","BN")),
         source = factor(type,levels=c("emp","sim"))
  )

nrow(dat_summary_m0)
nrow(dat_summary_m1)
nrow(dat_summary_m2)
nrow(dat_summary_m3)

# Use Aligned rank
art_model_mch_m0 <- art(mch ~ Dx*cond*foodType * source + (1|idx), data = dat_summary_m0)
art_model_mch_m1 <- art(mch ~ Dx*cond*foodType * source + (1|idx), data = dat_summary_m1)
art_model_mch_m2 <- art(mch ~ Dx*cond*foodType * source + (1|idx), data = dat_summary_m2)
art_model_mch_m3 <- art(mch ~ Dx*cond*foodType * source + (1|idx), data = dat_summary_m3)

art_model_mrt_m0 <- art(mrt ~ Dx*cond*foodType * source + (1|idx), data = dat_summary_m0)
art_model_mrt_m1 <- art(mrt ~ Dx*cond*foodType * source + (1|idx),  data = dat_summary_m1)
art_model_mrt_m2 <- art(mrt ~ Dx*cond*foodType * source + (1|idx), data = dat_summary_m2)
art_model_mrt_m3 <- art(mrt ~ Dx*cond*foodType * source + (1|idx), data = dat_summary_m3)

# Perform the ANOVA on aligned ranks

#M0
anova(art_model_mch_m0)
# source: Main effect and foodType
anova(art_model_mrt_m0)
# source: Main effect

anova(art_model_mch_m1)
# source: Main effect, Group and foodType
anova(art_model_mrt_m1)
# ns

anova(art_model_mch_m2)
# source: Main effect and foodType
anova(art_model_mrt_m2)
# ns

anova(art_model_mch_m3)
# source: Main effect
anova(art_model_mrt_m3)
# ns

# Figures (not-used)
dat_summary_m3 %>%
  group_by(idx,Dx,cond,type) %>%
  summarise(mchoice = mean(mch),
            mrt = mean(mrt)) %>%
  group_by(Dx,cond) %>%
  summarise(mch_diff = mean(mchoice[type=="emp"] - mchoice[type=="sim"]),
          mch_diff_pct = mch_diff / mean(mchoice[type=="emp"] ) )
  
  
  mutate(acc_diff = mean(mch[type == "emp"]) - mean(mch[type == "sim"]),
         rt_diff = median(mrt[type == "emp"]) - median(mrt[type == "sim"]),
         acc_diff_pct = ifelse(acc_diff == 0, acc_diff, acc_diff/ median(mch[type == "emp"]) ),
         rt_diff_pct = ifelse(rt_diff == 0, rt_diff, rt_diff/ median(mrt[type == "emp"]) )
         ) %>%
  ungroup() %>%
  summarise(avg_acc_diff = mean(acc_diff_pct,na.rm=T),
            avg_rt_diff = mean(rt_diff_pct,na.rm=T)
            )


dat_summary_m3 %>%
  mutate(type = factor(type, levels = c("emp","sim"),
                       labels = c("Empirical","Simulated"))) %>%
  ggplot(aes(x = foodType, y = mch, color = type, fill = type)) +
  theme_pubr(base_size = 18) +
  facet_grid(Dx~cond) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray", linewidth = 1.1) +
  stat_summary(geom = "bar",alpha=.25,
        # position=position_dodge(width = 1)
        ) +
  scale_color_brewer(type = "qual",palette = 6,direction = 1) +
  scale_fill_brewer(type = "qual",palette = 6,direction = 1) +
  labs(x = element_blank(),
       y = "Choice proportion\n(over reference item)",
       color = "Data type",
       fill = "Data type")
  
  
dat_summary_m3 %>%
  mutate(type = factor(type, levels = c("emp","sim"),
                       labels = c("Empirical","Simulated"))) %>%
  ggplot(aes(x = foodType, y = mrt, color = type, fill = type)) +
  theme_pubr(base_size = 18) +
  facet_grid(Dx~cond) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray", linewidth = 1.1) +
  stat_summary(geom = "bar",alpha=.25,
  ) +
  scale_color_brewer(type = "qual",palette = 6,direction = 1) +
  scale_fill_brewer(type = "qual",palette = 6,direction = 1) +
  labs(x = element_blank(),
       y = "Resspone time (s)",
       color = "Data type",
       fill = "Data type")

