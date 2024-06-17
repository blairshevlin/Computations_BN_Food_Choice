# fit_stddm_models.R - Code for fitting models from "Negative affect influences the computations underlying food choice in bulimia nervosa" 
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
# 06/17/24      Blair Shevlin                            Wrote code for original manuscript


# This version fits foodType WITHIN a condition

# Packages required
required_packages <- c(
    "DEoptim", 
    "Rcpp",
    "parallel", 
    "here", 
    "fs",
    "RcppParallel",
    "stats4",
    "pracma",
    "runjags",
    "tidyverse",
    "loo",
    "coda"
)

# Check and install missing packages
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
resFolder = path(base_path) / "results" / "stDDM"

Data<-read.csv(file=paste0(dataFolder, "/deidentified_ChoiceData.csv.csv"))

# M2 (Old M17): food-type -> Time

# M3 (Old M18): food-type -> Time, b1, b2

# M1 (Old M19): food-type -> b1, b2

# M0 (Old M31): Null

# Fit initial
for (m in c("M0","M1","M2","M3") ) {
  for (gg in c("bn","hc")){
    for (cc in c("Neutral","Negative")) {
      Data_partial <- 
        Data %>%
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
      
      fat = ifelse(Data_partial$foodType=="hf",2,1)
      
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
      
      model = file.path(scriptFolder / paste0("tssmHT_model_priors_",m,".txt"))
      
      results <- run.jags(model=model, 
                          monitor=monitor, data=dat, n.chains=3, inits=c(inits1,inits2, inits3), 
                          plots = TRUE, method="parallel", module="wiener", burnin=85000, sample=15000, thin=10)
      save(results,
        Data_partial, file= paste0(resFolder,"/params_HtSSM_FIT_",m,"_Dx-",gg,"_","Cond-",cc,".RData")) 
      
    }
  }
  
}


# Calculate model fit metrics
dic.df = NULL
models_of_interest= c("M0","M1","M2","M3")

nM = length(models_of_interest)
for (gg in c("hc","bn")) {
  for (cc in c("Negative","Neutral")) {
    for (m in models_of_interest) {
      load(file.path(resultFolder,paste("params_HtSSM_FIT_",m,"_Dx-",gg,"_","Cond-",cc,".RData",sep="")) )
      code.samples<- as.mcmc.list(results)
      mcmc.mat <- as.matrix(code.samples, chains = F)
      rm(code.samples)
      dic = round(mean(mcmc.mat[,"deviance"])  + .5 * var(mcmc.mat[,"deviance"]))
      N = nrow(Data_partial)
      chain <- rbind(results$mcmc[[1]], results$mcmc[[2]], results$mcmc[[3]])
      rm(results,mcmc.mat)
      loglik <- chain[,paste0("log_lik[",1:N,"]")]
      rm(chain)
      r_eff =  relative_eff(exp(loglik), chain_id = rep(1:3, each = 15000),
                              cores = 3)
      w = waic(loglik)
      l = loo(loglik,r_eff = r_eff)
      tmp.df = data.frame(Dx = gg,
                          cond = cc,
                          model = m,
                          dic = dic,
                           waic = w$estimates[3],
                         loo = l$estimates[3]
      )
      dic.df = rbind(dic.df,tmp.df)
    }
  }
}


dic.df %>%
  pivot_longer(cols = c(dic,waic,loo)) %>%
  group_by(Dx,cond,name) %>%
  summarise(model = model[value == min(value)])


dic.df %>%
  pivot_longer(cols = c(dic,waic,loo)) %>%
  group_by(Dx,name,model) %>%
  dplyr::summarise(val = sum(value) ) %>%
  group_by(Dx,name) %>%
  dplyr::summarise(model = model[val == min(val)])

dic.df %>%
  pivot_longer(cols = c(dic,waic,loo)) %>%
  filter(name != "dic") %>%
  group_by(name,Dx,model) %>%
  summarise(val = round(sum(value)))%>%
  mutate(model= factor(model,
                      levels = model_order),
         Dx = factor(Dx, levels = c("hc","bn"))
        ) %>%
  arrange(model,Dx) %>%
  as.data.frame()

dic.df %>%
  pivot_longer(cols = c(dic,waic,loo)) %>%
  filter(name != "dic") %>%
  group_by(name,Dx,model) %>%
  summarise(val = round(sum(value)))%>%
  mutate(Dx = factor(Dx, levels = c("hc","bn"))
  ) %>%
  filter(model == model[val == min(val)])


dic.df %>%
  filter(model %in% models_of_interest) %>%
  pivot_longer(cols = c(dic,waic,loo)) %>%
  group_by(name,model,Dx) %>%
  summarise(val = value[cond == "Negative"] + value[cond =="Neutral"])%>%
  group_by(name,Dx) %>%
  summarise(win = model[val == min(val)],
            v = val[val == min(val)])
  

dic.df %>%
  #filter(model %in% models_of_interest) %>%
  mutate(model = factor(model,
                        levels = c("M0","M1","M2","M3"))) %>%
  pivot_longer(cols = c(dic,waic,loo)) %>%
  group_by(model,Dx,name) %>%
  summarise(val = value[cond == "Negative"] + value[cond =="Neutral"]) %>%
  as.data.frame()



bn_neg = dic.df %>% filter()
#
bn_neg_waic = loo_compare(x = list(waic_m0,waic_m1,
                                   waic_m2,waic_m3))
# M20 wins
bn_neg_loo = loo_compare(x = list(loo_m16,loo_m17,loo_m18,
                                  loo_m19,loo_m20,loo_m21,loo_m22))


# Do loo compare

gg = "bn";
cc = "Negative";

m = "M16";
load(file.path(resultFolder,paste("params_HtSSM_FIT_",m,"_Dx-",gg,"_","Cond-",cc,"_rawRatings_v6B.RData",sep="")) )
code.samples<- as.mcmc.list(results)
mcmc.mat <- as.matrix(code.samples, chains = F)
dic = round(mean(mcmc.mat[,"deviance"])  + .5 * var(mcmc.mat[,"deviance"]))
N = nrow(Data_partial)
chain <- rbind(results$mcmc[[1]], results$mcmc[[2]], results$mcmc[[3]])
loglik <- chain[,paste0("log_lik[",1:N,"]")]
waic_m16 = waic(loglik)
r_eff =  relative_eff(exp(loglik), chain_id = rep(1:3, each = 15000),
                      cores = 3)
loo_m16 = loo(loglik,r_eff = r_eff,cores = 3)
m = "M17";
load(file.path(resultFolder,paste("params_HtSSM_FIT_",m,"_Dx-",gg,"_","Cond-",cc,"_rawRatings_v6B.RData",sep="")) )
code.samples<- as.mcmc.list(results)
mcmc.mat <- as.matrix(code.samples, chains = F)
dic = round(mean(mcmc.mat[,"deviance"])  + .5 * var(mcmc.mat[,"deviance"]))
N = nrow(Data_partial)
chain <- rbind(results$mcmc[[1]], results$mcmc[[2]], results$mcmc[[3]])
loglik <- chain[,paste0("log_lik[",1:N,"]")]
waic_m17 = waic(loglik)
r_eff =  relative_eff(exp(loglik), chain_id = rep(1:3, each = 15000),
                      cores = 3)
loo_m17 = loo(loglik,r_eff = r_eff,cores = 3)
m = "M18";
load(file.path(resultFolder,paste("params_HtSSM_FIT_",m,"_Dx-",gg,"_","Cond-",cc,"_rawRatings_v6B.RData",sep="")) )
code.samples<- as.mcmc.list(results)
mcmc.mat <- as.matrix(code.samples, chains = F)
dic = round(mean(mcmc.mat[,"deviance"])  + .5 * var(mcmc.mat[,"deviance"]))
N = nrow(Data_partial)
chain <- rbind(results$mcmc[[1]], results$mcmc[[2]], results$mcmc[[3]])
loglik <- chain[,paste0("log_lik[",1:N,"]")]
waic_m18 = waic(loglik)
r_eff =  relative_eff(exp(loglik), chain_id = rep(1:3, each = 15000),
                      cores = 3)
loo_m18 = loo(loglik,r_eff = r_eff,cores = 3)
m = "M19";
load(file.path(resultFolder,paste("params_HtSSM_FIT_",m,"_Dx-",gg,"_","Cond-",cc,"_rawRatings_v6B.RData",sep="")) )
code.samples<- as.mcmc.list(results)
mcmc.mat <- as.matrix(code.samples, chains = F)
dic = round(mean(mcmc.mat[,"deviance"])  + .5 * var(mcmc.mat[,"deviance"]))
N = nrow(Data_partial)
chain <- rbind(results$mcmc[[1]], results$mcmc[[2]], results$mcmc[[3]])
loglik <- chain[,paste0("log_lik[",1:N,"]")]
waic_m19 = waic(loglik)
r_eff =  relative_eff(exp(loglik), chain_id = rep(1:3, each = 15000),
                      cores = 3)
loo_m19 = loo(loglik,r_eff = r_eff,cores = 3)
m = "M20";
load(file.path(resultFolder,paste("params_HtSSM_FIT_",m,"_Dx-",gg,"_","Cond-",cc,"_rawRatings_v6B.RData",sep="")) )
code.samples<- as.mcmc.list(results)
mcmc.mat <- as.matrix(code.samples, chains = F)
dic = round(mean(mcmc.mat[,"deviance"])  + .5 * var(mcmc.mat[,"deviance"]))
N = nrow(Data_partial)
chain <- rbind(results$mcmc[[1]], results$mcmc[[2]], results$mcmc[[3]])
loglik <- chain[,paste0("log_lik[",1:N,"]")]
waic_m20 = waic(loglik)
r_eff =  relative_eff(exp(loglik), chain_id = rep(1:3, each = 15000),
                      cores = 3)
loo_m20 = loo(loglik,r_eff = r_eff,cores = 3)
m = "M21";
load(file.path(resultFolder,paste("params_HtSSM_FIT_",m,"_Dx-",gg,"_","Cond-",cc,"_rawRatings_v6B.RData",sep="")) )
code.samples<- as.mcmc.list(results)
mcmc.mat <- as.matrix(code.samples, chains = F)
dic = round(mean(mcmc.mat[,"deviance"])  + .5 * var(mcmc.mat[,"deviance"]))
N = nrow(Data_partial)
chain <- rbind(results$mcmc[[1]], results$mcmc[[2]], results$mcmc[[3]])
loglik <- chain[,paste0("log_lik[",1:N,"]")]
waic_m21 = waic(loglik)
r_eff =  relative_eff(exp(loglik), chain_id = rep(1:3, each = 15000),
                      cores = 3)
loo_m21 = loo(loglik,r_eff = r_eff,cores = 3)
m = "M22";
load(file.path(resultFolder,paste("params_HtSSM_FIT_",m,"_Dx-",gg,"_","Cond-",cc,"_rawRatings_v6B.RData",sep="")) )
code.samples<- as.mcmc.list(results)
mcmc.mat <- as.matrix(code.samples, chains = F)
dic = round(mean(mcmc.mat[,"deviance"])  + .5 * var(mcmc.mat[,"deviance"]))
N = nrow(Data_partial)
chain <- rbind(results$mcmc[[1]], results$mcmc[[2]], results$mcmc[[3]])
loglik <- chain[,paste0("log_lik[",1:N,"]")]
waic_m22 = waic(loglik)
r_eff =  relative_eff(exp(loglik), chain_id = rep(1:3, each = 15000),
                      cores = 3)
loo_m22 = loo(loglik,r_eff = r_eff,cores = 3)

bn_neg_waic = loo_compare(x = list(waic_m16,waic_m17,waic_m18,
                                   waic_m19,waic_m20,waic_m21,waic_m22))
# M20 wins
bn_neg_loo = loo_compare(x = list(loo_m16,loo_m17,loo_m18,
                                  loo_m19,loo_m20,loo_m21,loo_m22))
# M20 wins
bn_neg_waic_s = loo_compare(x = list(waic_m18,waic_m21))
bn_neg_loo_s = loo_compare(x = list(loo_m18,loo_m21))


gg = "bn";
cc = "Neutral";
m = "M16";
load(file.path(resultFolder,paste("params_HtSSM_FIT_",m,"_Dx-",gg,"_","Cond-",cc,"_rawRatings_v6B.RData",sep="")) )
code.samples<- as.mcmc.list(results)
mcmc.mat <- as.matrix(code.samples, chains = F)
dic = round(mean(mcmc.mat[,"deviance"])  + .5 * var(mcmc.mat[,"deviance"]))
N = nrow(Data_partial)
chain <- rbind(results$mcmc[[1]], results$mcmc[[2]], results$mcmc[[3]])
loglik <- chain[,paste0("log_lik[",1:N,"]")]
waic_m16 = waic(loglik)
r_eff =  relative_eff(exp(loglik), chain_id = rep(1:3, each = 15000),
                      cores = 3)
loo_m16 = loo(loglik,r_eff = r_eff,cores = 3)
m = "M17";
load(file.path(resultFolder,paste("params_HtSSM_FIT_",m,"_Dx-",gg,"_","Cond-",cc,"_rawRatings_v6B.RData",sep="")) )
code.samples<- as.mcmc.list(results)
mcmc.mat <- as.matrix(code.samples, chains = F)
dic = round(mean(mcmc.mat[,"deviance"])  + .5 * var(mcmc.mat[,"deviance"]))
N = nrow(Data_partial)
chain <- rbind(results$mcmc[[1]], results$mcmc[[2]], results$mcmc[[3]])
loglik <- chain[,paste0("log_lik[",1:N,"]")]
waic_m17 = waic(loglik)
r_eff =  relative_eff(exp(loglik), chain_id = rep(1:3, each = 15000),
                      cores = 3)
loo_m17 = loo(loglik,r_eff = r_eff,cores = 3)
m = "M18";
load(file.path(resultFolder,paste("params_HtSSM_FIT_",m,"_Dx-",gg,"_","Cond-",cc,"_rawRatings_v6B.RData",sep="")) )
code.samples<- as.mcmc.list(results)
mcmc.mat <- as.matrix(code.samples, chains = F)
dic = round(mean(mcmc.mat[,"deviance"])  + .5 * var(mcmc.mat[,"deviance"]))
N = nrow(Data_partial)
chain <- rbind(results$mcmc[[1]], results$mcmc[[2]], results$mcmc[[3]])
loglik <- chain[,paste0("log_lik[",1:N,"]")]
waic_m18 = waic(loglik)
r_eff =  relative_eff(exp(loglik), chain_id = rep(1:3, each = 15000),
                      cores = 3)
loo_m18 = loo(loglik,r_eff = r_eff,cores = 3)
m = "M19";
load(file.path(resultFolder,paste("params_HtSSM_FIT_",m,"_Dx-",gg,"_","Cond-",cc,"_rawRatings_v6B.RData",sep="")) )
code.samples<- as.mcmc.list(results)
mcmc.mat <- as.matrix(code.samples, chains = F)
dic = round(mean(mcmc.mat[,"deviance"])  + .5 * var(mcmc.mat[,"deviance"]))
N = nrow(Data_partial)
chain <- rbind(results$mcmc[[1]], results$mcmc[[2]], results$mcmc[[3]])
loglik <- chain[,paste0("log_lik[",1:N,"]")]
waic_m19 = waic(loglik)
r_eff =  relative_eff(exp(loglik), chain_id = rep(1:3, each = 15000),
                      cores = 3)
loo_m19 = loo(loglik,r_eff = r_eff,cores = 3)
m = "M20";
load(file.path(resultFolder,paste("params_HtSSM_FIT_",m,"_Dx-",gg,"_","Cond-",cc,"_rawRatings_v6B.RData",sep="")) )
code.samples<- as.mcmc.list(results)
mcmc.mat <- as.matrix(code.samples, chains = F)
dic = round(mean(mcmc.mat[,"deviance"])  + .5 * var(mcmc.mat[,"deviance"]))
N = nrow(Data_partial)
chain <- rbind(results$mcmc[[1]], results$mcmc[[2]], results$mcmc[[3]])
loglik <- chain[,paste0("log_lik[",1:N,"]")]
waic_m20 = waic(loglik)
r_eff =  relative_eff(exp(loglik), chain_id = rep(1:3, each = 15000),
                      cores = 3)
loo_m20 = loo(loglik,r_eff = r_eff,cores = 3)
m = "M21";
load(file.path(resultFolder,paste("params_HtSSM_FIT_",m,"_Dx-",gg,"_","Cond-",cc,"_rawRatings_v6B.RData",sep="")) )
code.samples<- as.mcmc.list(results)
mcmc.mat <- as.matrix(code.samples, chains = F)
dic = round(mean(mcmc.mat[,"deviance"])  + .5 * var(mcmc.mat[,"deviance"]))
N = nrow(Data_partial)
chain <- rbind(results$mcmc[[1]], results$mcmc[[2]], results$mcmc[[3]])
loglik <- chain[,paste0("log_lik[",1:N,"]")]
waic_m21 = waic(loglik)
r_eff =  relative_eff(exp(loglik), chain_id = rep(1:3, each = 15000),
                      cores = 3)
loo_m21 = loo(loglik,r_eff = r_eff,cores = 3)
m = "M22";
load(file.path(resultFolder,paste("params_HtSSM_FIT_",m,"_Dx-",gg,"_","Cond-",cc,"_rawRatings_v6B.RData",sep="")) )
code.samples<- as.mcmc.list(results)
mcmc.mat <- as.matrix(code.samples, chains = F)
dic = round(mean(mcmc.mat[,"deviance"])  + .5 * var(mcmc.mat[,"deviance"]))
N = nrow(Data_partial)
chain <- rbind(results$mcmc[[1]], results$mcmc[[2]], results$mcmc[[3]])
loglik <- chain[,paste0("log_lik[",1:N,"]")]
waic_m22 = waic(loglik)
r_eff =  relative_eff(exp(loglik), chain_id = rep(1:3, each = 15000),
                      cores = 3)
loo_m22 = loo(loglik,r_eff = r_eff,cores = 3)

bn_neu_waic = loo_compare(x = list(waic_m16,waic_m17,waic_m18,
                                   waic_m19,waic_m20,waic_m21,waic_m22))
# M15 wins (followed by M18)
bn_neu_loo = loo_compare(x = list(loo_m16,loo_m17,loo_m18,
                                  loo_m19,loo_m20,loo_m21,loo_m22))
# M17 Wins

bn_neu_waic_s = loo_compare(x = list(waic_m17,waic_m19,waic_m18))
bn_neu_loo_s = loo_compare(x = list(loo_m17,loo_m19,loo_m18))

gg = "hc";
cc = "Negative";

m = "M16";
load(file.path(resultFolder,paste("params_HtSSM_FIT_",m,"_Dx-",gg,"_","Cond-",cc,"_rawRatings_v6B.RData",sep="")) )
code.samples<- as.mcmc.list(results)
mcmc.mat <- as.matrix(code.samples, chains = F)
dic = round(mean(mcmc.mat[,"deviance"])  + .5 * var(mcmc.mat[,"deviance"]))
N = nrow(Data_partial)
chain <- rbind(results$mcmc[[1]], results$mcmc[[2]], results$mcmc[[3]])
loglik <- chain[,paste0("log_lik[",1:N,"]")]
waic_m16 = waic(loglik)
r_eff =  relative_eff(exp(loglik), chain_id = rep(1:3, each = 15000),
                      cores = 3)
loo_m16 = loo(loglik,r_eff = r_eff,cores = 3)
m = "M17";
load(file.path(resultFolder,paste("params_HtSSM_FIT_",m,"_Dx-",gg,"_","Cond-",cc,"_rawRatings_v6B.RData",sep="")) )
code.samples<- as.mcmc.list(results)
mcmc.mat <- as.matrix(code.samples, chains = F)
dic = round(mean(mcmc.mat[,"deviance"])  + .5 * var(mcmc.mat[,"deviance"]))
N = nrow(Data_partial)
chain <- rbind(results$mcmc[[1]], results$mcmc[[2]], results$mcmc[[3]])
loglik <- chain[,paste0("log_lik[",1:N,"]")]
waic_m17 = waic(loglik)
r_eff =  relative_eff(exp(loglik), chain_id = rep(1:3, each = 15000),
                      cores = 3)
loo_m17 = loo(loglik,r_eff = r_eff,cores = 3)
m = "M18";
load(file.path(resultFolder,paste("params_HtSSM_FIT_",m,"_Dx-",gg,"_","Cond-",cc,"_rawRatings_v6B.RData",sep="")) )
code.samples<- as.mcmc.list(results)
mcmc.mat <- as.matrix(code.samples, chains = F)
dic = round(mean(mcmc.mat[,"deviance"])  + .5 * var(mcmc.mat[,"deviance"]))
N = nrow(Data_partial)
chain <- rbind(results$mcmc[[1]], results$mcmc[[2]], results$mcmc[[3]])
loglik <- chain[,paste0("log_lik[",1:N,"]")]
waic_m18 = waic(loglik)
r_eff =  relative_eff(exp(loglik), chain_id = rep(1:3, each = 15000),
                      cores = 3)
loo_m18 = loo(loglik,r_eff = r_eff,cores = 3)
m = "M19";
load(file.path(resultFolder,paste("params_HtSSM_FIT_",m,"_Dx-",gg,"_","Cond-",cc,"_rawRatings_v6B.RData",sep="")) )
code.samples<- as.mcmc.list(results)
mcmc.mat <- as.matrix(code.samples, chains = F)
dic = round(mean(mcmc.mat[,"deviance"])  + .5 * var(mcmc.mat[,"deviance"]))
N = nrow(Data_partial)
chain <- rbind(results$mcmc[[1]], results$mcmc[[2]], results$mcmc[[3]])
loglik <- chain[,paste0("log_lik[",1:N,"]")]
waic_m19 = waic(loglik)
r_eff =  relative_eff(exp(loglik), chain_id = rep(1:3, each = 15000),
                      cores = 3)
loo_m19 = loo(loglik,r_eff = r_eff,cores = 3)
m = "M20";
load(file.path(resultFolder,paste("params_HtSSM_FIT_",m,"_Dx-",gg,"_","Cond-",cc,"_rawRatings_v6B.RData",sep="")) )
code.samples<- as.mcmc.list(results)
mcmc.mat <- as.matrix(code.samples, chains = F)
dic = round(mean(mcmc.mat[,"deviance"])  + .5 * var(mcmc.mat[,"deviance"]))
N = nrow(Data_partial)
chain <- rbind(results$mcmc[[1]], results$mcmc[[2]], results$mcmc[[3]])
loglik <- chain[,paste0("log_lik[",1:N,"]")]
waic_m20 = waic(loglik)
r_eff =  relative_eff(exp(loglik), chain_id = rep(1:3, each = 15000),
                      cores = 3)
loo_m20 = loo(loglik,r_eff = r_eff,cores = 3)
m = "M21";
load(file.path(resultFolder,paste("params_HtSSM_FIT_",m,"_Dx-",gg,"_","Cond-",cc,"_rawRatings_v6B.RData",sep="")) )
code.samples<- as.mcmc.list(results)
mcmc.mat <- as.matrix(code.samples, chains = F)
dic = round(mean(mcmc.mat[,"deviance"])  + .5 * var(mcmc.mat[,"deviance"]))
N = nrow(Data_partial)
chain <- rbind(results$mcmc[[1]], results$mcmc[[2]], results$mcmc[[3]])
loglik <- chain[,paste0("log_lik[",1:N,"]")]
waic_m21 = waic(loglik)
r_eff =  relative_eff(exp(loglik), chain_id = rep(1:3, each = 15000),
                      cores = 3)
loo_m21 = loo(loglik,r_eff = r_eff,cores = 3)
m = "M22";
load(file.path(resultFolder,paste("params_HtSSM_FIT_",m,"_Dx-",gg,"_","Cond-",cc,"_rawRatings_v6B.RData",sep="")) )
code.samples<- as.mcmc.list(results)
mcmc.mat <- as.matrix(code.samples, chains = F)
dic = round(mean(mcmc.mat[,"deviance"])  + .5 * var(mcmc.mat[,"deviance"]))
N = nrow(Data_partial)
chain <- rbind(results$mcmc[[1]], results$mcmc[[2]], results$mcmc[[3]])
loglik <- chain[,paste0("log_lik[",1:N,"]")]
waic_m22 = waic(loglik)
r_eff =  relative_eff(exp(loglik), chain_id = rep(1:3, each = 15000),
                      cores = 3)
loo_m22 = loo(loglik,r_eff = r_eff,cores = 3)

hc_neg_waic = loo_compare(x = list(waic_m16,waic_m17,waic_m18,
                                   waic_m19,waic_m20,waic_m21,waic_m22))
# M20 wins
hc_neg_loo = loo_compare(x = list(loo_m16,loo_m17,loo_m18,
                                  loo_m19,loo_m20,loo_m21,loo_m22))


hc_neg_waic_s = loo_compare(x = list(waic_m17,waic_m19,waic_m18))
hc_neg_loo_s = loo_compare(x = list(loo_m17,loo_m19,loo_m18))


gg = "hc";
cc = "Neutral";

m = "M16";
load(file.path(resultFolder,paste("params_HtSSM_FIT_",m,"_Dx-",gg,"_","Cond-",cc,"_rawRatings_v6B.RData",sep="")) )
code.samples<- as.mcmc.list(results)
mcmc.mat <- as.matrix(code.samples, chains = F)
dic = round(mean(mcmc.mat[,"deviance"])  + .5 * var(mcmc.mat[,"deviance"]))
N = nrow(Data_partial)
chain <- rbind(results$mcmc[[1]], results$mcmc[[2]], results$mcmc[[3]])
loglik <- chain[,paste0("log_lik[",1:N,"]")]
waic_m16 = waic(loglik)
r_eff =  relative_eff(exp(loglik), chain_id = rep(1:3, each = 15000),
                      cores = 3)
loo_m16 = loo(loglik,r_eff = r_eff,cores = 3)
m = "M17";
load(file.path(resultFolder,paste("params_HtSSM_FIT_",m,"_Dx-",gg,"_","Cond-",cc,"_rawRatings_v6B.RData",sep="")) )
code.samples<- as.mcmc.list(results)
mcmc.mat <- as.matrix(code.samples, chains = F)
dic = round(mean(mcmc.mat[,"deviance"])  + .5 * var(mcmc.mat[,"deviance"]))
N = nrow(Data_partial)
chain <- rbind(results$mcmc[[1]], results$mcmc[[2]], results$mcmc[[3]])
loglik <- chain[,paste0("log_lik[",1:N,"]")]
waic_m17 = waic(loglik)
r_eff =  relative_eff(exp(loglik), chain_id = rep(1:3, each = 15000),
                      cores = 3)
loo_m17 = loo(loglik,r_eff = r_eff,cores = 3)
m = "M18";
load(file.path(resultFolder,paste("params_HtSSM_FIT_",m,"_Dx-",gg,"_","Cond-",cc,"_rawRatings_v6B.RData",sep="")) )
code.samples<- as.mcmc.list(results)
mcmc.mat <- as.matrix(code.samples, chains = F)
dic = round(mean(mcmc.mat[,"deviance"])  + .5 * var(mcmc.mat[,"deviance"]))
N = nrow(Data_partial)
chain <- rbind(results$mcmc[[1]], results$mcmc[[2]], results$mcmc[[3]])
loglik <- chain[,paste0("log_lik[",1:N,"]")]
waic_m18 = waic(loglik)
r_eff =  relative_eff(exp(loglik), chain_id = rep(1:3, each = 15000),
                      cores = 3)
loo_m18 = loo(loglik,r_eff = r_eff,cores = 3)
m = "M19";
load(file.path(resultFolder,paste("params_HtSSM_FIT_",m,"_Dx-",gg,"_","Cond-",cc,"_rawRatings_v6B.RData",sep="")) )
code.samples<- as.mcmc.list(results)
mcmc.mat <- as.matrix(code.samples, chains = F)
dic = round(mean(mcmc.mat[,"deviance"])  + .5 * var(mcmc.mat[,"deviance"]))
N = nrow(Data_partial)
chain <- rbind(results$mcmc[[1]], results$mcmc[[2]], results$mcmc[[3]])
loglik <- chain[,paste0("log_lik[",1:N,"]")]
waic_m19 = waic(loglik)
r_eff =  relative_eff(exp(loglik), chain_id = rep(1:3, each = 15000),
                      cores = 3)
loo_m19 = loo(loglik,r_eff = r_eff,cores = 3)
m = "M20";
load(file.path(resultFolder,paste("params_HtSSM_FIT_",m,"_Dx-",gg,"_","Cond-",cc,"_rawRatings_v6B.RData",sep="")) )
code.samples<- as.mcmc.list(results)
mcmc.mat <- as.matrix(code.samples, chains = F)
dic = round(mean(mcmc.mat[,"deviance"])  + .5 * var(mcmc.mat[,"deviance"]))
N = nrow(Data_partial)
chain <- rbind(results$mcmc[[1]], results$mcmc[[2]], results$mcmc[[3]])
loglik <- chain[,paste0("log_lik[",1:N,"]")]
waic_m20 = waic(loglik)
r_eff =  relative_eff(exp(loglik), chain_id = rep(1:3, each = 15000),
                      cores = 3)
loo_m20 = loo(loglik,r_eff = r_eff,cores = 3)
m = "M21";
load(file.path(resultFolder,paste("params_HtSSM_FIT_",m,"_Dx-",gg,"_","Cond-",cc,"_rawRatings_v6B.RData",sep="")) )
code.samples<- as.mcmc.list(results)
mcmc.mat <- as.matrix(code.samples, chains = F)
dic = round(mean(mcmc.mat[,"deviance"])  + .5 * var(mcmc.mat[,"deviance"]))
N = nrow(Data_partial)
chain <- rbind(results$mcmc[[1]], results$mcmc[[2]], results$mcmc[[3]])
loglik <- chain[,paste0("log_lik[",1:N,"]")]
waic_m21 = waic(loglik)
r_eff =  relative_eff(exp(loglik), chain_id = rep(1:3, each = 15000),
                      cores = 3)
loo_m21 = loo(loglik,r_eff = r_eff,cores = 3)
m = "M22";
load(file.path(resultFolder,paste("params_HtSSM_FIT_",m,"_Dx-",gg,"_","Cond-",cc,"_rawRatings_v6B.RData",sep="")) )
code.samples<- as.mcmc.list(results)
mcmc.mat <- as.matrix(code.samples, chains = F)
dic = round(mean(mcmc.mat[,"deviance"])  + .5 * var(mcmc.mat[,"deviance"]))
N = nrow(Data_partial)
chain <- rbind(results$mcmc[[1]], results$mcmc[[2]], results$mcmc[[3]])
loglik <- chain[,paste0("log_lik[",1:N,"]")]
waic_m22 = waic(loglik)
r_eff =  relative_eff(exp(loglik), chain_id = rep(1:3, each = 15000),
                      cores = 3)
loo_m22 = loo(loglik,r_eff = r_eff,cores = 3)

hc_neu_waic = loo_compare(x = list(waic_m16,waic_m17,waic_m18,
                                   waic_m19,waic_m20,waic_m21,waic_m22))
# M18 wins
hc_neu_loo = loo_compare(x = list(loo_m16,loo_m17,loo_m18,
                                  loo_m19,loo_m20,loo_m21,loo_m22))
#m20

hc_neu_waic_s = loo_compare(x = list(waic_m17,waic_m19,waic_m18))
hc_neu_loo_s = loo_compare(x = list(loo_m17,loo_m19,loo_m18))


hc_neu_loo # M21 (weights, NDT, Time)
hc_neg_loo # M20 (weights, ndt)
bn_neu_loo # M18 (Time, b1, b2)
bn_neg_loo # M20 (weights, ndt)

hc_neu_waic # M18 (time, weights)
hc_neg_waic # M20 (weights, ndt)
bn_neu_waic # M18 (time, weights)
bn_neg_waic # M20 (weights, ndt)



# Just models with time, b1/b2
hc_neu_loo_s # M18 (weights)
hc_neg_loo_s # M19 (weights)
bn_neu_loo_s # M18 (time)
bn_neg_loo_s # M19 (weights)


hc_neu_waic_s
hc_neg_waic_s
bn_neu_waic_s
bn_neg_waic_s


# What we want to know is how food x condition affects time, weights
# So let's just see if NDT helps

