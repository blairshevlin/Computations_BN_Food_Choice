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
# 06/18/24      Blair Shevlin                         Wrote code for original manuscript


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
resFolder = path(base_path) / "results" / "stDDM" / "estimation"

Data<-read.csv(file=paste0(dataFolder, "/deidentified_ChoiceData.csv"))

# M2 (Old M17): food-type -> Time

# M3 (Old M18): food-type -> Time, b1, b2

# M1 (Old M19): food-type -> b1, b2

# M0 (Old M31): Null

# Fit initial
for (m in c("M0","M1","M2","M3") ) {
  for (gg in c("BN","HC")){
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
for (gg in c("HC","BN")) {
  for (cc in c("Negative","Neutral")) {
    for (m in models_of_interest) {
      load(file.path(resFolder,paste("params_HtSSM_FIT_",m,"_Dx-",gg,"_","Cond-",cc,"_rawRatings.RData",sep="")) )
      code.samples<- as.mcmc.list(results)
      mcmc.mat <- as.matrix(code.samples, chains = F)
      rm(code.samples)
      N = nrow(Data_partial)
      chain <- rbind(results$mcmc[[1]], results$mcmc[[2]], results$mcmc[[3]])
      rm(results,mcmc.mat)
      loglik <- chain[,paste0("log_lik[",1:N,"]")]
      rm(chain)
      w = waic(loglik)
      tmp.df = data.frame(Dx = gg,
                          cond = cc,
                          model = m,
                           waic = w$estimates[3]
      )
      dic.df = rbind(dic.df,tmp.df)
    }
  }
}

# Report WAIC scores
dic.df %>%
  pivot_longer(cols = c(dic,waic,loo)) %>%
  filter(name == "waic") %>%
  group_by(name,Dx,model) %>%
  summarise(val = round(sum(value)))%>%
  mutate(
         Dx = factor(Dx, levels = c("HC","BN"))
        ) %>%
  arrange(model,Dx) %>%
  as.data.frame()