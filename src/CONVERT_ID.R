

# Script for de-identifying participants

library(tidyverse)
library(readxl)
datPath <- '~/Research/Sinai-Food-Choice/data'
datstDDM <- "~/Research/Sinai-Food-Choice/stDDM/Data/"

datSave <- "~/Research/Sinai-Food-Choice/Computations_BN_Food_Choice/data/"

# Data used to fit models
fit_beh <- read.csv(file.path(datstDDM, "fc_data_JAGS_022824.csv"))
# Filter data so it only includes participants fit with model
good_ids <- unique(fit_beh$subject)

# Data with Self-report measures
beh <- read_xlsx(file.path(datPath,"bn_fct_23.5.18.xlsx"))

mh <- beh %>%
  filter(Subject %in% good_ids)  %>%
  # Add POMS Change
  mutate(POMS_Neg = ( `PostPOMSTotal NA_Neg` - `PrePOMSTotal NA_Neg`),
         POMS_Neu = ( `PostPOMSTotal NA_Neu` - `PrePOMSTotal NA_Neu` ),
         POMS_DEP_Neg = ( `PostPOMSDepression_Neg` - `PrePOMSDepression_Neg`),
         POMS_DEP_Neu = ( `PostPOMSDepression_Neu` -  `PrePOMSDepression_Neu` ),
         subject = Subject) %>%
  select(!Subject)


# Behavioral figures
beh.df = fit_beh %>%
  mutate(cond = factor(cond, levels = c("Neutral","Negative")),
         Dx = factor(Dx, levels = c("hc","bn"),
                     labels = c("HC","BN")),
         food = factor(foodType, levels = c("lf","hf"),
                       labels = c("Low fat", "High fat")),
         taste_z = tdz_group,
         health_z = hdz_group) %>%
  select(subject, Dx, cond, foodType, choice, td, hd, rt, taste_z, health_z)

# Just subj id and Dx
subj.df = beh.df %>% dplyr::select(subject,Dx) %>% distinct() %>%
  mutate(idx = paste0("SINAI_", order(subject)))

new_beh = merge(beh.df, subj.df) %>%
  select(!subject)

new_mh = merge(subj.df, mh, by = "subject") %>%
  select(!subject)

write.csv(new_beh, file = paste0(datSave,"deidentified_ChoiceData.csv"))
write.csv(new_mh, file = paste0(datSave,"deidentified_SelfReportData.csv"))




