# manuscript_analysis.R - Code for reproducing results from "Negative affect influences the computations underlying food choice in bulimia nervosa" 
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
# 06/17/24      Blair Shevlin                         Behavioral code for original manuscript
# 06/18/24      Blair Shevlin                         Modeling code for original manuscript
# 06/20/24      Blair Shevlin                         Symptom severity code for original manuscript
# 03/12/25      Blair Shevlin                         Affect change code for response to reviewers
# 03/18/25      Blair Shevlin                         Simple effects analysis for response to reviewers
# 03/21/25      Blair Shevlin                         Simple effects analyses of affect change
# 03/23/25      Blair Shevlin                         Alt. symptom severity for reviewers
# 05/27/25      Blair Shevlin                         Final edits for resubmission
# 11/07/25      Blair Shevlin                         Assessing marginal means for affect change analyses
# 11/13/25      Blair Shevlin                         Correlations between restriction and binge frequency 

# Packages required
required_packages <- c(
  "tidyverse",
  "sjPlot",
  "here",
  "fs",
  "broom",
  "MASS",
  "patchwork",
  "ggpubr",
  "cetcolor",
  "HDInterval",
  "coda",
  "doBy",
  "lmerTest",
  "faux",
  "ggeffects",
  "lsmeans",
  "glmmTMB")

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

# Standard error function
se = function(x) sd(x) / sqrt(length(x))

# Folders
datPath <- path(here()) / 'data' 
resPath <- path(here()) / 'results'
figPath <- path(here()) / 'results' / 'figures'

# Load behavioral data
beh.df <- read.csv(file = datPath / "deidentified_ChoiceData.csv") %>%
  mutate(cond = factor(cond, levels = c("Neutral","Negative")),
         food = factor(foodType, levels = c("Low Fat","High Fat"), labels = c("Low fat","High fat")),
         Dx = factor(Dx, levels = c("HC","BN")))

# Contrasts
contrasts(beh.df$cond) <- c(-1,1)
contrasts(beh.df$food) <- c(-1,1)
contrasts(beh.df$Dx) <- c(-1,1)

# Model 1: choice ~ food_type x affect x group
glm.1 <- glmer(data = beh.df,
               formula = choice ~
                 Dx * cond * food +
                 (1 + cond * food | idx),
               family=binomial,
               control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=20000)))
summary(glm.1)
# Table S.9 for supplements

# Marginal Means
lsmeans(glm.1, specs = ~ Dx , type = "response")
### Note --- subtract 1 from these probs to get reference item selection
lsmeans(glm.1, specs = ~ Dx * food , type = "response")
lsmeans(glm.1, specs = ~ Dx * food * cond, type = "response")

emm <- emmeans(glm.1, ~ Dx | food * cond)
pairs(emm)

# Model 2: affect x group + group x health + group x taste
glm.2 <- glmer(data = beh.df,
               formula = choice ~ Dx * cond +
                 taste_z + health_z + 
                 taste_z * Dx  + health_z * Dx  +
                 taste_z * cond  + health_z * cond  +
                 taste_z * Dx * cond + 
                 health_z * Dx * cond +
                 (1 + taste_z * cond + health_z * cond| idx),
               family=binomial,
               control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=20000)))
summary(glm.2)
# Table S10 for supplements

# Model 3: Self-control trials
data.sc <- beh.df %>%
  mutate(liked_item = ifelse(td>=5,1,0),
         disliked_item = ifelse(td <= 1, 1, 0),
         unhealthy_item = ifelse(hd <= 1,1,0),
         healthy_item = ifelse(hd >= 5, 1, 0),
         item_type = ifelse(liked_item == 1,
                            ifelse(healthy_item==1,"Liked Healthy", "Liked Unhealthy"),
                            ifelse(healthy_item==1,"Disliked Healthy", "Disliked Unhealthy"))) %>%
  # Only use trials with liked, unhealthy items or disliked, health items
  filter(item_type %in% c("Liked Unhealthy", "Disliked Healthy")) %>%
  mutate(sc = ifelse( (item_type == "Liked Unhealthy" & choice == 0) | (item_type == "Disliked Healthy" & choice == 1) ,1,0))

glm.3 <- glmer(data = data.sc,
               formula = sc ~
                 Dx * cond +
               (1 + cond|idx),
               family=binomial,
               control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=20000)))
summary(glm.3)
# Supplementary Table S11

# Model 4: Response times
beh.df$choice_c = factor(beh.df$choice, levels = c(0,1), labels = c("Reference item","Presented item"))
contrasts(beh.df$choice_c) <- c(-1,1)

lm.1 <- lmer(data = beh.df,
               formula = log(rt) ~   Dx * cond * food * choice_c +
                 (1 + cond * food | idx),
               REML = F,
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=20000)))
summary(lm.1)
# Supplementary Table S12


# Panels for Figure 2
choice.pred <- ggpredict(glm.1,terms = c("Dx","cond","food"))

fig2a <- 
  ggplot(choice.pred, aes(x = facet, y = predicted, 
                          color = x, fill = x)) +
  facet_wrap(.~group)+
  geom_hline(yintercept = 0.5,
             linetype = "dashed",
             color = "gray",
             linewidth = 1) +
  scale_colour_viridis_d(begin = 0,
                         end = .8,
                         direction = -1
                         ) +
  scale_fill_viridis_d(begin = 0,
                       end = .8,
                       direction = -1,
                       guide="none") +
  scale_linetype_manual(values = c("solid", "dashed")) +
  geom_bar(stat = "identity",alpha=.5,aes(linetype=x),
           position=position_dodge(width = 1)) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high), width=.2,
                position=position_dodge(1))+
  theme_pubr(base_size=18) +
  labs(color = element_blank(), #"Group",
       linetype = element_blank(), #"Group",
       fill = element_blank(), #"Group",
       y="Choice proportion\n(over reference item)",
       x = "Food type",
       title=element_blank()) +
    guides(linetype = guide_legend(override.aes = list(fill = c("#5ec962","#440154")))) +
    theme(legend.position=c(.45,.95))

taste_coeff_HC = coef(glm.2)$idx["taste_z"] %>% mutate(val = taste_z ,Dx = "HC", attribute = "taste") %>% dplyr::select(!c(taste_z))
taste_coeff_BN = (coef(glm.2)$idx["taste_z"] + coef(glm.2)$idx["Dx1:taste_z"]) %>% mutate(val = taste_z ,Dx = "BN", attribute = "taste") %>% dplyr::select(!c(taste_z))
health_coeff_HC = coef(glm.2)$idx["health_z"] %>% mutate(val = health_z ,Dx = "HC", attribute = "health") %>% dplyr::select(!c(health_z))
health_coeff_BN = (coef(glm.2)$idx["health_z"] + coef(glm.2)$idx["Dx1:health_z"]) %>% mutate(val = health_z ,Dx = "BN", attribute = "health") %>% dplyr::select(!c(health_z))

hr.beta = data.frame(Dx = c("HC","BN","HC","BN"),
                     attr = c("Healthiness","Healthiness","Tastiness","Tastiness"),
                     m = c(mean(health_coeff_HC$val),mean(health_coeff_BN$val),mean(taste_coeff_HC$val),mean(taste_coeff_BN$val)),
                     s = c(se(health_coeff_HC$val),se(health_coeff_BN$val),se(taste_coeff_HC$val),se(taste_coeff_BN$val))
) %>%
mutate(Dx = factor(Dx, levels = c("HC","BN")))


fig2b <- 
ggplot(hr.beta, aes(x = attr , y = m, color = Dx, fill = Dx, group = Dx)) +
  scale_colour_viridis_d(begin = 0,
                         end = .8,
                         direction = -1
  ) +
  scale_fill_viridis_d(begin = 0,
                       end = .8,
                       direction = -1,
                       guide="none") +
  scale_linetype_manual(values = c("solid", "dashed")) +
  geom_bar(stat = "identity",aes(linetype = Dx),position = position_dodge(width=1),alpha=0.5) +
  geom_errorbar(aes(ymin = m - s, ymax = m + s),linewidth = 1,width=.2,
    position = position_dodge(width=1)) +
  theme_pubr(base_size=18) +
  labs(color = element_blank(), #"Group",
       linetype = element_blank(), #"Group",
       fill = element_blank(), #"Group",
       y="Decision Bias",
       x =  element_blank(),
       title=element_blank()) +
  guides(fill = "none",
         color = "none",
         linetype = "none") 

figure2 = fig2a + fig2b  +
  plot_annotation(tag_levels = 'A')

ggsave(file = figPath / "figure2_final.tiff", plot = figure2, width = 12, height = 8)

##################################
# Computational modeling results #
##################################

# Load subject-level parameters
df.fit.full <- NULL

for (cc in c("Neutral","Negative")) {
  for (gg in c("BN","HC")) {
    if (gg == "BN" & cc == "Negative"){
      # Load in file that required re-estimation to ensure convergence
      load(file.path(resPath,paste("/stDDM/estimation/params_HtSSM_FIT_M3_Dx-",gg,"_","Cond-",cc,"_rawRatings_converted_rest.RData",sep="")) )
    } else {
        load(file.path(resPath,paste("/stDDM/estimation/params_HtSSM_FIT_M3_Dx-",gg,"_","Cond-",cc,"_rawRatings_converted.RData",sep="")) )
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
params <- df.fit.full %>%
  mutate(foodType = ifelse(grepl("_lf",params),"Low-Fat",
                           ifelse(grepl("_hf",params),"High-Fat", "NA")),
         params = recode(params,
                         "wt_lf" = "wt",
                         "wt_hf" = "wt",
                         "wh_lf" = "wh",
                         "wh_hf" = "wh",
                         "tHin_lf" = "tHin",
                         "tHin_hf" = "tHin")
  )

# Attribute timing
tHin.df <- params %>%
  filter(params == "tHin") %>%
  mutate( foodType = factor(foodType,
                            levels=c("Low-Fat","High-Fat"),
                            labels = c("Low-Fat","High-Fat")),
          cond = factor(cond,levels=c("Neutral","Negative")),
          Dx = factor(Dx,levels=c("HC","BN"),labels=c("Healthy Controls",
                                                      "Bulimia Nervosa")))                                                   
tHin.lm <- lmer (data=tHin.df,
                 formula = vals ~ Dx * cond * foodType + (1|idx),
                 REML=F,
                 control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=20000)))
summary(tHin.lm) 
# Supplementary Table S2

# Laura: If we decide not to use simple effects: we are just going to stay at the highest level of interaction.
# Then we just interpret three-way and figure. Then  the simple effects are not needed for any of these.

# Alternative: leave in simple effects, what makes sense to report?
# Hold hand to interpret each. Refer to plot: opposite patterns in neutral condition. 
# Then after


tHin.lm.ne <- lmer (data=tHin.df[tHin.df$cond == "Neutral",],
                 formula = vals ~ foodType + Dx + (1|idx),
                 REML=F,
                 control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=20000)))
summary(tHin.lm.ne) 

tHin.lm.na <- lmer (data=tHin.df[tHin.df$cond == "Negative",],
                 formula = vals ~ foodType + Dx + (1|idx),
                 REML=F,
                 control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=20000)))
summary(tHin.lm.na) 

tHin.lm.noft <- lmer (data=tHin.df,
                 formula = vals ~ Dx * cond + (1|idx),
                 REML=F,
                 control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=20000)))
summary(tHin.lm.noft) 

tHin.lm.twoways <- lmer (data=tHin.df,
                 formula = vals ~ foodType * cond + Dx * cond + (1|idx),
                 REML=F,
                 control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=20000)))
summary(tHin.lm.twoways) 

# 1. Simple effects analysis - Group differences within each condition/food type combination
emmeans(tHin.lm, pairwise ~ Dx | cond * foodType)
emmeans(tHin.lm, pairwise ~ Dx | cond)


# 2. Simple slopes - How each group changes from neutral to negative within each food type
emmeans(tHin.lm, pairwise ~ cond | Dx * foodType)

# 3. Difference in differences - How the Group x Food_Type interaction changes across affect conditions
emmeans(tHin.lm, pairwise ~ foodType | Dx * cond)

# 4. Three-way interaction contrasts
emmeans(tHin.lm, pairwise ~ Dx * cond * foodType)

# For Group by Food Type interaction for each condition
emm_ft_cond <- emmeans(tHin.lm, ~ Dx * foodType | cond)
pairs(emm_ft_cond, by = "Dx", adjust ="none")        # Looking at food type effect for each group in Neutral
# Confirmed that the difference between LF-HF in Neutral condition is significant for BN but not HC group

emm_dx_cond <- emmeans(tHin.lm, ~ Dx * cond | foodType)
pairs(emm_dx_cond, by = "foodType", adjust ="none")    

emm_dx_cond.collapse = emmeans(tHin.lm, ~ Dx * cond)
pairs(emm_dx_cond.collapse, adjust ="none")   

# Want difference in difference


# Analyze the effect separately for each food type
emm_by_foodtype <- emmeans(tHin.lm, ~ Dx * cond | foodType)

foodType_contrasts <- contrast(emm_by_foodtype, 
                              interaction = "pairwise", 
                              by = "foodType")
print(foodType_contrasts)
  # Showing the difference between HC and BN across conditions occured for LF but not HF foods
# Supplementary Table S3

# Laura's recommendation - effect of condition in each group (collapsing food type)
emmeans(tHin.lm,  pairwise ~ cond * Dx )


emm <- emmeans(tHin.lm, ~ Dx * cond)

# Define the difference-in-differences contrast
# (Dx1 conda - Dx1 condab) - (Dx2 conda - Dx2 condab)
contrast_list <- list(
  "DinD" = c(1, -1, -1, 1)  # Assuming factor level order: Dx1:conda, Dx1:condab, Dx2:conda, Dx2:condab
)

# Test the contrast
contrast_result <- contrast(emm, contrast_list)
summary(contrast_result)

# Attribute Weights

## Taste
taste.df <- params %>%
  filter(params == "wt") %>%
  mutate( foodType = factor(foodType,
                            levels=c("Low-Fat","High-Fat"),
                            labels = c("Low-Fat","High-Fat")),
          cond = factor(cond,levels=c("Neutral","Negative")),
          Dx = factor(Dx,levels=c("HC","BN"),labels=c("Healthy Controls",
                                                      "Bulimia Nervosa")))
taste.lm1 <- lmer (data=taste.df,
               formula = vals ~ Dx * cond * foodType+ (1|idx),
               REML=F,
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=20000)))
summary(taste.lm1)
# Supplementary Table S4

# Simple effects
emm_taste_group <- emmeans(taste.lm1, ~ Dx | cond)
emm_taste_group_cont = pairs(emm_taste_group)  %>% as.data.frame() # Group differences for each condition
emm_taste_group_cont
emm_taste_group_cont$p.value

emm_taste_cond <- emmeans(taste.lm1, ~ cond | Dx)
emm_taste_cond_cont = pairs(emm_taste_cond)
emm_taste_cond_cont

emm_taste_ft <- emmeans(taste.lm1, ~ foodType | Dx)
emm_taste_group_food_cont = pairs(emm_taste_ft) %>% as.data.frame()
emm_taste_group_food_cont
emm_taste_group_food_cont$p.value

# Supplementary Table S5

## Health
health.df <- params %>%
  filter(params == "wh") %>%
  mutate( foodType = factor(foodType,
                            levels=c("Low-Fat","High-Fat"),
                            labels = c("Low-Fat","High-Fat")),
          cond = factor(cond,levels=c("Neutral","Negative")),
          Dx = factor(Dx,levels=c("HC","BN"),labels=c("Healthy Controls",
                                                      "Bulimia Nervosa")))
health.lm1 <- lmer (data=health.df,
                    formula = vals ~ Dx * cond * foodType + (1|idx),
                    REML=F,
                    control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=20000)))
summary(health.lm1)
# Supplementary Table S6

emm_health_group <- emmeans(health.lm1, ~ Dx | cond)
emm_health_group_cont = pairs(emm_health_group)  %>% as.data.frame() # Group differences for each condition
emm_health_group_cont
emm_taste_group_cont$p.value

emm_health_cond <- emmeans(health.lm1, ~ cond | Dx)
emm_health_cond_cont = pairs(emm_health_cond)
emm_health_cond_cont
# Food type differences for each group
emm_health_group_food <- emmeans(health.lm1, ~ foodType | Dx)
emm_health_group_food_cont = pairs(emm_health_group_food) %>% as.data.frame()
emm_health_group_food_cont
emm_health_group_food_cont$p.value

# Supplementary Table S7

# Figure 4
fig4a <- 
params %>%
  filter(params == "tHin") %>%
  mutate(grouping = ifelse(Dx == "BN" & cond == "Negative",
                           ifelse(foodType == "High-Fat","Negative\nHigh-Fat",
                                  "Negative\nLow-Fat"),cond
  ),
  grouping = factor(grouping,
                    levels = c("Neutral","Negative",
                               "Negative\nLow-Fat",
                               "Negative\nHigh-Fat")),
  cond = factor(cond,levels = c("Neutral","Negative")),
  foodType = factor(foodType,
                    levels=c("NA","Low-Fat","High-Fat"),
                    labels = c("All Foods","Low Fat","High Fat")),
  Dx = factor(Dx,levels=c("HC","BN"),labels=c("HC",
                                              "BN"))
  ) %>%
  ggplot(aes(x = foodType, y = vals, color = Dx, 
             group = Dx)) +
  theme_pubr(base_size = 18) +
  facet_wrap(~factor(cond, levels = c("Neutral", "Negative"))) +  
  geom_hline(yintercept = 0,
             linetype = "dashed") +
  geom_point(position = position_dodge2(width = .75),
             alpha = .5) +
  stat_summary(geom = "line",
               position = position_dodge2(width = .75),
               size = 1.5,alpha=.75,
               linewidth = 1.5) +
  stat_summary(position = position_dodge2(width = .75),
               size = 1.5,
               linewidth = 1.5) +
                # Add text above and below the line for Negative facet only
 geom_text(data = data.frame(cond = "Negative", 
                             x = 1.75, 
                             y = c(0.03, -0.03),
                             label = c("Health sooner", "Taste sooner")),
            aes(x = x, y = y, label = label),
            inherit.aes = FALSE,
            size = 4, color = "black") +
  labs(
       y = "Attribute Onset",
       color = "Group",
      x= element_blank()
      ) +
  coord_cartesian(ylim = c(-0.75, 0.25)) +
  scale_color_brewer(type = "qual") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1, hjust=1)) +
  theme(legend.position = "right", legend.box = "vertical")  +
  theme(legend.position=c(.9,.6))

fig4b <- params %>%
  filter(params == "wt") %>%
  mutate(grouping = ifelse(Dx == "BN" & cond == "Negative",
                           ifelse(foodType == "High-Fat","Negative\nHigh-Fat",
                                  "Negative\nLow-Fat"),cond
  ),
  grouping = factor(grouping,
                    levels = c("Neutral","Negative",
                               "Negative\nLow-Fat",
                               "Negative\nHigh-Fat")),
  cond = factor(cond,levels = c("Neutral","Negative")),
  foodType = factor(foodType,
                    levels=c("NA","Low-Fat","High-Fat"),
                    labels = c("All Foods","Low Fat","High Fat")),
  Dx = factor(Dx,levels=c("HC","BN"),labels=c("HC",
                                              "BN"))
  ) %>%
  ggplot(aes(x = foodType, y = vals, color = Dx, 
             group = Dx)) +
  theme_pubr(base_size = 18) +
  facet_wrap(~cond) +
  geom_hline(yintercept = 0,
             linetype = "dashed") +
  geom_point(position = position_dodge2(width = .75),
             alpha = .5) +
  stat_summary(geom = "line",
               position = position_dodge2(width = .75),
               size = 1.5,alpha=.75,
               linewidth = 1.5) +
  stat_summary(position = position_dodge2(width = .75),
               size = 1.5,
               linewidth = 1.5) +
  labs(
    y = "Tastiness Weight",
    color = "Group",
    x = element_blank()
    ) +
  scale_color_brewer(type = "qual") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1, hjust=1)) +
  theme(legend.position = "right", legend.box = "vertical")  +
  theme(legend.position=c(.9,.6))

fig4c = params %>%
  filter(params == "wh") %>%
  mutate(grouping = ifelse(Dx == "BN" & cond == "Negative",
                           ifelse(foodType == "High-Fat","Negative\nHigh-Fat",
                                  "Negative\nLow-Fat"),cond
  ),
  grouping = factor(grouping,
                    levels = c("Neutral","Negative",
                               "Negative\nLow-Fat",
                               "Negative\nHigh-Fat")),
  cond = factor(cond,levels = c("Neutral","Negative")),
  foodType = factor(foodType,
                    levels=c("NA","Low-Fat","High-Fat"),
                    labels = c("All Foods","Low Fat","High Fat")),
  Dx = factor(Dx,levels=c("HC","BN"),labels=c("HC",
                                              "BN"))
  ) %>%
  ggplot(aes(x = foodType, y = vals, color = Dx, 
             group = Dx)) +
  theme_pubr(base_size = 18) +
  facet_wrap(~cond) +
  geom_hline(yintercept = 0,
             linetype = "dashed") +
  geom_point(position = position_dodge2(width = .75),
             alpha = .5) +
  stat_summary(geom = "line",
               position = position_dodge2(width = .75),
               size = 1.5,alpha=.75,
               linewidth = 1.5) +
  stat_summary(position = position_dodge2(width = .75),
               size = 1.5,
               linewidth = 1.5) +
  labs(
    y = "Healthiness Weight",
    color = "Group",
    x = element_blank()
    ) +
  scale_color_brewer(type = "qual") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1, hjust=1)) +
  theme(legend.position = "right", legend.box = "vertical")  +
  theme(legend.position=c(.9,.6))

figure4 = fig4a + fig4b + fig4c + 
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'A') & theme(legend.position = 'bottom')

ggsave(file = figPath / "figure4_final.tiff", plot = figure4, width = 12, height = 8)

###################
# Symtom Severity #
###################
params = params %>%
  mutate(Dx = factor(Dx, levels = c("HC","BN"), labels = c("HC","BN")),)

mh = read.csv(file=datPath / "deidentified_SelfReportData.csv")

params_mh <- merge(params,mh) 

LOC <- params_mh %>%
  filter(Dx == "BN") %>%
  dplyr::select(idx,cond,foodType,params,vals,
         EDE.OBE.Month.1..episodes.,
         EDE.SBE.Month.1..episodes.,
         EDE.OBE.Month.2..episodes.,
         EDE.SBE.Month.2..episodes.,
         EDE.OBE.Month.3..episodes.,
         EDE.SBE.Month.3..episodes.
         ) %>%
  mutate(OBE_M1 = EDE.OBE.Month.1..episodes.,
         OBE_M2 = EDE.OBE.Month.2..episodes.,
         OBE_M3 = EDE.OBE.Month.3..episodes.,
         OBE_SUM = OBE_M1 + OBE_M2 + OBE_M3,
         SBE_M1 = EDE.SBE.Month.1..episodes.,
         SBE_M2 = EDE.SBE.Month.2..episodes.,
         SBE_M3 = EDE.SBE.Month.3..episodes.,
         SBE_SUM = SBE_M1 + SBE_M2 + SBE_M3,
         LOC_M1 = OBE_M1 + SBE_M1,
         LOC_SUM = OBE_SUM + SBE_SUM
  ) %>%
  group_by(idx,params) %>%
  mutate(params = recode(params,
                        boundary = 'alpha',
                        wt = 'omega[taste]',
                        wh = 'omega[health]',
                        bias = 'z',
                        nDT = 'tau[DT]',
                        tHin = 'tau[s]'),
         foodType = factor(foodType,
                           levels = c("Low-Fat","High-Fat")),
         cond = factor(cond,
                       levels = c("Neutral","Negative"))
  ) %>%
  dplyr::select(!c(EDE.OBE.Month.1..episodes.,
                   EDE.OBE.Month.2..episodes.,
                   EDE.OBE.Month.3..episodes.,
                   EDE.SBE.Month.1..episodes.,
                   EDE.SBE.Month.2..episodes.,
                   EDE.SBE.Month.3..episodes.
                     )) 

# wt, wh, tHin, ndt
LOC.df  = LOC %>% 
  filter(foodType!="NA") %>%
  pivot_longer(cols = c("OBE_SUM",
                        "SBE_SUM",
                        "LOC_SUM"),
               names_to = "LOC",
               values_to = "SUM") %>% distinct()


fig5 = ggplot(LOC.df[LOC.df$params == "tau[s]" & LOC.df$LOC == "SBE_SUM",], 
       aes(x = vals, y = SUM, color = cond, fill = cond,
           linetype = foodType, shape = foodType)) +
  theme_pubr(base_size = 18) +
  geom_vline(data=filter(LOC.df, cond=="Neutral"), aes(xintercept=0), 
             linetype = "dashed",
             linewidth = 1.5,
             colour="gray")  +
  geom_point(alpha = .5,size=3) + 
  geom_smooth(aes(linetype = foodType),
              alpha = 0.5,
              size = 3,
              method = "glm.nb",
              se = F) +
  facet_wrap(~cond, scales = "free_x") + 
  scale_color_brewer(type = "qual",palette = 6,direction = -1) +
  scale_fill_brewer(type = "qual",palette = 6,direction = -1) +
  coord_cartesian(ylim = c(-10,200)) +
  labs(x = "Attribute Onset Estimate",
       y = "Subjective Binge Episodes\n(3 Month Total)",
       shape = "Food Type",
       linetype = "Food Type") +
  theme(legend.position="right",
        legend.key.size =  unit(0.5, "in"),
        axis.text.x = element_text(
                                   angle = 45,
                                   vjust=1, hjust=1),
        panel.spacing.x = unit(10, "mm"))+
  guides(color = "none",
         fill = "none",
         linetype = guide_legend(order=2,
                                 override.aes=list(color = "black")),
         shape = guide_legend(order=3,
                              override.aes=list(size = 5,
                                                color = "black",alpha=1))
  )


taste_text_neg <- data.frame(SUM = -10,vals = -.75,lab = "Taste Earlier",
                             cond = factor( "Negative",c("Neutral","Negative")),
                             foodType = c("Low-Fat","Low-Fat"))
health_text_neg <- data.frame(SUM = -10,vals = -.63,lab = "Health Earlier",
                             cond = factor( "Negative",c("Neutral","Negative")),
                             foodType = c("Low-Fat","Low-Fat"))
taste_text_neu <- data.frame(SUM = -10,vals = -.6,lab = "Taste Earlier",
                       cond = factor( "Neutral",c("Neutral","Negative")),
                       foodType = c("Low-Fat","Low-Fat"))
health_text_neu <- data.frame(SUM = -10,vals = .3,lab = "Health Earlier",
                         cond = factor( "Neutral",c("Neutral","Negative")),
                         foodType = c("Low-Fat","Low-Fat"))
figure5 = fig5 + 
  geom_text(data = taste_text_neg,label = "Taste Earlier",color="black",size=6) +
  geom_text(data = taste_text_neu,label = "Taste Earlier",color="black",size=6) +
  geom_text(data = health_text_neg,label = "Health Earlier",color="black",size=6) +
  geom_text(data = health_text_neu,label = "Health Earlier",color="black",size=6)

ggsave(file = figPath / "figure5_final.tiff", plot = figure5, width = 12, height = 8)

fig5b = ggplot(LOC.df[LOC.df$params == "tau[s]" & LOC.df$LOC == "OBE_SUM",], 
       aes(x = vals, y = SUM, color = cond, fill = cond,
           linetype = foodType, shape = foodType)) +
  theme_pubr(base_size = 18) +
  geom_vline(data=filter(LOC.df, cond=="Neutral"), aes(xintercept=0), 
             linetype = "dashed",
             linewidth = 1.5,
             colour="gray")  +
  geom_point(alpha = .5,size=3) + 
  geom_smooth(aes(linetype = foodType),
              alpha = 0.5,
              size = 3,
              method = "glm.nb",
              se = F) +
  facet_wrap(~cond, scales = "free_x") + 
  scale_color_brewer(type = "qual",palette = 6,direction = -1) +
  scale_fill_brewer(type = "qual",palette = 6,direction = -1) +
  coord_cartesian(ylim = c(-10,200)) +
  labs(x = "Attribute Onset Estimate",
       y = "Subjective Binge Episodes\n(3 Month Total)",
       shape = "Food Type",
       linetype = "Food Type") +
  theme(legend.position="right",
        legend.key.size =  unit(0.5, "in"),
        axis.text.x = element_text(
                                   angle = 45,
                                   vjust=1, hjust=1),
        panel.spacing.x = unit(10, "mm"))+
  guides(color = "none",
         fill = "none",
         linetype = guide_legend(order=2,
                                 override.aes=list(color = "black")),
         shape = guide_legend(order=3,
                              override.aes=list(size = 5,
                                                color = "black",alpha=1))
  )
ggsave(file = figPath / "figure5_obe_NOTUSED.tiff", plot = fig5b, width = 12, height = 8)

LOC.tHin.wide = 
LOC %>%
  filter(params == "tau[s]") %>%
  dplyr::select(idx,cond,foodType,vals,SBE_SUM,OBE_SUM,params) %>%
  pivot_wider(names_from = c(cond,foodType), values_from = vals) %>%
  dplyr::select(idx, SBE_SUM,OBE_SUM,`Neutral_Low-Fat`,`Negative_Low-Fat`,`Neutral_High-Fat`,`Negative_High-Fat`) %>%
  as.data.frame()

sbe.m.full = glmmTMB(SBE_SUM ~ `Neutral_Low-Fat` + `Neutral_High-Fat`  + `Negative_Low-Fat` +`Negative_High-Fat`,
                   data=LOC.tHin.wide,
                   ziformula=~1,
                   family=nbinom1)                          
summary(sbe.m.full) # Supplementary Table S8 (Top)


obe.m.full = glmmTMB(OBE_SUM ~ `Neutral_Low-Fat` + `Neutral_High-Fat`  + `Negative_Low-Fat` +`Negative_High-Fat`,
                   data=LOC.tHin.wide,
                  #  ziformula=~1,
                   family=nbinom1)

# Supplementary Table S8 (Bottom)
summary(obe.m.full)


############################
# Supplementary Parameters #

# Non-decision time
ndt.df <- params %>%
  filter(params == "nDT") %>%
  mutate( cond = factor(cond,levels=c("Neutral","Negative")),
          Dx = factor(Dx,levels=c("HC","BN"),labels=c("Healthy Controls",
                                                      "Bulimia Nervosa")))
ndt.lm <- lmer (data=ndt.df,
                    formula = vals ~ Dx * cond  + (1|idx),
                    REML=F,
                    control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=20000)))
summary(ndt.lm) # Table S23

bound.df <- params %>%
  filter(params == "boundary") %>%
  mutate( cond = factor(cond,levels=c("Neutral","Negative")),
          Dx = factor(Dx,levels=c("HC","BN"),labels=c("Healthy Controls",
                                                      "Bulimia Nervosa")))
bound.lm <- lmer (data=bound.df,
                    formula = vals ~ Dx * cond  + (1|idx),
                    REML=F,
                    control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=20000)))
summary(bound.lm) # Table 24

emm_bound <- emmeans(bound.lm, ~ cond | Dx)

emm_bound_cont = pairs(emm_bound)  %>% as.data.frame() # Group differences for each condition
emm_bound_cont # Table 25

bias.df <- params %>%
  filter(params == "bias") %>%
  mutate( cond = factor(cond,levels=c("Neutral","Negative")),
          Dx = factor(Dx,levels=c("HC","BN"),labels=c("Healthy Controls",
                                                      "Bulimia Nervosa")))
bias.lm <- lmer (data=bias.df,
                    formula = vals ~ Dx * cond  + (1|idx),
                    REML=F,
                    control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=20000)))
summary(bias.lm) # Table 26


#########################
# Response to reviewers #
#########################

# What emotions significantly changed?
POMS = mh %>%
  dplyr::select(idx,Dx,
    PostPOMSTotal.NA_Neg , PrePOMSTotal.NA_Neg,
    PostPOMSTotal.NA_Neu , PrePOMSTotal.NA_Neu,
    PostPOMSAnger_Neg , PrePOMSAnger_Neg,
    PostPOMSAnger_Neu , PrePOMSAnger_Neu,
    PostPOMSConfusion_Neg , PrePOMSConfusion_Neg,
    PostPOMSConfusion_Neu , PrePOMSConfusion_Neu,
    PostPOMSDepression_Neg , PrePOMSDepression_Neg,
    PostPOMSDepression_Neu , PrePOMSDepression_Neu,
    PostPOMSTension_Neg , PrePOMSTension_Neg,
    PostPOMSTension_Neu , PrePOMSTension_Neu,
    PostPOMSFatigue_Neg , PrePOMSFatigue_Neg,
    PostPOMSFatigue_Neu , PrePOMSFatigue_Neu,
    PostPOMSVigour_Neg , PrePOMSVigour_Neg,
    PostPOMSVigour_Neu , PrePOMSVigour_Neu,
    ) %>%
    pivot_longer(cols=c(PostPOMSTotal.NA_Neg , PrePOMSTotal.NA_Neg,
    PostPOMSTotal.NA_Neu , PrePOMSTotal.NA_Neu,
    PostPOMSAnger_Neg , PrePOMSAnger_Neg,
    PostPOMSAnger_Neu , PrePOMSAnger_Neu,
    PostPOMSConfusion_Neg , PrePOMSConfusion_Neg,
    PostPOMSConfusion_Neu , PrePOMSConfusion_Neu,
    PostPOMSDepression_Neg , PrePOMSDepression_Neg,
    PostPOMSDepression_Neu , PrePOMSDepression_Neu,
    PostPOMSTension_Neg , PrePOMSTension_Neg,
    PostPOMSTension_Neu , PrePOMSTension_Neu,
    PostPOMSFatigue_Neg , PrePOMSFatigue_Neg,
    PostPOMSFatigue_Neu , PrePOMSFatigue_Neu,
    PostPOMSVigour_Neg , PrePOMSVigour_Neg,
    PostPOMSVigour_Neu , PrePOMSVigour_Neu,
                )) %>%
    mutate(item = ifelse(grepl("Total",name), "Total", 
            ifelse(grepl("Anger",name), "Anger",
            ifelse(grepl("Confusion",name), "Confusion",
            ifelse(grepl("Depression",name), "Depression",
            ifelse(grepl("Tension",name), "Tension", 
            ifelse(grepl("Fatigue",name), "Fatigue", "Vigor")
            ))))),
            time = factor(ifelse(grepl("Post",name),"Post","Pre"),
                    levels = c("Pre","Post")),
            cond = factor(ifelse(grepl("Neg",name), "Negative", "Neutral"),
                    levels = c("Neutral","Negative")),
            Dx = factor(Dx, levels = c("HC","BN")),
            item = factor(item, 
            levels = c("Total","Anger", "Confusion", "Depression", "Fatigue","Tension","Vigor"))
            ) %>%
    dplyr::select(!name)

Total_lmer = lmer(data = POMS[POMS$item == "Total",], 
                  formula = value ~ Dx * cond * time + (1|idx), REML=F)
summary(Total_lmer) # Table S1

# Simple effects showing affect change was no different between groups
emm_by_cond <- emmeans(Total_lmer, ~ Dx * time | cond)
pairs(emm_by_cond)
cond_contrasts <- contrast(emm_by_cond, 
                              interaction = "pairwise", 
                              by = "cond")
print(cond_contrasts)

# Simple effects analyses of group differences within each time point
emm_by_Dx <- emmeans(Total_lmer, ~ cond * time | Dx)
pairs(emm_by_Dx)

# Coehen's D
BN_n = length(unique(POMS$idx[POMS$Dx == "BN"]))
BE_negpre_negpost = 4.837 / sqrt(BN_n)
BE_neupre_neupost = 1.701 / sqrt(BN_n)
BE_neg_neu = 4.134 / sqrt(BN_n)

dx_contrasts <- contrast(emm_by_Dx, 
                              interaction = "pairwise", 
                              by = "cond")
print(dx_contrasts)

# For the Supplementary Materials
Anger_lmer = lmer(data = POMS[POMS$item == "Anger",], formula = value ~  Dx * cond * time + (1|idx), REML=F)
Confusion_lmer = lmer(data = POMS[POMS$item == "Confusion",], formula = value ~ Dx * cond * time + (1|idx), REML=F)
Depression_lmer = lmer(data = POMS[POMS$item == "Depression",], formula = value ~ Dx * cond * time  + (1|idx), REML=F)
Tension_lmer = lmer(data = POMS[POMS$item == "Tension",], formula = value ~ Dx * cond * time  + (1|idx), REML=F)
Fatigue_lmer = lmer(data = POMS[POMS$item == "Fatigue",], formula = value ~ Dx * cond * time + (1|idx), REML=F)
Vigour_lmer = lmer(data = POMS[POMS$item == "Vigor",], formula = value ~ Dx * cond * time  + (1|idx), REML=F)

summary(Anger_lmer) # Greater in negative (Table S27)
summary(Confusion_lmer) # Greater in negative (Table S28)
summary(Depression_lmer) # Greater in negative (Table S29)
summary(Fatigue_lmer) # NS greater in negative (Table S30)
summary(Tension_lmer) # Greater in negative (Table S31)
summary(Vigour_lmer) # Less in negative (Table S32)

F.s2 = 
POMS %>%
  filter(item != "Total") %>%
  ggplot(aes(x = cond, y = value, 
      group = interaction(time,Dx),
      color = Dx, fill = Dx,
      alpha = time) ) +
    theme_pubr(base_size = 44) +
    facet_wrap(~item, scales = "free") +
    stat_summary(geom="col",stat= "identity",
      position = position_dodge2(width = .5)) +
    scale_color_brewer(type = "qual") +
    scale_fill_brewer(type = "qual") +
    scale_alpha_discrete(range = c(0.5, 0.9)) +
    labs(x = element_blank(), y = "POMS Score",
    color = "Group", fill = "Group", alpha = "Time",group = element_blank() ) +
    theme(axis.text.x = element_text(angle = 25, vjust = 1, hjust=1)) +
    theme(legend.position=c(.07,.94),legend.box = "horizontal")

png(filename=figPath / "figureS2.png",  
width = 3000, height = 2400)
plot(F.s2)
dev.off()

# Relationship with other self report measures

alt_selfreport = params_mh %>%
  filter(Dx == "BN", params %in% c("wt","wh","tHin")) %>%
  mutate(SBE_M1 = EDE.SBE.Month.1..episodes.,
         SBE_M2 = EDE.SBE.Month.2..episodes.,
         SBE_M3 = EDE.SBE.Month.3..episodes.,
         SBE_SUM = SBE_M1 + SBE_M2 + SBE_M3) %>%
  dplyr::select(idx,cond,foodType,params,vals,EDE.Q.Restraint,
  DERS.Total.Score, BDI, EDE.Q.Total.Score,PostNeg_STAI.S,
  PostNeu_STAI.S,STAI.T,UPPS.P.Negative.Urgency, SBE_SUM
         ) %>%
  group_by(idx,params) %>%
  mutate(params = recode(params,
                        wt = 'omega[taste]',
                        wh = 'omega[health]',
                        tHin = 'tau[s]'),
         foodType = factor(foodType,
                           levels = c("Low-Fat","High-Fat")),
         cond = factor(cond,
                       levels = c("Neutral","Negative"))
  ) %>% 
  pivot_longer(values_to = "score", cols = c(DERS.Total.Score, BDI, EDE.Q.Total.Score,PostNeg_STAI.S,
  PostNeu_STAI.S,STAI.T,UPPS.P.Negative.Urgency,EDE.Q.Restraint, SBE_SUM)) %>%
  pivot_wider(names_from = c(cond,foodType), values_from = vals)

# Negative urgency
tau.urg = lm(data = alt_selfreport[alt_selfreport$params == "tau[s]" & alt_selfreport$name == "UPPS.P.Negative.Urgency", ],
          formula = score ~  `Neutral_Low-Fat` + `Neutral_High-Fat` + `Negative_Low-Fat` + `Negative_High-Fat`
  )
taste.urg = lm(data = alt_selfreport[alt_selfreport$params == "omega[taste]" & alt_selfreport$name == "UPPS.P.Negative.Urgency", ],
          formula = score ~  `Neutral_Low-Fat` + `Neutral_High-Fat` + `Negative_Low-Fat` + `Negative_High-Fat`
  )
health.urg = lm(data = alt_selfreport[alt_selfreport$params == "omega[health]" & alt_selfreport$name == "UPPS.P.Negative.Urgency", ],
          formula = score ~  `Neutral_Low-Fat` + `Neutral_High-Fat` + `Negative_Low-Fat` + `Negative_High-Fat`
  )
summary(tau.urg) 
summary(taste.urg)
summary(health.urg)
# Table S33

# Restrainttau.urg = lm(data = alt_selfreport[alt_selfreport$params == "tau[s]" & alt_selfreport$name == "UPPS.P.Negative.Urgency", ],
tau.restraint = lm(data = alt_selfreport[alt_selfreport$params == "tau[s]" & alt_selfreport$name == "EDE.Q.Restraint", ],
          formula = score ~  `Neutral_Low-Fat` + `Neutral_High-Fat` + `Negative_Low-Fat` + `Negative_High-Fat`
  )
taste.restraint = lm(data = alt_selfreport[alt_selfreport$params == "omega[taste]" & alt_selfreport$name == "EDE.Q.Restraint", ],
          formula = score ~  `Neutral_Low-Fat` + `Neutral_High-Fat` + `Negative_Low-Fat` + `Negative_High-Fat`
  )
health.restraint = lm(data = alt_selfreport[alt_selfreport$params == "omega[health]" & alt_selfreport$name == "EDE.Q.Restraint", ],
          formula = score ~  `Neutral_Low-Fat` + `Neutral_High-Fat` + `Negative_Low-Fat` + `Negative_High-Fat`
  )
summary(tau.restraint) 
summary(taste.restraint)
summary(health.restraint)
# Table S34

# Evaluate whether binge frequency measures are correlated with Restrictive subscale
selfreport = mh %>%
  filter(Dx == "BN") %>%
  mutate(SBE_M1 = EDE.SBE.Month.1..episodes.,
         SBE_M2 = EDE.SBE.Month.2..episodes.,
         SBE_M3 = EDE.SBE.Month.3..episodes.,
         SBE_SUM = SBE_M1 + SBE_M2 + SBE_M3,
         OBE_M1 = EDE.OBE.Month.1..episodes.,
         OBE_M2 = EDE.OBE.Month.2..episodes.,
         OBE_M3 = EDE.OBE.Month.3..episodes.,
         OBE_SUM = OBE_M1 + OBE_M2 + OBE_M3,
         ) %>%
  dplyr::select(idx,EDE.Q.Restraint, SBE_SUM,OBE_SUM) 


cor.test(selfreport$EDE.Q.Restraint, selfreport$SBE_SUM, method = "spearman")
cor.test(selfreport$EDE.Q.Restraint, selfreport$OBE_SUM, method = "spearman")
