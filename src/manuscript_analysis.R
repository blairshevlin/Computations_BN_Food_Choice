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
# 06/17/24      Blair Shevlin                         Modeling code for original manuscript

# Packages required
required_packages <- c(
  "tidyverse",
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
  "glmmTMB",
  "bbmle",
  "DHARMa"
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

# Standard error function
se = function(x) sd(x) / sqrt(length(x))

# Folders
datPath <- path(here()) / 'data' 
resPath <- path(here()) / 'results'
figPath <- path(here()) / 'results' / 'figures'

# Load behavioral data
beh.df <- read.csv(file = datPath / "deidentified_ChoiceData.csv") %>%
  mutate(cond = factor(cond, levels = c("Neutral","Negative")),
         food = factor(foodType, levels = c("lf","hf"), labels = c("Low fat","High fat")),
         Dx = factor(Dx, levels = c("HC","BN")))

# Contrasts
contrasts(beh.df$cond) <- c(-1,1)
contrasts(beh.df$food) <- c(-1,1)
contrasts(beh.df$Dx) <- c(-1,1)

# Model 1: choice ~ food_type x affect x group
glm.1 <- glmer(data = beh.df,
               formula = choice ~
                 Dx * cond * food +
                 (1 + cond * food|idx),
               family=binomial,
               control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=20000)))
# summary(glm.1)
# Table S.8 for supplements

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
#summary(glm.2)
# Table S.9 for supplements

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
#summary(glm.3)
# Supplementary Table S10

# Model 4: Response times
beh.df$choice_c = factor(beh.df$choice, levels = c(0,1), labels = c("Reference item","Presented item"))
contrasts(beh.df$choice_c) <- c(-1,1)

lm.1 <- lmer(data = beh.df,
               formula = log(rt) ~   Dx * cond * food * choice_c +
                 (1 + cond * food | idx),
               REML = F,
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=20000)))
#summary(lm.1)
# Supplementary Table S.11


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

ggsave(file = figPath / "figure2.png", plot = figure2, width = 12, height = 8)

##################################
# Computational modeling results #
##################################

# Load subject-level parameters
df.fit.full <- NULL

for (cc in c("Neutral","Negative")) {
  for (gg in c("bn","hc")) {
    load(file.path(resPath / 'stDDM',paste("/params_HtSSM_FIT_M3_Dx-",gg,"_","Cond-",cc,"_rawRatings.RData",sep="")))
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
          Dx = factor(Dx,levels=c("hc","bn"),labels=c("Healthy Controls",
                                                      "Bulimia Nervosa")))
tHin.lm <- lmer (data=tHin.df,
                 formula = vals ~ Dx * cond * foodType+ (1|idx),
                 REML=F,
                 control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=20000)))
summary(tHin.lm) 
# Supplementary Table 1

#  Within group: cond x Food
tHin.lm.b1 <- lmer(data=tHin.df[tHin.df$Dx == "Bulimia Nervosa",],
                 formula = vals ~ cond / foodType + (1|idx),
                 REML=F,
                 control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=20000)))
tHin.lm.b2 <- lmer(data=tHin.df[tHin.df$Dx == "Healthy Controls",],
                 formula = vals ~ cond / foodType + (1|idx),
                 REML=F,
                 control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=20000)))
summary(tHin.lm.b2) 
summary(tHin.lm.b1);
# Supplementary Table 2

# Equivalence of cond term

# Clogg et al. (1995) formula as cited by Ray Paternoster et al. (1998)
b1 = summary(tHin.lm.b1)$coeff[2,1] # mean est of BN
s1 = summary(tHin.lm.b1)$coeff[2,2] # se of BN
b2 = summary(tHin.lm.b2)$coeff[2,1] # mean est of HC
s2 = summary(tHin.lm.b2)$coeff[2,2] # see of HC
v = (b1 - b2) / sqrt(s1^2 + s2^2)
data.frame(diff=b, zdiff=v, `p-value`=format(2*pnorm(-abs(v)), scientific=FALSE))

# Attribute Weights

## Taste
taste.df <- params %>%
  filter(params == "wt") %>%
  mutate( foodType = factor(foodType,
                            levels=c("Low-Fat","High-Fat"),
                            labels = c("Low-Fat","High-Fat")),
          cond = factor(cond,levels=c("Neutral","Negative")),
          Dx = factor(Dx,levels=c("hc","bn"),labels=c("Healthy Controls",
                                                      "Bulimia Nervosa")))
taste.lm1 <- lmer (data=taste.df,
               formula = vals ~ Dx * cond * foodType+ (1|idx),
               REML=F,
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=20000)))
summary(taste.lm1)
# Supplementary Table 3

# Within group: cond x Food
taste.lm.b1 <- lmer(data=taste.df[taste.df$Dx == "Bulimia Nervosa",],
                 formula = vals ~ cond * foodType + (1|idx),
                 REML=F,
                 control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=20000)))
taste.lm.b2 <- lmer(data=taste.df[taste.df$Dx == "Healthy Controls",],
                 formula = vals ~ cond * foodType + (1|idx),
                 REML=F,
                 control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=20000)))
summary(taste.lm.b1); 
summary(taste.lm.b2) 
# Supplementary Table 5

## Health
health.df <- params %>%
  filter(params == "wh") %>%
  mutate( foodType = factor(foodType,
                            levels=c("Low-Fat","High-Fat"),
                            labels = c("Low-Fat","High-Fat")),
          cond = factor(cond,levels=c("Neutral","Negative")),
          Dx = factor(Dx,levels=c("hc","bn"),labels=c("Healthy Controls",
                                                      "Bulimia Nervosa")))
health.lm1 <- lmer (data=health.df,
                    formula = vals ~ Dx * cond * foodType + (1|idx),
                    REML=F,
                    control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=20000)))
summary(health.lm1)
# Supplementary Table 4

# Within group: cond x Food
health.lm.b1 <- lmer(data=health.df[health.df$Dx == "Bulimia Nervosa",],
                 formula = vals ~ cond * foodType + (1|idx),
                 REML=F,
                 control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=20000)))
health.lm.b2 <- lmer(data=health.df[health.df$Dx == "Healthy Controls",],
                 formula = vals ~ cond * foodType + (1|idx),
                 REML=F,
                 control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=20000)))
summary(health.lm.b1); 
summary(health.lm.b2) 
# Supplementary Table 6

# Figure 4
fig4a <- params %>%
  filter(params == "tHin") %>%
  mutate(grouping = ifelse(Dx == "bn" & cond == "Negative",
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
  Dx = factor(Dx,levels=c("hc","bn"),labels=c("HC",
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
  mutate(grouping = ifelse(Dx == "bn" & cond == "Negative",
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
  Dx = factor(Dx,levels=c("hc","bn"),labels=c("HC",
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
  mutate(grouping = ifelse(Dx == "bn" & cond == "Negative",
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
  Dx = factor(Dx,levels=c("hc","bn"),labels=c("HC",
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

ggsave(file = figPath / "figure4.png", plot = figure4, width = 12, height = 8)

###################
# Symtom Severity #
###################
params = params %>%
  mutate(Dx = factor(Dx, levels = c("hc","bn"), labels = c("HC","BN")),)

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
health_text_neg <- data.frame(SUM = -10,vals = -.6,lab = "Health Earlier",
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

ggsave(file = figPath / "figure5.png", plot = figure5, width = 12, height = 8)


# Change scores
LOC.change = LOC %>%
filter(params == "tau[s]") %>%
  group_by(SBE_SUM,OBE_SUM,idx,cond,params) %>%
  summarise(change = vals[foodType=="High-Fat"] - vals[foodType == "Low-Fat"])

summary(glm.nb (data = LOC.change,
                formula = SBE_SUM ~ change * cond ) )
summary(glm.nb (data = LOC.change,
                formula = OBE_SUM ~ change * cond ) )
                
# Supplementary Table 7



