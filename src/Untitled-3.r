
fig4a_alt <- 
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
                    levels=c("NA",
                    "Low-Fat","High-Fat"),
                    labels = c("All Foods",
                    "Low Fat","High Fat")),
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
                             x = 1.35, 
                             y = c(0.03, -0.03),
                             label = c("Health sooner", "Taste sooner")),
            aes(x = x, y = y, label = label),
            inherit.aes = FALSE,
            size = 6, color = "black") +
  labs(
       y = "Attribute Onset",
       color = "Group",
      x= element_blank()
      ) +
  coord_cartesian(ylim = c(-0.75, 0.25)) +
  scale_color_brewer(type = "qual") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1, hjust=1)) +
  theme(legend.position = "right", legend.box = "vertical")  +
  theme(legend.position=c(.9,.58))


fig4b_alt <- params %>%
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


fig5_alt = 
      ggplot(LOC.df[LOC.df$params == "tau[s]" & LOC.df$LOC == "SBE_SUM" & LOC.df$cond == "Negative",], 
       aes(x = vals, y = SUM, color = foodType, fill = foodType,
           linetype = foodType, shape = foodType)) +
  theme_pubr(base_size = 18) +
  geom_point(alpha = .5,size=3) + 
  geom_smooth(              alpha = 0.5,
              size = 3,
              method = "glm.nb",
              se = F) + 
  scale_color_brewer(type = "qual",palette = 6,direction = -1) +
  scale_fill_brewer(type = "qual",palette = 6,direction = -1) +
  coord_cartesian(ylim = c(-10,200)) +
  labs(x = "Attribute Onset Estimate",
       y = "Subjective Binge Episodes\n(3 Month Total)",
       shape = "Food Type",
       linetype = "Food Type",
       color = "Food Type", fill = "Food Type",) +
  theme(legend.position="right",
        legend.key.size =  unit(0.5, "in"),
        axis.text.x = element_text(
                                   angle = 45,
                                   vjust=1, hjust=1),
        panel.spacing.x = unit(10, "mm"))

taste_text_neg <- data.frame(SUM = -10,vals = -.75,lab = "Taste Earlier",
                             cond = factor( "Negative",c("Neutral","Negative")),
                             foodType = c("Low-Fat","Low-Fat"))
health_text_neg <- data.frame(SUM = -10,vals = -.63,lab = "Health Earlier",
                             cond = factor( "Negative",c("Neutral","Negative")),
                             foodType = c("Low-Fat","Low-Fat"))

figure5_alt = fig5_alt + 
  geom_text(data = taste_text_neg,label = "Taste Earlier",color="black",size=6) +
  geom_text(data = health_text_neg,label = "Health Earlier",color="black",size=6) 


tmp = fig4a_alt + fig4b_alt + plot_layout(guides = "collect") & theme(legend.position = "bottom")
ggsave(file = figPath / "tmp.tiff", plot = tmp, width = 12, height = 8)

tmp2 = fig4a_alt + figure5_alt
ggsave(file = figPath / "tmp2.tiff", plot = tmp2, width = 12, height = 8)
