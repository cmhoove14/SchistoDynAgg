library(tidyverse)
library(geepack)
library(grid)
library(gridExtra)
library(patchwork)

devtools::load_all()

load(here::here("Data/Derived/abc_processed_results.Rdata"))

# Case1 PLot ----------------------
kap_W_case1_plot <- abc_fin_df_case_long %>%
  filter(Case=="Case1" & pop == "Comm") %>% 
  ggplot(aes(x = obsW.med,
             y = obsalphaW.med^-1, 
             weight = obsalphaW.med/(obsalphaW.hiq-obsalphaW.loq))) +
  geom_point(aes(size = obsalphaW.med/(obsalphaW.hiq-obsalphaW.loq)),
             col = "#058dd6",
             alpha = 0.3,
             show.legend = FALSE) +
  stat_smooth(col = "black") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.1, vjust = -1)) +
  scale_size(range = c(1,4)) +
  scale_x_continuous(trans = "log",
                     breaks = c(0.01, 0.1, 1, 10, 100),
                     labels = c("0.01", "0.1", "1", "10", "100"),
                     limits = c(min(abc_fin_df_case_long$obsW.med, na.rm = T),
                                max(abc_fin_df_case_long$obsW.med, na.rm = T))) +
  scale_y_continuous(trans = "log",
                     breaks = c(0.001,0.01, 0.1, 1, 10, 100),
                     labels = c("0.001","0.01", "0.1", "1", "10", "100"),
                     limits = c(0.001,1000)) +
  labs(y = expression(Worm~Dispersion~Parameter~(kappa[st]^W)),
       title = "A")

# Case2 PLot ----------------------
kap_W_case2_plot <- abc_fin_df_case_long %>%
  filter(Case=="Case2" & pop == "Comm") %>% 
  ggplot(aes(x = obsW.med,
             y = obsalphaW.med^-1, 
             weight = obsalphaW.med/(obsalphaW.hiq-obsalphaW.loq))) +
  geom_point(aes(size = obsalphaW.med/(obsalphaW.hiq-obsalphaW.loq)),
             col = "#cc4a49",
             alpha = 0.3,
             show.legend = FALSE) +
  stat_smooth(col = "black") +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.title = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 0.1, vjust = -1)) +
  scale_size(range = c(1,4)) +
  scale_x_continuous(trans = "log",
                     breaks = c(0.01, 0.1, 1, 10, 100),
                     labels = c("0.01", "0.1", "1", "10", "100"),
                     limits = c(min(abc_fin_df_case_long$obsW.med, na.rm = T),
                                max(abc_fin_df_case_long$obsW.med, na.rm = T))) +
  scale_y_continuous(trans = "log",
                     breaks = c(0.001,0.01, 0.1, 1, 10, 100),
                     labels = c("0.001","0.01", "0.1", "1", "10", "100"),
                     limits = c(0.001,1000)) +
  labs(title = "B")

# Case3 PLot ----------------------
kap_W_case3_plot <- abc_fin_df_case_long %>%
  filter(Case=="Case3" & pop == "Comm") %>% 
  ggplot(aes(x = obsW.med,
             y = obsalphaW.med^-1, 
             weight = obsalphaW.med/(obsalphaW.hiq-obsalphaW.loq))) +
  geom_point(aes(size = obsalphaW.med/(obsalphaW.hiq-obsalphaW.loq)),
             col = "#7b6c7c",
             alpha = 0.3,
             show.legend = FALSE) +
  stat_smooth(col = "black") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.1, vjust = -1)) +
  scale_size(range = c(1,4)) +
  scale_x_continuous(trans = "log",
                     breaks = c(0.01, 0.1, 1, 10, 100),
                     labels = c("0.01", "0.1", "1", "10", "100"),
                     limits = c(min(abc_fin_df_case_long$obsW.med, na.rm = T),
                                max(abc_fin_df_case_long$obsW.med, na.rm = T))) +
  scale_y_continuous(trans = "log",
                     breaks = c(0.001,0.01, 0.1, 1, 10, 100),
                     labels = c("0.001","0.01", "0.1", "1", "10", "100"),
                     limits = c(0.001,1000)) +
  labs(x = expression(Mean~worm~burden~(W[st])),
       y = expression(Worm~Dispersion~Parameter~(kappa[st]^W)),
       title = "C")

# Case4 PLot ----------------------
kap_W_case4_plot <- abc_fin_df_case_long %>%
  filter(Case=="Case4" & pop == "Comm") %>% 
  ggplot(aes(x = obsW.med,
             y = obsalphaW.med^-1, 
             weight = obsalphaW.med/(obsalphaW.hiq-obsalphaW.loq))) +
  geom_point(aes(size = obsalphaW.med/(obsalphaW.hiq-obsalphaW.loq)),
             col = "#64a860",
             alpha = 0.3,
             show.legend = FALSE) +
  stat_smooth(col = "black") +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 0.1, vjust = -1)) +
  scale_size(range = c(1,4)) +
  scale_x_continuous(trans = "log",
                     breaks = c(0.01, 0.1, 1, 10, 100),
                     labels = c("0.01", "0.1", "1", "10", "100"),
                     limits = c(min(abc_fin_df_case_long$obsW.med, na.rm = T),
                                max(abc_fin_df_case_long$obsW.med, na.rm = T))) +
  scale_y_continuous(trans = "log",
                     breaks = c(0.001,0.01, 0.1, 1, 10, 100),
                     labels = c("0.001","0.01", "0.1", "1", "10", "100"),
                     limits = c(0.001,1000)) +
  labs(x = expression(Mean~worm~burden~(W[st])),
       title = "D")

#Fig 2: ABC kappa by W estimates ------------
png(here::here("Figures/Fig2_ABC_Kap_by_W_Results.png"),
    height = 6, width = 6, units = "in", res = 300)

grid.arrange(arrangeGrob(kap_W_case1_plot + theme(axis.title = element_blank()), 
                         kap_W_case2_plot + theme(axis.title = element_blank()),
                         kap_W_case3_plot + theme(axis.title = element_blank()),
                         kap_W_case4_plot + theme(axis.title = element_blank()), 
                         nrow = 2, ncol = 2,
                         left = textGrob(expression(Worm~Dispersion~Parameter~(kappa[st]^W)), rot = 90, vjust = 1),
                         bottom = textGrob(expression(Mean~worm~burden~(W[st])), vjust = 0.2)))



#(kap_W_case1_plot | kap_W_case2_plot) / (kap_W_case3_plot | kap_W_case4_plot)

dev.off()

pdf(here::here("Figures/Fig2_ABC_Kap_by_W_Results.pdf"),
    height = 6, width = 6)

grid.arrange(arrangeGrob(kap_W_case1_plot + theme(axis.title = element_blank()), 
                         kap_W_case2_plot + theme(axis.title = element_blank()),
                         kap_W_case3_plot + theme(axis.title = element_blank()),
                         kap_W_case4_plot + theme(axis.title = element_blank()), 
                         nrow = 2, ncol = 2,
                         left = textGrob(expression(Worm~Dispersion~Parameter~(kappa[st]^W)), rot = 90, vjust = 1),
                         bottom = textGrob(expression(Mean~worm~burden~(W[st])), vjust = 0.2)))

#(kap_W_case1_plot | kap_W_case2_plot) / (kap_W_case3_plot | kap_W_case4_plot) 

dev.off()

