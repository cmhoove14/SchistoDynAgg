library(tidyverse)
library(geepack)
library(grid)
library(gridExtra)
library(patchwork)

devtools::load_all()

load(here::here("Data/Derived/abc_processed_results.Rdata"))


# Case 4 partition parameter summary ----------------
Case4_part_par_by_W <- abc_fin_df_case_long %>%
  filter(Case=="Case4" & pop == "Comm") %>% 
  ggplot(aes(x = obsW.med,
             y = partW.med, 
             weight = partW.med/(partW.hiq-partW.loq))) +
  geom_point(aes(size = partW.med/(partW.hiq-partW.loq)),
             col = "#64a860",
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
  ylim(c(0,1)) +
  labs(y = "Case 4 partition parameter",
       x = "W")

Case4_part_par_by_W

ggsave(here::here("Figures/ABC_Case4_Partition_Parameter_by_W.png"),
       height = 4, width = 5, units = "in")

# Case 3 exposure susceptibility correlation----------------
Case3_corSC_by_W <- abc_fin_df_case_long %>%
  filter(Case=="Case3" & pop == "Comm") %>%
  ggplot(aes(x      = obsW.med,
             y      = corSC.med, 
             weight = corSC.med/(corSC.hiq-corSC.loq))) +
  geom_point(aes(size        = corSC.med/(corSC.hiq-corSC.loq)),
                 col         = "#7b6c7c",
                 alpha       = 0.3,
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
  #ylim(c(0,1)) +
  labs(y = "Case 3 Exp/Susc Correlation",
       x = "W")

Case3_corSC_by_W

ggsave(here::here("Figures/ABC_Case3_Susc_Exp_Corr_by_W.png"),
       height = 4, width = 5, units = "in")

