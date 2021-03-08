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
       height = 4, width = 4, units = "in")

# Bayes factors (Comparative estimate of model performance) for each case across observed egg burden -------------
abc_bayesF_comp <- abc_fin_df_case_long %>% 
  filter(Case %in% c("IIItoI", "IIItoII", "IItoI",
                     "IVtoI", "IVtoII", "IVtoIII") & pop == "Comm") %>% 
  ggplot(aes(x = UF_mean, y = BayesF, col = Case)) +
  stat_smooth() +
  geom_hline(yintercept = 1) +
  geom_point(alpha = 0.5) +
  scale_y_continuous(trans = "log",
                     breaks = c(0.1, 1, 10, 100, 1000),
                     labels = c("0.1", "1", "10", "100", "1000"),
                     limits = c(0.1,1000)) +
  scale_x_continuous(trans = "log",
                     breaks = c(0.01, 0.1, 1, 10, 100),
                     labels = c("0.01", "0.1", "1", "10", "100")) +
  # scale_color_manual(values = c("darkblue", "darkred", "forestgreen"),
  #                    labels = c("Case 3 to Case 1", "Case 3 to Case 2","Case 2 to Case 1")) +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10)) +
  labs(x = expression(paste("Mean egg burden (\U2130"[st],")")),
       y = "Bayes Factor",
       col = "Comparison")

abc_bayesF_comp

ggsave(here::here("Figures/ABC_Bayes_Factors.png"),
       height = 5, width = 7, units = "in")

