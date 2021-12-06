library(grid)
library(gridExtra)
library(patchwork)

devtools::load_all()

load(here::here("Data/Derived/abc_processed_results.Rdata"))

# Plot dispersion by worm burden   
W_by_E_plot <- abc_fin_df_case_long %>%
  filter(Case %in% c("Case1", "Case2", "Case3", "Case4") & pop == "Comm") %>% 
  ggplot(aes(x = UF_mean,
             y = obsW.med, 
             #shape = Intervention,
             col = Case,
             weight = obsW.med/(obsW.hiq-obsW.loq))) +
  geom_point(aes(size = obsW.med/(obsW.hiq-obsW.loq)),
             alpha = 0.5,
             show.legend = FALSE) +
  stat_smooth() +
  theme_classic() +
  theme(legend.position = "bottom") +
  scale_x_continuous(trans = "log",
                     breaks = c(0.01, 0.1, 1, 10, 100),
                     labels = c("0.01", "0.1", "1", "10", "100")) +
  scale_y_continuous(trans = "log",
                     breaks = c(0.01, 0.1, 1, 10, 100),
                     labels = c("0.01", "0.1", "1", "10", "100")) +
  scale_color_manual(values = c("#058dd6", "#cc4a49", "#7b6c7c","#64a860")) +
  labs(y = expression(Estimated~Mean~Worm~Burden~(W[st])),
       x = expression(Observed~Mean~Egg~Burden~(E[st])),
       col = "")

# Plot dispersion by worm burden   
kW_by_E_plot <- abc_fin_df_case_long %>%
  filter(Case %in% c("Case1", "Case2", "Case3", "Case4") & pop == "Comm") %>% 
  ggplot(aes(x = UF_mean,
             y = obsalphaW.med^-1, 
             #shape = Intervention,
             col = Case,
             weight = obsalphaW.med/(obsalphaW.hiq-obsalphaW.loq))) +
  stat_smooth() +
  geom_point(aes(size = obsalphaW.med/(obsalphaW.hiq-obsalphaW.loq)),
             alpha = 0.5,
             show.legend = FALSE) +
  #facet_grid(Intervention~Isl) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_x_continuous(trans = "log",
                     breaks = c(0.01, 0.1, 1, 10, 100),
                     labels = c("0.01", "0.1", "1", "10", "100")) +
  scale_y_continuous(trans = "log",
                     breaks = c(0.001, 0.01, 0.1, 1),
                     labels = c("0.001", "0.01", "0.1", "1")) +
  scale_color_manual(values = c("#058dd6", "#cc4a49", "#7b6c7c", "#64a860")) +
  labs(y = expression(Estimated~Worm~Aggregation~Parameter~(k[st]^W)),
       x = expression(Observed~Mean~Egg~Burden~(E[st])))

(W_by_E_plot | kW_by_E_plot)

ggsave(here::here("Figures/ABC-Comp_W_to_obsE.png"),
       height = 5, width = 7, units = "in")
