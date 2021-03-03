load(here::here("Data/Derived/abc_processed_results.Rdata"))

# Case1 PLot ----------------------
kap_W_case1_plot <- abc_fin_df_case_long %>%
  filter(Case=="Case1" & Pop == "Child") %>% 
  ggplot(aes(x = obsW.med,
             y = obsalphaW.med^-1, 
             weight = obsalphaW.med/(obsalphaW.hiq-obsalphaW.loq))) +
  geom_point(aes(size = obsalphaW.med/(obsalphaW.hiq-obsalphaW.loq)),
             col = "#058dd6",
             alpha = 0.3,
             show.legend = FALSE) +
  stat_smooth(col = "#058dd6") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  scale_size(range = c(1,4)) +
  scale_x_continuous(trans = "log",
                     breaks = c(0.01, 0.1, 1, 10, 100),
                     labels = c("0.01", "0.1", "1", "10", "100"),
                     limits = c(min(abc_fin_df_case_long$obsW.med, na.rm = T),
                                max(abc_fin_df_case_long$obsW.med, na.rm = T))) +
  scale_y_continuous(trans = "log",
                     breaks = c(0.01, 0.1, 1, 10, 100),
                     labels = c("0.01", "0.1", "1", "10", "100"),
                     limits = c(min(abc_fin_df_case_long$obsalphaW.med^-1, na.rm = T),
                                max(abc_fin_df_case_long$obsalphaW.med^-1, na.rm = T))) +
  labs(y = expression(Worm~Dispersion~Parameter~(kappa[st]^W)),
       title = "Case 1",
       tag = "A")

# Case2 PLot ----------------------
kap_W_case2_plot <- abc_fin_df_case_long %>%
  filter(Case=="Case2" & Pop == "Child") %>% 
  ggplot(aes(x = obsW.med,
             y = obsalphaW.med^-1, 
             weight = obsalphaW.med/(obsalphaW.hiq-obsalphaW.loq))) +
  geom_point(aes(size = obsalphaW.med/(obsalphaW.hiq-obsalphaW.loq)),
             col = "#cc4a49",
             alpha = 0.3,
             show.legend = FALSE) +
  stat_smooth(col = "#cc4a49") +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  scale_size(range = c(1,4)) +
  scale_x_continuous(trans = "log",
                     breaks = c(0.01, 0.1, 1, 10, 100),
                     labels = c("0.01", "0.1", "1", "10", "100"),
                     limits = c(min(abc_fin_df_case_long$obsW.med, na.rm = T),
                                max(abc_fin_df_case_long$obsW.med, na.rm = T))) +
  scale_y_continuous(trans = "log",
                     breaks = c(0.01, 0.1, 1, 10, 100),
                     labels = c("0.01", "0.1", "1", "10", "100"),
                     limits = c(min(abc_fin_df_case_long$obsalphaW.med^-1, na.rm = T),
                                max(abc_fin_df_case_long$obsalphaW.med^-1, na.rm = T))) +
  labs(x = expression(Mean~worm~burden~(W[st])),
       title = "Case 2",
       tag = "B")

# Case3 PLot ----------------------
kap_W_case3_plot <- abc_fin_df_case_long %>%
  filter(Case=="Case3" & Pop == "Child") %>% 
  ggplot(aes(x = obsW.med,
             y = obsalphaW.med^-1, 
             weight = obsalphaW.med/(obsalphaW.hiq-obsalphaW.loq))) +
  geom_point(aes(size = obsalphaW.med/(obsalphaW.hiq-obsalphaW.loq)),
             col = "#7b6c7c",
             alpha = 0.3,
             show.legend = FALSE) +
  stat_smooth(col = "#7b6c7c") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  scale_size(range = c(1,4)) +
  scale_x_continuous(trans = "log",
                     breaks = c(0.01, 0.1, 1, 10, 100),
                     labels = c("0.01", "0.1", "1", "10", "100"),
                     limits = c(min(abc_fin_df_case_long$obsW.med, na.rm = T),
                                max(abc_fin_df_case_long$obsW.med, na.rm = T))) +
  scale_y_continuous(trans = "log",
                     breaks = c(0.01, 0.1, 1, 10, 100),
                     labels = c("0.01", "0.1", "1", "10", "100"),
                     limits = c(min(abc_fin_df_case_long$obsalphaW.med^-1, na.rm = T),
                                max(abc_fin_df_case_long$obsalphaW.med^-1, na.rm = T))) +
  labs(title = "Case 3",
       tag = "C")

# Bayes factors (omparative estimate of model performance) for each case across observed egg burden
abc_bayesF_comp <- abc_fin_df_case_long %>% 
  filter(Case %in% c("IIItoI", "IIItoII", "IItoI") & Pop == "Child") %>% 
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
  scale_color_manual(values = c("darkblue", "darkred", "forestgreen"),
                     labels = c("Case 3 to Case 1", "Case 3 to Case 2","Case 2 to Case 1")) +
  theme_classic() +
  theme(legend.position = c(0.15,0.8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10)) +
  labs(x = expression(paste("Mean egg burden (\U2130"[st],")")),
       y = "Bayes Factor",
       col = "Comparison",
       tag = "D")


png(here::here("Figures/Fig2_ABC_Kap_by_W_Results.png"),
    height = 7, width = 7, units = "in", res = 300)

(kap_W_case1_plot | kap_W_case2_plot | kap_W_case3_plot) / abc_bayesF_comp +plot_layout(heights = c(2,3))

dev.off()

pdf(here::here("Figures/Fig1_GEE_Results.pdf"),
    height = 7, width = 7)

(kap_W_case1_plot | kap_W_case2_plot | kap_W_case3_plot) / abc_bayesF_comp +plot_layout(heights = c(2,3))

dev.off()


test_ws <- exp_seq(1e-3,100,200)

phi_pred <- as.data.frame(cbind("W" = test_ws,
                                "Phi_k005" = sapply(test_ws, phi_Wk, k = 0.05),
                                "Phi_dynak" = mapply(phi_Wk, test_ws, zanz_kappa_W_gee_case1(test_ws)),
                                "Phi_case2k005" = sapply(test_ws, phi_wk_sep, k = 0.05),
                                "Phi_case2dynak" = mapply(phi_wk_sep, test_ws, zanz_kappa_W_gee_case1(test_ws))))
