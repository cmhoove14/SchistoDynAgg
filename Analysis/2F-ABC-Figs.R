library(tidyverse)
library(geepack)
library(patchwork)

devtools::load_all()

load(here::here("Data/Derived/abc_processed_results.Rdata"))

# Case1 PLot ----------------------
kap_W_case1_plot <- abc_fin_df_case_long %>%
  filter(Case=="Case1" & pop == "Child") %>% 
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
  filter(Case=="Case2" & pop == "Child") %>% 
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
  filter(Case=="Case3" & pop == "Child") %>% 
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

# Case4 PLot ----------------------
kap_W_case4_plot <- abc_fin_df_case_long %>%
  filter(Case=="Case4" & pop == "Child") %>% 
  ggplot(aes(x = obsW.med,
             y = obsalphaW.med^-1, 
             weight = obsalphaW.med/(obsalphaW.hiq-obsalphaW.loq))) +
  geom_point(aes(size = obsalphaW.med/(obsalphaW.hiq-obsalphaW.loq)),
             col = "#64a860",
             alpha = 0.3,
             show.legend = FALSE) +
  stat_smooth(col = "#64a860") +
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
  labs(title = "Case 4",
       tag = "D")


# Bayes factors (omparative estimate of model performance) for each case across observed egg burden
abc_bayesF_comp <- abc_fin_df_case_long %>% 
  filter(Case %in% c("IIItoI", "IIItoII", "IItoI") & pop == "Child") %>% 
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
       tag = "E")

#Fig 2: ABC kappa by W estimates ------------
png(here::here("Figures/Fig2_ABC_Kap_by_W_Results.png"),
    height = 7, width = 7, units = "in", res = 300)

(kap_W_case1_plot | kap_W_case2_plot) / (kap_W_case3_plot | kap_W_case4_plot) / abc_bayesF_comp #+ plot_layout(heights = c(2,3))

dev.off()

pdf(here::here("Figures/Fig2_ABC_Kap_by_W_Results.pdf"),
    height = 7, width = 7)

(kap_W_case1_plot | kap_W_case2_plot) / (kap_W_case3_plot | kap_W_case4_plot) / abc_bayesF_comp #+ plot_layout(heights = c(2,3))

dev.off()


# Fig 3: Mating probability -------------
test_ws <- exp_seq(1e-3,100,200)

# Mean worm burden by aggregation for case 1
alphaW_gee_case1 <- geeglm(log(obsalphaW.med) ~ log(obsW.med), id = as.factor(Shehia),
                           family = "gaussian", corstr = "unstructured",
                           weights = obsalphaW.med/(obsalphaW.hiq-obsalphaW.loq),
                           data = abc_fin_df_case_long %>% filter(Case == "Case1" & pop == "Child"))

geeW_coef1_case1 <- coef(alphaW_gee_case1)[1]
geeW_coef2_case1 <- coef(alphaW_gee_case1)[2]

zanz_kappa_W_gee_case1 <- function(W){
  exp(geeW_coef1_case1+geeW_coef2_case1*log(W))^-1
}


phi_pred <- as.data.frame(cbind("W" = test_ws,
                                "Phi_k005" = sapply(test_ws, phi_Wk, k = 0.05),
                                "Phi_dynak" = mapply(phi_Wk, test_ws, zanz_kappa_W_gee_case1(test_ws)),
                                "Phi_case2k005" = sapply(test_ws, phi_wk_sep, k = 0.05),
                                "Phi_case2dynak" = mapply(phi_wk_sep, test_ws, zanz_kappa_W_gee_case1(test_ws))))

# Case1 mate prob plot ----------------------
abc_case1_mate_prob <- abc_fin_df_case_long %>% 
  filter(Case=="Case1" & pop == "Child") %>% 
  ggplot(aes(x = obsW.med,
             y = obsPhiProbW.med, 
             size = obsPhiProbW.med/(obsPhiProbW.hiq-obsPhiProbW.loq))) +
  geom_point(col = "#058dd6",
             alpha = 0.3,
             show.legend = FALSE) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none") +
  scale_size(range = c(1,4)) +
  ylim(c(0,1)) +
  scale_x_continuous(trans = "log",
                     breaks = c(0.01, 0.1, 1, 10, 100),
                     labels = c("0.01", "0.1", "1", "10", "100"),
                     limits = c(min(abc_fin_df_case_long$obsW.med, na.rm = T),
                                max(abc_fin_df_case_long$obsW.med, na.rm = T))) +
  labs(title = "Case 1",
       tag = "A")

# Case2 mate prob plot ----------------------
abc_case2_mate_prob <- abc_fin_df_case_long %>% 
  filter(Case=="Case2" & pop == "Child") %>% 
  ggplot(aes(x = obsW.med,
             y = obsPhiProbW.med, 
             size = obsPhiProbW.med/(obsPhiProbW.hiq-obsPhiProbW.loq))) +
  geom_point(col = "#cc4a49",
             alpha = 0.3,
             show.legend = FALSE) +
  theme_classic() +
  theme(axis.title.y = element_text(size = 9.5),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none") +
  scale_size(range = c(1,4)) +
  ylim(c(0,1)) +
  scale_x_continuous(trans = "log",
                     breaks = c(0.01, 0.1, 1, 10, 100),
                     labels = c("0.01", "0.1", "1", "10", "100"),
                     limits = c(min(abc_fin_df_case_long$obsW.med, na.rm = T),
                                max(abc_fin_df_case_long$obsW.med, na.rm = T)))+
  labs(y = expression(Mating~probability~(Phi[st])),
       title = "Case 2",
       tag = "B")

# Case3 mate prob plot ----------------------
abc_case3_mate_prob <- abc_fin_df_case_long %>% 
  filter(Case=="Case3" & pop == "Child") %>% 
  ggplot(aes(x = obsW.med,
             y = obsPhiProbW.med, 
             size = obsPhiProbW.med/(obsPhiProbW.hiq-obsPhiProbW.loq))) +
  geom_point(col = "#7b6c7c",
             alpha = 0.3,
             show.legend = FALSE) +
  theme_classic() +
  theme(axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none") +
  scale_size(range = c(1,4)) +
  ylim(c(0,1)) +
  scale_x_continuous(trans = "log",
                     breaks = c(0.01, 0.1, 1, 10, 100),
                     labels = c("0.01", "0.1", "1", "10", "100"),
                     limits = c(min(abc_fin_df_case_long$obsW.med, na.rm = T),
                                max(abc_fin_df_case_long$obsW.med, na.rm = T)))+
  labs(x = expression(Mean~worm~burden~(W[st])),
       title = "Case 3",
       tag = "C")

# Case4 mate prob plot ----------------------
abc_case4_mate_prob <- abc_fin_df_case_long %>% 
  filter(Case=="Case4" & pop == "Child") %>% 
  ggplot(aes(x = obsW.med,
             y = obsPhiProbW.med, 
             size = obsPhiProbW.med/(obsPhiProbW.hiq-obsPhiProbW.loq))) +
  geom_point(col = "#64a860",
             alpha = 0.3,
             show.legend = FALSE) +
  theme_classic() +
  theme(axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none") +
  scale_size(range = c(1,4)) +
  ylim(c(0,1)) +
  scale_x_continuous(trans = "log",
                     breaks = c(0.01, 0.1, 1, 10, 100),
                     labels = c("0.01", "0.1", "1", "10", "100"),
                     limits = c(min(abc_fin_df_case_long$obsW.med, na.rm = T),
                                max(abc_fin_df_case_long$obsW.med, na.rm = T)))+
  labs(x = expression(Mean~worm~burden~(W[st])),
       title = "Case 4",
       tag = "D")

# Mating probabilities ------------
mate_prob_lines <- phi_pred %>% 
  pivot_longer(Phi_k005:Phi_case2dynak) %>% 
  mutate(kap_dyna = if_else(grepl("dynak", name), "Dynamic", "Constant"),
         Case1_2 = if_else(grepl("case2", name), "Case2", "Case1")) %>% 
  filter(Case1_2 == "Case1") %>% 
  ggplot() +
  geom_line(aes(x = W, y = value, lty = kap_dyna),
            size = 1.2) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_x_continuous(trans = "log",
                     breaks = c(0.01, 0.1, 1, 10, 100),
                     labels = c("0.01", "0.1", "1", "10", "100"),
                     limits = c(0.05, 200)) +
  ylim(c(0,1)) +
  labs(y = expression(Mating~Probability~(Phi(W, kappa))),
       x = expression(Mean~worm~burden~(W[st])),
       tag = "D")

abc_w_mate_prob_lines <- mate_prob_lines +
  stat_smooth(data = abc_fin_df_case_long %>% filter(Case %in% c("Case1", "Case2", "Case3", "Case4") & pop == "Child"),
              aes(x = obsW.med, y = obsPhiProbW.med, col = Case)) +
  scale_color_manual(values = c("#058dd6", "#cc4a49", "#7b6c7c", "#64a860")) +
  theme(plot.title = element_text(size = 14,hjust = 0.5))# + labs(title = "Observed and Estimated\nMating Probability")

png(here::here("Figures/Fig3_ABC_MateProb_Results.png"),
    height = 5.5, width = 7, units = "in", res = 300)

(abc_case1_mate_prob | abc_case2_mate_prob) / (abc_case3_mate_prob | abc_case4_mate_prob) | abc_w_mate_prob_lines + plot_layout(widths = c(3,2))

dev.off()

pdf(here::here("Figures/Fig3_ABC_MateProb_Results.pdf"),
    height = 5.5, width = 7)

(abc_case1_mate_prob | abc_case2_mate_prob) / (abc_case3_mate_prob | abc_case4_mate_prob) | abc_w_mate_prob_lines + plot_layout(widths = c(3,2))

dev.off()
