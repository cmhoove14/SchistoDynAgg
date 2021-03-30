library(tidyverse)
library(geepack)
library(grid)
library(gridExtra)
library(patchwork)

devtools::load_all()

yO <- readRDS(here::here("Data/Derived/ABC_yO_data.rds"))

abc_post_preds_mateProb <- readRDS(here::here("Data/Derived/abc_post_pred_checks.rds")) %>% 
  filter(SumStat== "prob_Phi") %>% 
  
  # Temp solution to remove duplicates
  filter(!Shehia %in% c("KINYASINI", "MBUZINI")) %>% 
  
  
  left_join(yO, 
            by = c(#"Isl"    = "Isl", 
              "Shehia" = "Shehia", 
              "Year"   = "Year", 
              "Pop"    = "pop")) 
  


load(here::here("Data/Derived/ABC_W_by_kap_GEEs.Rdata"))
load(here::here("Data/Derived/gee_results.Rdata"))

# Estimate Mating probability under different circumstances -------------
test_ws <- exp_seq(1e-3,100,200)


phi_pred <- as.data.frame(cbind("W" = test_ws,
                                "Phi_k01" = sapply(test_ws, phi_Wk, k = 0.1),
                                "Phi_dynak_E" = mapply(phi_Wk, test_ws, zanz_kappa_E_gee(test_ws)),
                                "Phi_dynak_W1" = mapply(phi_Wk, test_ws, case1_gee_fx(test_ws)),
                                "Phi_dynak_W3" = mapply(phi_Wk, test_ws, case3_gee_fx(test_ws)),
                                "Phi_dynak_W4" = mapply(phi_Wk, test_ws, case4_gee_fx(test_ws))))
#Analytic Mating probabilities ------------
mate_prob_lines <- phi_pred %>% 
  pivot_longer(Phi_k01:Phi_dynak_W4) %>% 
  # mutate(kap_dyna = if_else(grepl("dynak", name), "Dynamic", "Constant"),
  #        Case1_2 = if_else(grepl("case2", name), "Case2", "Case1")) %>% 
  # filter(Case1_2 == "Case1") %>% 
  ggplot() +
  geom_line(aes(x = W, y = value, col = name),
            size = 1.2) +
  theme_classic() +
  #theme(legend.position = "none") +
  scale_x_continuous(trans = "log",
                     breaks = c(0.01, 0.1, 1, 10, 100),
                     labels = c("0.01", "0.1", "1", "10", "100"),
                     limits = c(0.01, 200)) +
  ylim(c(0,1)) +
  labs(y = expression(Mating~Probability~(Phi(W, kappa))),
       x = expression(Mean~worm~burden~(W)))

mate_prob_lines








# Case1 mate prob plot ----------------------
abc_case1_mate_prob <- abc_post_preds_mateProb %>% 
  filter(Case=="case1" & Pop == "Comm") %>% 
  ggplot(aes(x = UF_mean,
             y = q5, 
             size = 1/(IQR))) +
  geom_point(col = "#3b46ca",
             alpha = 0.3,
             show.legend = FALSE) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        # axis.text.x = element_blank(),
        # axis.title.y = element_blank(),
        # plot.title = element_text(hjust = 0.5),
        legend.position = "none") +
  scale_size(range = c(1,4)) +
  ylim(c(0,1)) +
  scale_x_continuous(trans = "log",
                     breaks = c(0.01, 0.1, 1, 10, 100),
                     labels = c("0.01", "0.1", "1", "10", "100")) +
  labs(title = "Case 1",
       tag = "A")

# Case2 mate prob plot ----------------------
abc_case2_mate_prob <- abc_post_preds_mateProb %>% 
  filter(Case=="case2" & Pop == "Comm") %>% 
  ggplot(aes(x = UF_mean,
             y = q5, 
             size = 1/(IQR))) +
  geom_point(col = "#fc45b7",
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
                     labels = c("0.01", "0.1", "1", "10", "100"))+
  labs(y = expression(Mating~probability~(Phi[st])),
       title = "Case 2",
       tag = "B")

# Case3 mate prob plot ----------------------
abc_case3_mate_prob <- abc_post_preds_mateProb %>% 
  filter(Case=="case3" & Pop == "Comm") %>% 
  ggplot(aes(x = UF_mean,
             y = q5, 
             size = 1/(IQR))) +
  geom_point(col = "#bc86af",
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
                     labels = c("0.01", "0.1", "1", "10", "100"))+
  labs(x = expression(Mean~worm~burden~(W[st])),
       title = "Case 3",
       tag = "C")

# Case4 mate prob plot ----------------------
abc_case3_mate_prob <- abc_post_preds_mateProb %>% 
  filter(Case=="case4" & Pop == "Comm") %>% 
  ggplot(aes(x = UF_mean,
             y = q5, 
             size = 1/(IQR))) +
  geom_point(col = "#009d23",
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
                     labels = c("0.01", "0.1", "1", "10", "100"))+
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
  stat_smooth(data = abc_fin_df_case_long %>% filter(Case %in% c("Case1", "Case2", "Case3", "Case4") & pop == "Comm"),
              aes(x = obsW.med, y = obsPhiProbW.med, col = Case)) +
  scale_color_manual(values = c("#058dd6", "#cc4a49", "#7b6c7c", "#64a860")) +
  theme(plot.title = element_text(size = 14,hjust = 0.5))# + labs(title = "Observed and Estimated\nMating Probability")

png(here::here("Figures/Fig3_ABC_MateProb_Results.png"),
    height = 5.5, width = 7, units = "in", res = 300)

(abc_case1_mate_prob | abc_case2_mate_prob) / (abc_case3_mate_prob | abc_case4_mate_prob) | abc_w_mate_prob_lines# + plot_layout(widths = c(3,3,2))

dev.off()

pdf(here::here("Figures/Fig3_ABC_MateProb_Results.pdf"),
    height = 5.5, width = 7)

(abc_case1_mate_prob | abc_case2_mate_prob) / (abc_case3_mate_prob | abc_case4_mate_prob) | abc_w_mate_prob_lines# + plot_layout(widths = c(3,2))

dev.off()
