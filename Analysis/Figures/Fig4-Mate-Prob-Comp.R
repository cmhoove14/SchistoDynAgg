library(tidyverse)
library(geepack)
library(grid)
library(gridExtra)
library(patchwork)

devtools::load_all()

yO <- readRDS(here::here("Data/Derived/ABC_yO_data.rds"))

abc_post_preds_MateProb <- readRDS(here::here("Data/Derived/abc_post_pred_checks.rds")) %>% 
  filter(SumStat %in% c("E", "E_se","E_pos2n", "mean_W", "kap_W", "prob_Phi")) %>% 
  pivot_wider(names_from = SumStat,
              values_from = q025:IQR,
              names_sep = "_") %>% 
  left_join(yO, 
            by = c("Isl"    = "Isl", 
                   "Shehia" = "Shehia", 
                   "Year"   = "Year", 
                   "Pop"    = "pop")) %>% 
  mutate(UF_prev     = UF_pos/n_ppl,
         E_mse       = (UF_mean - q5_E)^2/UF_mean,
         E_se_mse    = (UF_se - q5_E_se)^2/UF_se,
         E_pos2n_mse = (UFpos2n - q5_E_pos2n)^2/UFpos2n,
         MSE_sum     = E_mse+E_se_mse+E_pos2n_mse)

load(here::here("data/derived/ABC_W_by_kap_GEEs.Rdata"))
load(here::here("data/derived/gee_results.Rdata"))

# Analytic Mating probability estimates-------------
test_ws <- exp_seq(1e-3,100,200)

phi_opts <- expand.grid(Case = paste0("case", c(1:4)),
                        kap  = c("const", "dynaE", "dynaW"))

phi_preds <- bind_rows(lapply(1:nrow(phi_opts), function(x){
  dat <- abc_post_preds_MateProb %>% 
    filter(Pop == "Comm" & Case == x[1])
  
  kap_fx <- if_else(x[2] == "dynaE")
  
  if(x[2] == "const"){
    out_phi <- sapply(test_ws,
                      phi_Wk,
                      k = median(dat %>% pull(q5_kap_W)))
  } else if(x[2] == "dynaE"){
    out_phi <- mapply(phi_Wk,
                      test_ws,
                      zanz_kappa_W_gee_case1(test_ws))
    
  } else if(x[2] == "dynaW"){
    out_phi <- mapply(phi_Wk,
                      test_ws,
                      case1_gee_fx(test_ws))
  }
  
}))

phi_pred <- as.data.frame(cbind("W" = test_ws,
                                "Phi_c1_constk" = sapply(test_ws, phi_Wk, k = median(abc_post_preds_MateProb %>% 
                                                                                       filter(Pop == "Comm" & Case == "case1") %>% 
                                                                                       pull(q5_kap_W))),
                                "Phi_c1_dynak_E" = mapply(phi_Wk, test_ws, zanz_kappa_E_gee(test_ws)),
                                "Phi_c1_dynak_W" = mapply(phi_Wk, test_ws, case1_gee_fx(test_ws)),
                                "Phi_c2_constk" = sapply(test_ws, phi_wk_sep, k = median(abc_post_preds_MateProb %>% 
                                                                                       filter(Pop == "Comm" & Case == "case2") %>% 
                                                                                       pull(q5_kap_W))),
                                "Phi_c2_dynak_E" = mapply(phi_wk_sep, test_ws, zanz_kappa_E_gee(test_ws)),
                                "Phi_c2_dynak_W" = mapply(phi_wk_sep, test_ws, case2_gee_fx(test_ws)),
                                "Phi_c3_constk" = sapply(test_ws, phi_Wk, k = median(abc_post_preds_MateProb %>% 
                                                                                       filter(Pop == "Comm" & Case == "case3") %>% 
                                                                                       pull(q5_kap_W))),
                                "Phi_c3_dynak_E" = mapply(phi_Wk, test_ws, zanz_kappa_E_gee(test_ws)),
                                "Phi_c3_dynak_W" = mapply(phi_Wk, test_ws, case3_gee_fx(test_ws)),
                                "Phi_c4_constk" = sapply(test_ws, phi_Wk, k = median(abc_post_preds_MateProb %>% 
                                                                                       filter(Pop == "Comm" & Case == "case4") %>% 
                                                                                       pull(q5_kap_W))),
                                "Phi_c4_dynak_E" = mapply(phi_Wk, test_ws, zanz_kappa_E_gee(test_ws)),
                                "Phi_c4_dynak_W" = mapply(phi_Wk, test_ws, case4_gee_fx(test_ws))
                                
))


# Mating probability plots for individual cases -------------
# Case 1
abc_post_preds_MateProb1_plot <- abc_post_preds_MateProb %>% 
  filter(Case == "case1") %>% 
  ggplot(aes(x = q5_mean_W,
             y = q5_prob_Phi)) +
    geom_point(aes(size = 1/MSE_sum),
               col = "#3b46ca",
               alpha = 0.5,
               show.legend = FALSE) +
    geom_line(data = phi_pred %>% 
                dplyr::select(W:Phi_c1_dynak_W) %>% 
                rename("Constant k" = Phi_c1_constk,
                       "k(E)" = Phi_c1_dynak_E,
                       "k(W)" = Phi_c1_dynak_W) %>% 
                pivot_longer(`Constant k`:`k(W)`,
                             names_to = "Aggregation Function",
                             values_to = "Phi"),
              aes(x = W, y = Phi, lty = `Aggregation Function`)) +
    theme_classic() +
    theme(legend.position = c(0.8,0.25),
          legend.title = element_text(size = 9),
          legend.text = element_text(size = 7),
          plot.title = element_text(hjust = 0.1, vjust = -1)) +
    scale_size(range = c(1,4)) +
    ylim(c(0,1)) +
    scale_x_continuous(trans = "log",
                       breaks = c(0.01, 0.1, 1, 10, 100),
                       labels = c("0.01", "0.1", "1", "10", "100"),
                       limits = c(0.01, max(abc_post_preds_MateProb %>% 
                                              filter(Pop == "Comm") %>% 
                                              pull(q5_mean_W)))) +
  labs(y = expression(Mating~probability~(Phi[st])),
       lty = "Aggregation\nFunction",
       title = "A.   Case 1")
  





# Case 2
abc_post_preds_MateProb2_plot <- abc_post_preds_MateProb %>% 
  filter(Case == "case2") %>% 
  ggplot(aes(x = q5_mean_W,
             y = q5_prob_Phi)) +
  geom_point(aes(size = 1/MSE_sum),
             col = "#fc45b7",
             alpha = 0.5,
             show.legend = FALSE) +
  geom_line(data = phi_pred %>% 
              dplyr::select(W, Phi_c2_constk:Phi_c2_dynak_W) %>% 
              pivot_longer(Phi_c2_constk:Phi_c2_dynak_W,
                           names_to = "Aggregation Function",
                           values_to = "Phi"),
            aes(x = W, y = Phi, lty = `Aggregation Function`)) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 0.1, vjust = -1)) +
  scale_size(range = c(1,4)) +
  ylim(c(0,1)) +
  scale_x_continuous(trans = "log",
                     breaks = c(0.01, 0.1, 1, 10, 100),
                     labels = c("0.01", "0.1", "1", "10", "100"),
                     limits = c(0.01, max(abc_post_preds_MateProb %>% 
                                            filter(Pop == "Comm") %>% 
                                            pull(q5_mean_W)))) +
  labs(title = "B.   Case 2")





# Case 3
abc_post_preds_MateProb3_plot <- abc_post_preds_MateProb %>% 
  filter(Case == "case3") %>% 
  ggplot(aes(x = q5_mean_W,
             y = q5_prob_Phi)) +
  geom_point(col = "#bc86af",
             alpha = 0.5,
             show.legend = FALSE) +
  geom_line(data = phi_pred %>% 
              dplyr::select(W, Phi_c3_constk:Phi_c3_dynak_W) %>% 
              pivot_longer(Phi_c3_constk:Phi_c3_dynak_W,
                           names_to = "Aggregation Function",
                           values_to = "Phi"),
            aes(x = W, y = Phi, lty = `Aggregation Function`)) +
  theme_classic() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.1, vjust = -1)) +
  scale_size(range = c(1,4)) +
  ylim(c(0,1)) +
  scale_x_continuous(trans = "log",
                     breaks = c(0.01, 0.1, 1, 10, 100),
                     labels = c("0.01", "0.1", "1", "10", "100"),
                     limits = c(0.01, max(abc_post_preds_MateProb %>% 
                                            filter(Pop == "Comm") %>% 
                                            pull(q5_mean_W)))) +
  labs(x = expression(Mean~worm~burden~(W[st])),
       title = "C.   Case 3")





#Case 4
abc_post_preds_MateProb4_plot <- abc_post_preds_MateProb %>% 
  filter(Case == "case4") %>% 
  ggplot(aes(x = q5_mean_W,
             y = q5_prob_Phi)) +
  geom_point(aes(size = 1/MSE_sum),
             col = "#009d23",
             alpha = 0.5,
             show.legend = FALSE) +
  geom_line(data = phi_pred %>% 
              dplyr::select(W, Phi_c4_constk:Phi_c4_dynak_W) %>% 
              pivot_longer(Phi_c4_constk:Phi_c4_dynak_W,
                           names_to = "Aggregation Function",
                           values_to = "Phi"),
            aes(x = W, y = Phi, lty = `Aggregation Function`)) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 0.1, vjust = -1)) +
  scale_size(range = c(1,4)) +
  ylim(c(0,1)) +
  scale_x_continuous(trans = "log",
                     breaks = c(0.01, 0.1, 1, 10, 100),
                     labels = c("0.01", "0.1", "1", "10", "100"),
                     limits = c(0.01, max(abc_post_preds_MateProb %>% 
                                            filter(Pop == "Comm") %>% 
                                            pull(q5_mean_W)))) +
  labs(x = expression(Mean~worm~burden~(W[st])),
       title = "D.   Case 4")







#Combine individual cases in final plot ----------

png(here::here("Figures/Fig4_Mate_Probs.png"),
    height = 7, width = 7, units = "in", res = 300)

grid.arrange(arrangeGrob(abc_post_preds_MateProb1_plot + theme(axis.title = element_blank()), 
                         abc_post_preds_MateProb2_plot + theme(axis.title = element_blank()),
                         abc_post_preds_MateProb3_plot + theme(axis.title = element_blank()),
                         abc_post_preds_MateProb4_plot + theme(axis.title = element_blank()), 
                         nrow = 2, ncol = 2,
                         left = textGrob(expression(Mating~probability~(Phi[st])), rot = 90, vjust = 1),
                         bottom = textGrob(expression(Mean~worm~burden~(W[st])), vjust = 0.2)))

dev.off()




pdf(here::here("Figures/Fig4_Mate_Probs.pdf"),
    height = 7, width = 7)

grid.arrange(arrangeGrob(abc_post_preds_MateProb1_plot + theme(axis.title = element_blank()), 
                         abc_post_preds_MateProb2_plot + theme(axis.title = element_blank()),
                         abc_post_preds_MateProb3_plot + theme(axis.title = element_blank()),
                         abc_post_preds_MateProb4_plot + theme(axis.title = element_blank()), 
                         nrow = 2, ncol = 2,
                         left = textGrob(expression(Mating~probability~(Phi[st])), rot = 90, vjust = 1),
                         bottom = textGrob(expression(Mean~worm~burden~(W[st])), vjust = 0.2)))

dev.off()


