library(tidyverse)
library(geepack)
library(grid)
library(gridExtra)
library(patchwork)

devtools::load_all()

yO <- readRDS(here::here("Data/Derived/ABC_yO_data.rds"))

abc_post_preds_kapW <- readRDS(here::here("Data/Derived/abc_post_pred_checks.rds")) %>% 
  filter(SumStat %in% c("E", "E_se","E_pos2n", "mean_W", "var_W","kap_W")) %>% 
  pivot_wider(names_from = SumStat,
              values_from = q025:IQR,
              names_sep = "_") %>% 
  left_join(yO, 
            by = c("Isl"    = "Isl", 
                   "Shehia" = "Shehia", 
                   "Year"   = "Year", 
                   "Pop"    = "pop")) %>% 
  mutate(UF_prev     = UF_pos/n_ppl,
         E_mse       = (UF_mean - q5_E)^2/n_ppl,
         E_se_mse    = (UF_se - q5_E_se)^2/n_ppl,
         E_pos2n_mse = (UFpos2n - q5_E_pos2n)^2/n_ppl,
         MSE_sum     = E_mse+E_se_mse+E_pos2n_mse)



# Case1 PLot ----------------------
case1_gee <- geeglm(log(q5_kap_W^-1) ~ log(q5_mean_W), id = as.factor(Shehia),
                    family = "gaussian", corstr = "unstructured",
                    weights = 1/MSE_sum,
                    data = abc_post_preds_kapW  %>% filter(Case=="case1" & Pop == "Comm"))

case1_gee_coef1 <- coef(case1_gee)[1]
case1_gee_coef2 <- coef(case1_gee)[2]

case1_gee_fx <- function(W){
  exp(case1_gee_coef1+case1_gee_coef2*log(W))^-1
}


kap_W_case1_plot <- abc_post_preds_kapW %>%
  filter(Case=="case1" & Pop == "Comm") %>% 
  ggplot(aes(x = q5_mean_W,
             y = q5_kap_W)) +
  geom_point(aes(size = 1/MSE_sum),
             col = "#3b46ca",
             alpha = 0.3,
             show.legend = FALSE) +
  #stat_smooth(col = "black") +
  stat_function(fun = case1_gee_fx, size = 1, lty = 2, col = "black", alpha = 0.5) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.1, vjust = -1)) +
  scale_size(range = c(1,4)) +
  scale_x_continuous(trans = "log",
                     breaks = c(0.01, 0.1, 1, 10, 100),
                     labels = c("0.01", "0.1", "1", "10", "100"),
                     limits = c(0.005, max(abc_post_preds_kapW$q5_mean_W))) +
  scale_y_continuous(trans = "log",
                     breaks = c(0.001,0.01, 0.1, 1, 10, 100),
                     labels = c("0.001","0.01", "0.1", "1", "10", "100"),
                     limits = c(0.001,1000)) +
  labs(y = expression(Worm~Dispersion~Parameter~(kappa[st]^W)),
       title = "A.   Case1")

# Case2 PLot ----------------------
case2_gee <- geeglm(log(q5_kap_W^-1) ~ log(q5_mean_W), id = as.factor(Shehia),
                    family = "gaussian", corstr = "unstructured",
                    weights = 1/MSE_sum,
                    data = abc_post_preds_kapW  %>% filter(Case=="case2" & Pop == "Comm"))

case2_gee_coef1 <- coef(case2_gee)[1]
case2_gee_coef2 <- coef(case2_gee)[2]

case2_gee_fx <- function(W){
  exp(case2_gee_coef1+case2_gee_coef2*log(W))^-1
}


kap_W_case2_plot <- abc_post_preds_kapW %>%
  filter(Case=="case2" & Pop == "Comm") %>% 
  ggplot(aes(x = q5_mean_W,
             y = q5_kap_W)) +
  geom_point(aes(size = 1/MSE_sum),
             col = "#fc45b7",
             alpha = 0.3,
             show.legend = FALSE) +
  #stat_smooth(col = "black") +
  stat_function(fun = case2_gee_fx, size = 1, lty = 2, col = "black", alpha = 0.5) +
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
                     limits = c(0.005, max(abc_post_preds_kapW$q5_mean_W))) +
  scale_y_continuous(trans = "log",
                     breaks = c(0.001,0.01, 0.1, 1, 10, 100),
                     labels = c("0.001","0.01", "0.1", "1", "10", "100"),
                     limits = c(0.001,1000)) +
  labs(title = "B.   Case 2")

# Case3 PLot ----------------------
case3_gee <- geeglm(log(q5_kap_W^-1) ~ log(q5_mean_W), id = as.factor(Shehia),
                    family = "gaussian", corstr = "unstructured",
                    weights = 1/MSE_sum,
                    data = abc_post_preds_kapW  %>% filter(Case=="case3" & Pop == "Comm" & q5_mean_W > 0))

case3_gee_coef1 <- coef(case3_gee)[1]
case3_gee_coef2 <- coef(case3_gee)[2]

case3_gee_fx <- function(W){
  exp(case3_gee_coef1+case3_gee_coef2*log(W))^-1
}

kap_W_case3_plot <- abc_post_preds_kapW %>%
  filter(Case=="case3" & Pop == "Comm") %>% 
  ggplot(aes(x = q5_mean_W,
             y = q5_kap_W)) +
  geom_point(aes(size = 1/MSE_sum),
             col = "#bc86af",
             alpha = 0.3,
             show.legend = FALSE) +
  #stat_smooth(col = "black") +
  stat_function(fun = case3_gee_fx, size = 1, lty = 2, col = "black", alpha = 0.5) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.1, vjust = -1)) +
  scale_size(range = c(1,4)) +
  scale_x_continuous(trans = "log",
                     breaks = c(0.01, 0.1, 1, 10, 100),
                     labels = c("0.01", "0.1", "1", "10", "100"),
                     limits = c(0.005, max(abc_post_preds_kapW$q5_mean_W))) +
  scale_y_continuous(trans = "log",
                     breaks = c(0.001,0.01, 0.1, 1, 10, 100),
                     labels = c("0.001","0.01", "0.1", "1", "10", "100"),
                     limits = c(0.001,1000)) +
  labs(x = expression(Mean~worm~burden~(W[st])),
       y = expression(Worm~Dispersion~Parameter~(kappa[st]^W)),
       title = "C.   Case 3")

# Case4 PLot ----------------------
case4_gee <- geeglm(log(q5_kap_W^-1) ~ log(q5_mean_W), id = as.factor(Shehia),
                    family = "gaussian", corstr = "unstructured",
                    weights = 1/MSE_sum,
                    data = abc_post_preds_kapW  %>% filter(Case=="case4" & Pop == "Comm"))

case4_gee_coef1 <- coef(case4_gee)[1]
case4_gee_coef2 <- coef(case4_gee)[2]

case4_gee_fx <- function(W){
  exp(case4_gee_coef1+case4_gee_coef2*log(W))^-1
}

kap_W_case4_plot <- abc_post_preds_kapW %>%
  filter(Case=="case4" & Pop == "Comm") %>% 
  ggplot(aes(x = q5_mean_W,
             y = q5_kap_W)) +
  geom_point(aes(size = 1/MSE_sum),
             col = "#009d23",
             alpha = 0.3,
             show.legend = FALSE) +
  #stat_smooth(col = "black") +
  stat_function(fun = case4_gee_fx, size = 1, lty = 2, col = "black", alpha = 0.5) +
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
                     limits = c(0.005, max(abc_post_preds_kapW$q5_mean_W))) +
  scale_y_continuous(trans = "log",
                     breaks = c(0.001,0.01, 0.1, 1, 10, 100),
                     labels = c("0.001","0.01", "0.1", "1", "10", "100"),
                     limits = c(0.001,1000)) +
  labs(x = expression(Mean~worm~burden~(W[st])),
       title = "D.   Case 4")

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

save(list = c("case1_gee_fx", "case1_gee_coef1", "case1_gee_coef2",
              "case2_gee_fx", "case2_gee_coef1", "case2_gee_coef2",
              "case3_gee_fx", "case3_gee_coef1", "case3_gee_coef2",
              "case4_gee_fx", "case4_gee_coef1", "case4_gee_coef2"),
     file = here::here("Data/Derived/ABC_W_by_kap_GEEs.Rdata"))




















# Miscellaneous checks -------------
  






abc_post_preds_kapW %>% 
  filter(Pop == "Comm") %>% 
  ggplot(aes(x = UF_mean, y = q5_kap_W, col = Case)) +
    geom_point() +
    theme_classic() +
    scale_x_continuous(trans = "log",
                       breaks = c(0.01,0.1,1,10,100)) +
  ylim(c(0,1)) +
  stat_smooth()








# Investigate cluster -------
# Plots all look like there are two distinct clusters, so want to look into what might be causing that
abc_post_preds_kapW_clust <- readRDS(here::here("Data/Derived/abc_post_pred_checks.rds")) %>% 
  filter(SumStat %in% c("E", "E_se","E_pos2n", "mean_W", "var_W","kap_W")) %>% 
  pivot_wider(names_from = SumStat,
              values_from = q025:IQR,
              names_sep = "_") %>% 
  left_join(yO, 
            by = c("Isl"    = "Isl", 
                   "Shehia" = "Shehia", 
                   "Year"   = "Year", 
                   "Pop"    = "pop")) %>% 
  mutate(UF_prev     = UF_pos/n_ppl,
         E_mse       = (UF_mean - q5_E)^2/n_ppl,
         E_se_mse    = (UF_se - q5_E_se)^2/n_ppl,
         E_pos2n_mse = (UFpos2n - q5_E_pos2n)^2/n_ppl,
         MSE_sum     = E_mse+E_se_mse+E_pos2n_mse) %>% 
  filter(Pop == "Comm" & Case == "case3") %>% 
  mutate(clust = if_else(q5_kap_W <0.11, 1, 0))



abc_post_preds_kapW_clust %>%
  ggplot(aes(x = q5_mean_W,
             y = q5_kap_W,
             col = clust)) +
  geom_point(aes(size = 1/MSE_sum),
             show.legend = FALSE) +
  #stat_smooth(col = "black") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.1, vjust = -1)) +
  scale_x_continuous(trans = "log",
                     breaks = c(0.01, 0.1, 1, 10, 100),
                     labels = c("0.01", "0.1", "1", "10", "100"),
                     limits = c(0.005, max(abc_post_preds_kapW$q5_mean_W))) +
  scale_y_continuous(trans = "log",
                     breaks = c(0.001,0.01, 0.1, 1, 10, 100),
                     labels = c("0.001","0.01", "0.1", "1", "10", "100"),
                     limits = c(0.001,1000)) +
  labs(y = expression(Worm~Dispersion~Parameter~(kappa[st]^W)),
       title = "A.   Case1")


boxplot(abc_post_preds_kapW_clust$UF_mean ~ abc_post_preds_kapW_clust$clust)
boxplot(abc_post_preds_kapW_clust$UF_prev ~ abc_post_preds_kapW_clust$clust)
boxplot(abc_post_preds_kapW_clust$UF_pos ~ abc_post_preds_kapW_clust$clust)
