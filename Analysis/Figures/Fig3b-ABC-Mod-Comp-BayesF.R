library(tidyverse)
library(geepack)
library(grid)
library(gridExtra)
library(patchwork)

devtools::load_all()

abc_bayesFs <- readRDS(here::here("Data/Derived/abc_mod_comps.rds"))
yO          <- readRDS(here::here("Data/Derived/ABC_yO_data.rds"))

abc_bayes_comps <- abc_bayesFs %>% 
  rename("Isl_Shehia" = Shehia) %>% 
  left_join(yO %>% mutate(Isl_Shehia = paste(Isl, Shehia,sep = "_")),
            by = c("Isl_Shehia" = "Isl_Shehia",
                   "Year" = "Year",
                   "Pop" = "pop")) %>% 
  mutate(UF_prev = UF_pos/n_ppl)

# Bayes factors (Comparative estimate of model performance) in reference to case 1 -------------
abc_bayesF_comp1 <- abc_bayes_comps %>% 
  filter(Pop == "Comm") %>% 
  dplyr::select(c("Isl","Shehia", "Year", "UF_mean", "UF_prev",
                  "IItoI_BayesF", "IIItoI_BayesF", "IVtoI_BayesF")) %>% 
  pivot_longer(cols      = c("IItoI_BayesF", "IIItoI_BayesF", "IVtoI_BayesF"),
               names_to  = "Case",
               values_to = "BayesF")

abc_bayesF_comp1_plot <- abc_bayesF_comp1 %>% 
  ggplot(aes(x = UF_mean, y = BayesF, col = Case)) +
    stat_smooth() +
    geom_hline(yintercept = 1) +
    geom_point(alpha = 0.7) +
    scale_y_continuous(trans = "log1p",
                       breaks = c(0.1, 1, 10, 100, 1000),
                       labels = c("0.1", "1", "10", "100", "1000"),
                       limits = c(0.1,max(abc_bayesF_comp1$BayesF))) +
    scale_x_continuous(trans = "log",
                       breaks = c(0.01, 0.1, 1, 10, 100),
                       labels = c("0.01", "0.1", "1", "10", "100")) +
    scale_color_manual(values = c("IItoI_BayesF" = "#fc45b7", 
                                  "IIItoI_BayesF" = "#bc86af", 
                                  "IVtoI_BayesF" = "#009d23"),
                       labels = c("Case 3 to Case 1", "Case 2 to Case 1","Case 4 to Case 1")) +
    theme_classic() +
    theme(legend.position = c(0.2, 0.85),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 10)) +
    labs(x = expression(paste("Mean egg burden (\U2130"[st],")")),
         y = "Bayes Factor",
         col = "Comparison")

abc_bayesF_comp1_plot

ggsave(here::here("Figures/ABC_Bayes_Factors.png"),
       height = 5, width = 7, units = "in")

