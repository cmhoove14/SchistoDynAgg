library(grid)
library(gridExtra)
library(patchwork)

devtools::load_all()

yO <- readRDS(here::here("Data/Derived/ABC_yO_data.rds"))

abc_post_preds_SumStats <- readRDS(here::here("Data/Derived/abc_post_pred_checks.rds")) %>% 
  filter(SumStat %in% c("E", "E_se","E_pos2n"))
  



  
abc_post_preds_SumStats_wide <- abc_post_preds_SumStats %>% 
  #mutate(row = row_number()) %>% 
  pivot_wider(names_from = SumStat,
              values_from = q025:IQR,
              names_sep = "_") %>% 
  left_join(yO, 
            by = c("Isl"    = "Isl", 
                   "Shehia" = "Shehia", 
                   "Year"   = "Year", 
                   "Pop"    = "pop")) %>% 
  mutate(E_mse = (UF_mean - q5_E)^2/n_ppl,
         E_se_mse = (UF_se - q5_E_se)^2/n_ppl,
         E_pos2n_mse = (UFpos2n - q5_E_pos2n)^2/n_ppl)







abc_obs_gen_comp <- abc_post_preds_SumStats_wide %>% 
  dplyr::select(Isl, Shehia, Year, Pop, Case, UF_mean, q5_E, q25_E, q75_E) %>% 
  rename("Obs" = UF_mean,
         "GenMed" = q5_E,
         "GenHiq" = q75_E,
         "GenLoq" = q25_E) %>% 
  mutate(SumStat = "Mean Egg Burden") %>% 
  bind_rows(abc_post_preds_SumStats_wide %>% 
              dplyr::select(Isl, Shehia, Year, Pop, Case, UF_se, q5_E_se, q25_E_se, q75_E_se) %>% 
              rename("Obs" = UF_se,
                     "GenMed" = q5_E_se,
                     "GenHiq" = q75_E_se,
                     "GenLoq" = q25_E_se) %>% 
              mutate(SumStat = "Std. Er Egg Burden")) %>% 
  bind_rows(abc_post_preds_SumStats_wide %>% 
              dplyr::select(Isl, Shehia, Year, Pop, Case, UFpos2n, q5_E_pos2n, q25_E_pos2n, q75_E_pos2n) %>% 
              rename("Obs" = UFpos2n,
                     "GenMed" = q5_E_pos2n,
                     "GenHiq" = q75_E_pos2n,
                     "GenLoq" = q25_E_pos2n) %>% 
              mutate(SumStat = "Adjusted Prevalence"))

abc_obs_gen_comp_plot <- abc_obs_gen_comp %>% 
  ggplot(aes(x = Obs, y = GenMed, col = Case,
             ymin = GenLoq, ymax = GenHiq)) +
  geom_point(alpha = 0.5) +
  #geom_errorbar(alpha = 0.5) +
  geom_abline(slope = 1, lty = 2) +
  scale_x_continuous(trans = "log1p",
                     breaks = c(0, 1, 10, 100),
                     limits = c(0,100),
                     labels = c("0", "1", "10", "100")) +
  scale_y_continuous(trans = "log1p",
                     breaks = c(0, 1, 10, 100),
                     limits = c(0,100),
                     labels = c("0", "1", "10", "100")) +
  scale_color_manual(values = c("#3b46ca", "#fc45b7", "#bc86af", "#009d23")) +
  facet_grid(Case~SumStat) +
  theme_classic() +
  theme(legend.position = "none",
        strip.background = element_rect(fill = NULL)) +
  labs(x = "Observed value (log+1 transformed)",
       y = "Generated value (log+1 transformed)")

abc_obs_gen_comp_plot

#"Comparison of generated to observed summary statistics used in approximate bayesian computation estimation of community parasite burdens. Colors indicate the data generating Case and the 1:1 line implying perfect agreement between observed and generated data is shown. Error bars correspond to interquartile ranges of the generated summary statistics from parameter sets included in the posterior distribution."


