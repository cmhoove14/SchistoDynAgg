library(grid)
library(gridExtra)
library(patchwork)

devtools::load_all()

yO <- readRDS(here::here("Data/Derived/ABC_yO_data.rds"))

abc_post_preds_SumStats <- readRDS(here::here("Data/Derived/abc_post_pred_checks.rds")) %>% 
  filter(SumStat %in% c("E", "E_se","E_pos2n"))
  
  
abc_post_preds_SumStats_wide <- abc_post_preds_SumStats %>% 
  pivot_wider(names_from = SumStat,
              values_from = q025:IQR,
              names_sep = "_") %>% 
  left_join(yO, 
            by = c("Isl"    = "Isl", 
                   "Shehia" = "Shehia", 
                   "Year"   = "Year", 
                   "Pop"    = "pop")) %>% 
  mutate(UF_prev     = UF_pos/n_ppl,
         prev_group  = cut(UF_prev,
                           c(0,0.01,0.025,0.05,0.1,1.0)),
         E_mse       = (UF_mean - q5_E)^2/UF_mean,
         E_se_mse    = (UF_se - q5_E_se)^2/UF_se,
         E_pos2n_mse = (UFpos2n - q5_E_pos2n)^2/UFpos2n)


abc_post_preds_SumStats_MSE_sum <- abc_post_preds_SumStats_wide %>% 
  pivot_longer(E_mse:E_pos2n_mse,
               names_to = "SumStat2",
               values_to = "MSE") %>% 
  mutate(SumStat = case_when(SumStat2 == "E_mse" ~ "Mean Egg Burden",
                              SumStat2 == "E_se_mse" ~ "Std. Er Egg Burden",
                              SumStat2 == "E_pos2n_mse" ~ "Adjusted Prevalence")) %>% 
  group_by(SumStat, Case) %>% 
  summarise(mean_mse = round(mean(MSE), 3),
            label = paste0("Avg. MSE = ", mean_mse))



# compare observed to generated summary statistics --------------
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
  filter(Pop == "Comm") %>% 
  ggplot(aes(x = Obs, y = GenMed)) +
  geom_point(aes(col = Case),
             alpha = 0.5) +
  #geom_errorbar(aes(ymin = GenLoq, ymax = GenHiq), alpha = 0.5) +
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
  geom_text(
    data    = abc_post_preds_SumStats_MSE_sum,
    mapping = aes(x = 3, y = 50, label = label),
    size = 2.5
  ) +
  theme_classic() +
  theme(legend.position = "none",
        strip.background = element_rect(fill = NULL)) +
  labs(x = "Observed value (log+1 transformed)",
       y = "Generated value (log+1 transformed)")

png(here::here("Figures/Fig3a_ABC_PostPred_Obs_Comp.png"),
    height = 6, width = 6, units = "in", res = 300)

abc_obs_gen_comp_plot

dev.off()

pdf(here::here("Figures/Fig3a_ABC_PostPred_Obs_Comp.pdf"),
    height = 6, width = 6)

abc_obs_gen_comp_plot

dev.off()





# Look at mean squared error by grouped prevalence -------------
abc_mse_by_prev_plot <- abc_post_preds_SumStats_wide %>% 
  filter(Pop == "Comm") %>% 
  pivot_longer(cols = E_mse:E_pos2n_mse,
               names_to = "MSE") %>% 
  ggplot(aes(x = prev_group, 
             y = value,
             fill = Case)) +
    geom_boxplot() +
    theme_classic() +
    scale_y_continuous(trans = "log1p") +
    facet_wrap(.~MSE, ncol = 3)
  
abc_mse_by_prev_plot
