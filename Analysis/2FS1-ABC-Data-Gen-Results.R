library(grid)
library(gridExtra)
library(patchwork)

devtools::load_all()

load(here::here("Data/Derived/abc_processed_results.Rdata"))

abc_obs_gen_comp <- abc_fin_df_case_long %>% 
  filter(Case %in% c("Case1", "Case2", "Case3", "Case4")) %>% 
  dplyr::select(Isl, Shehia, Year, pop, Case, UF_mean, genE.med, genE.hiq, genE.loq) %>% 
  rename("Obs" = UF_mean,
         "GenMed" = genE.med,
         "GenHiq" = genE.hiq,
         "GenLoq" = genE.loq) %>% 
  mutate(SumStat = "Mean Egg Burden") %>% 
  bind_rows(abc_fin_df_case_long %>% 
              filter(Case %in% c("Case1", "Case2", "Case3", "Case4")) %>% 
              dplyr::select(Isl, Shehia, Year, pop, Case, UF_se, genEse.med, genEse.hiq, genEse.loq) %>% 
              rename("Obs" = UF_se,
                     "GenMed" = genEse.med,
                     "GenHiq" = genEse.hiq,
                     "GenLoq" = genEse.loq) %>% 
              mutate(SumStat = "Std. Er Egg Burden")) %>% 
  bind_rows(abc_fin_df_case_long %>% 
              filter(Case %in% c("Case1", "Case2", "Case3", "Case4")) %>% 
              dplyr::select(Isl, Shehia, Year, pop, Case, UFpos2n, genEpos2n.med, genEpos2n.hiq, genEpos2n.loq) %>% 
              rename("Obs" = UFpos2n,
                     "GenMed" = genEpos2n.med,
                     "GenHiq" = genEpos2n.hiq,
                     "GenLoq" = genEpos2n.loq) %>% 
              mutate(SumStat = "Adjusted Prevalence"))

abc_obs_gen_comp_plot <- abc_obs_gen_comp %>% 
  ggplot(aes(x = Obs, y = GenMed, col = Case,
             ymin = GenLoq, ymax = GenHiq)) +
  geom_point(alpha = 0.5) +
  geom_errorbar(alpha = 0.5) +
  geom_abline(slope = 1, lty = 2) +
  scale_x_continuous(trans = "log",
                     breaks = c(0.01, 0.1, 1, 10, 100),
                     labels = c("0.01", "0.1", "1", "10", "100")) +
  scale_y_continuous(trans = "log",
                     breaks = c(0.01, 0.1, 1, 10, 100),
                     labels = c("0.01", "0.1", "1", "10", "100")) +
  scale_color_manual(values = c("#058dd6", "#cc4a49", "#7b6c7c", "#64a860")) +
  facet_wrap(.~SumStat, ncol = 3) +
  theme_bw() +
  labs(x = "Observed value",
       y = "Generated value")

abc_obs_gen_comp_plot

ggsave(here("Figures/Supp1-abc_obs_gen_validation.png"),
       height = 5, width = 8, units = "in")

#"Comparison of generated to observed summary statistics used in approximate bayesian computation estimation of community parasite burdens. Colors indicate the data generating Case and the 1:1 line implying perfect agreement between observed and generated data is shown. Error bars correspond to interquartile ranges of the generated summary statistics from parameter sets included in the posterior distribution."