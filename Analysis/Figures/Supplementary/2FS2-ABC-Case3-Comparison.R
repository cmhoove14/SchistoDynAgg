library(grid)
library(gridExtra)
library(patchwork)

devtools::load_all()

load(here::here("Data/Derived/abc_processed_results.Rdata"))


# Add first principle component then do terrible data cleaning
sum_stat_pca <- prcomp(abc_fin_df2 %>% dplyr::select(UF_mean, UF_se, UFpos2n))

abc_fin_df2_pca <- abc_fin_df2 %>% 
  mutate(Case3_pca1st = -sum_stat_pca$x[,1]+1)

abc_fin_df_long_pca <- abc_fin_df2_pca %>% 
  pivot_longer(cols = Case1_W.med:Case3_pca1st,
               names_to = c("Case", "Measure"),
               names_sep = "_")

abc_fin_df_case_long_pca <- abc_fin_df_long_pca %>% 
  pivot_wider(names_from = "Measure", 
              values_from = "value",
              values_fn = {first})

case3_comp12_df <- abc_fin_df_case_long_pca %>% filter(closer > 0) %>% 
  pivot_longer(cols = c("UF_mean", "UF_se", "UFpos2n", "pca1st"),
               names_to = "Sum Stat")

abc_case3_compcase12 <- case3_comp12_df %>% 
  dplyr::select(Isl,Shehia, Year, `Sum Stat`, value, closer) %>% 
  #filter(`Sum Stat` != "pca1st") %>% 
  mutate(`Sum Stat` = case_when(`Sum Stat` == "UF_mean" ~ "Mean Egg burden",
                                `Sum Stat` == "UF_se" ~ "Std. Error Egg burden",
                                `Sum Stat` == "UFpos2n" ~ "Adjusted Egg prevalence",
                                `Sum Stat` == "pca1st" ~ "Normalized 1st PC")) %>% 
  ggplot(aes(x = as.factor(`Sum Stat`), y = value, fill = as.factor(closer))) +
  geom_boxplot() +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_y_continuous(trans = "log", breaks = c(0.001, 0.01, 0.1, 1, 10, 100)) +
  scale_fill_manual(values = c("#058dd6", "#cc4a49"),
                    labels = c("Case 1", "Case 2")) +
  labs(x = "Summary Statistic", y = "Value",
       fill = "Case 3\nComparison")

abc_case3_compcase12

ggsave(here::here("Figures/Supp2-abc_case3_comp_to_case1and2_by_sumstats.png"),
       height = 5, width = 8, units = "in")

#Distribution of summary statistics used in approximate bayesian computation estimation of worm burdens from egg burdens and their first principle component, stratified by whether the Case 3 worm burden estimates were closer to the Case 1 (blue) or Case 2 (red) estimates. This demonstrates that Case 2 dynamics are more likely to be estimated at lower parasite burdens, prevalences, and standard errors--indicative of lower overall transmission--while Case 1 dynamics are recovered in higher transmission setting