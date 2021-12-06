library(tidyverse)

devtools::load_all()

kap_diffs <- readRDS(here::here("Data/Derived/abc_post_pred_checks.rds")) %>% 
  filter(SumStat == "kap_W") %>% 
  dplyr::select(-c(q025,q25,q75,q975, n_0s, IQR)) %>% 
  pivot_wider(names_from = Case, values_from = q5) %>% 
  mutate(case1_3_diff = case3 - case1) %>% 
  left_join(yO, 
            by = c("Isl"    = "Isl", 
                   "Shehia" = "Shehia", 
                   "Year"   = "Year", 
                   "Pop"    = "pop"))

kap_diff1_3_plot <- kap_diffs %>% 
  ggplot(aes(x = UF_mean, y = case1_3_diff)) + 
  geom_point(alpha = 0.5) +
  theme_classic() +
  ylim(c(-50,1)) +
  scale_x_continuous(trans = "log",
                     breaks = c(0.01, 0.1, 1, 10, 100),
                     labels = c("0.01", "0.1", "1", "10", "100"),
                     limits = c(0.005, max(kap_diffs$UF_mean))) +
  labs(x = "Mean community egg burden",
       y = expression(paste(kappa[st]^3, " - ", kappa[st]^1)))


png(here::here("Figures/Fig5_Kap_Diff_Case1_3.png"),
    height = 4, width = 4, units = "in", res = 300)

kap_diff1_3_plot

dev.off()

pdf(here::here("Figures/Fig5_Kap_Diff_Case1_3.pdf"),
    height = 4, width = 4)

kap_diff1_3_plot

dev.off()

