library(tidyverse)
library(geepack)
library(grid)
library(gridExtra)
library(patchwork)

devtools::load_all()

load(here::here("Data/Derived/abc_processed_results.Rdata"))


# Bayes factors (Comparative estimate of model performance) for each case across observed egg burden -------------
abc_bayesF_comp <- abc_fin_df_case_long %>% 
  filter(Case %in% c("IIItoI", "IIItoII", "IItoI",
                     "IVtoI", "IVtoII", "IVtoIII") & pop == "Child") %>% 
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
  # scale_color_manual(values = c("darkblue", "darkred", "forestgreen"),
  #                    labels = c("Case 3 to Case 1", "Case 3 to Case 2","Case 2 to Case 1")) +
  theme_classic() +
  theme(legend.position = c(0.15,0.8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10)) +
  labs(x = expression(paste("Mean egg burden (\U2130"[st],")")),
       y = "Bayes Factor",
       col = "Comparison",
       tag = "E")


abc_bayesF_comp

