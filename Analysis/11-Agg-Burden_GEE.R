# -----------------------------
# SchistoDynAgg Estimate aggregation by burden relationship with GEE
# Chris Hoover
# -----------------------------

# Setup ---------------
library(tidyverse)
devtools::load_all()


# Load data ---------------
# Filter out pilot observations of 50 children in year 2011 and only focus on study years 2012-2017
comm_sums <- readRDS(here::here("Data", "Derived", "adults&chldrn_shehia_sums_Unguja_Pemba_2012_2017.rds")) %>% 
  filter(Year != 2011)

adlt_sums <- readRDS(here::here("Data", "Derived", "adults_shehia_sums_Unguja_Pemba_2012_2017.rds")) %>% 
  filter(Year != 2011)

chld_sums <- readRDS(here::here("Data", "Derived", "chldrn_shehia_sums_Unguja_Pemba_2012_2017.rds")) %>% 
  filter(Year != 2011)


# Estimate agg-buirden relationship for all community observations ------------
alpha_gee <- geeglm(log(UF_alpha_mle) ~ log(UF_mean), id = as.factor(Shehia),
                    family = "gaussian", corstr = "unstructured",
                    weights = 1/UF_alpha_mle_se,
                    data = comm_sums %>% filter(!is.na(UF_alpha_mle_se)))

gee_coef1 <- coef(alpha_gee)[1]
gee_coef2 <- coef(alpha_gee)[2]

zanz_kappa_E_gee <- function(E){
  exp(gee_coef1+gee_coef2*log(E))^-1
}

kap_e_gee_plot <- comm_sums %>%
  ggplot(aes(x     = UF_mean,
             y     = UF_alpha_mle^-1, 
             ymin  = (UF_alpha_mle - UF_alpha_mle_se)^-1, 
             ymax  = (UF_alpha_mle + UF_alpha_mle_se)^-1,
             shape = Intervention,
             col   = Isl,
             size  = 1/UF_alpha_mle_se)) +
    geom_point(alpha = 0.6) +
    #geom_errorbar(alpha = 0.6) +
    #facet_grid(Intervention~Isl) +
    theme_classic() +
    scale_x_continuous(trans = "log",
                       breaks = c(0.01, 0.1, 1, 10, 100),
                       labels = c("0.01", "0.1", "1", "10", "100"),
                       limits = c(0.01, 100)) +
    scale_y_continuous(breaks = c(0,0.025,0.05,0.075),
                       limits = c(0,0.075)) +
    stat_function(fun = zanz_kappa_E_gee, size = 1.2, col = "black") +
    scale_color_manual(values = c("#d0885f", "#2d474c")) +
    scale_shape_manual(values = c(17,18,15)) +
    labs(y = expression(Egg~Aggregation~Parameter~(kappa[st]^E)),
         x = expression(Mean~egg~burden~(E[st])))

# kap_e_gee_plot




# Estimate agg-burden relationship for all community observations, with island and intervention effects ------------

alpha_gee_adj <- geeglm(log(UF_alpha_mle) ~ log(UF_mean) + Isl + Intervention, id = as.factor(Shehia),
                    family = "gaussian", corstr = "unstructured",
                    weights = 1/UF_alpha_mle_se,
                    data = comm_sums %>% filter(!is.na(UF_alpha_mle_se)) %>% 
                      mutate(Intervention = factor(Intervention, levels = c("MDA", "Behaviour", "Snail"))))

gee_adj_sum <- broom::tidy(alpha_gee_adj) %>% 
  mutate(est.lo = estimate+qnorm(0.025)*std.error,
         est.hi = estimate+qnorm(0.975)*std.error)

save(chld_sums, adlt_sums, comm_sums,
     alpha_gee, gee_coef1, gee_coef2, zanz_kappa_E_gee, kap_e_gee_plot, alpha_gee_adj, gee_adj_sum,
     file = here::here("Data/Derived/gee_results.Rdata"))


