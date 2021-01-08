# -----------------------------
# SchistoDynAgg Estimate aggregation by burden relationship with GEE
# Chris Hoover
# -----------------------------

# Setup ---------------
devtools::load_all()

nboot <- 1000   # Bootstrap samples for uncertainty intervals around GEE estimates
rerun <- FALSE  # Should bootstrap samples be rerun or loaded from previous run?


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
  ggplot(aes(x = UF_mean_mle,
             y = UF_alpha_mle^-1, 
             ymin = (UF_alpha_mle - UF_alpha_mle_se)^-1, ymax = (UF_alpha_mle + UF_alpha_mle_se)^-1,
             shape = Intervention,
             col = Isl,
             size = 1/UF_alpha_mle_se)) +
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
    labs(y = expression(Egg~Dispersion~Parameter~(kappa[st]^E)),
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

gee_adj_sum














# Bootstrap GEE Estimates for uncertainty ---------------------------------

# Data frame of all stratifications to estiamte over
boot_sweeps <- expand.grid(Island = c("Pemba", "Unguja", "All"),
                           Treatment = c("MDA", "Snail", "Behaviour", "All"),
                           Population = c("Children", "Adults", "All"))



# Run on biostat computing cluster, self-contained script and bash script kappa_change_bootstrap.R & kappa_change_bootstrap.sh in Analysis/Cluster_Scripts/GEE_Bootstrap
if(rerun) {
  cl <- makeCluster(parallel::detectCores() - 1)
  registerDoParallel(cl)
  
  boot_ests <- foreach(
    x = 1:nrow(boot_sweeps),
    .combine = cbind,
    .packages = c("tidyverse", "geepack"),
    .options.RNG = 430
  ) %dorng% {
    boot_df <-
      get_boot_df(boot_sweeps[x, 1], boot_sweeps[x, 2], boot_sweeps[x, 3])
    get_boot_ests(boot_df, nboot)
    
  }
  
  parallel::stopCluster(cl)
  
  boot_meds <- matrixStats::colMedians(boot_ests, na.rm = TRUE)
  boot_loqs <-
    matrixStats::colQuantiles(boot_ests, probs = 0.25, na.rm = TRUE)
  boot_hiqs <-
    matrixStats::colQuantiles(boot_ests, probs = 0.75, na.rm = TRUE)
  
  boot_out <- cbind(boot_sweeps, boot_meds, boot_loqs, boot_hiqs)
  
  saveRDS(boot_out,
          paste0("dispersion_change_IQR_boot", nboot, ".rds"))
  
} else {
  
  boot_out <- readRDS(here::here("Data", "Derived", paste0("dispersion_change_IQR_boot", nboot,".rds")))

}




boot_out$E_IQR <- apply(boot_out, 1, function(x){
  df <- get_boot_df(x[1], x[2], x[3])
  
  # 25th and 75th quantiles for IQR  
  egg_burden025 <- quantile(df %>% 
                              filter(!is.na(UF_alpha_mle_se)) %>% 
                              pull(UF_mean), 0.25)
  egg_burden075 <- quantile(df %>% 
                              filter(!is.na(UF_alpha_mle_se)) %>% 
                              pull(UF_mean), 0.75)

  return(as.numeric(egg_burden075-egg_burden025))
})

boot_all_est <- boot_out %>% 
  filter(Island == "All" & Treatment == "All" & Population == "All")

strat_kap_iqr <- boot_out %>% 
  filter(Treatment != "All" & Island != "All") %>% 
  ggplot() +
    theme_classic() +
    geom_point(aes(x = Population, y = boot_meds, 
                   col = Island, shape = Treatment),
               size = 2, position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(x = Population, 
                      col = Island, shape = Treatment,
                      ymin = boot_loqs, ymax = boot_hiqs),
                  width = 0.2, position = position_dodge(width = 0.5)) +
    geom_point(data = boot_all_est,
               aes(x = 0.5, y = boot_meds),
               col = "black", size = 2) +
    geom_errorbar(data = boot_all_est,
                  aes(x = 0.5, ymin = boot_loqs, ymax = boot_hiqs),
                  col = "black", width = 0.1) +
    scale_color_manual(values = c("#d0885f", "#2d474c")) +
    scale_shape_manual(values = c(17,18,15),
                       labels = c("MDA Only", "MDA+Snail Control", "MDA+Behavior")) +
    #scale_y_continuous(position = "right") +
    geom_hline(yintercept = 0, lty = 2) +
    labs(y = expression(Delta~kappa["IQR"]),
         shape = "Intervention")


strat_kap_iqr


# Combine panels for final figure -------------------------------------------

#extract legend as in https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

f1_legend<-g_legend(strat_kap_iqr)

png(here::here("Figures", "Fig1_GEE_Results.png"),
               height = 7, width = 7, units = "in", res = 300)

grid.arrange(arrangeGrob(kap_e_gee_plot + theme(legend.position="none"),
                         strat_kap_iqr + theme(legend.position="none"),
                         nrow=2),
             f1_legend, widths=c(4, 1))

dev.off()


pdf(here::here("Figures", "Fig1_GEE_Results.pdf"),
               height = 7, width = 7)

grid.arrange(arrangeGrob(kap_e_gee_plot + theme(legend.position="none"),
                         strat_kap_iqr + theme(legend.position="none"),
                         nrow=2),
             f1_legend, widths=c(4, 1))

dev.off()