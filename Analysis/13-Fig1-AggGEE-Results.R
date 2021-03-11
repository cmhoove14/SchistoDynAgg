library(tidyverse)
library(gridExtra)

load(here::here("Data/Derived/gee_results.Rdata"))
boot_out <- readRDS(here::here("Data/Derived/dispersion_change_IQR_boot10000.rds"))


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



# Combine panels for final figure -------------------------------------------

#extract legend as in https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

f1_legend<-g_legend(strat_kap_iqr)

png(here::here("Figures/Fig1_GEE_Results.png"),
               height = 7, width = 7, units = "in", res = 300)

grid.arrange(arrangeGrob(kap_e_gee_plot + theme(legend.position="none"),
                         strat_kap_iqr + theme(legend.position="none"),
                         nrow=2),
             f1_legend, widths=c(4, 1))

dev.off()


pdf(here::here("Figures/Fig1_GEE_Results.pdf"),
               height = 7, width = 7)

grid.arrange(arrangeGrob(kap_e_gee_plot + theme(legend.position="none"),
                         strat_kap_iqr + theme(legend.position="none"),
                         nrow=2),
             f1_legend, widths=c(4, 1))

dev.off()