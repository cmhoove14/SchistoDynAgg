library(tidyverse)
library(abc)

devtools::load_all()

abc_sims <- readRDS(here::here("Data/Derived/abc_fit_ridge_log_4cases_hrg_fixed_1e+05iterations.rds"))
yO <- readRDS(here::here("Data/Derived/ABC_yO_data.rds"))

n_post_pred_runs <- 1000   # Posterior draws
post_pred_quants <- c(0.025,0.25,0.5,0.75,0.975) #Quantiles for summaries

# Get posterior summaries ----------------------
abc_successes <- sapply(1:length(abc_sims), function(s){
  ifelse(length(abc_sims[[s]]) == 4, NA, s)
})

abc_sims2 <- abc_sims[na.omit(abc_successes)]

abc_post_pred_checks <- bind_rows(lapply(1:length(abc_sims2), function(sim){
  abc_run <- abc_sims2[[sim]]
     <- yO$n_ppl[sim]
  
  shehia <- as.character(abc_run[[1]])
  year   <- as.integer(abc_run[[2]])
  pop    <- as.character(abc_run[[3]])
  
  n_ppl <- yO %>% 
    filter(Shehia == shehia,
           Year == year,
           pop == pop) %>% 
    pull(n_ppl)
  
    # Posterior predictions for case 1 runs
      abc_post_pred_case1 <- post_pred_data_gen(pars        = abc_run[[4]]$adj.values, 
                                                fixed_pars  = c("h"=10, "r" = 1, "g" = 0.001),
                                                n_ppl       = n_ppl,
                                                weights     = abc_run[[4]]$weights,
                                                data_gen_fx = gen_case1_data,
                                                n_reps      = n_post_pred_runs)
      
      abc_post_pred_case1_sum <- c(matrixStats::colQuantiles(abc_post_pred_case1, probs = post_pred_quants))
      names(abc_post_pred_case1_sum) <- paste0("case1_", 
                                               rep(c("E", "E.se", "E.pos2n"), times = length(post_pred_quants)),
                                               rep(post_pred_quants, each = 3))
      
      # Posterior predictions for case 2 runs
      abc_post_pred_case2 <- post_pred_data_gen(pars        = abc_run[[5]]$adj.values, 
                                                fixed_pars  = c("h"=10, "r" = 1, "g" = 0.001),
                                                n_ppl       = n_ppl,
                                                weights     = abc_run[[5]]$weights,
                                                data_gen_fx = gen_case2_data,
                                                n_reps      = n_post_pred_runs)
      
      abc_post_pred_case2_sum <- c(matrixStats::colQuantiles(abc_post_pred_case2, probs = post_pred_quants))
      names(abc_post_pred_case2_sum) <- paste0("case2_", 
                                               rep(c("E", "E.se", "E.pos2n"), times = length(post_pred_quants)),
                                               rep(post_pred_quants, each = 3))
      
      # Posterior predictions for case 3 runs
      abc_post_pred_case3 <- post_pred_data_gen(pars        = abc_run[[6]]$adj.values, 
                                                fixed_pars  = c("susc_shape" = 1, "susc_rate" = 1/500,
                                                                "h"=10, "r" = 1, "g" = 0.001),
                                                n_ppl       = n_ppl,
                                                weights     = abc_run[[6]]$weights,
                                                data_gen_fx = gen_case3_data,
                                                n_reps      = n_post_pred_runs)
      
      abc_post_pred_case3_sum <- c(matrixStats::colQuantiles(abc_post_pred_case3, probs = post_pred_quants))
      names(abc_post_pred_case3_sum) <- paste0("case3_", 
                                               rep(c("E", "E.se", "E.pos2n"), times = length(post_pred_quants)),
                                               rep(post_pred_quants, each = 3))
      
      # Posterior predictions for case 4 runs
      abc_post_pred_case4 <- post_pred_data_gen(pars        = abc_run[[7]]$adj.values, 
                                                fixed_pars  = c("h"=10, "r" = 1, "g" = 0.001),
                                                n_ppl       = n_ppl,
                                                weights     = abc_run[[7]]$weights,
                                                data_gen_fx = gen_case4_data,
                                                n_reps      = n_post_pred_runs)
      
      abc_post_pred_case4_sum <- c(matrixStats::colQuantiles(abc_post_pred_case4, probs = post_pred_quants))
      names(abc_post_pred_case4_sum) <- paste0("case4_", 
                                               rep(c("E", "E.se", "E.pos2n"), times = length(post_pred_quants)),
                                               rep(post_pred_quants, each = 3))
      
      
      df.out <- cbind(data.frame("Shehia" = shehia, 
                                 "Year" = year, 
                                 "Pop" = pop),
                      as.data.frame(t(c(abc_post_pred_case1_sum,
                                        abc_post_pred_case2_sum,
                                        abc_post_pred_case3_sum,
                                        abc_post_pred_case4_sum))))
  
  return(df.out)
  
}))
