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


n_cores <- parallel::detectCores()
cl <- parallel::makeCluster(n_cores)
registerDoParallel(cl)

clusterEvalQ(cl, devtools::load_all())

abc_post_pred_checks <- foreach(
  x = 1:length(abc_sims2),
  .packages = c("tidyverse", "abc"),
  .verbose = TRUE,
  .options.RNG = 7491
) %dorng% {
  
  abc_post_pred(abc_sims2, yO, n_post_pred_runs, x)

}  

saveRDS(abc_post_pred_checks,
        here::here("Data/Derived/abc_post_pred_checks.rds"))