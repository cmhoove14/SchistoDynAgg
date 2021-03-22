library(tidyverse)
library(abc)

devtools::load_all()

# Get options passed from BASH -----------------
# Number of weighted draws from initial posterior to construct final posterior  
opts <- commandArgs(TRUE)

if(length(opts) > 0){
  post_samp_size <- as.numeric(opts[1])
} else {
  post_samp_size <- 1000
  
}

# Load data and additional setup ----------------

abc_sims <- readRDS(here::here("Data/Derived/abc_fit_rejection_4cases_hrg_fixed_1e+05iterations.rds"))
yO <- readRDS(here::here("Data/Derived/ABC_yO_data.rds"))

# First remove shehia-years with failed convergence. Generally because of no infection in community
  abc_successes <- sapply(1:length(abc_sims), function(s){
    ifelse(length(abc_sims[[s]]) == 4, NA, s)
  })
  
  abc_sims2 <- abc_sims[na.omit(abc_successes)]
  
# Get posterior summaries ----------------------
abc_posteriors <- bind_rows(lapply(abc_sims2, function(abc_run){
  shehia <- as.character(abc_run[[1]])
  year   <- as.integer(abc_run[[2]])
  pop    <- as.character(abc_run[[3]])

# Process case 1 runs  
  case1_post_dists <- abc_run[[4]]$dist[abc_run[[4]]$region]                        # Distance metrics of inputs reaching tolerance
  case1_weights    <- max(case1_post_dists)/(case1_post_dists+1/post_samp_size)     # Weights normalized based on largest distance, truncated to avoid Inf 
  case1_init_post  <- abc_run[[4]]$unadj.values                                     # Initial posterior as values reaching tolerance with rejection method
  case1_post_samps <- sample(                                                       # Weighted sample
    x       = 1:nrow(case1_init_post), 
    size    = post_samp_size, 
    prob    = case1_weights, 
    replace = TRUE
  )
  case1_wgtd_post  <- case1_init_post[case1_post_samps,]       # Final weighted posterior
  
# Process case 2 runs  
  case2_post_dists <- abc_run[[5]]$dist[abc_run[[5]]$region]     # Distance metrics of inputs reaching tolerance
  case2_weights    <- max(case2_post_dists)/(case2_post_dists+1/post_samp_size)     # Weights normalized based on largest distance
  case2_init_post  <- abc_run[[5]]$unadj.values                  # Initial posterior as values reaching tolerance with rejection method
  case2_post_samps <- sample(                                    # Weighted sample
    x       = 1:nrow(case2_init_post), 
    size    = post_samp_size, 
    prob    = case2_weights, 
    replace = TRUE
  )
  case2_wgtd_post  <- case2_init_post[case2_post_samps,]       # Final weighted posterior
  
# Process case 3 runs  
  case3_post_dists <- abc_run[[6]]$dist[abc_run[[6]]$region]     # Distance metrics of inputs reaching tolerance
  case3_weights    <- max(case3_post_dists)/(case3_post_dists+1/post_samp_size)     # Weights normalized based on largest distance
  case3_init_post  <- abc_run[[6]]$unadj.values                  # Initial posterior as values reaching tolerance with rejection method
  case3_post_samps <- sample(                                    # Weighted sample
    x       = 1:nrow(case3_init_post), 
    size    = post_samp_size, 
    prob    = case3_weights, 
    replace = TRUE
  )
  case3_wgtd_post  <- case3_init_post[case3_post_samps,]       # Final weighted posterior
  
# Process case 4 runs  
  case4_post_dists <- abc_run[[7]]$dist[abc_run[[7]]$region]     # Distance metrics of inputs reaching tolerance
  case4_weights    <- max(case4_post_dists)/(case4_post_dists+1/post_samp_size)     # Weights normalized based on largest distance
  case4_init_post  <- abc_run[[7]]$unadj.values                  # Initial posterior as values reaching tolerance with rejection method
  case4_post_samps <- sample(                                    # Weighted sample
    x       = 1:nrow(case4_init_post), 
    size    = post_samp_size, 
    prob    = case4_weights, 
    replace = TRUE
  )
  case4_wgtd_post  <- case4_init_post[case4_post_samps,]       # Final weighted posterior
  
  return(list("Shehia"     = list(shehia),
              "Year"       = list(year),
              "Pop"        = list(pop),
              "Case1_Pstr" = list(case1_wgtd_post),
              "Case2_Post" = list(case2_wgtd_post),
              "Case3_Post" = list(case3_wgtd_post),
              "Case4_Post" = list(case4_wgtd_post)))
  
}))

saveRDS(abc_posteriors,
        file = here::here("Data/Derived/ABC_Weighted_Posteriors.rds"))
