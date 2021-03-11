# -----------------------------
# SchistoDynAgg ABC worm burden estimation from egg burden
# Chris Hoover
# -----------------------------


# Setup ---------------
  devtools::load_all()
  n_iter <- 1e5


# Load data and process for input to ABC ---------------
# Filter out pilot observations of 50 children in year 2011 and only focus on study years 2012-2017
comm_sums <- readRDS(here::here("Data", "Derived", "adults&chldrn_shehia_sums_Unguja_Pemba_2012_2017.rds")) %>% 
  filter(Year != 2011)

adlt_sums <- readRDS(here::here("Data", "Derived", "adults_shehia_sums_Unguja_Pemba_2012_2017.rds")) %>% 
  filter(Year != 2011)

chld_sums <- readRDS(here::here("Data", "Derived", "chldrn_shehia_sums_Unguja_Pemba_2012_2017.rds")) %>% 
  filter(Year != 2011)

  
  yOm <- comm_sums %>%
    filter(UF_max > 1) %>%
    mutate(
      UF_se = sqrt(UF_var) / sqrt(n_ppl),
      UFpos2n = UF_pos ^ 2 / n_ppl,
      pop = "Comm"
    ) %>%
    ungroup() %>%
    dplyr::select(Isl, Shehia, Year, UF_mean, UF_se, UF_pos, n_ppl, UFpos2n, pop)
  
  
  yOc <- chld_sums %>%
      filter(UF_max > 1) %>%
      mutate(
        UF_se = sqrt(UF_var) / sqrt(n_ppl),
        UFpos2n = UF_pos ^ 2 / n_ppl,
        pop = "Child"
      ) %>%
      ungroup() %>%
      dplyr::select(Isl, Shehia, Year, UF_mean, UF_se, UF_pos, n_ppl, UFpos2n, pop)
  
  
  
  yOa <- adlt_sums %>%
    filter(UF_max > 1) %>%
    mutate(
      UF_se = sqrt(UF_var) / sqrt(n_ppl),
      UFpos2n = UF_pos ^ 2 / n_ppl,
      pop = "Adult"
    ) %>%
    ungroup() %>%
    dplyr::select(Isl, Shehia, Year, UF_mean, UF_se, UF_pos, n_ppl, UFpos2n, pop)
  
  
  
  yO <- bind_rows(yOm, yOc, yOa)
  
  saveRDS(yO, here::here("Data/Derived/ABC_yO_data.rds"))
# Declare priors and run in parallel. takes ~ 3 hours on 8 core Lenovo ThinkPad with AMD Ryzen processor or ~1.5 hours on 24 compute node on Biostat cluster with 24 cores ----------------------
  abc_priors <- list(
    "mean_W_lo"     = 1e-6,
    "mean_W_hi"     = 1e3,
    "disp_W_lo"     = 1e-5,
    "disp_W_hi"     = 1e5,
    "susc_shape_lo" = 1,
    "susc_shape_hi" = 1,
    "susc_rate_lo"  = 1 / 500,
    "susc_rate_hi"  = 1 / 500,
    "mean_C_lo"     = 1,
    "mean_C_hi"     = 1e6,
    "disp_C_lo"     = 1e-5,
    "disp_C_hi"     = 1e5,
    "h_lo"          = 10,
    "h_hi"          = 10,
    "r_lo"          = 1,
    "r_hi"          = 1,
    "g_lo"          = 0.001,
    "g_hi"          = 0.001
  )
  
#Setup for running jobs across parallel nodes in cluster
  
  start.time <- Sys.time()
  n_cores <- parallel::detectCores()
  cl <- parallel::makeCluster(n_cores)
  registerDoParallel(cl)
  
  clusterEvalQ(cl, devtools::load_all())
  
  abc_sims <- foreach(
    x = 1:nrow(yO),
    .packages = c("tidyverse", "abc"),
    .verbose = TRUE,
    .options.RNG = 7491
  ) %dorng% {
    # See R/abc_worm_estimation_functions for abc_fit and other functions
    abc_fit(obs_data = yO[x, ],
            priors = abc_priors,
            iterations = n_iter)
  }
  
  parallel::stopCluster(cl)
  
  end.time <- Sys.time()
  
  end.time - start.time
  
  saveRDS(
    abc_sims,
    file = here::here("Data/Derived", 
                      paste0(
                        "abc_fit_ridge_log_4cases_hrg_fixed_",
                        n_iter,
                        "iterations.rds"
                      ))
  )    
