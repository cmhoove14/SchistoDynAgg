library(tidyverse)
library(abc)

devtools::load_all()

abc_posts <- readRDS(here::here("Data/Derived/ABC_Weighted_Posteriors.rds"))
yO <- readRDS(here::here("Data/Derived/ABC_yO_data.rds"))

post_pred_quants <- c(0.025,0.25,0.5,0.75,0.975) #Quantiles for summaries

n_cores <- parallel::detectCores()
cl <- parallel::makeCluster(n_cores)
registerDoParallel(cl)

clusterEvalQ(cl, devtools::load_all())

abc_post_pred_checks <- foreach(
  x = 1:nrow(abc_posts),
  .packages = c("tidyverse", "abc"),
  .export = c("abc_posts", "post_pred_quants", "yO"),
  .combine = rbind,
  .verbose = TRUE,
  .options.RNG = 7491
) %dorng% {
  Isl_Shehia <- str_split(as.character(abc_posts[[1]][[x]]), "_", simplify = T)
  ISL        <- Isl_Shehia[1]
  SH         <- Isl_Shehia[2]
  YR         <- as.integer(abc_posts[[2]][[x]])
  POP        <- as.character(abc_posts[[3]][[x]])
  
  n_ppl <- yO %>% 
    filter(Isl    == ISL,
           Shehia == SH,
           Year   == YR,
           pop    == POP) %>% 
    pull(n_ppl)
  
  # Posterior predictions for case 1 runs --------------
  abc_post_pred_case1 <- post_pred_data_gen(posteriors  = abc_posts[[4]][[x]], 
                                            fixed_pars  = c("h"=10, "r" = 1, "g" = 0.001),
                                            n_ppl       = n_ppl,
                                            data_gen_fx = gen_case1_data)
  
  abc_post_pred_case1_sum <- as.data.frame(matrixStats::rowQuantiles(x     = abc_post_pred_case1, 
                                                                     probs = post_pred_quants,
                                                                     na.rm = T))
  case1_0s <- sum(abc_post_pred_case1[1,] == 0)
  
  abc_post_pred_case1_df <- abc_post_pred_case1_sum %>% 
    rownames_to_column() %>% 
    rename("SumStat" = rowname,
           "q025"    = `2.5%`,
           "q25"     = `25%`,
           "q5"      = `50%`,
           "q75"     = `75%`,
           "q975"    = `97.5%`) %>% 
    mutate("Case"    = "case1",
           "Isl"     = ISL,
           "Shehia"  = SH, 
           "Year"    = YR, 
           "Pop"     = POP,
           "n_0s"    = case1_0s)
  
  rm("abc_post_pred_case1") ; gc()
  
  # Posterior predictions for case 2 runs ------------
  abc_post_pred_case2 <- post_pred_data_gen(posteriors  = abc_posts[[5]][[x]], 
                                            fixed_pars  = c("h"=10, "r" = 1, "g" = 0.001),
                                            n_ppl       = n_ppl,
                                            data_gen_fx = gen_case2_data)
  
  abc_post_pred_case2_sum <- as.data.frame(matrixStats::rowQuantiles(x     = abc_post_pred_case2, 
                                                                     probs = post_pred_quants,
                                                                     na.rm = T))
  case2_0s <- sum(abc_post_pred_case2[1,] == 0)
  
  abc_post_pred_case2_df <- abc_post_pred_case2_sum %>% 
    rownames_to_column() %>% 
    rename("SumStat" = rowname,
           "q025"    = `2.5%`,
           "q25"     = `25%`,
           "q5"      = `50%`,
           "q75"     = `75%`,
           "q975"    = `97.5%`) %>% 
    mutate("Case"    = "case2",
           "Isl"     = ISL,
           "Shehia"  = SH, 
           "Year"    = YR, 
           "Pop"     = POP,
           "n_0s"    = case2_0s)
  
  rm("abc_post_pred_case2") ; gc()
  
  # Posterior predictions for case 3 runs ---------------
  abc_post_pred_case3 <- post_pred_data_gen(posteriors  = abc_posts[[6]][[x]], 
                                            fixed_pars  = c("susc_shape" = 1, "susc_rate" = 1/500,
                                                            "h"=10, "r" = 1, "g" = 0.001),
                                            n_ppl       = n_ppl,
                                            data_gen_fx = gen_case3_data)
  
  abc_post_pred_case3_sum <- as.data.frame(matrixStats::rowQuantiles(x     = abc_post_pred_case3, 
                                                                     probs = post_pred_quants,
                                                                     na.rm = T))
  
  case3_0s <- sum(abc_post_pred_case3[1,] == 0)
  
  abc_post_pred_case3_df <- abc_post_pred_case3_sum %>% 
    rownames_to_column() %>% 
    rename("SumStat" = rowname,
           "q025"    = `2.5%`,
           "q25"     = `25%`,
           "q5"      = `50%`,
           "q75"     = `75%`,
           "q975"    = `97.5%`) %>% 
    mutate("Case"    = "case3",
           "Isl"     = ISL,
           "Shehia"  = SH, 
           "Year"    = YR, 
           "Pop"     = POP,
           "n_0s"    = case3_0s)
  
  rm("abc_post_pred_case3") ; gc()
  
  # Posterior predictions for case 4 runs -----------------
  abc_post_pred_case4 <- post_pred_data_gen(posteriors  = abc_posts[[7]][[x]], 
                                            fixed_pars  = c("h"=10, "r" = 1, "g" = 0.001),
                                            n_ppl       = n_ppl,
                                            data_gen_fx = gen_case4_data)
  
  abc_post_pred_case4_sum <- as.data.frame(matrixStats::rowQuantiles(x     = abc_post_pred_case4, 
                                                                     probs = post_pred_quants,
                                                                     na.rm = T))
  
  case4_0s <- sum(abc_post_pred_case4[1,] == 0)
  
  abc_post_pred_case4_df <- abc_post_pred_case4_sum %>% 
    rownames_to_column() %>% 
    rename("SumStat" = rowname,
           "q025"    = `2.5%`,
           "q25"     = `25%`,
           "q5"      = `50%`,
           "q75"     = `75%`,
           "q975"    = `97.5%`) %>% 
    mutate("Case"    = "case4",
           "Isl"     = ISL,
           "Shehia"  = SH, 
           "Year"    = YR, 
           "Pop"     = POP,
           "n_0s"    = case4_0s)
  
  rm(list = c("abc_post_pred_case4")) ; gc()
  
  df.out <- rbind(abc_post_pred_case1_df,
                  abc_post_pred_case2_df,
                  abc_post_pred_case3_df,
                  abc_post_pred_case4_df) %>% 
    relocate(Shehia, Year, Pop, Case) %>% 
    mutate(IQR = q75 - q25)
  
  return(df.out)
  
}  

saveRDS(abc_post_pred_checks,
        here::here("Data/Derived/abc_post_pred_checks.rds"))
