library(tidyverse)
library(abc)

abc_sims <- readRDS(here::here("Data/Derived/abc_fit_ridge_log_4cases_hrg_fixed_1e+05iterations.rds"))

# Get posterior summaries ----------------------
abc_posterior_sums <- bind_rows(lapply(abc_sims, function(abc_run){
  shehia <- as.character(abc_run[[1]])
  year <- as.integer(abc_run[[2]])
  pop <- as.character(abc_run[[3]])
  
  if(length(abc_run) == 4){
    
  df.out <- data.frame("Shehia" = shehia,
                       "Year" = year,
                       "Pop" = pop,
                      # case 1 summary stats 
                       "Case1_W.med" = NA,
                       "Case1_W.loq" = NA,
                       "Case1_W.hiq" = NA,
                       "Case1_alphaW.med" = NA,
                       "Case1_alphaW.loq" = NA,
                       "Case1_alphaW.hiq" = NA,
                       "Case1_h.med" = NA,
                       "Case1_h.loq" = NA,
                       "Case1_h.hiq" = NA,
                       # case 2 summary stats 
                       "Case2_W.med" = NA,
                       "Case2_W.loq" = NA,
                       "Case2_W.hiq" = NA,
                       "Case2_alphaW.med" = NA,
                       "Case2_alphaW.loq" = NA,
                       "Case2_alphaW.hiq" = NA,
                       "Case2_h.med" = NA,
                       "Case2_h.loq" = NA,
                       "Case2_h.hiq" = NA,
                      # case 3 summary stats 
                       "Case3_W.med" = NA,
                       "Case3_W.loq" = NA,
                       "Case3_W.hiq" = NA,
                       "Case3_alphaW.med" = NA,
                       "Case3_alphaW.loq" = NA,
                       "Case3_alphaW.hiq" = NA,
                       "Case3_h.med" = NA,
                       "Case3_h.loq" = NA,
                       "Case3_h.hiq" = NA,
                      # case 4 summary stats 
                       "Case4_W.med" = NA,
                       "Case4_W.loq" = NA,
                       "Case4_W.hiq" = NA,
                       "Case4_alphaW.med" = NA,
                       "Case4_alphaW.loq" = NA,
                       "Case4_alphaW.hiq" = NA,
                       "Case4_h.med" = NA,
                       "Case4_h.loq" = NA,
                       "Case4_h.hiq" = NA)
  } else {
  
  df.out <- tryCatch({
    # Process case 1 runs  
  abc_case1_sum <- get_abc_med_iqr(abc_run[[4]]$adj.values, abc_run[[4]]$weights)

# Process case 2 runs  
  abc_case2_sum <- get_abc_med_iqr(abc_run[[5]]$adj.values, abc_run[[5]]$weights)
  
# Process case 3 runs  
  abc_case3_pars <- abc_run[[6]]$adj.values
    abc_case3_pars <- cbind(abc_case3_pars, 
                            est_case3_meanW(alphaS = 1, betaS = 1/500, meanC = abc_case3_pars[,"mean_C"]),
                            est_case3_varW(alphaS = 1, betaS = 1/500, meanC = abc_case3_pars[,"mean_C"], alphaC = abc_case3_pars[,"disp_C"]))
      colnames(abc_case3_pars)[4:5] <- c("mean_W", "var_W")
  
    abc_case3_pars <- cbind(abc_case3_pars, 
                            est_dispW(abc_case3_pars[,"mean_W"], abc_case3_pars[,"var_W"]))
  
      colnames(abc_case3_pars)[6] <- "disp_W"
      
  abc_case3_sum <- get_abc_med_iqr(abc_case3_pars, abc_run[[6]]$weights)
  
# Process case 4 runs  
  abc_case4_sum <- get_abc_med_iqr(abc_run[[7]]$adj.values, abc_run[[7]]$weights)
  
  
  data.frame("Shehia" = shehia,
             "Year" = year,
             "Pop" = pop,
             # case 1 summary stats 
             "Case1_W.med" = abc_case1_sum[2,"mean_W"],
             "Case1_W.loq" = abc_case1_sum[1, "mean_W"],
             "Case1_W.hiq" = abc_case1_sum[3, "mean_W"],
             "Case1_alphaW.med" = abc_case1_sum[2,"disp_W"],
             "Case1_alphaW.loq" = abc_case1_sum[1,"disp_W"],
             "Case1_alphaW.hiq" = abc_case1_sum[3,"disp_W"],
             "Case1_h.med" = abc_case1_sum[2,"h"],
             "Case1_h.loq" = abc_case1_sum[1,"h"],
             "Case1_h.hiq" = abc_case1_sum[3,"h"],
             # case 2 summary stats 
             "Case2_W.med" = abc_case2_sum[2,"mean_W"],
             "Case2_W.loq" = abc_case2_sum[1, "mean_W"],
             "Case2_W.hiq" = abc_case2_sum[3, "mean_W"],
             "Case2_alphaW.med" = abc_case2_sum[2,"disp_W"],
             "Case2_alphaW.loq" = abc_case2_sum[1,"disp_W"],
             "Case2_alphaW.hiq" = abc_case2_sum[3,"disp_W"],
             "Case2_h.med" = abc_case2_sum[2,"h"],
             "Case2_h.loq" = abc_case2_sum[1,"h"],
             "Case2_h.hiq" = abc_case2_sum[3,"h"],
             # case 3 summary stats 
             "Case3_W.med" = abc_case3_sum[2,"mean_W"],
             "Case3_W.loq" = abc_case3_sum[1, "mean_W"],
             "Case3_W.hiq" = abc_case3_sum[3, "mean_W"],
             "Case3_alphaW.med" = abc_case3_sum[2,"disp_W"],
             "Case3_alphaW.loq" = abc_case3_sum[1,"disp_W"],
             "Case3_alphaW.hiq" = abc_case3_sum[3,"disp_W"],
             "Case3_h.med" = abc_case3_sum[2,"h"],
             "Case3_h.loq" = abc_case3_sum[1,"h"],
             "Case3_h.hiq" = abc_case3_sum[3,"h"],
             # case 4 summary stats 
             "Case4_W.med" = abc_case4_sum[2,"mean_W"],
             "Case4_W.loq" = abc_case4_sum[1, "mean_W"],
             "Case4_W.hiq" = abc_case4_sum[3, "mean_W"],
             "Case4_alphaW.med" = abc_case4_sum[2,"disp_W"],
             "Case4_alphaW.loq" = abc_case4_sum[1,"disp_W"],
             "Case4_alphaW.hiq" = abc_case4_sum[3,"disp_W"],
             "Case4_h.med" = abc_case4_sum[2,"h"],
             "Case4_h.loq" = abc_case4_sum[1,"h"],
             "Case4_h.hiq" = abc_case4_sum[3,"h"])

  },
error = function(cond){
  data.frame("Shehia" = shehia,
             "Year" = year,
             "Pop" = pop,
             # case 1 summary stats 
             "Case1_W.med" = NA,
             "Case1_W.loq" = NA,
             "Case1_W.hiq" = NA,
             "Case1_alphaW.med" = NA,
             "Case1_alphaW.loq" = NA,
             "Case1_alphaW.hiq" = NA,
             "Case1_h.med" = NA,
             "Case1_h.loq" = NA,
             "Case1_h.hiq" = NA,
             # case 2 summary stats 
             "Case2_W.med" = NA,
             "Case2_W.loq" = NA,
             "Case2_W.hiq" = NA,
             "Case2_alphaW.med" = NA,
             "Case2_alphaW.loq" = NA,
             "Case2_alphaW.hiq" = NA,
             "Case2_h.med" = NA,
             "Case2_h.loq" = NA,
             "Case2_h.hiq" = NA,
             # case 3 summary stats 
             "Case3_W.med" = NA,
             "Case3_W.loq" = NA,
             "Case3_W.hiq" = NA,
             "Case3_alphaW.med" = NA,
             "Case3_alphaW.loq" = NA,
             "Case3_alphaW.hiq" = NA,
             "Case3_h.med" = NA,
             "Case3_h.loq" = NA,
             "Case3_h.hiq" = NA,
             # case 4 summary stats 
             "Case4_W.med" = NA,
             "Case4_W.loq" = NA,
             "Case4_W.hiq" = NA,
             "Case4_alphaW.med" = NA,
             "Case4_alphaW.loq" = NA,
             "Case4_alphaW.hiq" = NA,
             "Case4_h.med" = NA,
             "Case4_h.loq" = NA,
             "Case4_h.hiq" = NA)
  
})    

  }
  
  return(df.out)
  
}))

# Get summaries from actual generated data, weighted by posterior fits --------------------
abc_gendata_sums <- do.call(rbind, lapply(abc_sims, function(abc_run){
  shehia <- as.character(abc_run[[1]])
  year <- as.integer(abc_run[[2]])
  pop <- as.character(abc_run[[3]])
  
  if(length(abc_run) == 4){
    df.out <- data.frame("Shehia" = shehia,
                         "Year" = year,
                         "Pop" = pop,
                         # case 1 summary stats 
                         "Case1_obsW.med" = NA,
                         "Case1_obsW.loq" = NA,
                         "Case1_obsW.hiq" = NA,
                         "Case1_obsalphaW.med" = NA,
                         "Case1_obsalphaW.loq" = NA,
                         "Case1_obsalphaW.hiq" = NA,
                         "Case1_obsPhiProbW.med" = NA,
                         "Case1_obsPhiProbW.loq" = NA,
                         "Case1_obsPhiProbW.hiq" = NA,
                         "Case1_obsNonPhiW.med" = NA,
                         "Case1_obsNonPhiW.loq" = NA,
                         "Case1_obsNonPhiW.hiq" = NA,
                         "Case1_corSC.med" = NA,
                         "Case1_corSC.loq" = NA,
                         "Case1_corSC.hiq" = NA,
                         # case 2 summary stats 
                         "Case2_obsW.med" = NA,
                         "Case2_obsW.loq" = NA,
                         "Case2_obsW.hiq" = NA,
                         "Case2_obsalphaW.med" = NA,
                         "Case2_obsalphaW.loq" = NA,
                         "Case2_obsalphaW.hiq" = NA,
                         "Case2_obsPhiProbW.med" = NA,
                         "Case2_obsPhiProbW.loq" = NA,
                         "Case2_obsPhiProbW.hiq" = NA,
                         "Case2_obsNonPhiW.med" = NA,
                         "Case2_obsNonPhiW.loq" = NA,
                         "Case2_obsNonPhiW.hiq" = NA,
                         "Case2_corSC.med" = NA,
                         "Case2_corSC.loq" = NA,
                         "Case2_corSC.hiq" = NA,
                         # case 3 summary stats 
                         "Case3_obsW.med" = NA,
                         "Case3_obsW.loq" = NA,
                         "Case3_obsW.hiq" = NA,
                         "Case3_obsalphaW.med" = NA,
                         "Case3_obsalphaW.loq" = NA,
                         "Case3_obsalphaW.hiq" = NA,
                         "Case3_obsPhiProbW.med" = NA,
                         "Case3_obsPhiProbW.loq" = NA,
                         "Case3_obsPhiProbW.hiq" = NA,
                         "Case3_obsNonPhiW.med" = NA,
                         "Case3_obsNonPhiW.loq" = NA,
                         "Case3_obsNonPhiW.hiq" = NA,
                         "Case3_corSC.med" = NA,
                         "Case3_corSC.loq" = NA,
                         "Case3_corSC.hiq" = NA,
                         # case 4 summary stats 
                         "Case4_obsW.med" = NA,
                         "Case4_obsW.loq" = NA,
                         "Case4_obsW.hiq" = NA,
                         "Case4_obsalphaW.med" = NA,
                         "Case4_obsalphaW.loq" = NA,
                         "Case4_obsalphaW.hiq" = NA,
                         "Case4_obsPhiProbW.med" = NA,
                         "Case4_obsPhiProbW.loq" = NA,
                         "Case4_obsPhiProbW.hiq" = NA,
                         "Case4_obsNonPhiW.med" = NA,
                         "Case4_obsNonPhiW.loq" = NA,
                         "Case4_obsNonPhiW.hiq" = NA,
                         "Case4_corSC.med" = NA,
                         "Case4_corSC.loq" = NA,
                         "Case4_corSC.hiq" = NA)
  } else {

  df.out <- tryCatch({
# Process case 1 runs  
  abc_case1_obs <- abc_run[[9]]
    abc_case1_obs <- cbind(abc_case1_obs, 
                           est_dispW(abc_case1_obs[,"mean_W"], abc_case1_obs[,"var_W"]))
    colnames(abc_case1_obs)[ncol(abc_case1_obs)] <- "disp_W"
    
  abc_case1_obs_sum <- get_abc_med_iqr(abc_case1_obs, abc_run[[4]]$weights)
  
#Process case 2 runs    
  abc_case2_obs <- abc_run[[10]]
    abc_case2_obs <- cbind(abc_case2_obs, 
                           est_dispW(abc_case2_obs[,"mean_W"], abc_case2_obs[,"var_W"]))
    colnames(abc_case2_obs)[ncol(abc_case2_obs)] <- "disp_W"
  
  abc_case2_obs_sum <- get_abc_med_iqr(abc_case2_obs, abc_run[[5]]$weights)

#Process case 3 runs (requires more since exposure and susceptibility are separate)
  abc_case3_obs <- abc_run[[11]]
    abc_case3_obs <- cbind(abc_case3_obs, 
                           est_dispW(abc_case3_obs[,"mean_W"], abc_case3_obs[,"var_W"]))
    colnames(abc_case3_obs)[ncol(abc_case3_obs)] <- "disp_W"

    abc_case3_obs <- cbind(abc_case3_obs, 
                           abc_run[[12]][,"cor_SC"])
    colnames(abc_case3_obs)[ncol(abc_case3_obs)] <- "cor_SC"
          
  abc_case3_obs_sum <- get_abc_med_iqr(abc_case3_obs, abc_run[[6]]$weights)

#Process case 4 runs    
  abc_case4_obs <- abc_run[[9]]
    abc_case4_obs <- cbind(abc_case4_obs, 
                           est_dispW(abc_case4_obs[,"mean_W"], abc_case4_obs[,"var_W"]))
    colnames(abc_case4_obs)[ncol(abc_case4_obs)] <- "disp_W"
  
  abc_case4_obs_sum <- get_abc_med_iqr(abc_case4_obs, abc_run[[5]]$weights)
  
  
    data.frame("Shehia" = shehia,
               "Year" = year,
               "Pop" = pop,
               # case 1 summary stats 
               "Case1_obsW.med" = abc_case1_obs_sum[2,"mean_W"],
               "Case1_obsW.loq" = abc_case1_obs_sum[1, "mean_W"],
               "Case1_obsW.hiq" = abc_case1_obs_sum[3, "mean_W"],
               "Case1_obsalphaW.med" = abc_case1_obs_sum[2,"disp_W"],
               "Case1_obsalphaW.loq" = abc_case1_obs_sum[1,"disp_W"],
               "Case1_obsalphaW.hiq" = abc_case1_obs_sum[3,"disp_W"],
               "Case1_obsPhiProbW.med" = abc_case1_obs_sum[2,"prob_Phi"],
               "Case1_obsPhiProbW.loq" = abc_case1_obs_sum[1,"prob_Phi"],
               "Case1_obsPhiProbW.hiq" = abc_case1_obs_sum[3,"prob_Phi"],
               "Case1_obsNonPhiW.med" = abc_case1_obs_sum[2,"non_Phi"],
               "Case1_obsNonPhiW.loq" = abc_case1_obs_sum[1,"non_Phi"],
               "Case1_obsNonPhiW.hiq" = abc_case1_obs_sum[3,"non_Phi"],
               "Case1_corSC.med" = NA,
               "Case1_corSC.loq" = NA,
               "Case1_corSC.hiq" = NA,
               # case 2 summary stats 
               "Case2_obsW.med" = abc_case2_obs_sum[2,"mean_W"],
               "Case2_obsW.loq" = abc_case2_obs_sum[1, "mean_W"],
               "Case2_obsW.hiq" = abc_case2_obs_sum[3, "mean_W"],
               "Case2_obsalphaW.med" = abc_case2_obs_sum[2,"disp_W"],
               "Case2_obsalphaW.loq" = abc_case2_obs_sum[1,"disp_W"],
               "Case2_obsalphaW.hiq" = abc_case2_obs_sum[3,"disp_W"],
               "Case2_obsPhiProbW.med" = abc_case2_obs_sum[2,"prob_Phi"],
               "Case2_obsPhiProbW.loq" = abc_case2_obs_sum[1,"prob_Phi"],
               "Case2_obsPhiProbW.hiq" = abc_case2_obs_sum[3,"prob_Phi"],
               "Case2_obsNonPhiW.med" = abc_case2_obs_sum[2,"non_Phi"],
               "Case2_obsNonPhiW.loq" = abc_case2_obs_sum[1,"non_Phi"],
               "Case2_obsNonPhiW.hiq" = abc_case2_obs_sum[3,"non_Phi"],
               "Case2_corSC.med" = NA,
               "Case2_corSC.loq" = NA,
               "Case2_corSC.hiq" = NA,
               # case 3 summary stats 
               "Case3_obsW.med" = abc_case3_obs_sum[2,"mean_W"],
               "Case3_obsW.loq" = abc_case3_obs_sum[1, "mean_W"],
               "Case3_obsW.hiq" = abc_case3_obs_sum[3, "mean_W"],
               "Case3_obsalphaW.med" = abc_case3_obs_sum[2,"disp_W"],
               "Case3_obsalphaW.loq" = abc_case3_obs_sum[1,"disp_W"],
               "Case3_obsalphaW.hiq" = abc_case3_obs_sum[3,"disp_W"],
               "Case3_obsPhiProbW.med" = abc_case3_obs_sum[2,"prob_Phi"],
               "Case3_obsPhiProbW.loq" = abc_case3_obs_sum[1,"prob_Phi"],
               "Case3_obsPhiProbW.hiq" = abc_case3_obs_sum[3,"prob_Phi"],
               "Case3_obsNonPhiW.med" = abc_case3_obs_sum[2,"non_Phi"],
               "Case3_obsNonPhiW.loq" = abc_case3_obs_sum[1,"non_Phi"],
               "Case3_obsNonPhiW.hiq" = abc_case3_obs_sum[3,"non_Phi"],
               "Case3_corSC.med" = abc_case3_obs_sum[2,"cor_SC"],
               "Case3_corSC.loq" = abc_case3_obs_sum[1,"cor_SC"],
               "Case3_corSC.hiq" = abc_case3_obs_sum[3,"cor_SC"],
               # case 4 summary stats 
               "Case4_obsW.med" = abc_case4_obs_sum[2,"mean_W"],
               "Case4_obsW.loq" = abc_case4_obs_sum[1, "mean_W"],
               "Case4_obsW.hiq" = abc_case4_obs_sum[3, "mean_W"],
               "Case4_obsalphaW.med" = abc_case4_obs_sum[2,"disp_W"],
               "Case4_obsalphaW.loq" = abc_case4_obs_sum[1,"disp_W"],
               "Case4_obsalphaW.hiq" = abc_case4_obs_sum[3,"disp_W"],
               "Case4_obsPhiProbW.med" = abc_case4_obs_sum[2,"prob_Phi"],
               "Case4_obsPhiProbW.loq" = abc_case4_obs_sum[1,"prob_Phi"],
               "Case4_obsPhiProbW.hiq" = abc_case4_obs_sum[3,"prob_Phi"],
               "Case4_obsNonPhiW.med" = abc_case4_obs_sum[2,"non_Phi"],
               "Case4_obsNonPhiW.loq" = abc_case4_obs_sum[1,"non_Phi"],
               "Case4_obsNonPhiW.hiq" = abc_case4_obs_sum[3,"non_Phi"],
               "Case4_corSC.med" = abc_case4_obs_sum[2,"cor_SC"],
               "Case4_corSC.loq" = abc_case4_obs_sum[1,"cor_SC"],
               "Case4_corSC.hiq" = abc_case4_obs_sum[3,"cor_SC"])

  },
error = function(cond){
  data.frame("Shehia" = shehia,
             "Year" = year,
             "Pop" = pop,
             # case 1 summary stats 
             "Case1_obsW.med" = NA,
             "Case1_obsW.loq" = NA,
             "Case1_obsW.hiq" = NA,
             "Case1_obsalphaW.med" = NA,
             "Case1_obsalphaW.loq" = NA,
             "Case1_obsalphaW.hiq" = NA,
             "Case1_obsPhiProbW.med" = NA,
             "Case1_obsPhiProbW.loq" = NA,
             "Case1_obsPhiProbW.hiq" = NA,
             "Case1_obsNonPhiW.med" = NA,
             "Case1_obsNonPhiW.loq" = NA,
             "Case1_obsNonPhiW.hiq" = NA,
             "Case1_corSC.med" = NA,
             "Case1_corSC.loq" = NA,
             "Case1_corSC.hiq" = NA,
             # case 2 summary stats 
             "Case2_obsW.med" = NA,
             "Case2_obsW.loq" = NA,
             "Case2_obsW.hiq" = NA,
             "Case2_obsalphaW.med" = NA,
             "Case2_obsalphaW.loq" = NA,
             "Case2_obsalphaW.hiq" = NA,
             "Case2_obsPhiProbW.med" = NA,
             "Case2_obsPhiProbW.loq" = NA,
             "Case2_obsPhiProbW.hiq" = NA,
             "Case2_obsNonPhiW.med" = NA,
             "Case2_obsNonPhiW.loq" = NA,
             "Case2_obsNonPhiW.hiq" = NA,
             "Case2_corSC.med" = NA,
             "Case2_corSC.loq" = NA,
             "Case2_corSC.hiq" = NA,
             # case 3 summary stats 
             "Case3_obsW.med" = NA,
             "Case3_obsW.loq" = NA,
             "Case3_obsW.hiq" = NA,
             "Case3_obsalphaW.med" = NA,
             "Case3_obsalphaW.loq" = NA,
             "Case3_obsalphaW.hiq" = NA,
             "Case3_obsPhiProbW.med" = NA,
             "Case3_obsPhiProbW.loq" = NA,
             "Case3_obsPhiProbW.hiq" = NA,
             "Case3_obsNonPhiW.med" = NA,
             "Case3_obsNonPhiW.loq" = NA,
             "Case3_obsNonPhiW.hiq" = NA,
             "Case3_corSC.med" = NA,
             "Case3_corSC.loq" = NA,
             "Case3_corSC.hiq" = NA,
             # case 4 summary stats 
             "Case4_obsW.med" = NA,
             "Case4_obsW.loq" = NA,
             "Case4_obsW.hiq" = NA,
             "Case4_obsalphaW.med" = NA,
             "Case4_obsalphaW.loq" = NA,
             "Case4_obsalphaW.hiq" = NA,
             "Case4_obsPhiProbW.med" = NA,
             "Case4_obsPhiProbW.loq" = NA,
             "Case4_obsPhiProbW.hiq" = NA,
             "Case4_obsNonPhiW.med" = NA,
             "Case4_obsNonPhiW.loq" = NA,
             "Case4_obsNonPhiW.hiq" = NA,
             "Case4_corSC.med" = NA,
             "Case4_corSC.loq" = NA,
             "Case4_corSC.hiq" = NA)


})    
  
  }
  
  return(df.out)
}))

# Get model comparison summaries ------------------------
abc_fit_stats <- do.call(rbind, lapply(abc_sims, function(abc_run){
  shehia <- as.character(abc_run[[1]])
  year <- as.integer(abc_run[[2]])
  pop <- as.character(abc_run[[3]])

  if(length(abc_run) == 4){
    out.df <- data.frame("Shehia" = shehia,
                         "Year" = year,
                         "Pop" = pop,
                         "Case1_Post.Prob" = NA,
                         "Case2_Post.Prob" = NA,
                         "Case3_Post.Prob" = NA,
                         "Case4_Post.Prob" = NA,
                         "IIItoI_BayesF" = NA,
                         "IIItoII_BayesF" = NA,
                         "IItoI_BayesF" = NA,
                         "IVtoI_BayesF" = NA,
                         "IVtoII_BayesF" = NA,
                         "IVtoIII_BayesF" = NA)
  } else {
    
    abc_fit_sum <- quiet(abc:::summary.postpr(abc_run[[8]]))
    
    if(is.null(abc_fit_sum$mnlogistic)){
      
      out.df <- data.frame("Shehia" = shehia,
                           "Year" = year,
                           "Pop" = pop,
                           "Case1_Post.Prob" = NA,
                           "Case2_Post.Prob" = NA,
                           "Case3_Post.Prob" = NA,
                           "Case4_Post.Prob" = NA,
                           "IIItoI_BayesF" = NA,
                           "IIItoII_BayesF" = NA,
                           "IItoI_BayesF" = NA,
                           "IVtoI_BayesF" = NA,
                           "IVtoII_BayesF" = NA,
                           "IVtoIII_BayesF" = NA)
      
    } else {
      
    out.df <- data.frame("Shehia" = shehia,
                         "Year" = year,
                         "Pop" = pop,
                         "Case1_Post.Prob" = abc_fit_sum$mnlogistic$Prob[1],
                         "Case2_Post.Prob" = abc_fit_sum$mnlogistic$Prob[2],
                         "Case3_Post.Prob" = abc_fit_sum$mnlogistic$Prob[3],
                         "Case4_Post.Prob" = abc_fit_sum$mnlogistic$Prob[4],
                         "IIItoI_BayesF" = abc_fit_sum$mnlogistic$BayesF[3,1],
                         "IIItoII_BayesF" = abc_fit_sum$mnlogistic$BayesF[3,2],
                         "IItoI_BayesF" = abc_fit_sum$mnlogistic$BayesF[2,1],
                         "IVtoI_BayesF" =  abc_fit_sum$mnlogistic$BayesF[4,1],
                         "IVtoII_BayesF" =  abc_fit_sum$mnlogistic$BayesF[4,2],
                         "IVtoIII_BayesF" =  abc_fit_sum$mnlogistic$BayesF[4,3])
      
    }
    
  }
  
  return(out.df)

}))

# Get generated datasets to compare to observed for model assessment
abc_sumstats <- do.call(rbind, lapply(abc_sims, function(abc_run){
  shehia <- as.character(abc_run[[1]])
  year <- as.integer(abc_run[[2]])
  pop <- abc_run[[3]]
  
  if(length(abc_run) == 4){
  out.df <- data.frame("Shehia" = shehia,
                       "Year" = year,
                       "Pop" = pop,
                      # case 1 summary stats 
                       "Case1_genE.med" = NA,
                       "Case1_genE.loq" = NA,
                       "Case1_genE.hiq" = NA,
                       "Case1_genEse.med" = NA,
                       "Case1_genEse.loq" = NA,
                       "Case1_genEse.hiq" = NA,
                       "Case1_genEpos2n.med" = NA,
                       "Case1_genEpos2n.loq" = NA,
                       "Case1_genEpos2n.hiq" = NA,
                      
                      # case 2 summary stats 
                       "Case2_genE.med" = NA,
                       "Case2_genE.loq" = NA,
                       "Case2_genE.hiq" = NA,
                       "Case2_genEse.med" = NA,
                       "Case2_genEse.loq" = NA,
                       "Case2_genEse.hiq" = NA,
                       "Case2_genEpos2n.med" = NA,
                       "Case2_genEpos2n.loq" = NA,
                       "Case2_genEpos2n.hiq" = NA,

                      # case 3 summary stats 
                       "Case3_genE.med" = NA,
                       "Case3_genE.loq" = NA,
                       "Case3_genE.hiq" = NA,
                       "Case3_genEse.med" = NA,
                       "Case3_genEse.loq" = NA,
                       "Case3_genEse.hiq" = NA,
                       "Case3_genEpos2n.med" = NA,
                       "Case3_genEpos2n.loq" = NA,
                       "Case3_genEpos2n.hiq" = NA,
                      
                      # case 4 summary stats 
                       "Case4_genE.med" = NA,
                       "Case4_genE.loq" = NA,
                       "Case4_genE.hiq" = NA,
                       "Case4_genEse.med" = NA,
                       "Case4_genEse.loq" = NA,
                       "Case4_genEse.hiq" = NA,
                       "Case4_genEpos2n.med" = NA,
                       "Case4_genEpos2n.loq" = NA,
                       "Case4_genEpos2n.hiq" = NA)

  } else {
    
  abc_case1_ss <- abc_run[[4]]$ss 

  abc_case2_ss <- abc_run[[5]]$ss 

  abc_case3_ss <- abc_run[[6]]$ss 

  abc_case4_ss <- abc_run[[7]]$ss 
  
  out.df <- data.frame("Shehia" = shehia,
                       "Year" = year,
                       "Pop" = pop,
                      # case 1 summary stats 
                       "Case1_genE.med" = quantile(abc_case1_ss[,1], 0.5),
                       "Case1_genE.loq" = quantile(abc_case1_ss[,1], 0.25),
                       "Case1_genE.hiq" = quantile(abc_case1_ss[,1], 0.75),
                       "Case1_genEse.med" = quantile(abc_case1_ss[,2], 0.5),
                       "Case1_genEse.loq" = quantile(abc_case1_ss[,2], 0.25),
                       "Case1_genEse.hiq" = quantile(abc_case1_ss[,2], 0.75),
                       "Case1_genEpos2n.med" = quantile(abc_case1_ss[,3], 0.5),
                       "Case1_genEpos2n.loq" = quantile(abc_case1_ss[,3], 0.25),
                       "Case1_genEpos2n.hiq" = quantile(abc_case1_ss[,3], 0.75),
                      
                      # case 2 summary stats 
                       "Case2_genE.med" = quantile(abc_case2_ss[,1], 0.5),
                       "Case2_genE.loq" = quantile(abc_case2_ss[,1], 0.25),
                       "Case2_genE.hiq" = quantile(abc_case2_ss[,1], 0.75),
                       "Case2_genEse.med" = quantile(abc_case2_ss[,2], 0.5),
                       "Case2_genEse.loq" = quantile(abc_case2_ss[,2], 0.25),
                       "Case2_genEse.hiq" = quantile(abc_case2_ss[,2], 0.75),
                       "Case2_genEpos2n.med" = quantile(abc_case2_ss[,3], 0.5),
                       "Case2_genEpos2n.loq" = quantile(abc_case2_ss[,3], 0.25),
                       "Case2_genEpos2n.hiq" = quantile(abc_case2_ss[,3], 0.75),

                      # case 3 summary stats 
                       "Case3_genE.med" = quantile(abc_case3_ss[,1], 0.5),
                       "Case3_genE.loq" = quantile(abc_case3_ss[,1], 0.25),
                       "Case3_genE.hiq" = quantile(abc_case3_ss[,1], 0.75),
                       "Case3_genEse.med" = quantile(abc_case3_ss[,2], 0.5),
                       "Case3_genEse.loq" = quantile(abc_case3_ss[,2], 0.25),
                       "Case3_genEse.hiq" = quantile(abc_case3_ss[,2], 0.75),
                       "Case3_genEpos2n.med" = quantile(abc_case3_ss[,3], 0.5),
                       "Case3_genEpos2n.loq" = quantile(abc_case3_ss[,3], 0.25),
                       "Case3_genEpos2n.hiq" = quantile(abc_case3_ss[,3], 0.75),
    
                      # case 4 summary stats 
                       "Case4_genE.med" = quantile(abc_case4_ss[,1], 0.5),
                       "Case4_genE.loq" = quantile(abc_case4_ss[,1], 0.25),
                       "Case4_genE.hiq" = quantile(abc_case4_ss[,1], 0.75),
                       "Case4_genEse.med" = quantile(abc_case4_ss[,2], 0.5),
                       "Case4_genEse.loq" = quantile(abc_case4_ss[,2], 0.25),
                       "Case4_genEse.hiq" = quantile(abc_case4_ss[,2], 0.75),
                       "Case4_genEpos2n.med" = quantile(abc_case4_ss[,3], 0.5),
                       "Case4_genEpos2n.loq" = quantile(abc_case4_ss[,3], 0.25),
                       "Case4_genEpos2n.hiq" = quantile(abc_case4_ss[,3], 0.75))
  }
  
  return(out.df)
}))

# Final datasets for analysis -----------
abc_fin_df <- cbind(yO, 
                    abc_posterior_sums[,-c(1:3)], 
                    abc_gendata_sums[,-c(1:3)], 
                    abc_fit_stats [,-c(1:3)], 
                    abc_sumstats[,-c(1:3)])
  
abc_fin_df2 <- abc_fin_df %>% 
  # Determine for Case 3 if closer to Case 1 or Case2 estimates
  mutate(Case3_diff1 = (Case1_W.med-Case3_W.med)^2,
         Case3_diff2 = (Case2_W.med-Case3_W.med)^2,
         Case3_closer = if_else(Case3_diff1 > Case3_diff2, 2, 1))

abc_fin_df_long <- abc_fin_df2 %>% 
  pivot_longer(cols = Case1_W.med:Case3_closer,
               names_to = c("Case", "Measure"),
               names_sep = "_")

abc_fin_df_case_long <- abc_fin_df_long %>% 
  pivot_wider(names_from = "Measure", 
              values_from = "value",
              values_fn = {first})

save(abc_fin_df_case_long, abc_fin_df_long, abc_fin_df2, abc_fin_df,
     here::here("Data/Derived/abc_processed_results.Rdata"))

rm(abc_sims) ; quiet(gc())
