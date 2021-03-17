library(tidyverse)
library(abc)

devtools::load_all()

abc_sims <- readRDS(here::here("Data/Derived/abc_fit_rejection_4cases_hrg_fixed_1e+05iterations.rds"))
yO <- readRDS(here::here("Data/Derived/ABC_yO_data.rds"))

# Get posterior summaries ----------------------
# First remove shehia-years with failed convergence. Generally because of no infection in community
abc_successes <- sapply(1:length(abc_sims), function(s){
  ifelse(length(abc_sims[[s]]) == 4, NA, s)
})

abc_sims2 <- abc_sims[na.omit(abc_successes)]

# Get model comparison summaries ------------------------
abc_mod_comp <- do.call(rbind, lapply(abc_sims2, function(abc_run){
  shehia <- as.character(abc_run[[1]])
  year <- as.integer(abc_run[[2]])
  pop <- as.character(abc_run[[3]])
    
    abc_fit_sum <- quiet(abc:::summary.postpr(abc_run[[8]]))
    
    if("mnlogistic" %in% names(abc_fit_sum)){
      Probs  = abc_fit_sum$mnlogistic$Prob
      BayesF = abc_fit_sum$mnlogistic$BayesF
    } else {
      Probs  = abc_fit_sum$Prob
      BayesF = abc_fit_sum$BayesF
    }
    
    # if(is.null(abc_fit_sum$mnlogistic)){
    #   
    #   out.df <- data.frame("Shehia" = shehia,
    #                        "Year" = year,
    #                        "Pop" = pop,
    #                        "Case1_Post.Prob" = NA,
    #                        "Case2_Post.Prob" = NA,
    #                        "Case3_Post.Prob" = NA,
    #                        "Case4_Post.Prob" = NA,
    #                        "IIItoI_BayesF" = NA,
    #                        "IIItoII_BayesF" = NA,
    #                        "IItoI_BayesF" = NA,
    #                        "IVtoI_BayesF" = NA,
    #                        "IVtoII_BayesF" = NA,
    #                        "IVtoIII_BayesF" = NA)
    #   
    # } else {
      
      out.df <- data.frame("Shehia" = shehia,
                           "Year" = year,
                           "Pop" = pop,
                           "Case1_Post.Prob" = Probs[1],
                           "Case2_Post.Prob" = Probs[2],
                           "Case3_Post.Prob" = Probs[3],
                           "Case4_Post.Prob" = Probs[4],
                           "IIItoI_BayesF" = BayesF[3,1],
                           "IIItoII_BayesF" = BayesF[3,2],
                           "IItoI_BayesF" = BayesF[2,1],
                           "IVtoI_BayesF" =  BayesF[4,1],
                           "IVtoII_BayesF" =  BayesF[4,2],
                           "IVtoIII_BayesF" =  BayesF[4,3])
      
  #   }
  #   
  # }
  
  return(out.df)
  
}))

saveRDS(abc_mod_comp,
        file = here::here("Data/Derived/abc_mod_comps.rds"))
