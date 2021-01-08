#---------------------
# SchistoDynAgg GEE Bootstrapping functions
# Chris Hoover
#---------------------

#' @title Get data frame of observations for desired ZEST strata 
#' 
#' @description 
#' 
#' @param Island Character indicating Island of observations (either "Pemba" or "Unguja")
#' @param Intervention Character indicating intervention of observations (either "MDA" or "Snail" or "Behavior)
#' @param Population Character indicating population of observations (either "Children or "Adults")
#' 
#' @return Data frame of observations from full ZEST dataset meeting input stratifications
#' @export
#' 

get_boot_df <- function(Island, Treatment, Population) {
  if (Population == "Children") {
    df <- chld_sums
  } else if (Population == "Adults") {
    df <- adlt_sums
  } else {
    df <- comm_sums
  }
  
  if (Treatment == "All") {
    df1 <- df
  } else {
    df1 <- df %>% filter(Intervention == Treatment)
  }
  
  if (Island == "All") {
    df2 <- df1
  } else {
    df2 <- df1 %>% filter(Isl == Island)
  }
  
  return(df2)
}







#' @title Get bootstrapped GEE estimates
#' 
#' @description 
#' 
#' @param df dataframe of observations to sample and fit GEE
#' @param boot_samps number of bootstrapped samples to take
#' 
#' @return Vector of point estimates for bootstrapped dataframes
#' @export
#'

get_boot_ests <- function(df, boot_samps){
# 25th and 75th quantiles for IQR  
  egg_burden025 <- quantile(df %>% 
                              filter(!is.na(UF_alpha_mle_se)) %>% 
                              pull(UF_mean), 0.25)
  egg_burden075 <- quantile(df %>% 
                              filter(!is.na(UF_alpha_mle_se)) %>% 
                              pull(UF_mean), 0.75)

# Shehia names to sample  
  shehias <- df %>% 
    filter(!is.na(UF_alpha_mle_se)) %>% 
    pull(Shehia) %>% 
    unique()
  
boot_ests <- sapply(1:boot_samps, function(...){
  # Get sample  
    boot_samp <- sample(shehias, 
                        size = length(shehias), 
                        replace = T)
  
    boot_dat <- bind_rows(lapply(boot_samp, function(d){
      df %>% filter(Shehia == d)
    }))
  
  # Got an error "contrasts can be applied only to factors with >=2 levels" which I think is thrown by initial glm fit when trying to fit a model with factor variable that only has one level. This was probably a 1 in a million occurrence where only one shehia ended up in the bootstrap data.frame, but inputting this ifelse in order to circumvent it should it happen again    
    
  if(length(unique(boot_dat$Shehia)) < 2){
    est.out <- NA_real_
  } else {
  # Fit model
    gee_fit <- geeglm(log(UF_alpha_mle) ~ log(UF_mean), id = as.factor(Shehia),
                      weights = 1/UF_alpha_mle_se,
                      family = "gaussian", corstr = "unstructured",
                      data = boot_dat)
    
  # Get transformed estimate in terms of kappa  
    est <- exp(coef(gee_fit)[1]+coef(gee_fit)[2]*log(egg_burden075))^-1 - exp(coef(gee_fit)[1]+coef(gee_fit)[2]*log(egg_burden025))^-1
    
    est.out <- as.numeric(est)

  }  
    
  return(est.out)  
    
})

  return(boot_ests)  

}
