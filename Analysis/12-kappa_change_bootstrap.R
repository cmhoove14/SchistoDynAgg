library(matrixStats)
library(geepack)
library(tidyverse)
library(doParallel)
library(doRNG)

#Datasets
  comm_sums <- readRDS(here::here("Data/Derived/adults&chldrn_shehia_sums_Unguja_Pemba_2012_2017.rds")) %>% 
    filter(Year != 2011 & !is.na(UF_alpha_mle_se))
  
  adlt_sums <- readRDS(here::here("Data/Derived/adults_shehia_sums_Unguja_Pemba_2012_2017.rds")) %>% 
    filter(Year != 2011 & !is.na(UF_alpha_mle_se))
  
  chld_sums <- readRDS(here::here("Data/Derived/chldrn_shehia_sums_Unguja_Pemba_2012_2017.rds")) %>% 
    filter(Year != 2011 & !is.na(UF_alpha_mle_se))

#Functions  
boot_sweeps <- expand.grid(Island = c("Pemba", "Unguja", "All"),
                           Treatment = c("MDA", "Snail", "Behaviour", "All"),
                           Population = c("Children", "Adults", "All"))

nboot <- 10000

get_boot_df <- function(Island, Treatment, Population){
  if(Population == "Children"){
    df <- chld_sums
  } else if(Population == "Adults"){
    df <- adlt_sums
  } else {
    df <- comm_sums
  }
  
  if(Treatment == "All"){
    df1 <- df
  } else {
    df1 <- df %>% filter(Intervention == Treatment)
  }
  
  if(Island == "All"){
    df2 <- df1
  } else {
    df2 <- df1 %>% filter(Isl == Island)
  }

  return(df2)
}

get_boot_ests <- function(df, boot_samps){
# 25th and 75th quantiles for IQR  
  egg_burden025 <- quantile(df %>% 
                              pull(UF_mean), 0.25)
  egg_burden075 <- quantile(df %>% 
                              pull(UF_mean), 0.75)

# Shehia names to sample  
  shehias <- df %>% 
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
  
  # For some strata with limited observations (e.g. many samples of adults resulted in very few positives leading to very few non-NA shehia-level observations), bootstrap sampling results in a sample that throws "Error in `contrasts<-`(`*tmp*`, value = contr.funs[1 + isOF[nn]]) : contrasts can be applied only to factors with 2 or more levels" due to all groups in the sample having two or fewer observations, so threw in this ifelse to skip over these samples
    
  if(max(count(boot_dat[!duplicated(boot_dat),"Shehia"]) %>% pull(n)) <= 2){
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

cl <- makeCluster(parallel::detectCores()-1)
registerDoParallel(cl)

boot_ests <- foreach(x            =1:nrow(boot_sweeps), 
                     .combine     = cbind,
                     .packages    = c("tidyverse", "geepack"),
                     .export      = c("get_boot_df", "get_boot_ests", "nboot",
                                      "boot_sweeps", "comm_sums", "adlt_sums", "chld_sums"),
                     .options.RNG = 7491) %dorng% {
                       
                       boot_df <- get_boot_df(boot_sweeps[x,1], boot_sweeps[x,2], boot_sweeps[x,3])
                       get_boot_ests(boot_df, nboot)
                       
                     }

parallel::stopCluster(cl)

boot_meds <- matrixStats::colMedians(boot_ests, na.rm = TRUE)
boot_loqs <- matrixStats::colQuantiles(boot_ests, probs = 0.25, na.rm = TRUE)
boot_hiqs <- matrixStats::colQuantiles(boot_ests, probs = 0.75, na.rm = TRUE)

boot_out <- cbind(boot_sweeps, boot_meds, boot_loqs, boot_hiqs)

saveRDS(boot_out, here::here("Data/Derived", paste0("dispersion_change_IQR_boot", nboot,".rds")))