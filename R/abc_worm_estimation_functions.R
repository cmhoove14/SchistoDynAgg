# ----------------------------
# SchistoDynAgg ABC Data generation functions
# Chris Hoover
# ----------------------------


# Couple util functions --------------

#' @title Quitely run
#' 
#' @description Run without outputting messages/warnings
#' 
#' @param x command call
#' 
#' @return silent result of command call
#' @export
#'

quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 





#' @title Case 3 mean from parameters
#' 
#' @description Analytic estimate of mean as combination of neg. bin. cercarial exposures and gamma susceptibility
#' 
#' @param alphaS gamma shape parameter for susceptibility
#' @param betaS gamma rate parameter for susceptibility
#' @param meanC mean cercarial exposure
#' 
#' @return estimate of mean
#' @export
#'

est_case3_meanW <- function(alphaS, betaS, meanC){
  alphaS*betaS*meanC
}




#' @title Case 3 variance from parameters
#' 
#' @description Analytic estimate of variance as combination of neg. bin. cercarial exposures and gamma susceptibility
#' 
#' @param alphaS gamma shape parameter for susceptibility
#' @param betaS gamma rate parameter for susceptibility
#' @param meanC mean cercarial exposure
#' @param alphaC inverse dispersion parameter of cercarial exposure
#' 
#' @return estimate of variance
#' @export
#'

est_case3_varW <- function(alphaS, betaS, meanC, alphaC){
  meanS <- alphaS*betaS
  varS <- alphaS*betaS^2
  varC <- meanC^2*alphaC+meanC
  (varC+meanC^2)*(varS+meanS^2)-(meanS*meanC)^2
}






#' @title Moments estimator of dispersion parameter
#' 
#' @description 
#' 
#' @param meanW mean worm burden
#' @param varW variance of worm burden
#' 
#' @return estimate of dispersion parameter
#' @export
#'

est_dispW <- function(meanW, varW){
  (varW-meanW)/meanW^2
}











#' @title Median and IQR of abc posterior
#' 
#' @description Utility function to return median and 25th and 75th quantiles of posterior distribution from ABC sim
#' 
#' @param abc_pars vector of abc pars from fitting
#' @param abc_weights vector of weights returned from abc posterior correction
#' 
#' @return vector with median and iqr of posterior
#' @export
#'
get_abc_med_iqr <- function(abc_pars, abc_weights){
  myprobs <- c(0.25, 0.5, 0.75)
  wt <- abc_weights
  
  quants <- apply(abc_pars, 2, function(x) quantreg::rq(x ~ 1, tau = myprobs, 
                                                        weights = wt)$coef)

  return(quants)
}











#' @title Generate synthetic case 1 dataset 
#' 
#' @description Generate worm pairs and resulting egg burdens for individuals assuming males and females distributed together from same negative binomial distribution, return summary statistics
#' 
#' @param pars numeric vector with mean worm burden, inverse dispersion parameter, mean eggs per mated pair, and dispersion of daily egg release
#' @param n number of individual human hosts
#' 
#' @return vector of summary statistics to compare to observed in ABC setup
#' @export
#'
 
gen_case1_data <- function(pars, n){
  mean_W <- pars[1]    # Mean worm burden
  disp_W <- 1/pars[2]  # Aggregation of mean worm burden
  
# Simulate individual worm burdens  
  Ni <- rnbinom(n, mu = mean_W, size = disp_W) 
  
# Simulate pairing assuming male and female distributed together  
  Nf <- rpois(n, Ni*0.5)
  Nm <- Ni-Nf
  
# Get number ated pairs as minimum of individual male and female burdens  
  Xi <- matrixStats::rowMins(cbind(Nm,Nf))
  
  h <- pars[3]      # Mean daily egg release per mated female worm per 10mL urine
  r <- pars[4]      # Aggregation parameter of h assuming neg. binomial dist'n
  g <- pars[5]      # Dens. dependent fecundity parameter
  
# Simulate individual egg burdens by randomly generating egg release for each mated pair
  ddf <- exp(-g*Ni)
  
  Ei <- rnbinom(n, mu = h*Xi*ddf, size = r)
  
# Get adjusted prevalence as number people egg positive squared over number people  
  E_pos2n <- sum(Ei > 0)^2/n

# Return summary statistics
  return(c(E = mean(Ei),
           E_se = sqrt(var(Ei))/sqrt(n),
           E_pos2n = E_pos2n,
           mean_W = mean(Ni),
           var_W = var(Ni),
           mean_Phi = mean(Xi),
           var_Phi = var(Xi),
           prob_Phi = (sum(Xi)*2)/sum(Ni),
           non_Phi = sum(Ni)-sum(Xi)*2))
}









#' @title Generate synthetic case 2 dataset 
#' 
#' @description Generate worm pairs and resulting egg burdens for individuals assuming males and females distributed from separate negative binomial distributions, return summary statistics. All prior distributions are drawn from uniform distributions 
#' 
#' @param pars numeric vector with mean worm burden, inverse dispersion parameter, mean eggs per mated pair, and dispersion of daily egg release
#' @param n number of individual human hosts
#' 
#' @return vector of summary statistics to compare to observed in ABC setup
#' @export
#'

gen_case2_data <- function(pars, n){
  mean_W <- pars[1]    # Mean worm burden
  disp_W <- 1/pars[2]     # Aggregation of mean worm burden
  
# Simulate pairing assuming male and female distributed separately  
  Nf <- rnbinom(n, mu = mean_W/2, size = disp_W)
  Nm <- rnbinom(n, mu = mean_W/2, size = disp_W)
  Ni <- Nf+Nm
  
  Xi <- matrixStats::rowMins(cbind(Nm,Nf))
  
  h <- pars[3]      # Mean daily egg release per mated female worm per 10mL urine
  r <- pars[4]      # Aggregation parameter of h assuming neg. binomial dist'n
  g <- pars[5]      # Dens. dependent fecundity parameter
  
# Simulate individual egg burdens by randomly generating egg release for each mated pair
  ddf <- exp(-g*Ni)
  
  Ei <- rnbinom(n, mu = h*Xi*ddf, size = r)
  
# Get adjusted prevalence as number people egg positive squared over number people  
  E_pos2n <- sum(Ei > 0)^2/n

# Return summary statistics
  return(c(E = mean(Ei),
           E_se = sqrt(var(Ei))/sqrt(n),
           E_pos2n = E_pos2n,
           mean_W = mean(Ni),
           var_W = var(Ni),
           mean_Phi = mean(Xi),
           var_Phi = var(Xi),
           prob_Phi = (sum(Xi)*2)/sum(Ni),
           non_Phi = sum(Ni)-sum(Xi)*2))
}










#' @title Generate priors for case 1/2 dataset generation 
#' 
#' @description Generate priors for ABC data generation for cases 1 and 2
#' 
#' @param iterations numeric of number of simulations to run
#' @param mean_w_lo lower range of mean worm burden
#' @param mean_w_hi upper range of mean worm burden
#' @param disp_w_lo lower range of worm burden dispersion
#' @param disp_w_hi upper range of worm burden dispersion
#' @param h_lo lower range of daily mean egg release per mated pair
#' @param h_hi upper range of daily mean egg release per mated pair
#' @param r_lo lower range of daily egg release dispersion
#' @param r_hi upper range of daily egg release dispersion
#' @param g_lo lower range of density dependent fecundity parameter
#' @param g_hi upper range of density dependent fecundity parameter
#' 
#' @return matrix of priors
#' @export
#'

gen_case12_pars <- function(iterations,
                            mean_W_lo, mean_W_hi,
                            disp_W_lo, disp_W_hi, 
                            h_lo, h_hi, 
                            r_lo, r_hi,
                            g_lo, g_hi){
  
  mean_W <- exp(runif(iterations, log(mean_W_lo), log(mean_W_hi)))
  disp_W <- exp(runif(iterations, log(disp_W_lo), log(disp_W_hi)))
  h <- exp(runif(iterations, log(h_lo), log(h_hi)))
  r <- exp(runif(iterations, log(r_lo), log(r_hi)))
  g <- exp(runif(iterations, log(g_lo), log(g_hi)))
  
  return(cbind(mean_W, disp_W, h, r, g))
  
}







#' @title Generate synthetic case 3 dataset 
#' 
#' @description Generate worm pairs and resulting egg burdens for individuals assuming gamma distributed susceptibility and negative binomially distributed cercarial exposures 
#' 
#' @param pars numeric vector with gamma shape and scale parameters for susceptibility distribution, mean cercarial exposure, and dispersion of cercarial exposure, mean eggs per mated pair, and dispersion of daily egg release
#' @param n number of individual human hosts
#' 
#' @return vector of summary statistics to compare to observed in ABC setup
#' @export
#'

gen_case3_data <- function(pars, n){
  susc_shape <- pars[1]    # Susceptibility shape1
  susc_rate <- pars[2]     # Susceptibility shape2
  
  mean_C <- pars[3]       # mean cercarial exposure
  disp_C <- 1/pars[4]      # distribution of cercarial exposures
# Individual susceptibilities  
  Si <- rgamma(n, shape = susc_shape, scale = susc_rate)
  
# Cercarial exposures with no correlation to susceptibility
  Ci <- rnbinom(n, mu = mean_C, size = disp_C)

# Cercarial exposures resulting in adult worm 
  Ni <- rpois(n,Ci*Si)
  
# Cercarial exposures which are female assuming equal sex distribution
  Fi <- rpois(n, Ci*0.5)
  
# Female worms from hypergeometric dist'n
  Nf <- rhyper(n, Ci-Fi, Fi, Ni)
  
# Male worms as worms not female  
  Nm <- Ni - Nf
  
# Mated pairs  
  Xi <- matrixStats::rowMins(cbind(Nm,Nf))
  
  h <- pars[5]      # Mean daily egg release per mated female worm per 10mL urine
  r <- pars[6]      # Aggregation parameter of h assuming neg. binomial dist'n
  g <- pars[7]      # Dens. dependent fecundity parameter
  
# Simulate individual egg burdens by randomly generating egg release for each mated pair
  ddf <- exp(-g*Ni)
  
  Ei <- rnbinom(n, mu = h*Xi*ddf, size = r)
  
# Get adjusted prevalence as number people egg positive squared over number people  
  E_pos2n <- sum(Ei > 0)^2/n

# Return summary statistics
  return(c(E = mean(Ei),
           E_se = sqrt(var(Ei))/sqrt(n),
           E_pos2n = E_pos2n,
           mean_W = mean(Ni),
           var_W = var(Ni),
           mean_Phi = mean(Xi),
           var_Phi = var(Xi),
           prob_Phi = (sum(Xi)*2)/sum(Ni),
           non_Phi = sum(Ni)-sum(Xi)*2,
           mean_C = mean(Ci),
           mean_S = mean(Si),
           var_C = var(Ci),
           var_S = var(Si),
           cov_SC = cov(Si, Ci),
           cov_S2C2 = cov(Si^2, Ci^2),
           cor_SC = cor(Si, Ci)))
}







#' @title Generate priors for case 3 dataset generation 
#' 
#' @description Generate priors for ABC data generation for case 3
#' 
#' @param iterations numeric of number of simulations to run
#' @param susc_shape_lo lower range of susceptibility shape parameter
#' @param susc_shape_hi upper range of susceptibility shape parameter
#' @param susc_rate_lo lower range of susceptibility rate parameter
#' @param susc_rate_hi upper range of susceptibility rate parameter
#' @param mean_C_lo lower range of mean cercarial exposure
#' @param mean_C_hi upper range of mean cercarial exposure
#' @param disp_C_lo lower range of cercarial exposure dispersion
#' @param disp_C_hi upper range of cercarial exposure dispersion
#' @param h_lo lower range of daily mean egg release per mated pair
#' @param h_hi upper range of daily mean egg release per mated pair
#' @param r_lo lower range of daily egg release dispersion
#' @param r_hi upper range of daily egg release dispersion
#' @param g_lo lower range of density dependent fecundity parameter
#' @param g_hi upper range of density dependent fecundity parameter
#' 
#' @return matrix of priors
#' @export
#'

gen_case3_pars <- function(iterations,
                           susc_shape_lo, susc_shape_hi,
                           susc_rate_lo, susc_rate_hi,
                           mean_C_lo, mean_C_hi,
                           disp_C_lo, disp_C_hi, 
                           h_lo, h_hi, 
                           r_lo, r_hi,
                           g_lo, g_hi){
  
  alpha_S <- exp(runif(iterations, log(susc_shape_lo), log(susc_shape_hi)))
  beta_S <- exp(runif(iterations, log(susc_rate_lo), log(susc_rate_hi)))
  mean_C <- exp(runif(iterations, log(mean_C_lo), log(mean_C_hi)))
  disp_C <- exp(runif(iterations, log(disp_C_lo), log(disp_C_hi)))
  h <- exp(runif(iterations, log(h_lo), log(h_hi)))
  r <- exp(runif(iterations, log(r_lo), log(r_hi)))
  g <- exp(runif(iterations, log(g_lo), log(g_hi)))
  
  return(cbind(alpha_S, beta_S, mean_C, disp_C, h, r, g))
  
}









#' @title Generate synthetic case 4 dataset 
#' 
#' @description Generate worm pairs and resulting egg burdens for individuals assuming hybrid of case 1 and 2 distributions with partitioning parameter
#' 
#' @param pars numeric vector with mean worm burden, inverse dispersion parameter, partitioning parameter, mean eggs per mated pair, and dispersion of daily egg release
#' @param n number of individual human hosts
#' 
#' @return vector of summary statistics to compare to observed in ABC setup
#' @export
#'
gen_case4_data <- function(pars, n){
  mean_W <- pars[1]    # Mean worm burden
  disp_W <- 1/pars[2]  # Aggregation of mean worm burden
  part_W <- pars[3]    # Proportion of worms following case1 (together) dynamics
  
# Simulate individual worm burdens following case 1 dynamics  
  Ni1 <- rnbinom(n, mu = mean_W * part_W, size = disp_W) 
  
  Nf1 <- rpois(n, Ni1*0.5)
  Nm1 <- Ni1-Nf1
  
# Simulate individual worm burdens following case 2 dynamics    
  Nf2 <- rnbinom(n, mu = mean_W/2 * (1-part_W), size = disp_W)
  Nm2 <- rnbinom(n, mu = mean_W/2 * (1-part_W), size = disp_W)
  Ni2 <- Nf2+Nm2
  
# Get pairing from both case 1 and two worms
  Nf <- Nf1 + Nf2
  Nm <- Nm1 + Nm2
  Ni <- Ni1 + Ni2
  
  Xi <- matrixStats::rowMins(cbind(Nm,Nf))
  
  h <- pars[4]      # Mean daily egg release per mated female worm per 10mL urine
  r <- pars[5]      # Aggregation parameter of h assuming neg. binomial dist'n
  g <- pars[6]      # Dens. dependent fecundity parameter
  
# Simulate individual egg burdens by randomly generating egg release for each mated pair
  ddf <- exp(-g*Ni)
  
  Ei <- rnbinom(n, mu = h*Xi*ddf, size = r)
  
# Get adjusted prevalence as number people egg positive squared over number people  
  E_pos2n <- sum(Ei > 0)^2/n

# Return summary statistics
  return(c(E = mean(Ei),
           E_se = sqrt(var(Ei))/sqrt(n),
           E_pos2n = E_pos2n,
           mean_W = mean(Ni),
           var_W = var(Ni),
           mean_Phi = mean(Xi),
           var_Phi = var(Xi),
           prob_Phi = (sum(Xi)*2)/sum(Ni),
           non_Phi = sum(Ni)-sum(Xi)*2))
  
}












#' @title Generate priors for case 4 dataset generation 
#' 
#' @description Generate priors for ABC data generation for case 4. Same as data generation for cases 1 and 2 but adds partitioning parameter distributed between 0 and 1
#' 
#' @param iterations numeric of number of simulations to run
#' @param mean_w_lo lower range of mean worm burden
#' @param mean_w_hi upper range of mean worm burden
#' @param disp_w_lo lower range of worm burden dispersion
#' @param disp_w_hi upper range of worm burden dispersion
#' @param h_lo lower range of daily mean egg release per mated pair
#' @param h_hi upper range of daily mean egg release per mated pair
#' @param r_lo lower range of daily egg release dispersion
#' @param r_hi upper range of daily egg release dispersion
#' @param g_lo lower range of density dependent fecundity parameter
#' @param g_hi upper range of density dependent fecundity parameter
#' 
#' @return matrix of priors
#' @export
#'
gen_case4_pars <- function(iterations,
                           mean_W_lo, mean_W_hi,
                           disp_W_lo, disp_W_hi, 
                           h_lo, h_hi, 
                           r_lo, r_hi,
                           g_lo, g_hi){
  
  mean_W <- exp(runif(iterations, log(mean_W_lo), log(mean_W_hi)))
  disp_W <- exp(runif(iterations, log(disp_W_lo), log(disp_W_hi)))
  part_W <- runif(iterations)
  h <- exp(runif(iterations, log(h_lo), log(h_hi)))
  r <- exp(runif(iterations, log(r_lo), log(r_hi)))
  g <- exp(runif(iterations, log(g_lo), log(g_hi)))
  
  return(cbind(mean_W, disp_W, part_W, h, r, g))
  
}



#' @title ABC fitting
#' 
#' @description Conduct approximate bayesian computation using four proposed data-generated mechanisms to estimate worm burdens from observed egg burdens under different data-generating mechanisms
#' 
#' @param obs_data observations on which to base the estimation from ZEST data
#' @param priors vector of ranges for prior distributions 
#' @param iterations number of synthetic datasets to generate and compare to observed data
#' 
#' @return list with abc objects that contain fitting statistics and generated datasets for each population
#' @export
#'

# Function to return ABC fits for all four cases for a given dataset ----------
  abc_fit <- function(obs_data, priors, iterations){
  # Extract priors  
    mean_W_lo = priors[["mean_W_lo"]] 
    mean_W_hi = priors[["mean_W_hi"]]
    disp_W_lo = priors[["disp_W_lo"]]
    disp_W_hi = priors[["disp_W_hi"]]
    susc_shape_lo = priors[["susc_shape_lo"]]  
    susc_shape_hi = priors[["susc_shape_hi"]]
    susc_rate_lo = priors[["susc_rate_lo"]] 
    susc_rate_hi = priors[["susc_rate_hi"]]
    mean_C_lo = priors[["mean_C_lo"]]
    mean_C_hi = priors[["mean_C_hi"]]
    disp_C_lo = priors[["disp_C_lo"]]
    disp_C_hi = priors[["disp_C_hi"]]
    h_lo = priors[["h_lo"]]
    h_hi = priors[["h_hi"]]
    r_lo = priors[["r_lo"]]
    r_hi = priors[["r_hi"]]
    g_lo = priors[["g_lo"]]
    g_hi = priors[["g_hi"]]
  
  #Extact observedsummary statistics
    obs_sum <- c(obs_data[["UF_mean"]], obs_data[["UF_se"]], obs_data[["UFpos2n"]])
    
  # Generate parameter sets for cases 1 and 2  
    pars12 <- gen_case12_pars(iterations = iterations,
                              mean_W_lo = mean_W_lo, mean_W_hi = mean_W_hi,
                              disp_W_lo = disp_W_lo, disp_W_hi = disp_W_hi,
                              h_lo = h_lo, h_hi = h_hi,
                              r_lo = r_lo, r_hi = r_hi,
                              g_lo = g_lo, g_hi = g_hi)
    
  # Generate parameter sets for case 3 
    pars3 <- gen_case3_pars(iterations = iterations,
                            susc_shape_lo = susc_shape_lo, susc_shape_hi = susc_shape_hi,
                            susc_rate_lo = susc_rate_lo, susc_rate_hi = susc_rate_hi,
                            mean_C_lo = mean_C_lo, mean_C_hi = mean_C_hi,
                            disp_C_lo = disp_C_lo, disp_C_hi = disp_C_hi, 
                            h_lo = h_lo, h_hi = h_hi, 
                            r_lo = r_lo, r_hi = r_hi,
                            g_lo = g_lo, g_hi = g_hi)
    
  # Generate parameter sets for case 4 
    pars4 <- gen_case4_pars(iterations = iterations,
                            mean_W_lo = mean_W_lo, mean_W_hi = mean_W_hi,
                            disp_W_lo = disp_W_lo, disp_W_hi = disp_W_hi,
                            h_lo = h_lo, h_hi = h_hi,
                            r_lo = r_lo, r_hi = r_hi,
                            g_lo = g_lo, g_hi = g_hi)
    
    
  # Simulate data and get summary statistics
    dat1 <- t(apply(X = pars12, 1, gen_case1_data, n = obs_data[["n_ppl"]]))
    dat2 <- t(apply(X = pars12, 1, gen_case2_data, n = obs_data[["n_ppl"]]))
    dat3 <- t(apply(X = pars3, 1, gen_case3_data, n = obs_data[["n_ppl"]]))
    dat4 <- t(apply(X = pars4, 1, gen_case4_data, n = obs_data[["n_ppl"]]))

  out_init <- tryCatch({
  # Use abc to get posteriors with ridge regression correction
    abc1 <- abc(target = obs_sum,
                param = pars12[,1:2],  # Mean and dispersion parameter
                sumstat = dat1[,1:3],  # Three summary stats
                tol = 100/iterations,
                transf = "log",
                method = "ridge")

    abc2 <- abc(target = obs_sum,
                param = pars12[,1:2],  # Mean and dispersion parameter
                sumstat = dat2[,1:3],  # Three summary stats
                tol = 100/iterations,
                transf = "log",
                method = "ridge")
  
    abc3 <- abc(target = obs_sum,
                param = pars3[,3:4],   # Mean and dispersion of cercarial exposure
                sumstat = dat3[,1:3],  # Three summary stats
                tol = 100/iterations,
                transf = "log",
                method = "ridge")

    abc4 <- abc(target = obs_sum,
                param = pars4[,1:3],   # Mean and dispersion parameters and partitioning parameter
                sumstat = dat4[,1:3],  # Three summary stats
                tol = 100/iterations,
                transf = "log",
                method = "ridge")
    
        
    abc_postpr <- postpr(target = obs_sum,
                         index = rep(c("I", "II", "III", "IV"), each = iterations),
                         sumstat = rbind(dat1[,1:3], dat2[,1:3], dat3[,1:3], dat4[,1:3]),
                         tol = 100/iterations,
                         method = "mnlogistic")
    
    abc1_wormstats <- dat1[which(abc1$region == TRUE),4:9]
    abc2_wormstats <- dat2[which(abc2$region == TRUE),4:9]
    abc3_wormstats <- dat3[which(abc3$region == TRUE),4:9]
    abc4_wormstats <- dat4[which(abc4$region == TRUE),4:9]
    abc3_parstats <- dat3[which(abc3$region == TRUE),10:16]
    
    list(abc1, abc2, abc3, abc4, abc_postpr,
         abc1_wormstats, abc2_wormstats, abc3_wormstats, abc4_wormstats,
         abc3_parstats)

  },
  error = function(cond){
    list(cond)
  })  

    return(c(list(obs_data[["Shehia"]], obs_data[["Year"]]), obs_data[["pop"]],
             out_init))

  }
