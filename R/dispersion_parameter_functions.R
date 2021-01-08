#' Get aggregation parameter from observed data  
#' 
#' Estimates the aggregation parameter, `kappa`, from vector of observed values (egg or worm burden)
#' based on estimation of the mean and variance as: 
#' $$\kappa=\frac{\mu^2-\frac{\sigma^2}{n}}{\sigma^2-\mu}$$
#' 
#' @param egg_vec vector of individual egg or worm burden estimates
#' 
#' @return estimate of the aggregation parameter of the negative binomial dist'n
#' @export
#' 

mean_var_agg_crctd <- function(vec){
  
  (mean(vec, na.rm = T)^2-var(vec, na.rm = T)/sum(!is.na(vec))) / 
    (var(vec, na.rm = T)-mean(vec, na.rm = T))
  
}

#' Estimate aggregation parameter, kappa, as function of mean infection intensity and prevalence using uniroot
#'
#' Prevalence, mean intensity and the aggregation parameter are all related,
#' therefore estimation of 1 can be achieved if the other two are known
#'
#' @param W Mean worm burden or infection intensity
#' @param prev prevalence of infection
#'
#' @return Estimate of the aggregation parameter
#' @export

prev_W_get_k <- function(W, prev){
  uniroot(function(k){ 1-(1+W/k)^-k-prev},
          interval = c(0,10))$root
}

#' Estimate aggregation parameter from mean-variance relationship  
#' 
#' @param vec vector of individual egg or worm burden intensities  
#' 
#' @return estimate of the aggregation parameter
#' @export
#'  
mean_var_agg_par <- function(vec){
  mean(vec, na.rm = T)^2/(var(vec, na.rm = T) - mean(vec, na.rm = T))
}

#' Estimate MLE dispersion parameter using `fitdistrplus::fitdist`  
#' 
#' Uses fitdist function from fitdistrplus package to estimate the dispersion parameter
#' of vector of count data  
#'  
#' @param vec vector of individual egg or worm burden intensities  
#' 
#' @return estimate of the dispersion parameter
#' @export

mle_disp_par <- function(vec){
  tryCatch(expr = {as.numeric(fitdistrplus::fitdist(as.numeric(na.omit(vec)),
                                   distr = "nbinom")$estimate[1])},
           error = function(err_mess){
             return(NA_real_)
           })
}

#' Estimate MLE standard error of dispersion parameter using `fitdistrplus::fitdist`  
#' 
#' Uses fitdist function from fitdistrplus package to estimate the standard error of the dispersion parameter
#' of vector of count data  
#'  
#' @param vec vector of individual egg or worm burden intensities  
#' 
#' @return estimate of the standard error of the dispersion parameter
#' @export

mle_disp_par_sd <- function(vec){
  tryCatch(expr = {as.numeric(fitdistrplus::fitdist(as.numeric(na.omit(vec)),
                                   distr = "nbinom")$sd[1])},
           error = function(err_mess){
             return(NA_real_)
           })
}

#' Estimate MLE mean parameter using `fitdistrplus::fitdist`  
#' 
#' Uses fitdist function from fitdistrplus package to estimate the mean
#' of vector of count data  
#'  
#' @param vec vector of individual egg or worm burden intensities  
#' 
#' @return estimate of the mean parameter
#' @export

mle_mean_par <- function(vec){
  tryCatch(expr = {as.numeric(fitdistrplus::fitdist(as.numeric(na.omit(vec)),
                                   distr = "nbinom")$estimate[2])},
           error = function(err_mess){
             return(NA_real_)
           })
}

#' Estimate MLE standard error of mean parameter using `fitdistrplus::fitdist`  
#' 
#' Uses fitdist function from fitdistrplus package to estimate the standard error of the mean parameter
#' of vector of count data  
#'  
#' @param vec vector of individual egg or worm burden intensities  
#' 
#' @return estimate of the standard error of the dispersion parameter
#' @export

mle_mean_par_sd <- function(vec){
  tryCatch(expr = {as.numeric(fitdistrplus::fitdist(as.numeric(na.omit(vec)),
                                   distr = "nbinom")$sd[2])},
           error = function(err_mess){
             return(NA_real_)
           })
}

#' Custom max likelihood function for dispersion parameter 
#' 
#' Estimates the likelihood of a vector of egg counts with negative binomial distribution 
#' fixed mean and inverse of the dispersion parameter, kappa, frequently denoted alpha
#' can be used in conjunction with an opimizer function to find max likelihood of the
#' dispersion parameter, implemented here in the `mle_kap` function 
#' 
#' @param alpha value of the inverse of the dispersion parameter, kappa
#' @param vec vector of individual egg or worm burden intensities
#' 
#' @return estimate of the negative log likelihood
#' @export

alpha_likelihood <- function(alpha, vec){
  mu_hat <- mean(vec, na.rm = TRUE)
  
  L = dnbinom(vec, mu = mu_hat, size = 1/alpha)
  
  -sum(log(L), na.rm = T)
}

#' Finds max likelihood of the dispersion parameter using `alpha_likelihood`
#' 
#' Finds the maximum likelihood estimate of the inverse of the dispersion parameter,
#' kappa, often denoted alpha, given a vector of egg or worm burdens. Assumes max likelihood
#' estimate of the mean of the negative binomial distribution is the empirical mean of the vector
#' Returns the message from `optim` if convergence fails
#' 
#' @param vec vector of individual egg or worm burden intensities
#' 
#' @return max likelihood estimate of the dispersion parameter
#' @export

mle_kap <- function(vec){
  if(sum(vec != 0, na.rm = T) == 0){
    return(NA_real_)
  } else {
    op_est <- optim(par = list("alpha" = 0.05),
                    fn = alpha_likelihood, 
                    method = "Brent",
                    lower = 1e-16, upper = 1e6, 
                    hessian = T,
                    vec = vec)
    
    if(op_est$convergence == 0){
      
      kap_est <- 1/op_est$par
      
      return(kap_est)
      
    } else {
      return(op_est$message)
    }
  }
}

#' Finds max likelihood of the dispersion parameter using `alpha_likelihood`
#' 
#' Finds the standard error of the maximum likelihood estimate of the inverse of the 
#' dispersion parameter,kappa, often denoted alpha, given a vector of egg or worm burdens. 
#' Uses sqrt of the inverse of the hessian returned by `optim` to estimate standard error of alpha
#' then transforms to estimate the standard error of kappa.
#' Assumes max likelihood estimate of the mean of the negative binomial distribution is the 
#' empirical mean of the vector
#' Returns the message from `optim` if convergence fails
#' 
#' @param vec vector of individual egg or worm burden intensities
#' 
#' @return standard error of the max likelihood estimate of the dispersion parameter
#' @export

mle_kap_se <- function(vec){
  if(sum(vec != 0, na.rm = T) == 0){
    return(NA_real_)
  } else {
    op_est <- optim(par = list("alpha" = 0.05),
                    fn = alpha_likelihood, 
                    method = "Brent",
                    lower = 1e-16, upper = 1e6, 
                    hessian = T,
                    vec = vec)
    
    if(op_est$convergence == 0){
      
      alpha <- op_est$par
      kappa <- alpha^-1
      alpha_se <- sqrt(diag(op_est$hessian^-1))
      alpha_lo <- alpha - alpha_se
      kappa_up <- alpha_lo^-1
      kappa_se <- kappa_up - kappa
      
      return(kappa_se)
      
    } else {
      return(op_est$message)
    }
  }
}

#' Finds max likelihood of the inverse of the dispersion parameter, alpha, dispersion parameter using `alpha_likelihood`
#' 
#' Finds the maximum likelihood estimate of the inverse of the dispersion parameter,
#' kappa, often denoted alpha, given a vector of egg or worm burdens. Assumes max likelihood
#' estimate of the mean of the negative binomial distribution is the empirical mean of the vector
#' Returns the message from `optim` if convergence fails
#' 
#' @param vec vector of individual egg or worm burden intensities
#' 
#' @return max likelihood estimate of the dispersion parameter
#' @export

mle_alpha <- function(vec){
  if(sum(vec != 0, na.rm = T) == 0){
    return(NA_real_)
  } else {
    op_est <- optim(par = list("alpha" = 0.05),
                    fn = alpha_likelihood, 
                    method = "Brent",
                    lower = 1e-16, upper = 1e6, 
                    hessian = T,
                    vec = vec)
    
    if(op_est$convergence == 0){
      
      alpha_est <- op_est$par
      
      return(alpha_est)
      
    } else {
      return(op_est$message)
    }
  }
}

#' Finds max likelihood of the dispersion parameter using `alpha_likelihood`
#' 
#' Finds the standard error of the maximum likelihood estimate of the inverse of the 
#' dispersion parameter,kappa, often denoted alpha, given a vector of egg or worm burdens. 
#' Uses sqrt of the inverse of the hessian returned by `optim` to estimate standard error of alpha
#' then transforms to estimate the standard error of kappa.
#' Assumes max likelihood estimate of the mean of the negative binomial distribution is the 
#' empirical mean of the vector
#' Returns the message from `optim` if convergence fails
#' 
#' @param vec vector of individual egg or worm burden intensities
#' 
#' @return standard error of the max likelihood estimate of the dispersion parameter
#' @export

mle_alpha_se <- function(vec){
  if(sum(vec != 0, na.rm = T) == 0){
    return(NA_real_)
  } else {
    op_est <- optim(par = list("alpha" = 0.05),
                    fn = alpha_likelihood, 
                    method = "Brent",
                    lower = 1e-16, upper = 1e6, 
                    hessian = T,
                    vec = vec)
    
    if(op_est$convergence == 0){
      
      alpha_se <- sqrt(diag(op_est$hessian^-1))

      return(alpha_se)
      
    } else {
      return(op_est$message)
    }
  }
}
