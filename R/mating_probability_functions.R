#' @title Schistosomiasis Case 1 mating probability
#'
#' @description Estimate the mating probability from the population mean worm burden and the
#' clumping parameter of the negative binomial distribution assuming male and female
#' parasites distributed together from the same negative binomial dist'n
#'
#' @param W mean population worm burden
#' @param k clumping parameter of the negative binomial distribution
#'
#' @return probability a female worm is succefully mated (ranges from 0-1)
#' @export

phi_Wk <- function(W, k) {
  if(W == 0){
    return(0)
  } else {
    b <- ((1-(W/(W + k)))^(1+k))/(2*pi)
    i = integrate(f = function(x, W, k){(1-cos(x))/((1 + (W/(W + k))*cos(x))^(1+k))},
                  lower = 0,
                  upper = 2*pi,
                  stop.on.error = FALSE,
                  W = W, k = k)$value
    
    return(1-b*i)
  }
}

#' @title Function to estimate the probability of an individual host having j mated pairs
#' 
#' @description derived in Anderson 1977 and Bradley and Anderson 1978 for schistosomes of opposite 
#' sexes drawn from separate negative binomial distributions (Case II: Distributed Separate)
#' 
#' @param j number of mated pairs
#' @param w mean worm burden in human population
#' @param k aggregation parameter of negative binomially distributed worms
#' 
#' @return estimate of the probability of j mated pairs in a human host
#' @export

prob_j_pairs <- function(j, w, k){
  #First part of sum in which i=j and therefore theta_ij=1
  i_eq_j = dnbinom(j, mu = w/2, size = k)
  
  #rest of the summation where i does not = j and therefore theta_ij=2
  #Use 2^14 rather than infinity means result is approximate, but accurate to several decimal points
  i_jp1_inf = sum(2*dnbinom((j+1):2^12, mu = w/2, size = k))
  
  dnbinom(j, mu = w/2, size = k) * (i_eq_j + i_jp1_inf)
}

#' @title Function to estimate mating probability for sexes distributed separately
#' 
#' @description derived in Anderson 1977 and Bradley and Anderson 1978 for schistosomes of opposite 
#' sexes drawn from separate negative binomial distributions (Case II: Distributed Separate)
#' This function requires `prob_j_pairs` function
#' 
#' @param w mean worm burden in human population
#' @param k aggregation parameter of negative binomially distributed worms  
#' 
#' @return estimate of the Case II mating probability  
#' @export 

phi_wk_sep <- function(w, k){
  sum(sapply(0:2^12, function(j){2*j*prob_j_pairs(j = j, w = w, k = k)}))/w
}
