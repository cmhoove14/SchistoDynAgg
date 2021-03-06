% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/abc_worm_estimation_functions.R
\name{gen_case2_data}
\alias{gen_case2_data}
\title{Generate synthetic case 2 dataset}
\usage{
gen_case2_data(pars, n)
}
\arguments{
\item{pars}{numeric vector with mean worm burden, inverse dispersion parameter, mean eggs per mated pair, and dispersion of daily egg release}

\item{n}{number of individual human hosts}
}
\value{
vector of summary statistics to compare to observed in ABC setup
}
\description{
Generate worm pairs and resulting egg burdens for individuals assuming males and females distributed from separate negative binomial distributions, return summary statistics. All prior distributions are drawn from uniform distributions
}
