% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/abc_worm_estimation_functions.R
\name{post_pred_data_gen}
\alias{post_pred_data_gen}
\title{ABC Posterior Predictive Check}
\usage{
post_pred_data_gen(pars, fixed_pars, n_ppl, weights, data_gen_fx, n_reps)
}
\arguments{
\item{pars}{vector of parameter values from `abc` that reach tolerance. Can be adjusted or unadjusted}

\item{fixed_pars}{vector of parameters needed for data generation but not drawn from posterior distribution}

\item{n_ppl}{number of people in community to generate data for}

\item{weights}{vector weights for each parameter set from `abc`}

\item{data_gen_fx}{data generation function (e.g. gen_case1_data; see above)}

\item{n_reps}{number of datasets to generate}
}
\value{
matrix with dim n_reps by number of summary statistics (3 for this analysis)
}
\description{
Generate data from parameter posterior distributions and compare to observed
}
