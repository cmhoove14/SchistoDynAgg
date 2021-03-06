% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dispersion_parameter_functions.R
\name{mle_disp_par}
\alias{mle_disp_par}
\title{Estimate MLE dispersion parameter using `fitdistrplus::fitdist`}
\usage{
mle_disp_par(vec)
}
\arguments{
\item{vec}{vector of individual egg or worm burden intensities}
}
\value{
estimate of the dispersion parameter
}
\description{
Uses fitdist function from fitdistrplus package to estimate the dispersion parameter
of vector of count data
}
