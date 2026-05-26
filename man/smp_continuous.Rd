\name{smp_continuous}
\alias{smp_continuous}
\title{Sequential mixture prior analysis for a continuous endpoint}
\usage{
smp_continuous(nt, yt, nc, yc, nn, yn, nh, yh, delta,
  n_samples = 5000L, pi2_prior = c("MAP", "pool"), seed = NULL)
}
\arguments{
\item{nt, nc, nn, nh}{Sample sizes for treatment, concurrent control, nonconcurrent control, and historical control.}
\item{yt, yc, yn, yh}{Numeric endpoint vectors for treatment, concurrent control, nonconcurrent control, and historical control.}
\item{delta}{Clinically meaningful mean difference used in the borrowing weights.}
\item{n_samples}{Number of posterior samples to draw.}
\item{pi2_prior}{Prior construction for the full-borrowing component: \code{"MAP"} or \code{"pool"}.}
\item{seed}{Optional random seed for reproducible posterior sampling and MAP fitting.}
}
\value{
An object of class \code{smp_fit}, a list with posterior samples, borrowing weights \code{w0}, \code{w1}, and \code{w2}, fitted control prior/posterior mixture components, and input metadata.
}
\description{
Fits the sequential mixture prior for a continuous endpoint using treatment, concurrent control, nonconcurrent control, and historical control observations.
}
\examples{
set.seed(1)
yc <- rnorm(80, 0, 3)
yt <- rnorm(120, 0.6, 3)
yn <- rnorm(100, 0, 3)
yh <- rnorm(100, 0.2, 3)
fit <- smp_continuous(
  nt = length(yt), yt = yt, nc = length(yc), yc = yc,
  nn = length(yn), yn = yn, nh = length(yh), yh = yh,
  delta = 0.60, n_samples = 1000, pi2_prior = "pool", seed = 1
)
summary(fit)
fit$weights
}