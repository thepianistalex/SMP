\name{smp_binary}
\alias{smp_binary}
\title{Sequential mixture prior analysis for a binary endpoint}
\usage{
smp_binary(nt, xt, nc, xc, nn, xn, nh, xh, delta,
  n_samples = 5000L, pi2_prior = c("MAP", "pool"), seed = NULL)
}
\arguments{
\item{nt, nc, nn, nh}{Sample sizes for treatment, concurrent control, nonconcurrent control, and historical control.}
\item{xt, xc, xn, xh}{Numbers of responders in treatment, concurrent control, nonconcurrent control, and historical control.}
\item{delta}{Clinically meaningful response-rate difference used in the borrowing weights.}
\item{n_samples}{Number of posterior samples to draw.}
\item{pi2_prior}{Prior construction for the full-borrowing component: \code{"MAP"} or \code{"pool"}.}
\item{seed}{Optional random seed for reproducible posterior sampling and MAP fitting.}
}
\value{
An object of class \code{smp_fit}, a list with posterior samples, borrowing weights \code{w0}, \code{w1}, and \code{w2}, fitted prior/posterior objects for the control arm, and input metadata.
}
\description{
Fits the sequential mixture prior for a binary endpoint using treatment, concurrent control, nonconcurrent control, and historical control counts.
}
\examples{
fit <- smp_binary(
  nt = 120, xt = 65, nc = 80, xc = 32,
  nn = 100, xn = 40, nh = 100, xh = 41,
  delta = 0.10, n_samples = 1000, pi2_prior = "pool", seed = 1
)
summary(fit)
fit$weights
}