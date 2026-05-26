\name{summary.smp_fit}
\alias{summary.smp_fit}
\title{Summarize an SMP fit}
\usage{
\method{summary}{smp_fit}(object, cred_mass = 0.95, ...)
}
\arguments{
\item{object}{An object returned by \code{\link{smp_binary}} or \code{\link{smp_continuous}}.}
\item{cred_mass}{Credible interval mass. The default is 0.95.}
\item{...}{Additional arguments, currently unused.}
}
\value{
A list containing posterior means and credible intervals, the posterior probability of superiority, and the SMP borrowing weights.
}
\description{
Computes posterior summaries for the concurrent control, treatment arm, and treatment-control difference.
}