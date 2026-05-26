\name{SMP-package}
\alias{SMP-package}
\alias{SMP}
\docType{package}
\title{Sequential Mixture Prior Methods for Platform Trials}
\description{
The SMP package implements sequential mixture prior analyses for binary and continuous endpoints in platform trials with concurrent, nonconcurrent, and historical controls.
}
\details{
The main exported functions are \code{\link{smp_binary}} and \code{\link{smp_continuous}}. Both return posterior samples for the concurrent control and treatment arms, posterior summaries through \code{summary()}, and the data-adaptive borrowing weights \code{w0}, \code{w1}, and \code{w2}.
}