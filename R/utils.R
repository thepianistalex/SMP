.clamp_prob <- function(p) {
	pmin(pmax(p, 1e-8), 1 - 1e-8)
}

.log_binom_lik <- function(x, n, p) {
	stats::dbinom(x, size = n, prob = .clamp_prob(p), log = TRUE)
}

.as_smp_fit <- function(endpoint, control_samples, treatment_samples, weights, inputs,
						control_prior = NULL, control_posterior = NULL, weight_statistics = NULL) {
	result <- list(
		endpoint = endpoint,
		control_samples = control_samples,
		treatment_samples = treatment_samples,
		weights = weights,
		weight_statistics = weight_statistics,
		control_prior = control_prior,
		control_posterior = control_posterior,
		inputs = inputs
	)
	class(result) <- "smp_fit"
	result
}

summary.smp_fit <- function(object, cred_mass = 0.95, ...) {
	alpha <- (1 - cred_mass) / 2
	diff_samples <- object$treatment_samples - object$control_samples
	parameter_names <- if (identical(object$endpoint, "binary")) {
		c("theta_C", "theta_T", "theta_T - theta_C")
	} else {
		c("mu_C", "mu_T", "mu_T - mu_C")
	}
	posterior <- data.frame(
		parameter = parameter_names,
		mean = c(mean(object$control_samples), mean(object$treatment_samples), mean(diff_samples)),
		lower = c(
			stats::quantile(object$control_samples, alpha, names = FALSE),
			stats::quantile(object$treatment_samples, alpha, names = FALSE),
			stats::quantile(diff_samples, alpha, names = FALSE)
		),
		upper = c(
			stats::quantile(object$control_samples, 1 - alpha, names = FALSE),
			stats::quantile(object$treatment_samples, 1 - alpha, names = FALSE),
			stats::quantile(diff_samples, 1 - alpha, names = FALSE)
		),
		stringsAsFactors = FALSE
	)
	weights <- data.frame(
		component = names(object$weights),
		weight = as.numeric(object$weights),
		stringsAsFactors = FALSE
	)
	result <- list(
		endpoint = object$endpoint,
		posterior = posterior,
		prob_superiority = mean(diff_samples > 0),
		weights = weights,
		cred_mass = cred_mass
	)
	class(result) <- "summary.smp_fit"
	result
}

print.smp_fit <- function(x, ...) {
	cat("Sequential mixture prior fit (", x$endpoint, " endpoint)\n", sep = "")
	cat("Posterior draws:", length(x$control_samples), "\n")
	cat("Borrowing weights:\n")
	print(data.frame(component = names(x$weights), weight = as.numeric(x$weights)), row.names = FALSE)
	invisible(x)
}

print.summary.smp_fit <- function(x, ...) {
	cat("Sequential mixture prior posterior summary\n")
	print(x$posterior, row.names = FALSE)
	cat("\nPosterior probability of superiority:", format(x$prob_superiority, digits = 4), "\n")
	cat("\nBorrowing weights:\n")
	print(x$weights, row.names = FALSE)
	invisible(x)
}