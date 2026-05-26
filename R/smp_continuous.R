smp_continuous <- function(nt, yt, nc, yc, nn, yn, nh, yh, delta,
						 n_samples = 5000L, pi2_prior = c("MAP", "pool"), seed = NULL) {
	pi2_prior <- match.arg(pi2_prior)
	if (!is.null(seed)) {
		set.seed(seed)
	}

	sigma2_hat <- stats::var(yc)
	sigma_hat <- sqrt(sigma2_hat)

	log_normal_lik <- function(y, mean_value, sd_value) {
		sum(stats::dnorm(y, mean = mean_value, sd = sd_value, log = TRUE))
	}

	treatment_samples <- .cont_posterior_samples(
		y = yt,
		sigma2 = stats::var(yt),
		n_draws = n_samples
	)

	pi0_components <- .cont_vague_prior_components()
	pi1_params <- .cont_posterior_params(
		y = yn,
		sigma2 = stats::var(yn)
	)
	pi1_components <- data.frame(weight = 1, mean = pi1_params$mean, sd = sqrt(pi1_params$var), stringsAsFactors = FALSE)

	if (pi2_prior == "MAP") {
		external_summary <- .cont_external_summary(nn, yn, nh, yh)
		pi2_prior_fit <- .cont_fit_map_prior(external_summary, sigma_hat)
		pi2_components <- .cont_components_from_prior(pi2_prior_fit)
	} else {
		pi2_params <- .cont_posterior_params(
			y = c(yn, yh),
			sigma2 = stats::var(c(yn, yh))
		)
		pi2_components <- data.frame(weight = 1, mean = pi2_params$mean, sd = sqrt(pi2_params$var), stringsAsFactors = FALSE)
	}

	theta_c_hat <- mean(yc)
	shifted <- theta_c_hat + c(-delta, delta)
	log_r1_num <- log_normal_lik(yn, theta_c_hat, sigma_hat)
	log_r1_den <- max(vapply(shifted, function(mu) log_normal_lik(yn, mu, sigma_hat), numeric(1)))
	log_r1 <- log_r1_num - log_r1_den

	shift_grid <- as.matrix(expand.grid(theta_n = shifted, theta_h = shifted))
	log_r2_num <- log_normal_lik(yn, theta_c_hat, sigma_hat) + log_normal_lik(yh, theta_c_hat, sigma_hat)
	log_r2_den <- max(apply(shift_grid, 1, function(mu_pair) {
		log_normal_lik(yn, mu_pair[1], sigma_hat) + log_normal_lik(yh, mu_pair[2], sigma_hat)
	}))
	log_r2 <- log_r2_num - log_r2_den

	log_weight_terms <- c(w0 = 0, w1 = log_r1, w2 = log_r2)
	log_weight_terms <- log_weight_terms - max(log_weight_terms)
	weights <- exp(log_weight_terms)
	weights <- weights / sum(weights)

	pi0_prior <- pi0_components
	pi1_prior <- pi1_components
	pi2_prior_weighted <- pi2_components
	pi0_prior$weight <- pi0_prior$weight * weights[["w0"]]
	pi1_prior$weight <- pi1_prior$weight * weights[["w1"]]
	pi2_prior_weighted$weight <- pi2_prior_weighted$weight * weights[["w2"]]

	control_prior <- rbind(pi0_prior, pi1_prior, pi2_prior_weighted)
	control_prior$weight <- control_prior$weight / sum(control_prior$weight)
	control_posterior <- .cont_update_mixture(control_prior, yc, sigma2_hat)
	control_samples <- .cont_sample_mixture(control_posterior, n_samples)

	.as_smp_fit(
		endpoint = "continuous",
		control_samples = control_samples,
		treatment_samples = treatment_samples,
		weights = weights,
		control_prior = control_prior,
		control_posterior = control_posterior,
		weight_statistics = c(r1 = exp(log_r1), r2 = exp(log_r2)),
		inputs = list(
			nt = nt, nc = nc, nn = nn, nh = nh, delta = delta,
			n_samples = n_samples, pi2_prior = pi2_prior, seed = seed
		)
	)
}
