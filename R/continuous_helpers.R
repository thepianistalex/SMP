.cont_vague_prior_spec <- function() {
	list(mean = 0, var = 1e6)
}

.cont_posterior_params <- function(y, sigma2, prior_mean = .cont_vague_prior_spec()$mean,
							   prior_var = .cont_vague_prior_spec()$var) {
	n_obs <- length(y)
	y_bar <- mean(y)
	prior_precision <- 1 / prior_var
	data_precision <- n_obs / sigma2
	posterior_var <- 1 / (prior_precision + data_precision)
	posterior_mean <- posterior_var * (prior_precision * prior_mean + data_precision * y_bar)
	list(mean = posterior_mean, var = posterior_var)
}

.cont_posterior_samples <- function(y, sigma2, n_draws, prior_mean = .cont_vague_prior_spec()$mean,
								prior_var = .cont_vague_prior_spec()$var) {
	params <- .cont_posterior_params(y, sigma2, prior_mean, prior_var)
	stats::rnorm(n_draws, mean = params$mean, sd = sqrt(params$var))
}

.cont_components_from_prior <- function(prior) {
	prior_matrix <- as.matrix(prior)
	if (nrow(prior_matrix) < 3L) {
		stop("Continuous normal prior must contain weight, mean, and sd rows.", call. = FALSE)
	}
	data.frame(
		weight = as.numeric(prior_matrix[1, ]),
		mean = as.numeric(prior_matrix[2, ]),
		sd = as.numeric(prior_matrix[3, ]),
		stringsAsFactors = FALSE
	)
}

.cont_update_mixture <- function(components, y, sigma2) {
	n_obs <- length(y)
	y_bar <- mean(y)
	data_precision <- n_obs / sigma2
	data_var <- 1 / data_precision
	prior_var <- components$sd ^ 2
	prior_precision <- 1 / prior_var
	posterior_var <- 1 / (prior_precision + data_precision)
	posterior_mean <- posterior_var * (prior_precision * components$mean + data_precision * y_bar)
	log_weight <- log(pmax(components$weight, .Machine$double.xmin)) +
		stats::dnorm(y_bar, mean = components$mean, sd = sqrt(prior_var + data_var), log = TRUE)
	log_weight <- log_weight - max(log_weight)
	posterior_weight <- exp(log_weight)
	posterior_weight <- posterior_weight / sum(posterior_weight)
	data.frame(weight = posterior_weight, mean = posterior_mean, sd = sqrt(posterior_var), stringsAsFactors = FALSE)
}

.cont_sample_mixture <- function(components, n_draws) {
	component_id <- sample.int(
		nrow(components),
		size = n_draws,
		replace = TRUE,
		prob = components$weight
	)
	stats::rnorm(n_draws, mean = components$mean[component_id], sd = components$sd[component_id])
}

.cont_vague_prior_components <- function(prior_mean = .cont_vague_prior_spec()$mean,
									prior_var = .cont_vague_prior_spec()$var) {
	data.frame(weight = 1, mean = prior_mean, sd = sqrt(prior_var), stringsAsFactors = FALSE)
}

.cont_external_summary <- function(nn, yn, nh, yh) {
	data.frame(
		study = c(1, 2),
		n = c(nn, nh),
		mean = c(mean(yn), mean(yh)),
		se = c(stats::sd(yn) / sqrt(nn), stats::sd(yh) / sqrt(nh)),
		stringsAsFactors = FALSE
	)
}

.cont_fit_map_prior <- function(external_summary, sigma_hat) {
	old_mc_control <- getOption("RBesT.MC.control")
	new_mc_control <- if (is.null(old_mc_control)) {
		list(adapt_delta = 0.999)
	} else {
		utils::modifyList(old_mc_control, list(adapt_delta = 0.999))
	}
	options(RBesT.MC.control = new_mc_control)
	on.exit(options(RBesT.MC.control = old_mc_control), add = TRUE)

	suppressMessages(
		RBesT::automixfit(
			RBesT::gMAP(
				cbind(mean, se) ~ 1 | study,
				weights = external_summary$n,
				data = external_summary,
				family = stats::gaussian(),
				beta.prior = cbind(0, sigma_hat),
				tau.dist = "HalfNormal",
				tau.prior = cbind(0, sigma_hat / 2)
			)
		)
	)
}