smp_binary <- function(nt, xt, nc, xc, nn, xn, nh, xh, delta,
					   n_samples = 5000L, pi2_prior = c("MAP", "pool"), seed = NULL) {
	pi2_prior <- match.arg(pi2_prior)
	if (!is.null(seed)) {
		set.seed(seed)
	}

	treatment_samples <- stats::rbeta(n_samples, xt + 1, nt - xt + 1)
	pi1_prior <- c(1, xn + 1, nn - xn + 1)

	if (pi2_prior == "MAP") {
		old_mc_control <- getOption("RBesT.MC.control")
		new_mc_control <- if (is.null(old_mc_control)) {
			list(adapt_delta = 0.999)
		} else {
			utils::modifyList(old_mc_control, list(adapt_delta = 0.999))
		}
		options(RBesT.MC.control = new_mc_control)
		on.exit(options(RBesT.MC.control = old_mc_control), add = TRUE)

		pi2_prior_fit <- suppressMessages(
			RBesT::automixfit(
				RBesT::gMAP(
					cbind(r, n - r) ~ 1 | study,
					family = stats::binomial(),
					data = data.frame(study = c(1, 2), n = c(nn, nh), r = c(xn, xh)),
					tau.dist = "HalfNormal",
					tau.prior = 1,
					beta.prior = 2
				)
			)
		)
	} else {
		pi2_prior_fit <- RBesT::mixbeta(c(1, xn + xh + 1, nn + nh - xn - xh + 1))
	}

	theta_c_hat <- xc / nc
	shifted <- theta_c_hat + c(-delta, delta)
	log_r1_num <- .log_binom_lik(xn, nn, theta_c_hat)
	log_r1_den <- max(vapply(shifted, function(p) .log_binom_lik(xn, nn, p), numeric(1)))
	r1 <- exp(log_r1_num - log_r1_den)

	shift_grid <- as.matrix(expand.grid(theta_n = shifted, theta_h = shifted))
	log_r2_num <- .log_binom_lik(xn, nn, theta_c_hat) + .log_binom_lik(xh, nh, theta_c_hat)
	log_r2_den <- max(apply(shift_grid, 1, function(p) {
		.log_binom_lik(xn, nn, p[1]) + .log_binom_lik(xh, nh, p[2])
	}))
	r2 <- exp(log_r2_num - log_r2_den)

	weights <- c(
		w0 = 1 / (1 + r1 + r2),
		w1 = r1 / (1 + r1 + r2),
		w2 = r2 / (1 + r1 + r2)
	)

	if (pi2_prior == "MAP") {
		control_prior <- do.call(
			RBesT::mixbeta,
			c(
				list(c(weights[["w0"]], 1, 1), c(weights[["w1"]], pi1_prior[2], pi1_prior[3])),
				lapply(seq_len(ncol(pi2_prior_fit)), function(j) {
					c(weights[["w2"]] * pi2_prior_fit["w", j], pi2_prior_fit["a", j], pi2_prior_fit["b", j])
				})
			)
		)
	} else {
		control_prior <- RBesT::mixbeta(
			c(weights[["w0"]], 1, 1),
			c(weights[["w1"]], pi1_prior[2], pi1_prior[3]),
			c(weights[["w2"]], xn + xh + 1, nn + nh - xn - xh + 1)
		)
	}

	control_posterior <- RBesT::postmix(control_prior, n = nc, r = xc)
	control_samples <- RBesT::rmix(control_posterior, n_samples)

	.as_smp_fit(
		endpoint = "binary",
		control_samples = control_samples,
		treatment_samples = treatment_samples,
		weights = weights,
		control_prior = control_prior,
		control_posterior = control_posterior,
		weight_statistics = c(r1 = r1, r2 = r2),
		inputs = list(
			nt = nt, xt = xt, nc = nc, xc = xc, nn = nn, xn = xn,
			nh = nh, xh = xh, delta = delta, n_samples = n_samples,
			pi2_prior = pi2_prior, seed = seed
		)
	)
}
