# SMP

`SMP` is an R package for sequential mixture prior analyses in platform trials with concurrent, nonconcurrent, and historical controls.

The package exports two main functions:

- `smp_binary()` for binary endpoints.
- `smp_continuous()` for continuous endpoints.

Both functions return posterior samples for the treatment and concurrent control arms, posterior summaries through `summary()`, and the SMP borrowing weights `w0`, `w1`, and `w2`.

See `vignettes/smp-demo.Rmd` for binary and continuous simulated examples with visualizations.