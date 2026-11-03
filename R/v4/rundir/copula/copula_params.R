# copula_params.R -- parameters of the Gaussian-copula family (S7e), chosen by the PI on 2026-10-23
# (docs/decisions/D-53.md) and pre-registered in docs/steps/S7e.md.
COPULA_RHO   <- 0.6   # latent (normal-score) correlation rho^|t_j - t_k| between visits of a subject
COPULA_FLOOR <- 0     # no subject-level constant component
