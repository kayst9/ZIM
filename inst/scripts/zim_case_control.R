drop_intercept <- function(X) {
  has_intercept <- colnames(X) == "(Intercept)"
  if (any(has_intercept)) X <- X[, !has_intercept, drop = FALSE]
  X
}

rhs_model_matrix <- function(formula, data, keep_intercept = FALSE) {
  X <- stats::model.matrix(formula, data = data)
  if (!keep_intercept) X <- drop_intercept(X)
  storage.mode(X) <- "double"
  X
}

check_binary_y <- function(y) {
  if (is.factor(y)) y <- as.integer(y) - 1L
  y <- as.numeric(y)
  if (!all(y %in% c(0, 1))) stop("Y must be coded as 0/1, or a two-level factor.")
  y
}

ZIM_DEFAULT_RIDGE <- 1e-8
ZIM_DEFAULT_EL_MASS_TOL <- 1e-8
ZIM_DEFAULT_COMPACT_CENTERING_TOL <- 1e-8
ZIM_DEFAULT_PSD_TOL <- 1e-8
ZIM_DEFAULT_BREAD_RANK_TOL <- 1e-10
ZIM_DEFAULT_CONDITION_WARN <- 1e12

validate_ridge <- function(ridge, argument = "ridge") {
  if (!is.character(argument) || length(argument) != 1L ||
      is.na(argument) || !nzchar(argument)) {
    stop("argument must be a nonempty character scalar.")
  }
  if (!is.numeric(ridge) || length(ridge) != 1L ||
      is.na(ridge) || !is.finite(ridge) || ridge < 0) {
    stop(argument, " must be a nonnegative finite numeric scalar.")
  }
  as.numeric(ridge)
}

validate_shrinkage_controls <- function(inverse_method, shrink_lambda,
                                        beta_method = NULL,
                                        shrink_structured = TRUE) {
  if (!identical(inverse_method, "shrinkage")) {
    if (!is.null(shrink_lambda)) {
      stop("shrink_lambda is used only with inverse_method = 'shrinkage'.")
    }
    return(NULL)
  }
  if (!is.null(beta_method) && !identical(beta_method, "GMM")) {
    stop("inverse_method = 'shrinkage' is available only for beta_method = 'GMM'; WKM is unchanged.")
  }
  if (!is.numeric(shrink_lambda) || length(shrink_lambda) != 1L ||
      is.na(shrink_lambda) || !is.finite(shrink_lambda) ||
      shrink_lambda < 0 || shrink_lambda > 1) {
    stop("inverse_method = 'shrinkage' requires shrink_lambda to be a finite numeric scalar in [0, 1].")
  }
  if (!isTRUE(shrink_structured)) {
    stop("inverse_method = 'shrinkage' always preserves the target component; shrink_structured must be TRUE.")
  }
  as.numeric(shrink_lambda)
}

validate_positive_tolerance <- function(x, argument) {
  if (!is.character(argument) || length(argument) != 1L ||
      is.na(argument) || !nzchar(argument)) {
    stop("argument must be a nonempty character scalar.")
  }
  if (!is.numeric(x) || length(x) != 1L || is.na(x) ||
      !is.finite(x) || x <= 0) {
    stop(argument, " must be a positive finite numeric scalar.")
  }
  as.numeric(x)
}

matrix_rank_diagnostics <- function(A,
                                    rank_tol = ZIM_DEFAULT_BREAD_RANK_TOL) {
  rank_tol <- validate_positive_tolerance(rank_tol, "rank_tol")
  A <- as.matrix(A)
  if (nrow(A) != ncol(A)) stop("A must be square for rank diagnostics.")
  if (any(!is.finite(A))) {
    return(list(
      rank = 0L,
      dimension = ncol(A),
      rank_tol = rank_tol,
      rank_cutoff = NA_real_,
      min_singular_value = NA_real_,
      max_singular_value = NA_real_,
      condition_number = Inf,
      singular_values = rep(NA_real_, ncol(A)),
      full_rank = FALSE
    ))
  }
  singular_values <- svd(A, nu = 0, nv = 0)$d
  max_sv <- if (length(singular_values) > 0L) max(singular_values) else 0
  cutoff <- rank_tol * max(max_sv, .Machine$double.eps)
  rank <- sum(singular_values > cutoff)
  min_sv <- if (length(singular_values) > 0L) min(singular_values) else 0
  condition_number <- if (rank == ncol(A) && min_sv > 0) max_sv / min_sv else Inf
  list(
    rank = rank,
    dimension = ncol(A),
    rank_tol = rank_tol,
    rank_cutoff = cutoff,
    min_singular_value = min_sv,
    max_singular_value = max_sv,
    condition_number = condition_number,
    singular_values = singular_values,
    full_rank = rank == ncol(A)
  )
}

warn_on_condition_number <- function(diagnostics,
                                     threshold = ZIM_DEFAULT_CONDITION_WARN,
                                     context = "Matrix") {
  threshold <- validate_positive_tolerance(threshold, "condition_warn")
  if (!is.list(diagnostics) || is.null(diagnostics$condition_number)) {
    stop("diagnostics must be returned by matrix_rank_diagnostics().")
  }
  kappa_value <- diagnostics$condition_number
  if (is.finite(kappa_value) && kappa_value > threshold) {
    warning(
      context, " is poorly conditioned: condition number = ",
      signif(kappa_value, 6), ", warning threshold = ",
      signif(threshold, 6), "."
    )
  }
  invisible(diagnostics)
}

matrix_psd_diagnostics <- function(A, tol = ZIM_DEFAULT_PSD_TOL) {
  tol <- validate_positive_tolerance(tol, "tol")
  A <- symmetrize(as.matrix(A))
  if (nrow(A) != ncol(A)) stop("A must be square for PSD diagnostics.")
  if (any(!is.finite(A))) {
    return(list(
      minimum_eigenvalue = NA_real_,
      maximum_eigenvalue = NA_real_,
      psd_tolerance = tol,
      psd_cutoff = NA_real_,
      is_psd = FALSE,
      eigenvalues = rep(NA_real_, nrow(A))
    ))
  }
  values <- eigen(A, symmetric = TRUE, only.values = TRUE)$values
  scale <- max(1, max(abs(values)))
  cutoff <- tol * scale
  list(
    minimum_eigenvalue = min(values),
    maximum_eigenvalue = max(values),
    psd_tolerance = tol,
    psd_cutoff = cutoff,
    is_psd = min(values) >= -cutoff,
    eigenvalues = values
  )
}

validate_joint_covariance_blocks <- function(H_theta_theta,
                                             H_theta_beta,
                                             H_beta_beta,
                                             tol = ZIM_DEFAULT_PSD_TOL,
                                             context = "Joint covariance") {
  H_theta_theta <- symmetrize(as.matrix(H_theta_theta))
  H_theta_beta <- as.matrix(H_theta_beta)
  H_beta_beta <- symmetrize(as.matrix(H_beta_beta))
  if (nrow(H_theta_theta) != ncol(H_theta_theta)) stop("H_theta_theta must be square.")
  if (nrow(H_beta_beta) != ncol(H_beta_beta)) stop("H_beta_beta must be square.")
  if (nrow(H_theta_beta) != nrow(H_theta_theta) ||
      ncol(H_theta_beta) != nrow(H_beta_beta)) {
    stop("H_theta_beta dimensions do not match the diagonal covariance blocks.")
  }
  H_joint <- rbind(
    cbind(H_theta_theta, H_theta_beta),
    cbind(t(H_theta_beta), H_beta_beta)
  )
  diagnostics <- matrix_psd_diagnostics(H_joint, tol = tol)
  if (!isTRUE(diagnostics$is_psd)) {
    stop(
      context, " is not positive semidefinite within tolerance: minimum eigenvalue = ",
      signif(diagnostics$minimum_eigenvalue, 6),
      ", allowed negative cutoff = -", signif(diagnostics$psd_cutoff, 6), "."
    )
  }
  diagnostics
}

validate_el_masses <- function(p_hat,
                               Delta_hat,
                               mass_tol = ZIM_DEFAULT_EL_MASS_TOL,
                               action = c("stop", "warn", "none"),
                               context = "EL reconstruction") {
  mass_tol <- validate_positive_tolerance(mass_tol, "mass_tol")
  action <- match.arg(action)
  p_hat <- as.numeric(p_hat)
  Delta_hat <- as.numeric(Delta_hat)
  if (length(p_hat) != length(Delta_hat)) stop("p_hat and Delta_hat must have the same length.")
  if (any(!is.finite(p_hat)) || any(p_hat < 0) ||
      any(!is.finite(Delta_hat)) || any(Delta_hat <= 0)) {
    stop("EL masses require finite nonnegative p_hat and finite positive Delta_hat.")
  }
  mass0 <- sum(p_hat)
  mass1 <- sum(p_hat * Delta_hat)
  max_abs_mass_error <- max(abs(c(mass0 - 1, mass1 - 1)))
  valid <- is.finite(mass0) && mass0 > 0 &&
    is.finite(mass1) && mass1 > 0 &&
    is.finite(max_abs_mass_error) && max_abs_mass_error <= mass_tol
  if (!valid && !identical(action, "none")) {
    msg <- paste0(
      context,
      ": reconstructed control/case masses are not sufficiently close to one; ",
      "mass0 = ", signif(mass0, 8),
      ", mass1 = ", signif(mass1, 8),
      ", max error = ", signif(max_abs_mass_error, 8),
      ", tolerance = ", signif(mass_tol, 8), "."
    )
    if (identical(action, "stop")) stop(msg) else warning(msg)
  }
  list(
    mass0 = mass0,
    mass1 = mass1,
    max_abs_mass_error = max_abs_mass_error,
    mass_tol = mass_tol,
    action = action,
    valid = valid
  )
}

ZIM_PCA_FAMILY_METHODS <- c("pca", "hbb_pca", "sigma_pca", "m_factor_diag", "hbb_factor_diag")

is_pca_family_method <- function(inverse_method) {
  is.character(inverse_method) && length(inverse_method) == 1L &&
    !is.na(inverse_method) && inverse_method %in% ZIM_PCA_FAMILY_METHODS
}

validate_pca_controls <- function(
    pca_rule = c("variance", "condition", "condition_or_variance", "fixed"),
    pc_var = 0.99,
    n_pc = NULL,
    pc_max_condition = 1e4,
    pca_scale = c("covariance", "correlation"),
    pca_rule_was_missing = FALSE,
    active = TRUE) {
  pca_rule <- match.arg(pca_rule)
  pca_scale <- match.arg(pca_scale)

  inferred_from_n_pc <- isTRUE(pca_rule_was_missing) && !is.null(n_pc)
  if (inferred_from_n_pc) pca_rule <- "fixed"

  if (!isTRUE(active)) {
    return(list(
      pca_rule = pca_rule,
      pc_var = pc_var,
      n_pc = n_pc,
      pc_max_condition = pc_max_condition,
      pca_scale = pca_scale,
      inferred_from_n_pc = inferred_from_n_pc
    ))
  }

  uses_pc_var <- pca_rule %in% c("variance", "condition_or_variance")
  uses_pc_max_condition <- pca_rule %in% c("condition", "condition_or_variance")

  if (uses_pc_var) {
    if (!is.numeric(pc_var) || length(pc_var) != 1L || is.na(pc_var) ||
        !is.finite(pc_var) || pc_var <= 0 || pc_var > 1) {
      stop(
        "pc_var must be a finite scalar in (0, 1] when pca_rule uses ",
        "the variance rule."
      )
    }
    pc_var <- as.numeric(pc_var)
  }
  if (uses_pc_max_condition) {
    if (!is.numeric(pc_max_condition) || length(pc_max_condition) != 1L ||
        is.na(pc_max_condition) || !is.finite(pc_max_condition) ||
        pc_max_condition < 1) {
      stop(
        "pc_max_condition must be a finite scalar >= 1 when pca_rule uses ",
        "the condition rule."
      )
    }
    pc_max_condition <- as.numeric(pc_max_condition)
  }

  if (pca_rule %in% c("variance", "condition", "condition_or_variance")) {
    if (!is.null(n_pc)) {
      stop("n_pc is used only when pca_rule = 'fixed'.")
    }
  } else {
    if (is.null(n_pc)) {
      stop("n_pc must be supplied when pca_rule = 'fixed'.")
    }
    if (!is.numeric(n_pc) || length(n_pc) != 1L || is.na(n_pc) ||
        !is.finite(n_pc) || n_pc < 1 ||
        abs(n_pc - round(n_pc)) > sqrt(.Machine$double.eps) ||
        n_pc > .Machine$integer.max) {
      stop("n_pc must be a positive integer when pca_rule = 'fixed'.")
    }
    n_pc <- as.integer(round(n_pc))
  }

  list(
    pca_rule = pca_rule,
    pc_var = pc_var,
    n_pc = n_pc,
    pc_max_condition = pc_max_condition,
    pca_scale = pca_scale,
    inferred_from_n_pc = inferred_from_n_pc
  )
}

safe_exp <- function(eta, limit = 35) {
  if (length(limit) != 1L || !is.finite(limit) || limit <= 0) {
    stop("limit must be a positive finite scalar.")
  }
  exp(pmin(pmax(eta, -limit), limit))
}

safe_exp_derivative <- function(eta, limit = 35) {

  eta <- as.numeric(eta)
  out <- safe_exp(eta, limit = limit)
  out[!is.finite(eta) | abs(eta) >= limit] <- 0
  out
}

log_mean_exp <- function(x) {
  x <- as.numeric(x)
  if (length(x) == 0L || any(!is.finite(x))) {
    stop("log_mean_exp requires a nonempty finite numeric vector.")
  }
  x_max <- max(x)
  x_max + log(mean(exp(x - x_max)))
}

check_exp_range <- function(eta,
                            context,
                            limit = 35,
                            action = c("none", "warn", "stop")) {
  action <- match.arg(action)
  if (length(context) != 1L || is.na(context) || !nzchar(context)) {
    stop("context must be a nonempty character scalar.")
  }
  if (length(limit) != 1L || !is.finite(limit) || limit <= 0) {
    stop("limit must be a positive finite scalar.")
  }

  eta <- as.numeric(eta)
  finite <- is.finite(eta)
  clipped_low <- finite & eta <= -limit
  clipped_high <- finite & eta >= limit
  n_nonfinite <- sum(!finite)
  n_clipped_low <- sum(clipped_low)
  n_clipped_high <- sum(clipped_high)
  n_clipped <- n_clipped_low + n_clipped_high
  max_abs_eta <- if (any(finite)) max(abs(eta[finite])) else Inf
  clipping_active <- n_nonfinite > 0L || n_clipped > 0L

  diagnostic <- list(
    context = context,
    limit = limit,
    n_total = length(eta),
    n_nonfinite = n_nonfinite,
    n_clipped_low = n_clipped_low,
    n_clipped_high = n_clipped_high,
    n_clipped = n_clipped,
    max_abs_eta = max_abs_eta,
    clipping_active = clipping_active
  )

  if (clipping_active && !identical(action, "none")) {
    msg <- paste0(
      context,
      ": exponential clipping or a non-finite linear predictor was detected; ",
      "n_clipped = ", n_clipped,
      ", n_nonfinite = ", n_nonfinite,
      ", max(abs(eta)) = ", signif(max_abs_eta, 6),
      ", clipping limit = ", limit,
      ". The untruncated exponential-model formulas and analytic Jacobian are ",
      "not valid at a clipped final estimate."
    )
    if (identical(action, "stop")) stop(msg) else warning(msg)
  }

  diagnostic
}

symmetrize <- function(A) {
  (A + t(A)) / 2
}

regularize_negative_bread <- function(A, ridge) {
  ridge <- validate_ridge(ridge)

  A <- as.matrix(A)
  if (nrow(A) != ncol(A)) stop("A must be square.")
  A <- symmetrize(A)
  if (ridge <= 0) return(A)
  scale <- max(mean(abs(diag(A))), 1)
  A - diag(ridge * scale, nrow(A))
}

right_solve <- function(A, X, ridge = 0) {
  ridge <- validate_ridge(ridge)

  A <- as.matrix(A)
  X <- as.matrix(X)
  if (nrow(A) != ncol(A)) stop("A must be square.")
  if (ncol(X) != nrow(A)) stop("ncol(X) must equal nrow(A).")
  if (ridge > 0) A <- regularize_jacobian(A, ridge = ridge)
  t(solve(t(A), t(X)))
}

right_solve_transpose <- function(A, X, ridge = 0) {
  ridge <- validate_ridge(ridge)

  A <- as.matrix(A)
  X <- as.matrix(X)
  if (nrow(A) != ncol(A)) stop("A must be square.")
  if (ncol(X) != nrow(A)) stop("ncol(X) must equal nrow(A).")
  if (ridge > 0) A <- regularize_jacobian(A, ridge = ridge)
  t(solve(A, t(X)))
}

regularize_matrix <- function(A, ridge) {
  ridge <- validate_ridge(ridge)

  A <- symmetrize(as.matrix(A))
  if (ridge <= 0) return(A)
  scale <- max(mean(abs(diag(A))), 1)
  A + diag(ridge * scale, nrow(A))
}

retrospective_intercept_index <- function(coef_names) {
  idx <- which(coef_names == "(Intercept)" | grepl(":\\(Intercept\\)$", coef_names))
  if (length(idx) == 0L) integer(0) else idx[1]
}

to_retrospective_from_glm <- function(coef_glm, rho_sample, intercept_index = 1L) {
  if (is.null(rho_sample) || length(rho_sample) != 1L ||
      !is.finite(rho_sample) || rho_sample <= 0) {
    stop("rho_sample must be a positive finite scalar.")
  }
  gamma_ret <- coef_glm
  if (!is.null(intercept_index) && length(intercept_index) > 0L &&
      !is.na(intercept_index[1])) {
    gamma_ret[intercept_index[1]] <- coef_glm[intercept_index[1]] - log(rho_sample)
  }
  gamma_ret
}

to_external_prospective <- function(gamma_ret, rho_E, intercept_index = 1L) {
  if (is.null(rho_E) || length(rho_E) != 1L || !is.finite(rho_E) || rho_E <= 0) {
    stop("rho_E must be a positive finite scalar.")
  }
  gamma_pros <- gamma_ret
  if (!is.null(intercept_index) && length(intercept_index) > 0L &&
      !is.na(intercept_index[1])) {
    gamma_pros[intercept_index[1]] <- gamma_ret[intercept_index[1]] + log(rho_E)
  }
  gamma_pros
}

regularize_jacobian <- function(A, ridge) {
  ridge <- validate_ridge(ridge)

  A <- as.matrix(A)
  if (nrow(A) != ncol(A)) stop("A must be square.")
  if (ridge <= 0) return(A)
  scale <- max(mean(abs(diag(A))), 1)
  A - diag(ridge * scale, nrow(A))
}

zim_null_default <- function(x, y) {
  if (is.null(x)) y else x
}

block_diag_dense <- function(blocks) {
  if (!is.list(blocks) || length(blocks) == 0L) {
    return(matrix(0, 0, 0))
  }
  nr <- sum(vapply(blocks, nrow, integer(1)))
  nc <- sum(vapply(blocks, ncol, integer(1)))
  out <- matrix(0, nr, nc)
  r0 <- 1L
  c0 <- 1L
  for (b in blocks) {
    rr <- r0:(r0 + nrow(b) - 1L)
    cc <- c0:(c0 + ncol(b) - 1L)
    out[rr, cc] <- b
    r0 <- r0 + nrow(b)
    c0 <- c0 + ncol(b)
  }
  out
}

compute_internal_el_weights <- function(Y,
                                        Delta_hat,
                                        rho_I,
                                        mass_tol = ZIM_DEFAULT_EL_MASS_TOL,
                                        mass_action = c("warn", "stop", "none")) {
  mass_action <- match.arg(mass_action)
  mass_tol <- validate_positive_tolerance(mass_tol, "mass_tol")
  Y <- check_binary_y(Y)
  if (length(Delta_hat) != length(Y)) stop("Delta_hat must have one value per internal subject.")
  if (any(!is.finite(Delta_hat)) || any(Delta_hat <= 0)) stop("Delta_hat must be finite and positive.")
  if (is.null(rho_I) || length(rho_I) != 1L || !is.finite(rho_I) || rho_I <= 0) {
    stop("rho_I must be a positive finite scalar.")
  }
  n0 <- sum(Y == 0)
  n1 <- sum(Y == 1)
  if (n0 <= 0 || n1 <= 0) stop("Internal data must contain both cases and controls.")
  p_hat <- 1 / n0 / (1 + rho_I * Delta_hat)
  mass_diagnostics <- validate_el_masses(
    p_hat = p_hat,
    Delta_hat = Delta_hat,
    mass_tol = mass_tol,
    action = mass_action,
    context = "Internal target EL reconstruction"
  )
  list(
    p_hat = as.numeric(p_hat),
    diagnostics = list(
      n_control_internal = n0,
      n_case_internal = n1,
      rho_I = rho_I,
      E0_mass = mass_diagnostics$mass0,
      E1_mass = mass_diagnostics$mass1,
      max_abs_mass_error = mass_diagnostics$max_abs_mass_error,
      mass_tol = mass_diagnostics$mass_tol,
      mass_action = mass_diagnostics$action,
      mass_valid = mass_diagnostics$valid,
      pooled_mass = sum(p_hat * (1 + rho_I * Delta_hat))
    )
  )
}

validate_external_case_control_counts <- function(external_n_case, external_n_control) {
  if (is.null(external_n_case) || is.null(external_n_control)) {
    stop("Provide both external_n_case and external_n_control.")
  }
  if (length(external_n_case) != 1L || length(external_n_control) != 1L ||
      !is.finite(external_n_case) || !is.finite(external_n_control) ||
      external_n_case <= 0 || external_n_control <= 0) {
    stop("external_n_case and external_n_control must be positive finite scalars.")
  }
  list(
    external_n_case = as.numeric(external_n_case),
    external_n_control = as.numeric(external_n_control),
    n_external = as.numeric(external_n_case + external_n_control),
    rho_E = as.numeric(external_n_case / external_n_control)
  )
}

has_complete_names <- function(x) {
  nm <- names(x)
  !is.null(nm) && length(nm) == length(x) &&
    all(!is.na(nm)) && all(nzchar(nm))
}

strip_beta_prefix <- function(beta_names) {
  sub("^beta[0-9]+:", "", beta_names)
}

match_summary_names <- function(input_names, beta_names, what = "beta_tilde") {
  if (length(input_names) != length(beta_names)) {
    stop(what, " has names, but its length does not match the required beta dimension.")
  }
  if (anyDuplicated(input_names)) {
    stop(what, " has duplicated names; cannot safely align summary components.")
  }

  missing_exact <- setdiff(beta_names, input_names)
  if (length(missing_exact) == 0L) {
    return(match(beta_names, input_names))
  }

  stripped <- strip_beta_prefix(beta_names)
  if (!anyDuplicated(stripped)) {
    missing_stripped <- setdiff(stripped, input_names)
    if (length(missing_stripped) == 0L) {
      return(match(stripped, input_names))
    }
  }

  stop(
    what, " names do not match the working-model beta names. Missing examples: ",
    paste(utils::head(missing_exact, 10L), collapse = ", "),
    if (length(missing_exact) > 10L) ", ..." else "",
    ". Supply ", what, " in the same order as working_models, or name it using the returned beta_names."
  )
}

coerce_beta_tilde_vector <- function(beta_tilde) {

  if (is.null(beta_tilde)) stop("beta_tilde must be supplied.")

  if (is.data.frame(beta_tilde)) {
    if (all(c("term", "estimate") %in% names(beta_tilde))) {
      out <- beta_tilde[["estimate"]]
      names(out) <- as.character(beta_tilde[["term"]])
      return(out)
    }
    if (ncol(beta_tilde) == 1L) {
      out <- beta_tilde[[1L]]
      rn <- rownames(beta_tilde)
      if (!is.null(rn) && length(rn) == length(out) && all(nzchar(rn))) names(out) <- rn
      return(out)
    }
    stop("beta_tilde as a data.frame must have one column or columns named 'term' and 'estimate'.")
  }

  if (is.matrix(beta_tilde)) {
    if (ncol(beta_tilde) == 1L) {
      out <- as.vector(beta_tilde[, 1L])
      rn <- rownames(beta_tilde)
      if (!is.null(rn) && length(rn) == length(out) && all(nzchar(rn))) names(out) <- rn
      return(out)
    }
    if (nrow(beta_tilde) == 1L) {
      out <- as.vector(beta_tilde[1L, ])
      cn <- colnames(beta_tilde)
      if (!is.null(cn) && length(cn) == length(out) && all(nzchar(cn))) names(out) <- cn
      return(out)
    }
    stop("beta_tilde as a matrix must be one row or one column.")
  }

  beta_tilde
}

align_beta_tilde_to_working <- function(beta_tilde, beta_names) {

  if (anyDuplicated(beta_names)) {
    stop("Internal working beta names are duplicated; cannot safely align beta_tilde.")
  }

  beta_tilde <- coerce_beta_tilde_vector(beta_tilde)
  if (has_complete_names(beta_tilde)) {
    idx <- match_summary_names(names(beta_tilde), beta_names, what = "beta_tilde")
    beta_tilde <- beta_tilde[idx]
  } else {
    if (length(beta_tilde) != length(beta_names)) {
      stop("length(beta_tilde) must equal total reported beta dimension across working models.")
    }
  }

  beta_tilde <- as.numeric(beta_tilde)
  if (any(!is.finite(beta_tilde))) stop("beta_tilde must contain only finite numeric values.")
  names(beta_tilde) <- beta_names
  beta_tilde
}

align_sigma_betabeta_to_working <- function(Sigma_beta_beta, beta_names) {

  Sigma_beta_beta <- as.matrix(Sigma_beta_beta)
  p_beta <- length(beta_names)
  if (nrow(Sigma_beta_beta) != ncol(Sigma_beta_beta)) {
    stop("Sigma_beta_beta must be a square matrix.")
  }

  rn <- rownames(Sigma_beta_beta)
  cn <- colnames(Sigma_beta_beta)
  has_matrix_names <- !is.null(rn) && !is.null(cn) &&
    length(rn) == nrow(Sigma_beta_beta) && length(cn) == ncol(Sigma_beta_beta) &&
    all(!is.na(rn)) && all(!is.na(cn)) && all(nzchar(rn)) && all(nzchar(cn))

  if (has_matrix_names) {
    row_idx <- match_summary_names(rn, beta_names, what = "rownames(Sigma_beta_beta)")
    col_idx <- match_summary_names(cn, beta_names, what = "colnames(Sigma_beta_beta)")
    Sigma_beta_beta <- Sigma_beta_beta[row_idx, col_idx, drop = FALSE]
  } else {
    if (!all(dim(Sigma_beta_beta) == c(p_beta, p_beta))) {
      stop("Sigma_beta_beta must be a p_beta by p_beta matrix on the sqrt(n_internal) scale.")
    }
  }

  if (any(!is.finite(Sigma_beta_beta))) {
    stop("Sigma_beta_beta must contain only finite numeric values.")
  }
  dimnames(Sigma_beta_beta) <- list(beta_names, beta_names)
  symmetrize(Sigma_beta_beta)
}

require_working_intercept <- function(D_k, k = NA_integer_) {
  intercept_idx <- retrospective_intercept_index(colnames(D_k))
  if (length(intercept_idx) == 0L) {
    stop(
      "Working model ", k,
      " has no intercept in its alpha design. Case-control WKM requires a free intercept ",
      "so the prospective glm intercept can be converted to the retrospective scale."
    )
  }
  intercept_idx
}

working_score_pieces <- function(D_k, gamma_k, rho) {
  if (ncol(D_k) != length(gamma_k)) {
    stop("gamma_k length does not match ncol(D_k).")
  }
  delta <- safe_exp(drop(D_k %*% gamma_k))
  denom <- 1 + rho * delta
  list(
    delta = delta,
    phi0 = D_k * as.numeric(-rho * delta / denom),
    phi1 = D_k * as.numeric(1 / denom)
  )
}

shrink_lambda_diag_target <- function(R, Y, fixed_case_control = TRUE) {
  Z <- if (fixed_case_control) {
    center_by_case_control_stratum(R, Y)
  } else {
    sweep(R, 2, colMeans(R))
  }
  n <- nrow(Z)
  if (n < 3) return(1)
  SS  <- crossprod(Z)
  Zsq <- Z * Z
  Q   <- crossprod(Zsq)
  S   <- SS / (n - 1)
  Var_s <- (n / (n - 1)^3) * (Q - (SS^2) / n)
  off <- row(S) != col(S)
  num <- sum(Var_s[off])
  den <- sum((S[off])^2)
  lambda <- if (den > 0) num / den else 1
  min(1, max(0, lambda))
}

select_pc_components <- function(eig,
                                 eig_tol,
                                 pca_rule,
                                 pc_var,
                                 n_pc,
                                 pc_max_condition,
                                 context = "PCA") {
  if (!is.numeric(eig_tol) || length(eig_tol) != 1L || is.na(eig_tol) ||
      !is.finite(eig_tol) || eig_tol < 0 || eig_tol >= 1) {
    stop("eig_tol must be a finite scalar in [0, 1) for ", context, ".")
  }
  if (!is.list(eig) || is.null(eig$values) || is.null(eig$vectors)) {
    stop("eig must be the result of eigen(..., symmetric = TRUE).")
  }

  values_raw <- as.numeric(eig$values)
  values <- pmax(values_raw, 0)
  eigen_max <- max(values)
  if (!is.finite(eigen_max) || eigen_max <= 0) {
    stop("No positive eigenvalues available for ", context, ".")
  }

  positive <- values > eig_tol * eigen_max
  n_positive <- sum(positive)
  if (n_positive == 0L) {
    stop("No positive eigenvalues available for ", context, ".")
  }

  values_positive <- values[seq_len(n_positive)]
  numerical_cutoff <- eig_tol * values_positive[1L]
  rule_cutoff <- NA_real_

  if (identical(pca_rule, "variance")) {
    cumulative <- cumsum(values_positive) / sum(values_positive)
    n_pc_use <- which(cumulative >= pc_var)[1L]
    rule_cutoff <- values_positive[n_pc_use]
  } else if (identical(pca_rule, "condition")) {
    rule_cutoff <- values_positive[1L] / pc_max_condition
    eligible <- which(values_positive >= rule_cutoff)
    n_pc_use <- max(eligible)
  } else if (identical(pca_rule, "condition_or_variance")) {
    cumulative <- cumsum(values_positive) / sum(values_positive)
    n_pc_variance <- which(cumulative >= pc_var)[1L]
    condition_cutoff <- values_positive[1L] / pc_max_condition
    eligible <- which(values_positive >= condition_cutoff)
    n_pc_condition <- max(eligible)
    n_pc_use <- min(n_pc_variance, n_pc_condition)
    rule_cutoff <- values_positive[n_pc_use]
  } else {
    n_pc_use <- n_pc
    if (n_pc_use > n_positive) {
      stop(
        "n_pc = ", n_pc_use,
        " exceeds the numerical positive rank ", n_positive,
        " determined by eig_tol = ", signif(eig_tol, 6),
        " for ", context, "."
      )
    }
    rule_cutoff <- values_positive[n_pc_use]
  }

  keep <- seq_len(n_pc_use)
  U <- eig$vectors[, keep, drop = FALSE]
  lambda <- values_positive[keep]
  first_discarded <- if (n_pc_use < n_positive) {
    values_positive[n_pc_use + 1L]
  } else {
    NA_real_
  }
  raw_scale <- max(max(abs(values_raw)), .Machine$double.eps)
  material_negative <- values_raw < -max(eig_tol, 100 * .Machine$double.eps) * raw_scale

  list(
    U = U,
    lambda = lambda,
    values_raw = values_raw,
    values_positive = values_positive,
    n_pc_used = n_pc_use,
    n_positive = n_positive,
    numerical_cutoff = numerical_cutoff,
    rule_cutoff = rule_cutoff,
    effective_cutoff = max(numerical_cutoff, rule_cutoff),
    pc_var_used = sum(lambda) / sum(values_positive),
    full_condition = values_positive[1L] / values_positive[n_positive],
    retained_condition = lambda[1L] / lambda[length(lambda)],
    first_discarded = first_discarded,
    raw_min_eigenvalue = min(values_raw),
    raw_max_eigenvalue = max(values_raw),
    n_material_negative = sum(material_negative)
  )
}

build_pc_reduction_basis <- function(A_raw,
                                     eig_tol,
                                     pca_rule,
                                     pc_var,
                                     n_pc,
                                     pc_max_condition,
                                     pca_scale = c("covariance", "correlation"),
                                     context = "PC reduction basis") {
  pca_scale <- match.arg(pca_scale)
  A_raw <- symmetrize(as.matrix(A_raw))
  if (nrow(A_raw) != ncol(A_raw)) stop(context, " matrix must be square.")
  if (any(!is.finite(A_raw))) stop(context, " matrix must be finite.")
  p <- nrow(A_raw)

  if (identical(pca_scale, "correlation")) {
    marginal_variance <- diag(A_raw)
    bad_variance <- !is.finite(marginal_variance) | marginal_variance <= 0
    if (any(bad_variance)) {
      bad_names <- names(marginal_variance)[bad_variance]
      if (is.null(bad_names) || length(bad_names) != sum(bad_variance) ||
          any(is.na(bad_names)) || any(!nzchar(bad_names))) {
        bad_names <- which(bad_variance)
      }
      stop(
        context, " with pca_scale = 'correlation' requires every diagonal ",
        "entry of its spectral basis to be positive and finite. Problematic ",
        "beta components: ", paste(utils::head(bad_names, 10L), collapse = ", "),
        if (sum(bad_variance) > 10L) ", ..." else "", "."
      )
    }
    scale_sd <- sqrt(marginal_variance)
    A_spectral <- sweep(
      sweep(A_raw, 1L, scale_sd, "/"),
      2L, scale_sd, "/"
    )
  } else {
    scale_sd <- rep(1, p)
    names(scale_sd) <- colnames(A_raw)
    A_spectral <- A_raw
  }
  A_spectral <- symmetrize(A_spectral)

  eig <- eigen(A_spectral, symmetric = TRUE)
  selected <- select_pc_components(
    eig = eig,
    eig_tol = eig_tol,
    pca_rule = pca_rule,
    pc_var = pc_var,
    n_pc = n_pc,
    pc_max_condition = pc_max_condition,
    context = context
  )

  if (selected$n_material_negative > 0L) {
    warning(
      context, " has ", selected$n_material_negative,
      " materially negative eigenvalue(s). Only the positive spectral ",
      "directions are eligible for PC dimension reduction."
    )
  }

  B <- sweep(selected$U, 1L, scale_sd, "/")
  rownames(B) <- colnames(A_raw)
  colnames(B) <- paste0("PC", seq_len(ncol(B)))

  list(
    reduction_matrix = B,
    U_spectral = selected$U,
    scale_sd = scale_sd,
    selected = selected,
    spectral_matrix = A_spectral,
    pca_scale = pca_scale
  )
}

build_factor_diag_covariance <- function(A_raw,
                                         eig_tol,
                                         pca_rule,
                                         pc_var,
                                         n_pc,
                                         pc_max_condition,
                                         pca_scale = c("covariance", "correlation"),
                                         context = "factor-plus-diagonal covariance") {
  pca_scale <- match.arg(pca_scale)
  A_raw <- symmetrize(as.matrix(A_raw))
  if (nrow(A_raw) != ncol(A_raw)) stop(context, " matrix must be square.")
  if (any(!is.finite(A_raw))) stop(context, " matrix must be finite.")
  p <- nrow(A_raw)

  if (identical(pca_scale, "correlation")) {
    marginal_variance <- diag(A_raw)
    bad_variance <- !is.finite(marginal_variance) | marginal_variance <= 0
    if (any(bad_variance)) {
      bad_names <- names(marginal_variance)[bad_variance]
      if (is.null(bad_names) || length(bad_names) != sum(bad_variance) ||
          any(is.na(bad_names)) || any(!nzchar(bad_names))) {
        bad_names <- which(bad_variance)
      }
      stop(
        context, " with pca_scale = 'correlation' requires every diagonal ",
        "entry of its spectral basis to be positive and finite. Problematic ",
        "beta components: ", paste(utils::head(bad_names, 10L), collapse = ", "),
        if (sum(bad_variance) > 10L) ", ..." else "", "."
      )
    }
    scale_sd <- sqrt(marginal_variance)
    A_spectral <- sweep(
      sweep(A_raw, 1L, scale_sd, "/"),
      2L, scale_sd, "/"
    )
  } else {
    scale_sd <- rep(1, p)
    names(scale_sd) <- colnames(A_raw)
    A_spectral <- A_raw
  }
  A_spectral <- symmetrize(A_spectral)

  eig <- eigen(A_spectral, symmetric = TRUE)
  selected <- select_pc_components(
    eig = eig,
    eig_tol = eig_tol,
    pca_rule = pca_rule,
    pc_var = pc_var,
    n_pc = n_pc,
    pc_max_condition = pc_max_condition,
    context = context
  )

  if (selected$n_material_negative > 0L) {
    warning(
      context, " has ", selected$n_material_negative,
      " materially negative eigenvalue(s). The factor-plus-diagonal ",
      "reconstruction uses the positive spectral part and restores the ",
      "original marginal variances on the diagonal."
    )
  }

  U <- selected$U
  lambda <- selected$lambda
  U_scaled <- sweep(U, 2L, sqrt(lambda), "*")
  low_rank_spectral <- tcrossprod(U_scaled)

  target_diag_spectral <- diag(A_spectral)
  residual_diag_spectral <- target_diag_spectral - diag(low_rank_spectral)
  residual_scale <- max(1, max(abs(target_diag_spectral)))
  residual_tol <- max(10 * eig_tol, 100 * .Machine$double.eps) * residual_scale
  if (any(residual_diag_spectral < -residual_tol)) {
    bad <- which(residual_diag_spectral < -residual_tol)
    stop(
      context, " cannot preserve the original diagonal because the retained ",
      "factor part exceeds the marginal variance for component(s): ",
      paste(utils::head(bad, 10L), collapse = ", "),
      if (length(bad) > 10L) ", ..." else "",
      ". The spectral basis may be materially indefinite or the requested ",
      "rank may be incompatible with a factor-plus-diagonal reconstruction."
    )
  }
  residual_diag_spectral[residual_diag_spectral < 0] <- 0

  factor_spectral <- low_rank_spectral +
    diag(residual_diag_spectral, nrow = p, ncol = p)
  factor_original <- sweep(
    sweep(factor_spectral, 1L, scale_sd, "*"),
    2L, scale_sd, "*"
  )
  factor_original <- symmetrize(factor_original)
  dimnames(factor_original) <- dimnames(A_raw)

  residual_diag_original <- residual_diag_spectral * scale_sd^2
  diagonal_error <- max(abs(diag(factor_original) - diag(A_raw)))
  diagonal_tolerance <- max(1, max(abs(diag(A_raw)))) *
    max(10 * eig_tol, 100 * .Machine$double.eps)
  factor_eigenvalues <- eigen(
    factor_original,
    symmetric = TRUE,
    only.values = TRUE
  )$values

  difference <- A_raw - factor_original
  diag(difference) <- 0

  list(
    matrix = factor_original,
    projection_spectral = tcrossprod(U),
    scale_sd = scale_sd,
    selected = selected,
    residual_diag_original = residual_diag_original,
    diagonal_error = diagonal_error,
    diagonal_tolerance = diagonal_tolerance,
    diagonal_preserved = diagonal_error <= diagonal_tolerance,
    offdiag_removed_frobenius = sqrt(sum(difference * difference)),
    factor_min_eigenvalue = min(factor_eigenvalues),
    factor_max_eigenvalue = max(factor_eigenvalues),
    pca_scale = pca_scale
  )
}

project_crosscov_to_factor_space <- function(H_theta_beta, factor_object) {
  H_theta_beta <- as.matrix(H_theta_beta)
  scale_sd <- factor_object$scale_sd
  if (ncol(H_theta_beta) != length(scale_sd)) {
    stop("H_theta_beta has the wrong number of columns for factor projection.")
  }
  H_standardized <- sweep(H_theta_beta, 2L, scale_sd, "/")
  H_projected_standardized <- H_standardized %*%
    factor_object$projection_spectral
  H_projected <- sweep(H_projected_standardized, 2L, scale_sd, "*")
  dimnames(H_projected) <- dimnames(H_theta_beta)
  H_projected
}

compute_zim_weight <- function(H_theta_beta, H_beta_beta, Sigma_beta_beta,
                               ridge,
                               inverse_method = c(
                                 "solve", "pca", "hbb_pca", "sigma_pca",
                                 "m_factor_diag", "hbb_factor_diag", "shrink", "shrinkage"
                               ),
                               pc_var = 0.99,
                               n_pc = NULL,
                               eig_tol = 1e-10,
                               pca_rule = c("variance", "condition", "condition_or_variance", "fixed"),
                               pc_max_condition = 1e4,
                               pca_scale = c("covariance", "correlation"),

                               H_theta_theta = NULL,
                               IF_theta = NULL,
                               IF_beta = NULL,
                               Y = NULL,
                               fixed_case_control = TRUE,
                               shrink_structured = TRUE,
                               shrink_lambda = NULL) {
  ridge <- validate_ridge(ridge)
  pca_rule_was_missing <- missing(pca_rule)
  inverse_method <- match.arg(inverse_method)
  shrink_lambda <- validate_shrinkage_controls(
    inverse_method, shrink_lambda, shrink_structured = shrink_structured
  )
  pca_family_active <- is_pca_family_method(inverse_method)
  pca_controls <- validate_pca_controls(
    pca_rule = pca_rule,
    pc_var = pc_var,
    n_pc = n_pc,
    pc_max_condition = pc_max_condition,
    pca_scale = pca_scale,
    pca_rule_was_missing = pca_rule_was_missing,
    active = pca_family_active
  )
  pca_rule <- pca_controls$pca_rule
  pc_var <- pca_controls$pc_var
  n_pc <- pca_controls$n_pc
  pc_max_condition <- pca_controls$pc_max_condition
  pca_scale <- pca_controls$pca_scale

  H_theta_beta <- as.matrix(H_theta_beta)
  H_beta_beta <- symmetrize(as.matrix(H_beta_beta))
  Sigma_beta_beta <- symmetrize(as.matrix(Sigma_beta_beta))
  if (nrow(H_beta_beta) != ncol(H_beta_beta) ||
      nrow(Sigma_beta_beta) != ncol(Sigma_beta_beta) ||
      !all(dim(H_beta_beta) == dim(Sigma_beta_beta))) {
    stop("H_beta_beta and Sigma_beta_beta must be square matrices with matching dimensions.")
  }
  if (ncol(H_theta_beta) != nrow(H_beta_beta)) {
    stop("ncol(H_theta_beta) must equal nrow(H_beta_beta).")
  }

  M_beta_raw <- symmetrize(H_beta_beta + Sigma_beta_beta)

  M_beta <- regularize_matrix(M_beta_raw, ridge = ridge)
  p_beta <- ncol(M_beta)

  empty_diag <- function(method) list(
    inverse_method = method,
    ridge = ridge,
    p_beta = p_beta,
    pca_rule = NA_character_,
    pca_scale = NA_character_,
    pca_rule_inferred_from_n_pc = FALSE,
    n_pc_requested = NA_integer_,
    n_pc_used = NA_integer_,
    pc_var_requested = NA_real_,
    pc_var_used = NA_real_,
    pc_max_condition_requested = NA_real_,
    eig_tol_requested = NA_real_,
    n_positive_eigenvalues = NA_integer_,
    numerical_eigen_cutoff = NA_real_,
    rule_eigen_cutoff = NA_real_,
    effective_eigen_cutoff = NA_real_,
    retained_eigen_min = NA_real_,
    first_discarded_eigenvalue = NA_real_,
    eigen_min = NA_real_,
    eigen_max = NA_real_,
    condition_number = NA_real_,
    retained_condition_number = NA_real_,
    covariance_condition_number = NA_real_,
    beta_scale_sd_min = NA_real_,
    beta_scale_sd_max = NA_real_,
    pca_basis = NA_character_,
    pca_basis_condition_number = NA_real_,
    pca_basis_raw_min_eigenvalue = NA_real_,
    pca_basis_material_negative_eigenvalues = NA_integer_,
    pca_reduced_dimension = NA_integer_,
    pca_reduction_matrix_rank = NA_integer_,
    pca_reduced_total_min_eigenvalue = NA_real_,
    pca_reduced_total_condition_number = NA_real_,
    factor_diag_basis = NA_character_,
    factor_diag_diagonal_target = NA_character_,
    factor_diag_diagonal_preserved = NA,
    factor_diag_max_abs_diagonal_error = NA_real_,
    factor_diag_total_M_max_abs_diagonal_error = NA_real_,
    factor_diag_residual_variance_min = NA_real_,
    factor_diag_residual_variance_max = NA_real_,
    factor_diag_offdiag_removed_frobenius = NA_real_,
    factor_diag_basis_raw_min_eigenvalue = NA_real_,
    factor_diag_basis_material_negative_eigenvalues = NA_integer_,
    factor_diag_matrix_min_eigenvalue = NA_real_,
    factor_diag_matrix_condition_number = NA_real_,
    factor_diag_H_theta_beta_projected = NA,
    factor_diag_H_theta_beta_projection_max_abs_change = NA_real_,
    shrink_lambda = NA_real_,
    shrink_structured = NA
  )

  if (inverse_method == "solve") {
    W_hat <- right_solve(M_beta, H_theta_beta)
    kappa_est <- tryCatch(kappa(M_beta, exact = FALSE), error = function(e) NA_real_)
    diagnostics <- empty_diag("solve")
    diagnostics$n_pc_used <- p_beta
    diagnostics$pc_var_used <- 1
    diagnostics$condition_number <- kappa_est
    diagnostics$retained_condition_number <- kappa_est
    diagnostics$covariance_condition_number <- kappa_est
    return(list(
      W_hat = W_hat,
      M_beta = M_beta,
      M_beta_weight_model = M_beta,
      M_beta_weight_model_raw = M_beta_raw,
      H_beta_beta_weight_model = H_beta_beta,
      H_theta_beta_weight_model = H_theta_beta,
      diagnostics = diagnostics
    ))
  }

  if (inverse_method == "shrinkage") {
    if (is.null(H_theta_theta)) {
      stop("inverse_method = 'shrinkage' requires H_theta_theta.")
    }
    Htt <- as.matrix(H_theta_theta)
    if (!is.numeric(Htt) ||
        !identical(dim(Htt), c(nrow(H_theta_beta), nrow(H_theta_beta))) ||
        any(!is.finite(Htt))) {
      stop("inverse_method = 'shrinkage' requires finite H_theta_theta matching the rows of H_theta_beta.")
    }
    Htt <- symmetrize(Htt)
    tryCatch(chol(Htt), error = function(e) {
      stop("inverse_method = 'shrinkage' requires positive-definite H_theta_theta.")
    })
    validate_joint_covariance_blocks(
      Htt, H_theta_beta, H_beta_beta,
      context = "Target-preserving shrinkage joint covariance"
    )
    if (any(!is.finite(Sigma_beta_beta))) {
      stop("inverse_method = 'shrinkage' requires finite Sigma_beta_beta.")
    }

    K_target <- symmetrize(t(H_theta_beta) %*% solve(Htt, H_theta_beta))
    R_residual <- symmetrize(H_beta_beta - K_target)
    residual_psd <- matrix_psd_diagnostics(R_residual)
    if (!isTRUE(residual_psd$is_psd)) {
      stop("Target-preserving shrinkage residual covariance is not positive semidefinite within tolerance.")
    }
    R_weight <- (1 - shrink_lambda) * R_residual
    diag(R_weight) <- diag(R_residual)
    Hbb_weight <- symmetrize(K_target + R_weight)

    if (shrink_lambda == 0) Hbb_weight <- H_beta_beta
    diag(Hbb_weight) <- diag(H_beta_beta)
    dimnames(Hbb_weight) <- dimnames(H_beta_beta)
    M_weight_raw <- symmetrize(Hbb_weight + Sigma_beta_beta)
    M_weight <- regularize_matrix(M_weight_raw, ridge = ridge)
    W_hat <- right_solve(M_weight, H_theta_beta)

    diagnostics <- empty_diag("shrinkage")
    diagnostics$n_pc_used <- p_beta
    diagnostics$pc_var_used <- 1
    diagnostics$shrink_lambda <- shrink_lambda
    diagnostics$shrink_structured <- TRUE
    diagnostics$shrink_target <- "diag(Hbb - Hbt solve(Htt, Htb))"
    diagnostics$shrink_residual_min_eigenvalue <- residual_psd$minimum_eigenvalue
    diagnostics$shrink_diagonal_max_abs_error <- max(abs(diag(Hbb_weight) - diag(H_beta_beta)))
    diagnostics$shrink_offdiag_removed_frobenius <- norm(Hbb_weight - H_beta_beta, "F")
    diagnostics$condition_number <- tryCatch(
      kappa(M_weight, exact = FALSE), error = function(e) NA_real_
    )
    diagnostics$retained_condition_number <- diagnostics$condition_number
    diagnostics$covariance_condition_number <- tryCatch(
      kappa(M_beta, exact = FALSE), error = function(e) NA_real_
    )
    return(list(
      W_hat = W_hat,
      M_beta = M_beta,
      M_beta_weight_model = M_weight,
      M_beta_weight_model_raw = M_weight_raw,
      H_beta_beta_weight_model = Hbb_weight,
      H_theta_beta_weight_model = H_theta_beta,
      diagnostics = diagnostics
    ))
  }

  if (inverse_method == "pca") {
    if (!is.numeric(eig_tol) || length(eig_tol) != 1L || is.na(eig_tol) ||
        !is.finite(eig_tol) || eig_tol < 0 || eig_tol >= 1) {
      stop("eig_tol must be a finite scalar in [0, 1) for the PCA inverse.")
    }

    if (identical(pca_scale, "correlation")) {
      beta_variance <- diag(M_beta_raw)
      bad_variance <- !is.finite(beta_variance) | beta_variance <= 0
      if (any(bad_variance)) {
        bad_names <- names(beta_variance)[bad_variance]
        if (is.null(bad_names) || length(bad_names) != sum(bad_variance) ||
            any(is.na(bad_names)) || any(!nzchar(bad_names))) {
          bad_names <- which(bad_variance)
        }
        stop(
          "pca_scale = 'correlation' requires every diagonal entry of ",
          "H_beta_beta + Sigma_beta_beta to be positive and finite. ",
          "Problematic beta components: ",
          paste(utils::head(bad_names, 10L), collapse = ", "),
          if (sum(bad_variance) > 10L) ", ..." else "",
          ". Harmonize/remove degenerate summaries or use pca_scale = 'covariance'."
        )
      }
      beta_scale_sd <- sqrt(beta_variance)
      M_pca_raw <- sweep(
        sweep(M_beta_raw, 1L, beta_scale_sd, "/"),
        2L, beta_scale_sd, "/"
      )
      M_pca <- regularize_matrix(M_pca_raw, ridge = ridge)
      H_theta_beta_pca <- sweep(H_theta_beta, 2L, beta_scale_sd, "/")
    } else {
      beta_scale_sd <- rep(1, p_beta)
      names(beta_scale_sd) <- colnames(M_beta)
      M_pca <- M_beta
      H_theta_beta_pca <- H_theta_beta
    }

    eig <- eigen(M_pca, symmetric = TRUE)
    selected <- select_pc_components(
      eig = eig,
      eig_tol = eig_tol,
      pca_rule = pca_rule,
      pc_var = pc_var,
      n_pc = n_pc,
      pc_max_condition = pc_max_condition,
      context = "the PCA inverse"
    )
    U <- selected$U
    lambda <- selected$lambda

    W_pca <- sweep(H_theta_beta_pca %*% U, 2, lambda, "/") %*% t(U)

    W_hat <- sweep(W_pca, 2L, beta_scale_sd, "/")

    diagnostics <- empty_diag("pca")
    diagnostics$pca_rule <- pca_rule
    diagnostics$pca_scale <- pca_scale
    diagnostics$pca_rule_inferred_from_n_pc <- pca_controls$inferred_from_n_pc
    diagnostics$n_pc_requested <- if (identical(pca_rule, "fixed")) n_pc else NA_integer_
    diagnostics$n_pc_used <- selected$n_pc_used
    diagnostics$pc_var_requested <- if (pca_rule %in% c("variance", "condition_or_variance")) {
      pc_var
    } else {
      NA_real_
    }
    diagnostics$pc_var_used <- selected$pc_var_used
    diagnostics$pc_max_condition_requested <- if (pca_rule %in% c("condition", "condition_or_variance")) {
      pc_max_condition
    } else {
      NA_real_
    }
    diagnostics$eig_tol_requested <- eig_tol
    diagnostics$n_positive_eigenvalues <- selected$n_positive
    diagnostics$numerical_eigen_cutoff <- selected$numerical_cutoff
    diagnostics$rule_eigen_cutoff <- selected$rule_cutoff
    diagnostics$effective_eigen_cutoff <- selected$effective_cutoff
    diagnostics$retained_eigen_min <- min(lambda)
    diagnostics$first_discarded_eigenvalue <- selected$first_discarded
    diagnostics$eigen_min <- min(selected$values_positive)
    diagnostics$eigen_max <- max(selected$values_positive)
    diagnostics$condition_number <- selected$full_condition
    diagnostics$retained_condition_number <- selected$retained_condition
    diagnostics$covariance_condition_number <- tryCatch(
      kappa(M_beta, exact = FALSE),
      error = function(e) NA_real_
    )
    diagnostics$beta_scale_sd_min <- min(beta_scale_sd)
    diagnostics$beta_scale_sd_max <- max(beta_scale_sd)
    diagnostics$pca_basis <- "M_beta = H_beta_beta + Sigma_beta_beta"
    diagnostics$pca_basis_condition_number <- selected$full_condition
    diagnostics$pca_basis_raw_min_eigenvalue <- selected$raw_min_eigenvalue
    diagnostics$pca_basis_material_negative_eigenvalues <- selected$n_material_negative
    diagnostics$pca_reduced_dimension <- selected$n_pc_used
    diagnostics$pca_reduction_matrix_rank <- selected$n_pc_used
    diagnostics$pca_reduced_total_min_eigenvalue <- min(lambda)
    diagnostics$pca_reduced_total_condition_number <- selected$retained_condition
    return(list(
      W_hat = W_hat,
      M_beta = M_beta,
      M_beta_weight_model = NULL,
      M_beta_weight_model_raw = M_beta_raw,
      H_beta_beta_weight_model = NULL,
      H_theta_beta_weight_model = H_theta_beta,
      diagnostics = diagnostics
    ))
  }

  if (inverse_method %in% c("hbb_pca", "sigma_pca")) {

    basis_name <- if (identical(inverse_method, "hbb_pca")) {
      "H_beta_beta"
    } else {
      "Sigma_beta_beta"
    }
    basis_matrix <- if (identical(inverse_method, "hbb_pca")) {
      H_beta_beta
    } else {
      Sigma_beta_beta
    }
    basis_label <- if (identical(inverse_method, "hbb_pca")) {
      "internal H_beta_beta"
    } else {
      "external Sigma_beta_beta"
    }

    basis <- build_pc_reduction_basis(
      A_raw = basis_matrix,
      eig_tol = eig_tol,
      pca_rule = pca_rule,
      pc_var = pc_var,
      n_pc = n_pc,
      pc_max_condition = pc_max_condition,
      pca_scale = pca_scale,
      context = paste0(inverse_method, " ", basis_name, " basis")
    )
    B <- basis$reduction_matrix
    selected <- basis$selected

    M_reduced_raw <- symmetrize(crossprod(B, M_beta_raw %*% B))
    M_reduced <- symmetrize(crossprod(B, M_beta %*% B))
    H_theta_beta_reduced <- H_theta_beta %*% B

    W_reduced <- tryCatch(
      right_solve(M_reduced, H_theta_beta_reduced),
      error = function(e) {
        stop(
          "The ", inverse_method,
          " projected total covariance could not be inverted: ",
          conditionMessage(e),
          ". Increase ridge, retain fewer components, or inspect the selected ",
          basis_label, " subspace."
        )
      }
    )
    W_hat <- W_reduced %*% t(B)

    reduced_values <- eigen(
      M_reduced,
      symmetric = TRUE,
      only.values = TRUE
    )$values
    reduced_min <- min(reduced_values)
    reduced_max <- max(reduced_values)
    reduced_condition <- if (is.finite(reduced_min) && reduced_min > 0) {
      reduced_max / reduced_min
    } else {
      Inf
    }
    original_condition <- tryCatch(
      kappa(M_beta, exact = FALSE),
      error = function(e) NA_real_
    )

    diagnostics <- empty_diag(inverse_method)
    diagnostics$pca_rule <- pca_rule
    diagnostics$pca_scale <- pca_scale
    diagnostics$pca_rule_inferred_from_n_pc <- pca_controls$inferred_from_n_pc
    diagnostics$n_pc_requested <- if (identical(pca_rule, "fixed")) n_pc else NA_integer_
    diagnostics$n_pc_used <- selected$n_pc_used
    diagnostics$pc_var_requested <- if (pca_rule %in% c("variance", "condition_or_variance")) {
      pc_var
    } else {
      NA_real_
    }
    diagnostics$pc_var_used <- selected$pc_var_used
    diagnostics$pc_max_condition_requested <- if (pca_rule %in% c("condition", "condition_or_variance")) {
      pc_max_condition
    } else {
      NA_real_
    }
    diagnostics$eig_tol_requested <- eig_tol
    diagnostics$n_positive_eigenvalues <- selected$n_positive
    diagnostics$numerical_eigen_cutoff <- selected$numerical_cutoff
    diagnostics$rule_eigen_cutoff <- selected$rule_cutoff
    diagnostics$effective_eigen_cutoff <- selected$effective_cutoff
    diagnostics$retained_eigen_min <- min(selected$lambda)
    diagnostics$first_discarded_eigenvalue <- selected$first_discarded
    diagnostics$eigen_min <- min(selected$values_positive)
    diagnostics$eigen_max <- max(selected$values_positive)

    diagnostics$condition_number <- reduced_condition
    diagnostics$retained_condition_number <- selected$retained_condition
    diagnostics$covariance_condition_number <- original_condition
    diagnostics$beta_scale_sd_min <- min(basis$scale_sd)
    diagnostics$beta_scale_sd_max <- max(basis$scale_sd)
    diagnostics$pca_basis <- basis_name
    diagnostics$pca_basis_condition_number <- selected$full_condition
    diagnostics$pca_basis_raw_min_eigenvalue <- selected$raw_min_eigenvalue
    diagnostics$pca_basis_material_negative_eigenvalues <- selected$n_material_negative
    diagnostics$pca_reduced_dimension <- selected$n_pc_used
    diagnostics$pca_reduction_matrix_rank <- qr(B)$rank
    diagnostics$pca_reduced_total_min_eigenvalue <- reduced_min
    diagnostics$pca_reduced_total_condition_number <- reduced_condition

    return(list(
      W_hat = W_hat,
      M_beta = M_beta,
      M_beta_weight_model = NULL,
      M_beta_weight_model_raw = NULL,
      H_beta_beta_weight_model = NULL,
      H_theta_beta_weight_model = H_theta_beta,
      M_beta_weight_reduced = M_reduced,
      M_beta_weight_reduced_raw = M_reduced_raw,
      beta_reduction_matrix = B,
      H_theta_beta_weight_reduced = H_theta_beta_reduced,
      diagnostics = diagnostics
    ))
  }

  if (inverse_method %in% c("m_factor_diag", "hbb_factor_diag")) {
    factor_basis_name <- if (identical(inverse_method, "m_factor_diag")) {
      "M_beta = H_beta_beta + Sigma_beta_beta"
    } else {
      "H_beta_beta"
    }
    factor_basis <- if (identical(inverse_method, "m_factor_diag")) {
      M_beta_raw
    } else {
      H_beta_beta
    }

    factor <- build_factor_diag_covariance(
      A_raw = factor_basis,
      eig_tol = eig_tol,
      pca_rule = pca_rule,
      pc_var = pc_var,
      n_pc = n_pc,
      pc_max_condition = pc_max_condition,
      pca_scale = pca_scale,
      context = paste0(inverse_method, " spectral basis")
    )

    if (identical(inverse_method, "m_factor_diag")) {
      H_beta_beta_model <- NULL
      H_theta_beta_model <- H_theta_beta
      M_weight_raw <- factor$matrix
      numerator_projected <- FALSE
      numerator_change <- 0
      diagonal_target <- "diag(M_beta)"
    } else {
      H_beta_beta_model <- factor$matrix
      H_theta_beta_model <- project_crosscov_to_factor_space(
        H_theta_beta,
        factor
      )
      M_weight_raw <- symmetrize(H_beta_beta_model + Sigma_beta_beta)
      numerator_projected <- TRUE
      numerator_change <- max(abs(H_theta_beta_model - H_theta_beta))
      diagonal_target <- "diag(H_beta_beta), hence diag(M_beta)"
    }

    M_weight <- regularize_matrix(M_weight_raw, ridge = ridge)
    W_hat <- right_solve(M_weight, H_theta_beta_model)
    actual_condition <- tryCatch(
      kappa(M_weight, exact = FALSE),
      error = function(e) NA_real_
    )
    original_condition <- tryCatch(
      kappa(M_beta, exact = FALSE),
      error = function(e) NA_real_
    )
    total_diagonal_error <- max(abs(diag(M_weight_raw) - diag(M_beta_raw)))
    selected <- factor$selected

    diagnostics <- empty_diag(inverse_method)
    diagnostics$pca_rule <- pca_rule
    diagnostics$pca_scale <- pca_scale
    diagnostics$pca_rule_inferred_from_n_pc <- pca_controls$inferred_from_n_pc
    diagnostics$n_pc_requested <- if (identical(pca_rule, "fixed")) n_pc else NA_integer_
    diagnostics$n_pc_used <- selected$n_pc_used
    diagnostics$pc_var_requested <- if (pca_rule %in% c("variance", "condition_or_variance")) {
      pc_var
    } else {
      NA_real_
    }
    diagnostics$pc_var_used <- selected$pc_var_used
    diagnostics$pc_max_condition_requested <- if (pca_rule %in% c("condition", "condition_or_variance")) {
      pc_max_condition
    } else {
      NA_real_
    }
    diagnostics$eig_tol_requested <- eig_tol
    diagnostics$n_positive_eigenvalues <- selected$n_positive
    diagnostics$numerical_eigen_cutoff <- selected$numerical_cutoff
    diagnostics$rule_eigen_cutoff <- selected$rule_cutoff
    diagnostics$effective_eigen_cutoff <- selected$effective_cutoff
    diagnostics$retained_eigen_min <- min(selected$lambda)
    diagnostics$first_discarded_eigenvalue <- selected$first_discarded
    diagnostics$eigen_min <- min(selected$values_positive)
    diagnostics$eigen_max <- max(selected$values_positive)
    diagnostics$condition_number <- actual_condition
    diagnostics$retained_condition_number <- selected$retained_condition
    diagnostics$covariance_condition_number <- original_condition
    diagnostics$beta_scale_sd_min <- min(factor$scale_sd)
    diagnostics$beta_scale_sd_max <- max(factor$scale_sd)
    diagnostics$pca_basis <- factor_basis_name
    diagnostics$pca_basis_condition_number <- selected$full_condition
    diagnostics$pca_basis_raw_min_eigenvalue <- selected$raw_min_eigenvalue
    diagnostics$pca_basis_material_negative_eigenvalues <- selected$n_material_negative
    diagnostics$pca_reduced_dimension <- selected$n_pc_used
    diagnostics$pca_reduction_matrix_rank <- selected$n_pc_used
    diagnostics$factor_diag_basis <- factor_basis_name
    diagnostics$factor_diag_diagonal_target <- diagonal_target
    diagnostics$factor_diag_diagonal_preserved <- isTRUE(factor$diagonal_preserved) &&
      total_diagonal_error <= factor$diagonal_tolerance
    diagnostics$factor_diag_max_abs_diagonal_error <- factor$diagonal_error
    diagnostics$factor_diag_total_M_max_abs_diagonal_error <- total_diagonal_error
    diagnostics$factor_diag_residual_variance_min <- min(factor$residual_diag_original)
    diagnostics$factor_diag_residual_variance_max <- max(factor$residual_diag_original)
    diagnostics$factor_diag_offdiag_removed_frobenius <- factor$offdiag_removed_frobenius
    diagnostics$factor_diag_basis_raw_min_eigenvalue <- selected$raw_min_eigenvalue
    diagnostics$factor_diag_basis_material_negative_eigenvalues <- selected$n_material_negative
    diagnostics$factor_diag_matrix_min_eigenvalue <- tryCatch(
      min(eigen(M_weight, symmetric = TRUE, only.values = TRUE)$values),
      error = function(e) NA_real_
    )
    diagnostics$factor_diag_matrix_condition_number <- actual_condition
    diagnostics$factor_diag_H_theta_beta_projected <- numerator_projected
    diagnostics$factor_diag_H_theta_beta_projection_max_abs_change <- numerator_change

    return(list(
      W_hat = W_hat,
      M_beta = M_beta,
      M_beta_weight_model = M_weight,
      M_beta_weight_model_raw = M_weight_raw,
      H_beta_beta_weight_model = H_beta_beta_model,
      H_theta_beta_weight_model = H_theta_beta_model,
      diagnostics = diagnostics
    ))
  }

  if (is.null(IF_beta) || is.null(IF_theta) || is.null(H_theta_theta) || is.null(Y)) {
    stop("inverse_method = 'shrink' requires IF_beta, IF_theta, H_theta_theta, and Y.")
  }
  if (ncol(IF_beta) != p_beta) {
    stop("IF_beta has the wrong number of columns for shrink.")
  }
  d <- ncol(IF_theta)

  A_t <- solve(H_theta_theta, H_theta_beta)
  A_hat <- t(A_t)

  if (shrink_structured) {

    R <- IF_beta - IF_theta %*% A_t
  } else {
    R <- IF_beta
  }

  lambda <- shrink_lambda_diag_target(R, Y, fixed_case_control)
  U_hat <- symmetrize(case_control_cov(R, Y, fixed_case_control))
  U_shrink <- (1 - lambda) * U_hat
  diag(U_shrink) <- diag(U_hat)

  if (shrink_structured) {
    M0 <- regularize_matrix(U_shrink + Sigma_beta_beta, ridge = ridge)
    M0inv <- solve(M0)
    mid <- solve(solve(H_theta_theta) + t(A_hat) %*% M0inv %*% A_hat)
    Minv <- M0inv - M0inv %*% A_hat %*% mid %*% t(A_hat) %*% M0inv
    M_inverted_for_diag <- M0
  } else {
    M_shrink <- regularize_matrix(U_shrink + Sigma_beta_beta, ridge = ridge)
    Minv <- solve(M_shrink)
    M_inverted_for_diag <- M_shrink
  }

  W_hat <- H_theta_beta %*% Minv

  kappa_est <- tryCatch(
    kappa(M_inverted_for_diag, exact = FALSE),
    error = function(e) NA_real_
  )
  diagnostics <- empty_diag("shrink")
  diagnostics$n_pc_used <- p_beta
  diagnostics$condition_number <- kappa_est
  diagnostics$retained_condition_number <- kappa_est
  diagnostics$covariance_condition_number <- tryCatch(
    kappa(M_beta, exact = FALSE),
    error = function(e) NA_real_
  )
  diagnostics$shrink_lambda <- lambda
  diagnostics$shrink_structured <- shrink_structured

  list(
    W_hat = W_hat,
    M_beta = M_beta,
    M_beta_weight_model = M_inverted_for_diag,
    M_beta_weight_model_raw = NULL,
    H_beta_beta_weight_model = NULL,
    H_theta_beta_weight_model = H_theta_beta,
    diagnostics = diagnostics
  )
}

center_by_case_control_stratum <- function(IF, Y) {

  out <- IF
  for (y in c(0, 1)) {
    idx <- which(Y == y)
    if (length(idx) > 0) {
      out[idx, ] <- sweep(out[idx, , drop = FALSE], 2, colMeans(out[idx, , drop = FALSE]))
    }
  }
  out
}

case_control_cov <- function(IF, Y, fixed_case_control = TRUE) {
  if (fixed_case_control) {
    IF <- center_by_case_control_stratum(IF, Y)
  } else {
    IF <- sweep(IF, 2, colMeans(IF), "-")
  }
  crossprod(IF) / nrow(IF)
}

case_control_crosscov <- function(IF_left, IF_right, Y, fixed_case_control = TRUE) {
  if (nrow(IF_left) != nrow(IF_right)) {
    stop("IF_left and IF_right must have the same number of rows.")
  }
  if (fixed_case_control) {
    IF_left <- center_by_case_control_stratum(IF_left, Y)
    IF_right <- center_by_case_control_stratum(IF_right, Y)
  } else {
    IF_left <- sweep(IF_left, 2, colMeans(IF_left), "-")
    IF_right <- sweep(IF_right, 2, colMeans(IF_right), "-")
  }
  crossprod(IF_left, IF_right) / nrow(IF_left)
}

el_reconstructed_crosscov <- function(A0, A1, B0, B1,
                                      p_hat, Delta_hat, rho,
                                      centering = c("separate", "compact"),
                                      compact_tol = ZIM_DEFAULT_COMPACT_CENTERING_TOL) {
  centering <- match.arg(centering)
  compact_tol <- validate_positive_tolerance(compact_tol, "compact_tol")
  A0 <- as.matrix(A0)
  A1 <- as.matrix(A1)
  B0 <- as.matrix(B0)
  B1 <- as.matrix(B1)
  if (length(rho) != 1L || !is.finite(rho) || rho <= 0) stop("rho must be a positive finite scalar.")
  p0 <- as.numeric(p_hat)
  p1 <- as.numeric(p_hat * Delta_hat)
  if (nrow(A0) != length(p0) || nrow(A1) != length(p0) ||
      nrow(B0) != length(p0) || nrow(B1) != length(p0)) {
    stop("EL reconstructed rows and weights must have the same number of rows.")
  }
  mass0 <- sum(p0)
  mass1 <- sum(p1)
  if (!is.finite(mass0) || mass0 <= 0 || !is.finite(mass1) || mass1 <= 0) {
    stop("EL reconstructed masses must be positive and finite.")
  }
  p0n <- p0 / mass0
  p1n <- p1 / mass1
  w0 <- p0n / (1 + rho)
  w1 <- rho * p1n / (1 + rho)
  muA0 <- colSums(A0 * p0n)
  muB0 <- colSums(B0 * p0n)
  muA1 <- colSums(A1 * p1n)
  muB1 <- colSums(B1 * p1n)
  if (identical(centering, "separate")) {
    A0c <- sweep(A0, 2, muA0, "-")
    B0c <- sweep(B0, 2, muB0, "-")
    A1c <- sweep(A1, 2, muA1, "-")
    B1c <- sweep(B1, 2, muB1, "-")
    return(crossprod(A0c, B0c * w0) + crossprod(A1c, B1c * w1))
  }
  center_A <- muA0 + rho * muA1
  center_B <- muB0 + rho * muB1
  scale_A <- max(1, max(abs(muA0)), rho * max(abs(muA1)))
  scale_B <- max(1, max(abs(muB0)), rho * max(abs(muB1)))
  resid_A <- max(abs(center_A))
  resid_B <- max(abs(center_B))
  if (!is.finite(resid_A) || !is.finite(resid_B) ||
      resid_A > compact_tol * scale_A || resid_B > compact_tol * scale_B) {
    stop(
      "Compact EL centering identity is not satisfied: max residual for A = ",
      signif(resid_A, 6), ", max residual for B = ", signif(resid_B, 6),
      ", relative tolerance = ", signif(compact_tol, 6),
      ". Use centering = 'separate' or solve the relevant identities more tightly."
    )
  }
  crossprod(A0, B0 * w0) + crossprod(A1, B1 * w1) - tcrossprod(muA0, muB0) / rho
}

el_reconstructed_cov <- function(A0, A1, p_hat, Delta_hat, rho,
                                 centering = c("separate", "compact"),
                                 compact_tol = ZIM_DEFAULT_COMPACT_CENTERING_TOL) {
  out <- el_reconstructed_crosscov(
    A0, A1, A0, A1,
    p_hat = p_hat,
    Delta_hat = Delta_hat,
    rho = rho,
    centering = centering,
    compact_tol = compact_tol
  )
  symmetrize(out)
}

formula6_beta_cov_case_control <- function(Y, working, gamma_tilde_list, rho_k,
                                           target_delta, ridge) {
  ridge <- validate_ridge(ridge)
  Y <- check_binary_y(Y)
  if (length(target_delta) != length(Y)) {
    stop("target_delta must have one value per external subject.")
  }
  if (any(!is.finite(target_delta)) || any(target_delta <= 0)) {
    stop("target_delta must be finite and positive.")
  }
  if (max(abs(rho_k - rho_k[1])) > 1e-8) {
    stop("formula6 external covariance currently requires a common rho_k.")
  }

  idx0 <- which(Y == 0)
  if (length(idx0) < 2L) stop("formula6 external covariance needs at least two controls.")
  rho <- rho_k[1]
  delta_target0 <- target_delta[idx0]

  beta_widths <- vapply(working, function(wk) wk$p_beta_k, integer(1))
  beta_starts <- cumsum(c(1L, head(beta_widths, -1L)))
  p_beta <- sum(beta_widths)
  U0_beta <- matrix(NA_real_, nrow = length(idx0), ncol = p_beta)
  U1_beta <- matrix(NA_real_, nrow = length(idx0), ncol = p_beta)
  beta_names <- unlist(lapply(working, function(wk) wk$beta_names), use.names = FALSE)
  colnames(U0_beta) <- beta_names
  colnames(U1_beta) <- beta_names

  for (k in seq_along(working)) {
    wk <- working[[k]]
    D0 <- wk$D_k[idx0, , drop = FALSE]
    gamma_k <- gamma_tilde_list[[k]]
    eta0 <- drop(D0 %*% gamma_k)
    delta0 <- safe_exp(eta0)
    denom0 <- 1 + rho * delta0

    phi0 <- D0 * as.numeric(-rho * delta0 / denom0)
    phi1 <- D0 * as.numeric(1 / denom0)

    deriv_weight <- -rho * delta0 * (1 + rho * delta_target0) / (denom0^2)
    A_k <- crossprod(D0, D0 * as.numeric(deriv_weight)) /
      (length(idx0) * (rho + 1))

    beta_idx <- wk$q_alpha_k + seq_len(wk$p_beta_k)
    out_idx <- beta_starts[k] + seq_len(wk$p_beta_k) - 1L
    U0_beta[, out_idx] <- right_solve_transpose(A_k, phi0, ridge = ridge)[, beta_idx, drop = FALSE]
    U1_beta[, out_idx] <- right_solve_transpose(A_k, phi1, ridge = ridge)[, beta_idx, drop = FALSE]
  }

  w0 <- 1 / (rho + 1)
  w1 <- rho * delta_target0 / (rho + 1)
  mean_u0 <- colMeans(U0_beta)

  H_beta <- w0 * crossprod(U0_beta) / length(idx0) +
    crossprod(U1_beta, U1_beta * as.numeric(w1)) / length(idx0) -
    w0 * tcrossprod(mean_u0)
  symmetrize(H_beta)
}

estimate_sigma_betabeta_one_external_cc <- function(
    working,
    gamma_plugin_list,
    Delta_hat,
    p_hat,
    n_internal,
    external_n_case,
    external_n_control,
    ridge,
    B_centering = c("separate", "compact"),
    mass_tol = ZIM_DEFAULT_EL_MASS_TOL,
    mass_action = c("stop", "warn", "none"),
    compact_tol = ZIM_DEFAULT_COMPACT_CENTERING_TOL,
    rank_tol = ZIM_DEFAULT_BREAD_RANK_TOL,
    condition_warn = ZIM_DEFAULT_CONDITION_WARN,
    psd_tol = ZIM_DEFAULT_PSD_TOL) {
  ridge <- validate_ridge(ridge)
  B_centering <- match.arg(B_centering)
  mass_action <- match.arg(mass_action)
  mass_tol <- validate_positive_tolerance(mass_tol, "mass_tol")
  compact_tol <- validate_positive_tolerance(compact_tol, "compact_tol")
  rank_tol <- validate_positive_tolerance(rank_tol, "rank_tol")
  condition_warn <- validate_positive_tolerance(condition_warn, "condition_warn")
  psd_tol <- validate_positive_tolerance(psd_tol, "psd_tol")
  ext <- validate_external_case_control_counts(external_n_case, external_n_control)
  rho_E <- ext$rho_E
  n_external <- ext$n_external
  n <- length(Delta_hat)
  if (length(p_hat) != n) stop("p_hat and Delta_hat must have the same length.")
  if (length(working) != length(gamma_plugin_list)) {
    stop("gamma_plugin_list must have one full retrospective gamma vector per working model.")
  }
  if (length(n_internal) != 1L || !is.finite(n_internal) || n_internal <= 0) {
    stop("n_internal must be a positive finite scalar.")
  }
  mass_diagnostics <- validate_el_masses(
    p_hat = p_hat, Delta_hat = Delta_hat,
    mass_tol = mass_tol, action = mass_action,
    context = "Automatic one-external-study covariance"
  )
  beta_widths <- vapply(working, function(wk) wk$p_beta_k, integer(1))
  p_beta <- sum(beta_widths)
  gamma_widths <- vapply(working, function(wk) ncol(wk$D_k), integer(1))
  p_gamma <- sum(gamma_widths)
  IF0_beta_E <- matrix(0, nrow = n, ncol = p_beta)
  IF1_beta_E <- matrix(0, nrow = n, ncol = p_beta)
  A_E_blocks <- vector("list", length(working))
  A_E_rank_diagnostics <- vector("list", length(working))
  beta_index_gamma <- integer(0)
  gamma_names <- character(0)
  beta_names <- unlist(lapply(working, `[[`, "beta_names"), use.names = FALSE)
  colnames(IF0_beta_E) <- colnames(IF1_beta_E) <- beta_names
  a0_E <- numeric(p_gamma)
  a1_E <- numeric(p_gamma)
  gamma_col_start <- 1L
  beta_col_start <- 1L
  for (k in seq_along(working)) {
    wk <- working[[k]]
    D_k <- wk$D_k
    gamma_k <- gamma_plugin_list[[k]]
    if (length(gamma_k) != ncol(D_k)) stop("gamma_plugin_list[[", k, "]] has the wrong length.")
    pieces <- working_score_pieces(D_k, gamma_k, rho = rho_E)
    deriv_scalar <- -rho_E * pieces$delta / (1 + rho_E * pieces$delta)^2
    A_E_k <- crossprod(
      D_k,
      D_k * as.numeric(p_hat * (1 + rho_E * Delta_hat) / (1 + rho_E) * deriv_scalar)
    )
    A_E_blocks[[k]] <- A_E_k
    A_E_diag <- matrix_rank_diagnostics(A_E_k, rank_tol = rank_tol)
    A_E_rank_diagnostics[[k]] <- A_E_diag
    if (!isTRUE(A_E_diag$full_rank)) {
      stop("External working bread for model ", k, " is rank deficient: rank ",
           A_E_diag$rank, " of ", A_E_diag$dimension, ".")
    }
    warn_on_condition_number(A_E_diag, threshold = condition_warn,
                             context = paste0("External working bread for model ", k))
    beta_local <- wk$q_alpha_k + seq_len(wk$p_beta_k)
    P_beta_k <- matrix(0, nrow = wk$p_beta_k, ncol = ncol(D_k))
    P_beta_k[cbind(seq_len(wk$p_beta_k), beta_local)] <- 1
    A_E_k_reg <- regularize_negative_bread(A_E_k, ridge = ridge)
    AinvT_PbetaT <- solve(t(A_E_k_reg), t(P_beta_k))
    beta_cols <- beta_col_start:(beta_col_start + wk$p_beta_k - 1L)
    IF0_beta_E[, beta_cols] <- -pieces$phi0 %*% AinvT_PbetaT
    IF1_beta_E[, beta_cols] <- -pieces$phi1 %*% AinvT_PbetaT
    gamma_cols <- gamma_col_start:(gamma_col_start + ncol(D_k) - 1L)
    a0_E[gamma_cols] <- colSums(pieces$phi0 * as.numeric(p_hat))
    a1_E[gamma_cols] <- colSums(pieces$phi1 * as.numeric(p_hat * Delta_hat))
    beta_index_gamma <- c(beta_index_gamma, gamma_cols[beta_local])
    gamma_names <- c(gamma_names, wk$gamma_names)
    gamma_col_start <- gamma_col_start + ncol(D_k)
    beta_col_start <- beta_col_start + wk$p_beta_k
  }
  names(a0_E) <- names(a1_E) <- gamma_names
  names(A_E_blocks) <- names(A_E_rank_diagnostics) <- paste0("working_", seq_along(working))
  mass0_E <- mass_diagnostics$mass0
  mass1_E <- mass_diagnostics$mass1
  mu0_E <- a0_E / mass0_E
  mu1_E <- a1_E / mass1_E
  center_resid_E <- a0_E + rho_E * a1_E
  p0n_E <- as.numeric(p_hat) / mass0_E
  p1n_E <- as.numeric(p_hat * Delta_hat) / mass1_E
  a0_beta_IF_E <- colSums(IF0_beta_E * as.numeric(p_hat))
  a1_beta_IF_E <- colSums(IF1_beta_E * as.numeric(p_hat * Delta_hat))
  center_resid_beta_IF_E <- a0_beta_IF_E + rho_E * a1_beta_IF_E
  a0_beta_IF_E_normalized <- colSums(IF0_beta_E * p0n_E)
  a1_beta_IF_E_normalized <- colSums(IF1_beta_E * p1n_E)
  center_resid_beta_IF_E_normalized <- a0_beta_IF_E_normalized + rho_E * a1_beta_IF_E_normalized
  H_ext_beta_beta <- el_reconstructed_cov(
    A0 = IF0_beta_E, A1 = IF1_beta_E,
    p_hat = p_hat, Delta_hat = Delta_hat, rho = rho_E,
    centering = B_centering, compact_tol = compact_tol
  )
  H_ext_psd_diagnostics <- matrix_psd_diagnostics(H_ext_beta_beta, tol = psd_tol)
  if (!isTRUE(H_ext_psd_diagnostics$is_psd)) {
    stop("Automatic external beta covariance is not positive semidefinite: minimum eigenvalue = ",
         signif(H_ext_psd_diagnostics$minimum_eigenvalue, 6), ".")
  }
  B_centering_used <- if (identical(B_centering, "separate")) {
    "separate_EL_external_stratum_centering_beta_IF_normalized_masses"
  } else "validated_compact_fixed_stratum_beta_IF_normalized_masses"
  dimnames(H_ext_beta_beta) <- list(beta_names, beta_names)
  Sigma_beta_beta <- symmetrize((n_internal / n_external) * H_ext_beta_beta)
  dimnames(Sigma_beta_beta) <- list(beta_names, beta_names)
  list(
    Sigma_beta_beta = Sigma_beta_beta,
    H_ext_beta_beta = H_ext_beta_beta,
    Sigma_gamma_E = NULL,
    A_E = NULL,
    B_E = NULL,
    A_E_blocks = A_E_blocks,
    A_E_rank_diagnostics = A_E_rank_diagnostics,
    H_ext_psd_diagnostics = H_ext_psd_diagnostics,
    B_centering = B_centering_used,
    B_centering_requested = B_centering,
    mu0_E = mu0_E,
    mu1_E = mu1_E,
    a0_E = a0_E,
    a1_E = a1_E,
    center_resid_E = center_resid_E,
    max_abs_center_resid_E = max(abs(center_resid_E)),
    a0_beta_IF_E = a0_beta_IF_E,
    a1_beta_IF_E = a1_beta_IF_E,
    center_resid_beta_IF_E = center_resid_beta_IF_E,
    max_abs_center_resid_beta_IF_E = max(abs(center_resid_beta_IF_E)),
    a0_beta_IF_E_normalized = a0_beta_IF_E_normalized,
    a1_beta_IF_E_normalized = a1_beta_IF_E_normalized,
    center_resid_beta_IF_E_normalized = center_resid_beta_IF_E_normalized,
    max_abs_center_resid_beta_IF_E_normalized = max(abs(center_resid_beta_IF_E_normalized)),
    EL_mass0_E = mass0_E,
    EL_mass1_E = mass1_E,
    EL_mass_diagnostics = mass_diagnostics,
    rho_E = rho_E,
    n_external = n_external,
    external_n_case = ext$external_n_case,
    external_n_control = ext$external_n_control,
    beta_index_gamma = beta_index_gamma,
    gamma_names = gamma_names,
    beta_names = beta_names,
    p_gamma = p_gamma,
    p_beta = p_beta,
    shortcut_used = TRUE,
    method = "retrospective_external_fixed_stratum_one_external_study_beta_IF_shortcut"
  )
}

align_gamma_start_to_working <- function(gamma_start, wk, what = "gamma_start") {
  input_names <- names(gamma_start)
  gamma_start <- as.numeric(gamma_start)
  if (!is.null(input_names) && length(input_names) == length(gamma_start)) {
    names(gamma_start) <- input_names
  }
  if (!is.null(input_names) && length(input_names) == length(gamma_start) &&
      all(!is.na(input_names)) && all(nzchar(input_names))) {
    if (all(wk$gamma_names %in% input_names)) {
      gamma_start <- gamma_start[match(wk$gamma_names, input_names)]
    } else {
      stop(what, " names do not match the reference working-model gamma names.")
    }
  } else {
    if (length(gamma_start) != ncol(wk$D_k)) {
      stop(what, " has the wrong length for the reference working model.")
    }
  }
  if (any(!is.finite(gamma_start))) {
    warning(what, " contains non-finite starting values; replacing them with zero.")
    gamma_start[!is.finite(gamma_start)] <- 0
  }
  names(gamma_start) <- wk$gamma_names
  gamma_start
}

calibrate_external_working_intercept_start <- function(wk_reference,
                                                       gamma_start,
                                                       Delta_ref,
                                                       rho_E,
                                                       bracket_width = 20,
                                                       max_expand = 6L) {
  gamma_start <- align_gamma_start_to_working(
    gamma_start = gamma_start,
    wk = wk_reference,
    what = "gamma_start"
  )
  if (length(Delta_ref) != nrow(wk_reference$D_k)) {
    stop("Delta_ref must have one value per reference control.")
  }
  if (any(!is.finite(Delta_ref)) || any(Delta_ref <= 0)) {
    stop("Delta_ref must be finite and positive.")
  }
  if (length(rho_E) != 1L || !is.finite(rho_E) || rho_E <= 0) {
    stop("rho_E must be a positive finite scalar.")
  }

  intercept_idx <- retrospective_intercept_index(wk_reference$gamma_names)
  if (length(intercept_idx) == 0L) return(gamma_start)

  intercept_col <- wk_reference$D_k[, intercept_idx[1], drop = TRUE]
  if (max(abs(intercept_col - 1)) > 1e-8) {
    stop("The external reference working intercept column must be a column of ones.")
  }

  D_ref <- wk_reference$D_k

  L <- drop(D_ref %*% gamma_start) - gamma_start[intercept_idx[1]]
  alpha_norm <- -log_mean_exp(L)

  score_intercept <- function(alpha0) {
    delta <- safe_exp(alpha0 + L)
    sum((Delta_ref - delta) / (1 + rho_E * delta))
  }

  alpha0 <- alpha_norm
  lo <- alpha_norm - bracket_width
  hi <- alpha_norm + bracket_width
  f_lo <- score_intercept(lo)
  f_hi <- score_intercept(hi)

  expand_iter <- 0L
  while (is.finite(f_lo) && is.finite(f_hi) && f_lo * f_hi > 0 &&
         expand_iter < as.integer(max_expand)) {
    bracket_width <- bracket_width * 2
    lo <- alpha_norm - bracket_width
    hi <- alpha_norm + bracket_width
    f_lo <- score_intercept(lo)
    f_hi <- score_intercept(hi)
    expand_iter <- expand_iter + 1L
  }

  if (is.finite(f_lo) && is.finite(f_hi) && f_lo * f_hi <= 0) {
    alpha0 <- stats::uniroot(score_intercept, lower = lo, upper = hi)$root
  } else {
    warning(
      "Could not bracket the external working intercept score; using normalized intercept start."
    )
  }

  gamma_start[intercept_idx[1]] <- alpha0
  names(gamma_start) <- wk_reference$gamma_names
  gamma_start
}

reference_reconstructed_crosscov <- function(
    A0, A1, B0, B1, Delta_ref, rho_E,
    centering = c("separate", "compact"),
    compact_tol = ZIM_DEFAULT_COMPACT_CENTERING_TOL) {
  centering <- match.arg(centering)
  compact_tol <- validate_positive_tolerance(compact_tol, "compact_tol")
  A0 <- as.matrix(A0)
  A1 <- as.matrix(A1)
  B0 <- as.matrix(B0)
  B1 <- as.matrix(B1)
  Delta_ref <- as.numeric(Delta_ref)
  m_ref <- length(Delta_ref)
  if (m_ref <= 1L) stop("Reference covariance requires at least two controls.")
  if (length(rho_E) != 1L || !is.finite(rho_E) || rho_E <= 0) {
    stop("rho_E must be a positive finite scalar.")
  }
  if (nrow(A0) != m_ref || nrow(A1) != m_ref ||
      nrow(B0) != m_ref || nrow(B1) != m_ref) {
    stop("Reference reconstructed rows and Delta_ref must have matching row counts.")
  }
  if (any(!is.finite(Delta_ref)) || any(Delta_ref <= 0)) {
    stop("Delta_ref must be finite and positive.")
  }

  if (identical(centering, "separate")) {
    p0 <- rep(1 / m_ref, m_ref)
    p1 <- Delta_ref / sum(Delta_ref)
    muA0 <- colSums(A0 * p0)
    muB0 <- colSums(B0 * p0)
    muA1 <- colSums(A1 * p1)
    muB1 <- colSums(B1 * p1)
    A0c <- sweep(A0, 2, muA0, "-")
    B0c <- sweep(B0, 2, muB0, "-")
    A1c <- sweep(A1, 2, muA1, "-")
    B1c <- sweep(B1, 2, muB1, "-")
    w0 <- p0 / (1 + rho_E)
    w1 <- rho_E * p1 / (1 + rho_E)
    return(
      crossprod(A0c, B0c * as.numeric(w0)) +
        crossprod(A1c, B1c * as.numeric(w1))
    )
  }

  mass1 <- mean(Delta_ref)
  aA0 <- colMeans(A0)
  aB0 <- colMeans(B0)
  aA1 <- colMeans(A1 * as.numeric(Delta_ref))
  aB1 <- colMeans(B1 * as.numeric(Delta_ref))
  residA <- max(abs(aA0 + rho_E * aA1))
  residB <- max(abs(aB0 + rho_E * aB1))
  scaleA <- max(1, max(abs(aA0)), rho_E * max(abs(aA1)))
  scaleB <- max(1, max(abs(aB0)), rho_E * max(abs(aB1)))
  if (abs(mass1 - 1) > compact_tol ||
      !is.finite(residA) || residA > compact_tol * scaleA ||
      !is.finite(residB) || residB > compact_tol * scaleB) {
    stop(
      "Compact reference centering identity is not satisfied: mean(Delta_ref) = ",
      signif(mass1, 8), ", residual A = ", signif(residA, 8),
      ", residual B = ", signif(residB, 8),
      ", tolerance = ", signif(compact_tol, 8),
      ". Use centering = 'separate'."
    )
  }
  w0 <- rep(1 / ((1 + rho_E) * m_ref), m_ref)
  w1 <- rho_E * Delta_ref / ((1 + rho_E) * m_ref)
  crossprod(A0, B0 * as.numeric(w0)) +
    crossprod(A1, B1 * as.numeric(w1)) -
    tcrossprod(aA0, aB0) / rho_E
}

reference_reconstructed_cov <- function(
    A0, A1, Delta_ref, rho_E,
    centering = c("separate", "compact"),
    compact_tol = ZIM_DEFAULT_COMPACT_CENTERING_TOL) {
  out <- reference_reconstructed_crosscov(
    A0 = A0, A1 = A1, B0 = A0, B1 = A1,
    Delta_ref = Delta_ref, rho_E = rho_E,
    centering = centering, compact_tol = compact_tol
  )
  symmetrize(out)
}

estimate_sigma_betabeta_reference_cc <- function(
    working_reference,
    gamma_plugin_list,
    Delta_ref,
    n_internal,
    external_n_case,
    external_n_control,
    ridge,
    B_centering = c("separate", "compact"),
    compact_tol = ZIM_DEFAULT_COMPACT_CENTERING_TOL,
    rank_tol = ZIM_DEFAULT_BREAD_RANK_TOL,
    condition_warn = ZIM_DEFAULT_CONDITION_WARN,
    psd_tol = ZIM_DEFAULT_PSD_TOL) {
  ridge <- validate_ridge(ridge)
  B_centering <- match.arg(B_centering)
  compact_tol <- validate_positive_tolerance(compact_tol, "compact_tol")
  rank_tol <- validate_positive_tolerance(rank_tol, "rank_tol")
  condition_warn <- validate_positive_tolerance(condition_warn, "condition_warn")
  psd_tol <- validate_positive_tolerance(psd_tol, "psd_tol")
  ext <- validate_external_case_control_counts(
    external_n_case, external_n_control
  )
  rho_E <- ext$rho_E
  n_external <- ext$n_external

  m_ref <- length(Delta_ref)
  if (m_ref <= 1L) stop("reference_data must contain at least two controls.")
  if (any(!is.finite(Delta_ref)) || any(Delta_ref <= 0)) {
    stop("Delta_ref must be finite and positive.")
  }
  if (length(working_reference) != length(gamma_plugin_list)) {
    stop(
      "gamma_plugin_list must have one full retrospective gamma vector per ",
      "working model."
    )
  }
  if (length(n_internal) != 1L || !is.finite(n_internal) || n_internal <= 0) {
    stop("n_internal must be a positive finite scalar.")
  }

  beta_widths <- vapply(
    working_reference, function(wk) wk$p_beta_k, integer(1)
  )
  gamma_widths <- vapply(
    working_reference, function(wk) ncol(wk$D_k), integer(1)
  )
  p_beta <- sum(beta_widths)
  p_gamma <- sum(gamma_widths)
  beta_names <- unlist(
    lapply(working_reference, `[[`, "beta_names"), use.names = FALSE
  )
  gamma_names <- unlist(
    lapply(working_reference, `[[`, "gamma_names"), use.names = FALSE
  )

  IF0_beta_external <- matrix(0, nrow = m_ref, ncol = p_beta)
  IF1_beta_external <- matrix(0, nrow = m_ref, ncol = p_beta)
  colnames(IF0_beta_external) <- colnames(IF1_beta_external) <- beta_names

  A_E_blocks <- vector("list", length(working_reference))
  A_E_reg_blocks <- vector("list", length(working_reference))
  A_E_rank_diagnostics <- vector("list", length(working_reference))
  A_E_regularized_rank_diagnostics <- vector(
    "list", length(working_reference)
  )
  working_exp_diagnostics <- vector("list", length(working_reference))
  beta_index_gamma <- integer(0)
  a0_E <- numeric(p_gamma)
  a1_E <- numeric(p_gamma)

  beta_col_start <- 1L
  gamma_col_start <- 1L
  for (k in seq_along(working_reference)) {
    wk <- working_reference[[k]]
    Qk <- wk$D_k
    gamma_k <- align_gamma_start_to_working(
      gamma_start = gamma_plugin_list[[k]],
      wk = wk,
      what = paste0("gamma_plugin_list[[", k, "]]" )
    )

    eta_k <- drop(Qk %*% gamma_k)
    working_exp_diagnostics[[k]] <- check_exp_range(
      eta_k,
      context = paste0("Final reconstructed external working model ", k),
      action = "stop"
    )
    delta_k <- exp(eta_k)
    denom_k <- 1 + rho_E * delta_k

    phi0_k <- Qk * as.numeric(-rho_E * delta_k / denom_k)
    phi1_k <- Qk * as.numeric(1 / denom_k)

    deriv_weight <- -rho_E * delta_k * (1 + rho_E * Delta_ref) /
      ((1 + rho_E) * denom_k^2)
    A_E_k <- crossprod(Qk, Qk * as.numeric(deriv_weight)) / m_ref
    A_E_blocks[[k]] <- A_E_k

    A_diag <- matrix_rank_diagnostics(A_E_k, rank_tol = rank_tol)
    A_E_rank_diagnostics[[k]] <- A_diag
    if (!isTRUE(A_diag$full_rank)) {
      stop(
        "Reference external working bread for model ", k,
        " is rank deficient: rank ", A_diag$rank,
        " of ", A_diag$dimension, "."
      )
    }
    warn_on_condition_number(
      A_diag,
      threshold = condition_warn,
      context = paste0("Reference external working bread for model ", k)
    )

    A_E_k_reg <- regularize_negative_bread(A_E_k, ridge = ridge)
    A_E_reg_blocks[[k]] <- A_E_k_reg
    A_E_regularized_rank_diagnostics[[k]] <- matrix_rank_diagnostics(
      A_E_k_reg, rank_tol = rank_tol
    )

    beta_local <- wk$q_alpha_k + seq_len(wk$p_beta_k)
    S_k <- matrix(0, nrow = wk$p_beta_k, ncol = ncol(Qk))
    S_k[cbind(seq_len(wk$p_beta_k), beta_local)] <- 1
    dimnames(S_k) <- list(wk$beta_names, wk$gamma_names)
    L_beta_k <- right_solve(A_E_k_reg, S_k)

    beta_cols <- beta_col_start:(beta_col_start + wk$p_beta_k - 1L)
    IF0_beta_external[, beta_cols] <- -phi0_k %*% t(L_beta_k)
    IF1_beta_external[, beta_cols] <- -phi1_k %*% t(L_beta_k)

    gamma_cols <- gamma_col_start:(gamma_col_start + ncol(Qk) - 1L)
    a0_E[gamma_cols] <- colMeans(phi0_k)
    a1_E[gamma_cols] <- colMeans(phi1_k * as.numeric(Delta_ref))
    beta_index_gamma <- c(beta_index_gamma, gamma_cols[beta_local])

    beta_col_start <- beta_col_start + wk$p_beta_k
    gamma_col_start <- gamma_col_start + ncol(Qk)
  }

  names(A_E_blocks) <- names(A_E_reg_blocks) <-
    names(A_E_rank_diagnostics) <-
    names(A_E_regularized_rank_diagnostics) <-
    names(working_exp_diagnostics) <- paste0(
      "working_", seq_along(working_reference)
    )
  names(a0_E) <- names(a1_E) <- gamma_names

  mass0_E <- 1
  mass1_E <- mean(Delta_ref)
  if (!is.finite(mass1_E) || mass1_E <= 0) {
    stop("The reconstructed external case mass must be positive and finite.")
  }
  mu0_E <- a0_E / mass0_E
  mu1_E <- a1_E / mass1_E
  center_resid_E <- a0_E + rho_E * a1_E

  p0n_E <- rep(1 / m_ref, m_ref)
  p1n_E <- Delta_ref / sum(Delta_ref)
  a0_beta_IF_E <- colMeans(IF0_beta_external)
  a1_beta_IF_E <- colMeans(
    IF1_beta_external * as.numeric(Delta_ref)
  )
  center_resid_beta_IF_E <- a0_beta_IF_E + rho_E * a1_beta_IF_E
  a0_beta_IF_E_normalized <- colSums(
    IF0_beta_external * as.numeric(p0n_E)
  )
  a1_beta_IF_E_normalized <- colSums(
    IF1_beta_external * as.numeric(p1n_E)
  )
  center_resid_beta_IF_E_normalized <-
    a0_beta_IF_E_normalized + rho_E * a1_beta_IF_E_normalized

  H_ext_beta_beta <- reference_reconstructed_cov(
    A0 = IF0_beta_external,
    A1 = IF1_beta_external,
    Delta_ref = Delta_ref,
    rho_E = rho_E,
    centering = B_centering,
    compact_tol = compact_tol
  )
  dimnames(H_ext_beta_beta) <- list(beta_names, beta_names)
  H_ext_psd_diagnostics <- matrix_psd_diagnostics(
    H_ext_beta_beta, tol = psd_tol
  )
  if (!isTRUE(H_ext_psd_diagnostics$is_psd)) {
    stop(
      "Reference external beta influence covariance is not positive ",
      "semidefinite: minimum eigenvalue = ",
      signif(H_ext_psd_diagnostics$minimum_eigenvalue, 6), "."
    )
  }

  Sigma_beta_beta <- symmetrize(
    (n_internal / n_external) * H_ext_beta_beta
  )
  dimnames(Sigma_beta_beta) <- list(beta_names, beta_names)
  Sigma_psd_diagnostics <- matrix_psd_diagnostics(
    Sigma_beta_beta, tol = psd_tol
  )
  if (!isTRUE(Sigma_psd_diagnostics$is_psd)) {
    stop(
      "Reference external beta covariance is not positive semidefinite: ",
      "minimum eigenvalue = ",
      signif(Sigma_psd_diagnostics$minimum_eigenvalue, 6), "."
    )
  }

  B_centering_used <- if (identical(B_centering, "separate")) {
    "separate_reference_external_stratum_centering_beta_IF"
  } else {
    "validated_compact_reference_external_beta_IF"
  }

  list(
    Sigma_beta_beta = Sigma_beta_beta,
    H_ext_beta_beta = H_ext_beta_beta,

    Sigma_gamma_E = NULL,
    A_E = NULL,
    A_E_reg = NULL,
    B_E = NULL,
    A_E_blocks = A_E_blocks,
    A_E_reg_blocks = A_E_reg_blocks,
    A_E_rank_diagnostics = A_E_rank_diagnostics,
    A_E_regularized_rank_diagnostics =
      A_E_regularized_rank_diagnostics,
    working_exp_diagnostics = working_exp_diagnostics,
    H_ext_psd_diagnostics = H_ext_psd_diagnostics,
    Sigma_psd_diagnostics = Sigma_psd_diagnostics,
    IF0_beta_external = IF0_beta_external,
    IF1_beta_external = IF1_beta_external,
    B_centering = B_centering_used,
    B_centering_requested = B_centering,
    mu0_E = mu0_E,
    mu1_E = mu1_E,
    a0_E = a0_E,
    a1_E = a1_E,
    center_resid_E = center_resid_E,
    max_abs_center_resid_E = max(abs(center_resid_E)),
    a0_beta_IF_E = a0_beta_IF_E,
    a1_beta_IF_E = a1_beta_IF_E,
    center_resid_beta_IF_E = center_resid_beta_IF_E,
    max_abs_center_resid_beta_IF_E = max(
      abs(center_resid_beta_IF_E)
    ),
    a0_beta_IF_E_normalized = a0_beta_IF_E_normalized,
    a1_beta_IF_E_normalized = a1_beta_IF_E_normalized,
    center_resid_beta_IF_E_normalized =
      center_resid_beta_IF_E_normalized,
    max_abs_center_resid_beta_IF_E_normalized = max(
      abs(center_resid_beta_IF_E_normalized)
    ),
    reference_mass0_E = mass0_E,
    reference_mass1_E = mass1_E,
    Delta_ref_mean = mean(Delta_ref),
    Delta_ref_min = min(Delta_ref),
    Delta_ref_max = max(Delta_ref),
    rho_E = rho_E,
    n_external = n_external,
    external_n_case = ext$external_n_case,
    external_n_control = ext$external_n_control,
    n_reference = m_ref,
    beta_index_gamma = beta_index_gamma,
    gamma_names = gamma_names,
    beta_names = beta_names,
    p_gamma = p_gamma,
    p_beta = p_beta,
    blockwise_beta_IF = TRUE,
    full_gamma_sandwich_formed = FALSE,
    method = paste0(
      "retrospective_external_fixed_stratum_reference_control_",
      "blockwise_beta_IF"
    )
  )
}

make_target_el_pseudo_if <- function(target,
                                     Delta_hat,
                                     rho_I,
                                     A_eta_eta_reg) {
  Z <- as.matrix(target$Z)
  d_eta <- ncol(Z)
  theta_index_eta <- if (d_eta >= 2L) 2L:d_eta else integer(0)
  if (length(theta_index_eta) == 0L) stop("The target model must contain at least one theta parameter.")
  if (length(Delta_hat) != nrow(Z)) stop("Delta_hat must have one value per internal subject.")
  if (nrow(A_eta_eta_reg) != d_eta || ncol(A_eta_eta_reg) != d_eta) {
    stop("A_eta_eta_reg has the wrong dimension.")
  }
  denom_target <- 1 + rho_I * Delta_hat
  s0 <- Z * as.numeric(-rho_I * Delta_hat / denom_target)
  s1 <- Z * as.numeric(1 / denom_target)
  P_theta <- matrix(0, nrow = length(theta_index_eta), ncol = d_eta)
  P_theta[cbind(seq_along(theta_index_eta), theta_index_eta)] <- 1
  K_theta <- right_solve(A_eta_eta_reg, P_theta)
  IF0_theta <- -s0 %*% t(K_theta)
  IF1_theta <- -s1 %*% t(K_theta)
  theta_names <- names(target$theta_hat)
  colnames(IF0_theta) <- colnames(IF1_theta) <- theta_names
  list(
    s0 = s0,
    s1 = s1,
    P_theta = P_theta,
    K_theta = K_theta,
    theta_index_eta = theta_index_eta,
    IF0_theta = IF0_theta,
    IF1_theta = IF1_theta
  )
}

estimate_H_onepop_gmm_el_cc <- function(target,
                                        working,
                                        gamma_hat_list,
                                        rho_I,
                                        p_hat,
                                        Delta_hat,
                                        rho_working,
                                        ridge,
                                        B_centering = c("separate", "compact"),
                                        return_el_pseudo_if = FALSE,
                                        mass_tol = ZIM_DEFAULT_EL_MASS_TOL,
                                        mass_action = c("stop", "warn", "none"),
                                        compact_tol = ZIM_DEFAULT_COMPACT_CENTERING_TOL,
                                        rank_tol = ZIM_DEFAULT_BREAD_RANK_TOL,
                                        condition_warn = ZIM_DEFAULT_CONDITION_WARN,
                                        joint_psd_tol = ZIM_DEFAULT_PSD_TOL) {
  ridge <- validate_ridge(ridge)
  B_centering <- match.arg(B_centering)
  mass_action <- match.arg(mass_action)
  mass_tol <- validate_positive_tolerance(mass_tol, "mass_tol")
  compact_tol <- validate_positive_tolerance(compact_tol, "compact_tol")
  rank_tol <- validate_positive_tolerance(rank_tol, "rank_tol")
  condition_warn <- validate_positive_tolerance(condition_warn, "condition_warn")
  joint_psd_tol <- validate_positive_tolerance(joint_psd_tol, "joint_psd_tol")
  if (length(rho_I) != 1L || !is.finite(rho_I) || rho_I <= 0) stop("rho_I must be a positive finite scalar.")
  if (length(rho_working) != 1L || !is.finite(rho_working) || rho_working <= 0) {
    stop("rho_working must be a positive finite scalar.")
  }
  Z <- target$Z
  n <- nrow(Z)
  d_eta <- ncol(Z)
  if (length(Delta_hat) != n || length(p_hat) != n) {
    stop("Delta_hat and p_hat must have one value per internal subject.")
  }
  if (length(working) != length(gamma_hat_list)) {
    stop("gamma_hat_list must have one full retrospective gamma vector per working model.")
  }
  mass_diagnostics <- validate_el_masses(
    p_hat = p_hat,
    Delta_hat = Delta_hat,
    mass_tol = mass_tol,
    action = mass_action,
    context = "One-population GMM EL covariance"
  )
  eta_names <- c("tau", names(target$theta_hat))
  if (length(eta_names) != d_eta) eta_names <- colnames(Z)
  beta_names <- unlist(lapply(working, `[[`, "beta_names"), use.names = FALSE)
  denom_target <- 1 + rho_I * Delta_hat
  deriv_eta <- -rho_I * Delta_hat / denom_target^2
  A_eta_eta <- crossprod(
    Z,
    Z * as.numeric(p_hat * (1 + rho_I * Delta_hat) / (1 + rho_I) * deriv_eta)
  )
  dimnames(A_eta_eta) <- list(eta_names, eta_names)
  A_rank_diagnostics <- matrix_rank_diagnostics(A_eta_eta, rank_tol = rank_tol)
  if (!isTRUE(A_rank_diagnostics$full_rank)) {
    stop("GMM EL target bread is rank deficient: rank ", A_rank_diagnostics$rank,
         " of ", A_rank_diagnostics$dimension, ".")
  }
  warn_on_condition_number(A_rank_diagnostics, threshold = condition_warn,
                           context = "GMM EL target bread")
  A_eta_eta_reg <- regularize_negative_bread(A_eta_eta, ridge = ridge)
  A_inv_reg <- solve(A_eta_eta_reg)
  target_pseudo <- make_target_el_pseudo_if(
    target = target,
    Delta_hat = Delta_hat,
    rho_I = rho_I,
    A_eta_eta_reg = A_eta_eta_reg
  )
  s0 <- target_pseudo$s0
  s1 <- target_pseudo$s1
  theta_index_eta <- target_pseudo$theta_index_eta
  IF0_theta <- target_pseudo$IF0_theta
  IF1_theta <- target_pseudo$IF1_theta
  H_theta_theta_shortcut <- symmetrize(
    -A_inv_reg[theta_index_eta, theta_index_eta, drop = FALSE]
  )
  K <- length(working)
  H_theta_beta_shortcut_parts <- vector("list", K)
  IF0_beta_parts <- vector("list", K)
  IF1_beta_parts <- vector("list", K)
  IF0_beta_unresidualized_parts <- vector("list", K)
  IF1_beta_unresidualized_parts <- vector("list", K)
  C_blocks <- vector("list", K)
  D_blocks <- vector("list", K)
  D_rank_diagnostics <- vector("list", K)
  gamma_names <- character(0)
  beta_index_full <- integer(0)
  beta_index_gamma <- integer(0)
  gamma_widths <- vapply(working, function(wk) ncol(wk$D_k), integer(1))
  gamma_offsets <- cumsum(c(0L, head(gamma_widths, -1L)))
  p0n <- as.numeric(p_hat) / mass_diagnostics$mass0
  p1n <- as.numeric(p_hat * Delta_hat) / mass_diagnostics$mass1
  a_target <- colSums(s0 * p0n)
  a1_target <- colSums(s1 * p1n)
  b_list <- vector("list", K)
  a1_gamma_list <- vector("list", K)
  u_beta_list <- vector("list", K)
  v_beta_list <- vector("list", K)
  theoretical_correction_available <- isTRUE(mass_diagnostics$valid)
  A_inv_raw <- if (theoretical_correction_available) {
    tryCatch(solve(A_eta_eta), error = function(e) NULL)
  } else NULL
  if (is.null(A_inv_raw)) theoretical_correction_available <- FALSE

  for (k in seq_along(working)) {
    wk <- working[[k]]
    D_design <- wk$D_k
    gamma_k <- gamma_hat_list[[k]]
    if (length(gamma_k) != ncol(D_design)) stop("gamma_hat_list[[", k, "]] has the wrong length.")
    pieces <- working_score_pieces(D_design, gamma_k, rho = rho_working)
    delta_k <- pieces$delta
    h_scalar <- rho_working * (Delta_hat - delta_k) /
      ((1 + rho_working * delta_k) * (1 + rho_I * Delta_hat))
    h_k <- D_design * as.numeric(h_scalar)
    C_k <- crossprod(
      D_design,
      Z * as.numeric(
        p_hat * rho_working * Delta_hat * (1 + rho_I * delta_k) /
          ((1 + rho_working * delta_k) * (1 + rho_I) * (1 + rho_I * Delta_hat))
      )
    )
    deriv_gamma <- -rho_working * delta_k / (1 + rho_working * delta_k)^2
    D_k_bread <- crossprod(
      D_design,
      D_design * as.numeric(
        p_hat * (1 + rho_working * Delta_hat) / (1 + rho_I) * deriv_gamma
      )
    )
    dimnames(C_k) <- list(wk$gamma_names, eta_names)
    dimnames(D_k_bread) <- list(wk$gamma_names, wk$gamma_names)
    D_diag <- matrix_rank_diagnostics(D_k_bread, rank_tol = rank_tol)
    D_rank_diagnostics[[k]] <- D_diag
    if (!isTRUE(D_diag$full_rank)) {
      stop("GMM EL working bread for model ", k, " is rank deficient: rank ",
           D_diag$rank, " of ", D_diag$dimension, ".")
    }
    warn_on_condition_number(D_diag, threshold = condition_warn,
                             context = paste0("GMM EL working bread for model ", k))
    D_k_reg <- regularize_negative_bread(D_k_bread, ridge = ridge)
    beta_local <- wk$q_alpha_k + seq_len(wk$p_beta_k)
    P_beta_k <- matrix(0, nrow = wk$p_beta_k, ncol = ncol(D_design))
    P_beta_k[cbind(seq_len(wk$p_beta_k), beta_local)] <- 1
    AinvT_Ct <- solve(t(A_eta_eta_reg), t(C_k))
    DinvT_PbetaT <- solve(t(D_k_reg), t(P_beta_k))
    H_eta_beta_shortcut_k <- AinvT_Ct %*% DinvT_PbetaT
    H_theta_beta_shortcut_parts[[k]] <- H_eta_beta_shortcut_k[theta_index_eta, , drop = FALSE]
    T_k <- right_solve(A_eta_eta_reg, C_k)
    R0_k <- h_k - s0 %*% t(T_k)
    R1_k <- h_k - s1 %*% t(T_k)
    L_k <- right_solve(D_k_reg, P_beta_k)
    IF0_beta_parts[[k]] <- -R0_k %*% t(L_k)
    IF1_beta_parts[[k]] <- -R1_k %*% t(L_k)
    IF0_beta_unresidualized_parts[[k]] <- -h_k %*% t(L_k)
    IF1_beta_unresidualized_parts[[k]] <- -h_k %*% t(L_k)
    colnames(IF0_beta_parts[[k]]) <- colnames(IF1_beta_parts[[k]]) <- wk$beta_names
    colnames(IF0_beta_unresidualized_parts[[k]]) <-
      colnames(IF1_beta_unresidualized_parts[[k]]) <- wk$beta_names
    b_k <- colSums(h_k * p0n)
    a1_gamma_k <- colSums(h_k * p1n)
    b_list[[k]] <- b_k
    a1_gamma_list[[k]] <- a1_gamma_k
    if (theoretical_correction_available) {
      D_inv_raw <- tryCatch(solve(D_k_bread), error = function(e) NULL)
      if (is.null(D_inv_raw)) {
        theoretical_correction_available <- FALSE
      } else {
        c_k_raw <- as.numeric(C_k %*% (A_inv_raw %*% a_target))
        u_beta_list[[k]] <- as.numeric((D_inv_raw %*% c_k_raw)[beta_local])
        v_beta_list[[k]] <- as.numeric((D_inv_raw %*% b_k)[beta_local])
      }
    }
    C_blocks[[k]] <- C_k
    D_blocks[[k]] <- D_k_bread
    gamma_names <- c(gamma_names, wk$gamma_names)
    beta_index_gamma <- c(beta_index_gamma, gamma_offsets[k] + beta_local)
    beta_index_full <- c(beta_index_full, d_eta + gamma_offsets[k] + beta_local)
  }
  H_theta_beta_shortcut <- do.call(cbind, H_theta_beta_shortcut_parts)
  IF0_beta <- do.call(cbind, IF0_beta_parts)
  IF1_beta <- do.call(cbind, IF1_beta_parts)
  IF0_beta_unresidualized <- do.call(cbind, IF0_beta_unresidualized_parts)
  IF1_beta_unresidualized <- do.call(cbind, IF1_beta_unresidualized_parts)
  colnames(IF0_beta) <- colnames(IF1_beta) <- beta_names
  colnames(IF0_beta_unresidualized) <- colnames(IF1_beta_unresidualized) <- beta_names
  H_theta_theta <- el_reconstructed_cov(
    A0 = IF0_theta, A1 = IF1_theta,
    p_hat = p_hat, Delta_hat = Delta_hat, rho = rho_I,
    centering = B_centering, compact_tol = compact_tol
  )
  H_theta_beta <- el_reconstructed_crosscov(
    A0 = IF0_theta, A1 = IF1_theta,
    B0 = IF0_beta, B1 = IF1_beta,
    p_hat = p_hat, Delta_hat = Delta_hat, rho = rho_I,
    centering = B_centering, compact_tol = compact_tol
  )
  H_beta_beta <- el_reconstructed_cov(
    A0 = IF0_beta, A1 = IF1_beta,
    p_hat = p_hat, Delta_hat = Delta_hat, rho = rho_I,
    centering = B_centering, compact_tol = compact_tol
  )
  H_beta_beta_unresidualized <- el_reconstructed_cov(
    A0 = IF0_beta_unresidualized, A1 = IF1_beta_unresidualized,
    p_hat = p_hat, Delta_hat = Delta_hat, rho = rho_I,
    centering = B_centering, compact_tol = compact_tol
  )
  correction_beta_max_abs <- NA_real_
  if (theoretical_correction_available &&
      length(u_beta_list) == K && length(v_beta_list) == K &&
      all(vapply(u_beta_list, function(x) !is.null(x), logical(1))) &&
      all(vapply(v_beta_list, function(x) !is.null(x), logical(1)))) {
    correction_beta_max_abs <- 0
    for (k in seq_along(u_beta_list)) {
      for (l in seq_along(u_beta_list)) {
        corr_kl <- (
          tcrossprod(u_beta_list[[k]], v_beta_list[[l]]) +
            tcrossprod(v_beta_list[[k]], u_beta_list[[l]]) -
            tcrossprod(u_beta_list[[k]], u_beta_list[[l]])
        ) / rho_I
        correction_beta_max_abs <- max(correction_beta_max_abs, max(abs(corr_kl)))
      }
    }
  } else theoretical_correction_available <- FALSE
  dimnames(H_theta_theta) <- list(names(target$theta_hat), names(target$theta_hat))
  dimnames(H_theta_beta) <- list(names(target$theta_hat), beta_names)
  dimnames(H_beta_beta) <- list(beta_names, beta_names)
  dimnames(H_theta_theta_shortcut) <- list(names(target$theta_hat), names(target$theta_hat))
  dimnames(H_theta_beta_shortcut) <- list(names(target$theta_hat), beta_names)
  joint_psd_diagnostics <- validate_joint_covariance_blocks(
    H_theta_theta, H_theta_beta, H_beta_beta,
    tol = joint_psd_tol,
    context = "GMM EL joint internal covariance"
  )
  center0 <- c(a_target, unlist(b_list, use.names = FALSE))
  center1 <- c(a1_target, unlist(a1_gamma_list, use.names = FALSE))
  center_resid <- center0 + rho_I * center1
  EL_IF0_theta <- IF0_theta
  EL_IF1_theta <- IF1_theta
  EL_IF0_beta <- IF0_beta
  EL_IF1_beta <- IF1_beta
  if (!isTRUE(return_el_pseudo_if)) {
    EL_IF0_theta <- EL_IF1_theta <- EL_IF0_beta <- EL_IF1_beta <- NULL
  }
  list(
    H_theta_theta = symmetrize(H_theta_theta),
    H_theta_beta = H_theta_beta,
    H_beta_beta = symmetrize(H_beta_beta),
    Omega = NULL,
    J = NULL,
    B = NULL,
    J_blocks = list(A_eta_eta = A_eta_eta, C_gamma_eta = C_blocks, D_gamma_gamma = D_blocks),
    bread_rank_diagnostics = list(A_eta_eta = A_rank_diagnostics,
                                  D_gamma_gamma = D_rank_diagnostics),
    B_centering = paste0(B_centering,
      "_EL_stratum_centering_GMM_target_and_residualized_beta_pseudo_IF"),
    B_centering_requested = B_centering,
    B_first_term_target_working_zero = TRUE,
    H_theta_theta_source = "GMM_EL_target_pseudo_IF_fixed_stratum_covariance",
    H_theta_beta_source = "GMM_EL_target_beta_pseudo_IF_crosscovariance",
    H_beta_beta_source = "GMM_EL_residualized_beta_pseudo_IF_fixed_stratum_covariance",
    H_theta_theta_shortcut = H_theta_theta_shortcut,
    H_theta_beta_shortcut = H_theta_beta_shortcut,
    H_theta_theta_max_abs_diff_from_shortcut = max(abs(H_theta_theta - H_theta_theta_shortcut)),
    H_theta_beta_max_abs_diff_from_shortcut = max(abs(H_theta_beta - H_theta_beta_shortcut)),
    H_theta_theta_sandwich = H_theta_theta,
    H_theta_beta_sandwich = H_theta_beta,
    H_theta_theta_max_abs_diff_from_sandwich = 0,
    H_theta_beta_max_abs_diff_from_sandwich = 0,
    gmm_Hbb_note_formula_correction_beta_max_abs = correction_beta_max_abs,
    gmm_Hbb_theoretical_correction_unregularized_max_abs = correction_beta_max_abs,
    gmm_Hbb_theoretical_correction_available = theoretical_correction_available,
    gmm_Hbb_theoretical_correction_uses_unregularized_breads = TRUE,
    gmm_Hbb_residualization_effect_beta_max_abs =
      max(abs(H_beta_beta - H_beta_beta_unresidualized)),
    H_beta_beta_unresidualized = H_beta_beta_unresidualized,
    H_joint_psd_diagnostics = joint_psd_diagnostics,
    H_joint_min_eigenvalue = joint_psd_diagnostics$minimum_eigenvalue,
    H_joint_max_eigenvalue = joint_psd_diagnostics$maximum_eigenvalue,
    H_joint_psd_cutoff = joint_psd_diagnostics$psd_cutoff,
    H_joint_is_psd = joint_psd_diagnostics$is_psd,
    a0 = center0,
    a1 = center1,
    mu0 = center0,
    mu1 = center1,
    center_resid = center_resid,
    max_abs_center_resid = max(abs(center_resid)),
    EL_mass0 = mass_diagnostics$mass0,
    EL_mass1 = mass_diagnostics$mass1,
    EL_mass_diagnostics = mass_diagnostics,
    theta_index_full = theta_index_eta,
    beta_index_full = beta_index_full,
    beta_index_gamma = beta_index_gamma,
    param_names = c(eta_names, gamma_names),
    beta_names = beta_names,
    rho_I = rho_I,
    rho_working = rho_working,
    method = paste0("retrospective_internal_fixed_stratum_GMM_EL_pseudo_IF_working_ratio_",
                    ifelse(abs(rho_working - rho_I) < 1e-12, "internal", "external")),
    EL_IF0_theta = EL_IF0_theta,
    EL_IF1_theta = EL_IF1_theta,
    EL_IF0_beta = EL_IF0_beta,
    EL_IF1_beta = EL_IF1_beta,
    EL_p_hat = if (isTRUE(return_el_pseudo_if)) p_hat else NULL,
    EL_Delta_hat = if (isTRUE(return_el_pseudo_if)) Delta_hat else NULL,
    EL_pseudo_if_scale = if (isTRUE(return_el_pseudo_if)) {
      "target_and_residualized_GMM_beta_pseudo_IF_uncentered_rows"
    } else NULL,
    IF_theta = NULL,
    IF_beta = NULL
  )
}

estimate_H_onepop_gmm_shortcut_cc <- function(...) {
  estimate_H_onepop_gmm_el_cc(...)
}

estimate_H_onepop_wkm_el_cc <- function(target,
                                        working,
                                        gamma_hat_list,
                                        rho_I,
                                        p_hat,
                                        Delta_hat,
                                        rho_working,
                                        ridge,
                                        B_centering = c("separate", "compact"),
                                        return_el_pseudo_if = FALSE,
                                        mass_tol = ZIM_DEFAULT_EL_MASS_TOL,
                                        mass_action = c("stop", "warn", "none"),
                                        compact_tol = ZIM_DEFAULT_COMPACT_CENTERING_TOL,
                                        rank_tol = ZIM_DEFAULT_BREAD_RANK_TOL,
                                        condition_warn = ZIM_DEFAULT_CONDITION_WARN,
                                        joint_psd_tol = ZIM_DEFAULT_PSD_TOL) {
  ridge <- validate_ridge(ridge)
  B_centering <- match.arg(B_centering)
  mass_action <- match.arg(mass_action)
  mass_tol <- validate_positive_tolerance(mass_tol, "mass_tol")
  compact_tol <- validate_positive_tolerance(compact_tol, "compact_tol")
  rank_tol <- validate_positive_tolerance(rank_tol, "rank_tol")
  condition_warn <- validate_positive_tolerance(condition_warn, "condition_warn")
  joint_psd_tol <- validate_positive_tolerance(joint_psd_tol, "joint_psd_tol")
  if (length(rho_I) != 1L || !is.finite(rho_I) || rho_I <= 0) stop("rho_I must be a positive finite scalar.")
  if (length(rho_working) != 1L || !is.finite(rho_working) || rho_working <= 0) {
    stop("rho_working must be a positive finite scalar.")
  }
  Z <- target$Z
  n <- nrow(Z)
  d_eta <- ncol(Z)
  if (length(Delta_hat) != n || length(p_hat) != n) {
    stop("Delta_hat and p_hat must have one value per internal subject.")
  }
  if (length(working) != length(gamma_hat_list)) {
    stop("gamma_hat_list must have one full retrospective gamma vector per working model.")
  }
  mass_diagnostics <- validate_el_masses(
    p_hat = p_hat,
    Delta_hat = Delta_hat,
    mass_tol = mass_tol,
    action = mass_action,
    context = "One-population WKM EL covariance"
  )
  eta_names <- c("tau", names(target$theta_hat))
  if (length(eta_names) != d_eta) eta_names <- colnames(Z)
  beta_names <- unlist(lapply(working, `[[`, "beta_names"), use.names = FALSE)
  denom_target <- 1 + rho_I * Delta_hat
  deriv_eta <- -rho_I * Delta_hat / denom_target^2
  A_eta_eta <- crossprod(
    Z,
    Z * as.numeric(p_hat * (1 + rho_I * Delta_hat) / (1 + rho_I) * deriv_eta)
  )
  dimnames(A_eta_eta) <- list(eta_names, eta_names)
  A_rank_diagnostics <- matrix_rank_diagnostics(A_eta_eta, rank_tol = rank_tol)
  if (!isTRUE(A_rank_diagnostics$full_rank)) {
    stop("WKM EL target bread is rank deficient: rank ", A_rank_diagnostics$rank,
         " of ", A_rank_diagnostics$dimension, ".")
  }
  warn_on_condition_number(A_rank_diagnostics, threshold = condition_warn,
                           context = "WKM EL target bread")
  A_eta_eta_reg <- regularize_negative_bread(A_eta_eta, ridge = ridge)
  A_inv_reg <- solve(A_eta_eta_reg)
  target_pseudo <- make_target_el_pseudo_if(
    target = target, Delta_hat = Delta_hat, rho_I = rho_I,
    A_eta_eta_reg = A_eta_eta_reg
  )
  s0 <- target_pseudo$s0
  s1 <- target_pseudo$s1
  theta_index_eta <- target_pseudo$theta_index_eta
  IF0_theta <- target_pseudo$IF0_theta
  IF1_theta <- target_pseudo$IF1_theta
  H_theta_theta_shortcut <- symmetrize(
    -A_inv_reg[theta_index_eta, theta_index_eta, drop = FALSE]
  )
  K <- length(working)
  H_theta_beta_shortcut_parts <- vector("list", K)
  IF0_beta_parts <- vector("list", K)
  IF1_beta_parts <- vector("list", K)
  C_blocks <- vector("list", K)
  D_blocks <- vector("list", K)
  D_rank_diagnostics <- vector("list", K)
  gamma_names <- character(0)
  beta_index_full <- integer(0)
  gamma_widths <- vapply(working, function(wk) ncol(wk$D_k), integer(1))
  gamma_offsets <- cumsum(c(0L, head(gamma_widths, -1L)))
  case_multiplier_wkm <- rho_working / rho_I
  p0n <- as.numeric(p_hat) / mass_diagnostics$mass0
  p1n <- as.numeric(p_hat * Delta_hat) / mass_diagnostics$mass1
  a_target <- colSums(s0 * p0n)
  a1_target <- colSums(s1 * p1n)
  b0_list <- vector("list", K)
  b1_list <- vector("list", K)

  for (k in seq_along(working)) {
    wk <- working[[k]]
    D_design <- wk$D_k
    gamma_k <- gamma_hat_list[[k]]
    if (length(gamma_k) != ncol(D_design)) stop("gamma_hat_list[[", k, "]] has the wrong length.")
    pieces <- working_score_pieces(D_design, gamma_k, rho = rho_working)
    delta_k <- pieces$delta
    m0_k <- pieces$phi0
    m1_k <- case_multiplier_wkm * pieces$phi1
    C_k <- crossprod(
      D_design,
      Z * as.numeric(
        p_hat * rho_working * Delta_hat * (1 + rho_I * delta_k) /
          ((1 + rho_working * delta_k) * (1 + rho_I) * (1 + rho_I * Delta_hat))
      )
    )
    deriv_gamma <- -rho_working * delta_k / (1 + rho_working * delta_k)^2
    D_k_bread <- crossprod(
      D_design,
      D_design * as.numeric(
        p_hat * (1 + rho_working * Delta_hat) / (1 + rho_I) * deriv_gamma
      )
    )
    dimnames(C_k) <- list(wk$gamma_names, eta_names)
    dimnames(D_k_bread) <- list(wk$gamma_names, wk$gamma_names)
    D_diag <- matrix_rank_diagnostics(D_k_bread, rank_tol = rank_tol)
    D_rank_diagnostics[[k]] <- D_diag
    if (!isTRUE(D_diag$full_rank)) {
      stop("WKM EL working bread for model ", k, " is rank deficient: rank ",
           D_diag$rank, " of ", D_diag$dimension, ".")
    }
    warn_on_condition_number(D_diag, threshold = condition_warn,
                             context = paste0("WKM EL working bread for model ", k))
    D_k_reg <- regularize_negative_bread(D_k_bread, ridge = ridge)
    beta_local <- wk$q_alpha_k + seq_len(wk$p_beta_k)
    P_beta_k <- matrix(0, nrow = wk$p_beta_k, ncol = ncol(D_design))
    P_beta_k[cbind(seq_len(wk$p_beta_k), beta_local)] <- 1
    AinvT_Ct <- solve(t(A_eta_eta_reg), t(C_k))
    DinvT_PbetaT <- solve(t(D_k_reg), t(P_beta_k))
    H_eta_beta_shortcut_k <- AinvT_Ct %*% DinvT_PbetaT
    H_theta_beta_shortcut_parts[[k]] <- H_eta_beta_shortcut_k[theta_index_eta, , drop = FALSE]
    L_k <- right_solve(D_k_reg, P_beta_k)
    IF0_beta_parts[[k]] <- -m0_k %*% t(L_k)
    IF1_beta_parts[[k]] <- -m1_k %*% t(L_k)
    colnames(IF0_beta_parts[[k]]) <- colnames(IF1_beta_parts[[k]]) <- wk$beta_names
    b0_list[[k]] <- colSums(m0_k * p0n)
    b1_list[[k]] <- colSums(m1_k * p1n)
    C_blocks[[k]] <- C_k
    D_blocks[[k]] <- D_k_bread
    gamma_names <- c(gamma_names, wk$gamma_names)
    beta_index_full <- c(beta_index_full, d_eta + gamma_offsets[k] + beta_local)
  }
  H_theta_beta_shortcut <- do.call(cbind, H_theta_beta_shortcut_parts)
  IF0_beta <- do.call(cbind, IF0_beta_parts)
  IF1_beta <- do.call(cbind, IF1_beta_parts)
  colnames(IF0_beta) <- colnames(IF1_beta) <- beta_names
  H_theta_theta <- el_reconstructed_cov(
    A0 = IF0_theta, A1 = IF1_theta,
    p_hat = p_hat, Delta_hat = Delta_hat, rho = rho_I,
    centering = B_centering, compact_tol = compact_tol
  )
  H_theta_beta <- el_reconstructed_crosscov(
    A0 = IF0_theta, A1 = IF1_theta,
    B0 = IF0_beta, B1 = IF1_beta,
    p_hat = p_hat, Delta_hat = Delta_hat, rho = rho_I,
    centering = B_centering, compact_tol = compact_tol
  )
  H_beta_beta <- el_reconstructed_cov(
    A0 = IF0_beta, A1 = IF1_beta,
    p_hat = p_hat, Delta_hat = Delta_hat, rho = rho_I,
    centering = B_centering, compact_tol = compact_tol
  )
  dimnames(H_theta_theta) <- list(names(target$theta_hat), names(target$theta_hat))
  dimnames(H_theta_beta) <- list(names(target$theta_hat), beta_names)
  dimnames(H_beta_beta) <- list(beta_names, beta_names)
  dimnames(H_theta_theta_shortcut) <- list(names(target$theta_hat), names(target$theta_hat))
  dimnames(H_theta_beta_shortcut) <- list(names(target$theta_hat), beta_names)
  joint_psd_diagnostics <- validate_joint_covariance_blocks(
    H_theta_theta, H_theta_beta, H_beta_beta,
    tol = joint_psd_tol,
    context = "WKM EL joint internal covariance"
  )
  center0 <- c(a_target, unlist(b0_list, use.names = FALSE))
  center1 <- c(a1_target, unlist(b1_list, use.names = FALSE))
  center_resid <- center0 + rho_I * center1
  EL_IF0_theta <- IF0_theta
  EL_IF1_theta <- IF1_theta
  EL_IF0_beta <- IF0_beta
  EL_IF1_beta <- IF1_beta
  if (!isTRUE(return_el_pseudo_if)) {
    EL_IF0_theta <- EL_IF1_theta <- EL_IF0_beta <- EL_IF1_beta <- NULL
  }
  list(
    H_theta_theta = symmetrize(H_theta_theta),
    H_theta_beta = H_theta_beta,
    H_beta_beta = symmetrize(H_beta_beta),
    Omega = NULL,
    J = NULL,
    B = NULL,
    J_blocks = list(A_eta_eta = A_eta_eta, C_gamma_eta = C_blocks, D_gamma_gamma = D_blocks),
    bread_rank_diagnostics = list(A_eta_eta = A_rank_diagnostics,
                                  D_gamma_gamma = D_rank_diagnostics),
    B_centering = paste0(B_centering,
      "_EL_stratum_centering_WKM_target_and_beta_pseudo_IF"),
    B_centering_requested = B_centering,
    B_first_term_target_working_zero = FALSE,
    H_theta_theta_source = "WKM_EL_target_pseudo_IF_fixed_stratum_covariance",
    H_theta_beta_source = "WKM_EL_target_beta_pseudo_IF_crosscovariance",
    H_beta_beta_source = "WKM_EL_beta_pseudo_IF_fixed_stratum_covariance",
    H_theta_theta_shortcut = H_theta_theta_shortcut,
    H_theta_beta_shortcut = H_theta_beta_shortcut,
    H_theta_theta_max_abs_diff_from_shortcut = max(abs(H_theta_theta - H_theta_theta_shortcut)),
    H_theta_beta_max_abs_diff_from_shortcut = max(abs(H_theta_beta - H_theta_beta_shortcut)),
    H_theta_theta_sandwich = H_theta_theta,
    H_theta_beta_sandwich = H_theta_beta,
    H_theta_theta_max_abs_diff_from_sandwich = 0,
    H_theta_beta_max_abs_diff_from_sandwich = 0,
    H_joint_psd_diagnostics = joint_psd_diagnostics,
    H_joint_min_eigenvalue = joint_psd_diagnostics$minimum_eigenvalue,
    H_joint_max_eigenvalue = joint_psd_diagnostics$maximum_eigenvalue,
    H_joint_psd_cutoff = joint_psd_diagnostics$psd_cutoff,
    H_joint_is_psd = joint_psd_diagnostics$is_psd,
    a0 = center0,
    a1 = center1,
    mu0 = center0,
    mu1 = center1,
    center_resid = center_resid,
    max_abs_center_resid = max(abs(center_resid)),
    EL_mass0 = mass_diagnostics$mass0,
    EL_mass1 = mass_diagnostics$mass1,
    EL_mass_diagnostics = mass_diagnostics,
    theta_index_full = theta_index_eta,
    beta_index_full = beta_index_full,
    param_names = c(eta_names, gamma_names),
    beta_names = beta_names,
    rho_I = rho_I,
    rho_working = rho_working,
    method = paste0("retrospective_internal_fixed_stratum_WKM_EL_pseudo_IF_working_ratio_",
                    ifelse(abs(rho_working - rho_I) < 1e-12, "internal", "external")),
    EL_IF0_theta = EL_IF0_theta,
    EL_IF1_theta = EL_IF1_theta,
    EL_IF0_beta = EL_IF0_beta,
    EL_IF1_beta = EL_IF1_beta,
    EL_p_hat = if (isTRUE(return_el_pseudo_if)) p_hat else NULL,
    EL_Delta_hat = if (isTRUE(return_el_pseudo_if)) Delta_hat else NULL,
    EL_pseudo_if_scale = if (isTRUE(return_el_pseudo_if)) {
      "target_and_WKM_beta_pseudo_IF_uncentered_rows"
    } else NULL,
    IF_theta = NULL,
    IF_beta = NULL
  )
}

estimate_H_onepop_wkm_el_shortcut_cc <- function(...) {
  estimate_H_onepop_wkm_el_cc(...)
}

estimate_H_onepop_cc <- function(method = c("GMM", "WKM"),
                                 target,
                                 working,
                                 gamma_hat_list,
                                 rho_I,
                                 p_hat,
                                 Delta_hat,
                                 rho_working = NULL,
                                 ridge,
                                 B_centering = c("separate", "compact"),
                                 return_el_pseudo_if = FALSE,
                                 mass_tol = ZIM_DEFAULT_EL_MASS_TOL,
                                 mass_action = c("stop", "warn", "none"),
                                 compact_tol = ZIM_DEFAULT_COMPACT_CENTERING_TOL,
                                 rank_tol = ZIM_DEFAULT_BREAD_RANK_TOL,
                                 condition_warn = ZIM_DEFAULT_CONDITION_WARN,
                                 joint_psd_tol = ZIM_DEFAULT_PSD_TOL) {
  ridge <- validate_ridge(ridge)
  method <- match.arg(method)
  B_centering <- match.arg(B_centering)
  mass_action <- match.arg(mass_action)
  mass_tol <- validate_positive_tolerance(mass_tol, "mass_tol")
  compact_tol <- validate_positive_tolerance(compact_tol, "compact_tol")
  rank_tol <- validate_positive_tolerance(rank_tol, "rank_tol")
  condition_warn <- validate_positive_tolerance(condition_warn, "condition_warn")
  joint_psd_tol <- validate_positive_tolerance(joint_psd_tol, "joint_psd_tol")
  if (is.null(rho_working)) rho_working <- rho_I
  if (length(rho_I) != 1L || !is.finite(rho_I) || rho_I <= 0) {
    stop("rho_I must be a positive finite scalar.")
  }
  if (length(rho_working) != 1L || !is.finite(rho_working) || rho_working <= 0) {
    stop("rho_working must be a positive finite scalar.")
  }

  if (identical(method, "GMM")) {
    return(estimate_H_onepop_gmm_el_cc(
      target = target,
      working = working,
      gamma_hat_list = gamma_hat_list,
      rho_I = rho_I,
      p_hat = p_hat,
      Delta_hat = Delta_hat,
      rho_working = rho_working,
      ridge = ridge,
      B_centering = B_centering,
      return_el_pseudo_if = return_el_pseudo_if,
      mass_tol = mass_tol,
      mass_action = mass_action,
      compact_tol = compact_tol,
      rank_tol = rank_tol,
      condition_warn = condition_warn,
      joint_psd_tol = joint_psd_tol
    ))
  }

  if (identical(method, "WKM")) {
    return(estimate_H_onepop_wkm_el_cc(
      target = target,
      working = working,
      gamma_hat_list = gamma_hat_list,
      rho_I = rho_I,
      p_hat = p_hat,
      Delta_hat = Delta_hat,
      rho_working = rho_working,
      ridge = ridge,
      B_centering = B_centering,
      return_el_pseudo_if = return_el_pseudo_if,
      mass_tol = mass_tol,
      mass_action = mass_action,
      compact_tol = compact_tol,
      rank_tol = rank_tol,
      condition_warn = condition_warn,
      joint_psd_tol = joint_psd_tol
    ))
  }

  Z <- target$Z
  n <- nrow(Z)
  d_eta <- ncol(Z)
  if (length(Delta_hat) != n || length(p_hat) != n) {
    stop("Delta_hat and p_hat must have one value per internal subject.")
  }
  if (length(working) != length(gamma_hat_list)) {
    stop("gamma_hat_list must have one full retrospective gamma vector per working model.")
  }

  gamma_widths <- vapply(working, function(wk) ncol(wk$D_k), integer(1))
  p_gamma_total <- sum(gamma_widths)
  p_total <- d_eta + p_gamma_total

  eta_names <- c("tau", names(target$theta_hat))
  if (length(eta_names) != d_eta) eta_names <- colnames(Z)
  gamma_names <- unlist(lapply(working, `[[`, "gamma_names"), use.names = FALSE)
  param_names <- c(eta_names, gamma_names)
  beta_names <- unlist(lapply(working, `[[`, "beta_names"), use.names = FALSE)

  m0 <- matrix(0, nrow = n, ncol = p_total)
  m1 <- matrix(0, nrow = n, ncol = p_total)
  J <- matrix(0, nrow = p_total, ncol = p_total)
  colnames(m0) <- colnames(m1) <- param_names
  dimnames(J) <- list(param_names, param_names)

  eta_cols <- seq_len(d_eta)
  denom_target <- 1 + rho_I * Delta_hat
  s0 <- Z * as.numeric(-rho_I * Delta_hat / denom_target)
  s1 <- Z * as.numeric(1 / denom_target)
  m0[, eta_cols] <- s0
  m1[, eta_cols] <- s1

  deriv_eta <- -rho_I * Delta_hat / denom_target^2
  J[eta_cols, eta_cols] <- crossprod(
    Z,
    Z * as.numeric(
      p_hat * (1 + rho_I * Delta_hat) / (1 + rho_I) * deriv_eta
    )
  )

  beta_index_full <- integer(0)
  col_start <- d_eta + 1L
  case_multiplier_wkm <- rho_working / rho_I

  for (k in seq_along(working)) {
    wk <- working[[k]]
    D_k <- wk$D_k
    gamma_k <- gamma_hat_list[[k]]
    if (length(gamma_k) != ncol(D_k)) {
      stop("gamma_hat_list[[", k, "]] has the wrong length.")
    }

    cols <- col_start:(col_start + ncol(D_k) - 1L)
    pieces <- working_score_pieces(D_k, gamma_k, rho = rho_working)

    if (identical(method, "WKM")) {

      m0[, cols] <- pieces$phi0
      m1[, cols] <- case_multiplier_wkm * pieces$phi1
    } else {

      h_scalar <- rho_working * (Delta_hat - pieces$delta) /
        ((1 + rho_working * pieces$delta) * (1 + rho_I * Delta_hat))
      h_k <- D_k * as.numeric(h_scalar)
      m0[, cols] <- h_k
      m1[, cols] <- h_k

      J[cols, eta_cols] <- crossprod(
        D_k,
        Z * as.numeric(
          p_hat * rho_working * Delta_hat * (1 + rho_I * pieces$delta) /
            ((1 + rho_working * pieces$delta) *
               (1 + rho_I) * (1 + rho_I * Delta_hat))
        )
      )
    }

    deriv_gamma <- -rho_working * pieces$delta /
      (1 + rho_working * pieces$delta)^2
    J[cols, cols] <- crossprod(
      D_k,
      D_k * as.numeric(
        p_hat * (1 + rho_working * Delta_hat) / (1 + rho_I) * deriv_gamma
      )
    )

    beta_local <- wk$q_alpha_k + seq_len(wk$p_beta_k)
    beta_index_full <- c(beta_index_full, cols[beta_local])
    col_start <- col_start + ncol(D_k)
  }

  w0 <- p_hat / (1 + rho_I)
  w1 <- p_hat * rho_I * Delta_hat / (1 + rho_I)

  mass0 <- sum(p_hat)
  mass1 <- sum(p_hat * Delta_hat)
  if (!is.finite(mass0) || mass0 <= 0 ||
      !is.finite(mass1) || mass1 <= 0) {
    stop("EL masses used for internal H centering must be positive and finite.")
  }
  a0 <- colSums(m0 * as.numeric(p_hat))
  a1 <- colSums(m1 * as.numeric(p_hat * Delta_hat))
  mu0 <- a0 / mass0
  mu1 <- a1 / mass1
  center_resid <- a0 + rho_I * a1

  if (identical(method, "GMM")) {

    gamma_cols <- if (p_gamma_total > 0L) d_eta + seq_len(p_gamma_total) else integer(0)
    B_uncentered <- matrix(0, nrow = p_total, ncol = p_total)
    dimnames(B_uncentered) <- list(param_names, param_names)

    B_uncentered[eta_cols, eta_cols] <-
      crossprod(m0[, eta_cols, drop = FALSE],
                m0[, eta_cols, drop = FALSE] * as.numeric(w0)) +
      crossprod(m1[, eta_cols, drop = FALSE],
                m1[, eta_cols, drop = FALSE] * as.numeric(w1))

    if (length(gamma_cols) > 0L) {
      h_G <- m0[, gamma_cols, drop = FALSE]
      B_uncentered[gamma_cols, gamma_cols] <-
        crossprod(h_G, h_G * as.numeric(w0 + w1))
    }

    if (identical(B_centering, "separate")) {
      B <- B_uncentered -
        tcrossprod(a0) / (mass0 * (1 + rho_I)) -
        rho_I * tcrossprod(a1) / (mass1 * (1 + rho_I))
      B_centering_used <- "separate_EL_stratum_centering_GMM_zero_target_working_first_term"
    } else {
      B <- B_uncentered - tcrossprod(a0) / rho_I
      B_centering_used <- "compact_fixed_stratum_correction_GMM_zero_target_working_first_term"
    }
    B_first_term_target_working_zero <- TRUE
  } else {
    if (identical(B_centering, "separate")) {
      m0c <- sweep(m0, 2, mu0, "-")
      m1c <- sweep(m1, 2, mu1, "-")
      B <- crossprod(m0c, m0c * as.numeric(w0)) +
        crossprod(m1c, m1c * as.numeric(w1))
      B_centering_used <- "separate_EL_stratum_centering"
    } else {
      B <- crossprod(m0, m0 * as.numeric(w0)) +
        crossprod(m1, m1 * as.numeric(w1)) -
        tcrossprod(a0) / rho_I
      B_centering_used <- "compact_fixed_stratum_correction"
    }
    B_first_term_target_working_zero <- FALSE
  }
  B <- symmetrize(B)
  dimnames(B) <- list(param_names, param_names)

  J_reg <- regularize_jacobian(J, ridge = ridge)
  J_inv <- solve(J_reg)
  Omega <- symmetrize(J_inv %*% B %*% t(J_inv))
  dimnames(Omega) <- list(param_names, param_names)

  theta_index_full <- if (d_eta >= 2L) 2L:d_eta else integer(0)
  if (length(theta_index_full) == 0L) stop("The target model must contain at least one theta parameter.")

  H_theta_theta_sandwich <- Omega[theta_index_full, theta_index_full, drop = FALSE]
  H_theta_beta_sandwich <- Omega[theta_index_full, beta_index_full, drop = FALSE]
  H_theta_theta_source <- "full_Jinv_B_JinvT_sandwich"
  H_theta_beta_source <- "full_Jinv_B_JinvT_sandwich"
  H_theta_theta_max_abs_diff_from_sandwich <- 0
  H_theta_beta_max_abs_diff_from_sandwich <- 0

  if (identical(method, "GMM")) {

    gamma_cols <- if (p_gamma_total > 0L) d_eta + seq_len(p_gamma_total) else integer(0)
    beta_index_gamma <- beta_index_full - d_eta
    A_eta_eta <- J_reg[eta_cols, eta_cols, drop = FALSE]
    C_gamma_eta <- J_reg[gamma_cols, eta_cols, drop = FALSE]
    D_gamma_gamma <- J_reg[gamma_cols, gamma_cols, drop = FALSE]

    A_inv <- solve(A_eta_eta)
    H_theta_theta <- -A_inv[theta_index_full, theta_index_full, drop = FALSE]

    P_beta_t <- matrix(0, nrow = length(gamma_cols), ncol = length(beta_index_gamma))
    P_beta_t[cbind(beta_index_gamma, seq_along(beta_index_gamma))] <- 1
    AinvT_Ct <- solve(t(A_eta_eta), t(C_gamma_eta))
    DinvT_PbetaT <- solve(t(D_gamma_gamma), P_beta_t)
    H_eta_beta <- AinvT_Ct %*% DinvT_PbetaT
    H_theta_beta <- H_eta_beta[theta_index_full, , drop = FALSE]

    H_theta_theta_source <- "GMM_J_only_minus_J_eta_eta_inverse"
    H_theta_beta_source <- "GMM_J_only_Jetaeta_invT_Jgammaeta_T_Jgammagamma_invT"
    H_theta_theta_max_abs_diff_from_sandwich <- max(abs(H_theta_theta - H_theta_theta_sandwich))
    H_theta_beta_max_abs_diff_from_sandwich <- max(abs(H_theta_beta - H_theta_beta_sandwich))
  } else {
    H_theta_theta <- H_theta_theta_sandwich
    H_theta_beta <- H_theta_beta_sandwich
  }

  H_beta_beta <- Omega[beta_index_full, beta_index_full, drop = FALSE]

  dimnames(H_theta_theta) <- list(names(target$theta_hat), names(target$theta_hat))
  dimnames(H_theta_beta) <- list(names(target$theta_hat), beta_names)
  dimnames(H_beta_beta) <- list(beta_names, beta_names)
  dimnames(H_theta_theta_sandwich) <- list(names(target$theta_hat), names(target$theta_hat))
  dimnames(H_theta_beta_sandwich) <- list(names(target$theta_hat), beta_names)

  EL_IF0_theta <- EL_IF1_theta <- EL_IF0_beta <- EL_IF1_beta <- NULL
  if (isTRUE(return_el_pseudo_if) && identical(method, "GMM")) {
    selected_cols <- c(theta_index_full, beta_index_full)
    EL_IF0_selected <- m0 %*% t(J_inv[selected_cols, , drop = FALSE])
    EL_IF1_selected <- m1 %*% t(J_inv[selected_cols, , drop = FALSE])
    d_theta <- length(theta_index_full)
    theta_sel <- seq_len(d_theta)
    beta_sel <- d_theta + seq_along(beta_index_full)
    EL_IF0_theta <- EL_IF0_selected[, theta_sel, drop = FALSE]
    EL_IF1_theta <- EL_IF1_selected[, theta_sel, drop = FALSE]
    EL_IF0_beta <- EL_IF0_selected[, beta_sel, drop = FALSE]
    EL_IF1_beta <- EL_IF1_selected[, beta_sel, drop = FALSE]
    colnames(EL_IF0_theta) <- colnames(EL_IF1_theta) <- names(target$theta_hat)
    colnames(EL_IF0_beta) <- colnames(EL_IF1_beta) <- beta_names
  }

  list(
    H_theta_theta = symmetrize(H_theta_theta),
    H_theta_beta = H_theta_beta,
    H_beta_beta = symmetrize(H_beta_beta),
    Omega = Omega,
    J = J,
    B = B,
    B_centering = B_centering_used,
    B_centering_requested = B_centering,
    B_first_term_target_working_zero = B_first_term_target_working_zero,
    H_theta_theta_source = H_theta_theta_source,
    H_theta_beta_source = H_theta_beta_source,
    H_theta_theta_sandwich = H_theta_theta_sandwich,
    H_theta_beta_sandwich = H_theta_beta_sandwich,
    H_theta_theta_max_abs_diff_from_sandwich = H_theta_theta_max_abs_diff_from_sandwich,
    H_theta_beta_max_abs_diff_from_sandwich = H_theta_beta_max_abs_diff_from_sandwich,
    a0 = a0,
    a1 = a1,
    mu0 = mu0,
    mu1 = mu1,
    center_resid = center_resid,
    max_abs_center_resid = max(abs(center_resid)),
    EL_mass0 = mass0,
    EL_mass1 = mass1,
    theta_index_full = theta_index_full,
    beta_index_full = beta_index_full,
    param_names = param_names,
    beta_names = beta_names,
    rho_I = rho_I,
    rho_working = rho_working,
    method = paste0("retrospective_internal_fixed_stratum_", method,
                    "_working_ratio_", ifelse(abs(rho_working - rho_I) < 1e-12,
                                                "internal", "external")),
    EL_IF0_theta = EL_IF0_theta,
    EL_IF1_theta = EL_IF1_theta,
    EL_IF0_beta = EL_IF0_beta,
    EL_IF1_beta = EL_IF1_beta,
    EL_p_hat = if (isTRUE(return_el_pseudo_if) && identical(method, "GMM")) p_hat else NULL,
    EL_Delta_hat = if (isTRUE(return_el_pseudo_if) && identical(method, "GMM")) Delta_hat else NULL,
    EL_pseudo_if_scale = if (isTRUE(return_el_pseudo_if) && identical(method, "GMM")) "m_times_Jinv_transpose_uncentered" else NULL,
    IF_theta = NULL,
    IF_beta = NULL
  )
}

pack_working_models <- function(working_models, data) {
  if (!is.list(working_models) || length(working_models) == 0) {
    stop("working_models must be a non-empty list.")
  }

  lapply(seq_along(working_models), function(k) {
    wk <- working_models[[k]]
    if (inherits(wk, "formula")) {
      wk <- list(alpha = ~ 1, beta = wk)
    }
    if (is.null(wk$alpha)) wk$alpha <- ~ 1
    if (is.null(wk$beta)) stop("Each working model needs a beta formula.")

    A_k <- rhs_model_matrix(wk$alpha, data, keep_intercept = TRUE)
    B_k <- rhs_model_matrix(wk$beta, data, keep_intercept = FALSE)
    if (ncol(B_k) == 0) stop("Each working beta formula must contain at least one term.")

    colnames(A_k) <- paste0("alpha", k, ":", colnames(A_k))
    colnames(B_k) <- paste0("beta", k, ":", colnames(B_k))
    D_k <- cbind(A_k, B_k)

    list(
      alpha_formula = wk$alpha,
      beta_formula = wk$beta,
      A_k = A_k,
      B_k = B_k,
      D_k = D_k,
      q_alpha_k = ncol(A_k),
      p_beta_k = ncol(B_k),
      gamma_names = colnames(D_k),
      beta_names = colnames(B_k)
    )
  })
}

fit_profile_logistic <- function(Y, X_no_intercept, rho,
                                 control = stats::glm.control(
                                   epsilon = 1e-10,
                                   maxit = 100,
                                   trace = FALSE
                                 ),
                                 score_tol = 1e-8,
                                 rank_tol = ZIM_DEFAULT_BREAD_RANK_TOL,
                                 condition_warn = ZIM_DEFAULT_CONDITION_WARN) {
  Y <- check_binary_y(Y)
  rank_tol <- validate_positive_tolerance(rank_tol, "rank_tol")
  condition_warn <- validate_positive_tolerance(condition_warn, "condition_warn")
  if (length(rho) != 1L || !is.finite(rho) || rho <= 0) {
    stop("rho must be a positive finite scalar.")
  }
  if (length(score_tol) != 1L || !is.finite(score_tol) || score_tol <= 0) {
    stop("score_tol must be a positive finite scalar.")
  }
  Z <- cbind("(Intercept)" = 1, X_no_intercept)
  fit <- stats::glm.fit(
    x = Z,
    y = Y,
    family = stats::binomial(),
    control = control
  )
  gamma_sample <- stats::coef(fit)
  if (!isTRUE(fit$converged)) stop("Target logistic fit did not converge.")
  if (anyNA(gamma_sample) || any(!is.finite(gamma_sample))) {
    stop("Target logistic fit produced non-finite or aliased coefficients.")
  }

  tau_hat <- unname(gamma_sample[1] - log(rho))
  theta_hat <- unname(gamma_sample[-1])
  names(theta_hat) <- colnames(X_no_intercept)

  eta_sample <- drop(Z %*% gamma_sample)
  q_hat <- stats::binomial()$linkinv(eta_sample)
  log_Delta_hat <- tau_hat + drop(X_no_intercept %*% theta_hat)
  target_exp_diagnostic <- check_exp_range(
    eta = log_Delta_hat,
    context = "Final target density-ratio model",
    action = "stop"
  )
  Delta_hat <- exp(log_Delta_hat)
  if (any(!is.finite(Delta_hat)) || any(Delta_hat <= 0)) {
    stop("Target density-ratio estimates are not finite and positive.")
  }

  target_score <- colSums(Z * as.numeric(Y - q_hat))
  target_score_scaled <- target_score / length(Y)
  target_score_max_abs_scaled <- max(abs(target_score_scaled))
  if (!is.finite(target_score_max_abs_scaled) ||
      target_score_max_abs_scaled > score_tol) {
    warning(
      "Target logistic score is not solved tightly: max absolute mean score = ",
      signif(target_score_max_abs_scaled, 6), "."
    )
  }

  information <- crossprod(Z, Z * as.numeric(q_hat * (1 - q_hat))) / length(Y)
  information_diagnostics <- matrix_rank_diagnostics(information, rank_tol = rank_tol)
  if (!isTRUE(information_diagnostics$full_rank)) {
    stop(
      "Target logistic information is rank deficient: rank ",
      information_diagnostics$rank, " of ", information_diagnostics$dimension, "."
    )
  }
  warn_on_condition_number(
    information_diagnostics,
    threshold = condition_warn,
    context = "Target logistic information"
  )

  q_from_delta <- rho * Delta_hat / (1 + rho * Delta_hat)
  max_abs_probability_identity <- max(abs(q_hat - q_from_delta))
  if (!is.finite(max_abs_probability_identity) ||
      max_abs_probability_identity > 1e-8) {
    warning(
      "Target probability/density-ratio identity is not solved tightly: max difference = ",
      signif(max_abs_probability_identity, 6), "."
    )
  }

  list(
    fit = fit,
    Z = Z,
    tau_hat = tau_hat,
    theta_hat = theta_hat,
    q_hat = q_hat,
    Delta_hat = Delta_hat,
    diagnostics = list(
      target_score = target_score,
      target_score_scaled = target_score_scaled,
      target_score_max_abs_scaled = target_score_max_abs_scaled,
      score_tol = score_tol,
      information_rank = information_diagnostics$rank,
      information_dimension = information_diagnostics$dimension,
      information_singular_values = information_diagnostics$singular_values,
      information_rank_cutoff = information_diagnostics$rank_cutoff,
      information_condition_number = information_diagnostics$condition_number,
      information_min_singular_value = information_diagnostics$min_singular_value,
      information_max_singular_value = information_diagnostics$max_singular_value,
      information_rank_diagnostics = information_diagnostics,
      rank_tol = rank_tol,
      condition_warn = condition_warn,
      exp_diagnostic = target_exp_diagnostic,
      max_abs_probability_identity = max_abs_probability_identity,
      glm_iterations = if (!is.null(fit$iter)) fit$iter else NA_integer_
    )
  )
}

target_glm_H_theta_theta <- function(target, ridge) {
  ridge <- validate_ridge(ridge)
  Z <- target$Z
  q_hat <- target$q_hat
  info <- crossprod(Z, Z * as.numeric(q_hat * (1 - q_hat)))
  info_inv <- tryCatch(
    solve(info),
    error = function(e) solve(regularize_matrix(info, ridge = ridge))
  )
  H_eta_eta <- nrow(Z) * symmetrize(info_inv)
  theta_idx <- if (ncol(Z) >= 2L) 2L:ncol(Z) else integer(0)
  if (length(theta_idx) == 0L) stop("The target model must contain at least one theta parameter.")

  H_theta_theta <- H_eta_eta[theta_idx, theta_idx, drop = FALSE]
  dimnames(H_theta_theta) <- list(names(target$theta_hat), names(target$theta_hat))
  symmetrize(H_theta_theta)
}

working_start_from_glm <- function(Y, wk, rho_for_intercept,
                                   beta_tilde_k = NULL,
                                   weights = NULL,
                                   family = stats::binomial(),
                                   control = stats::glm.control(epsilon = 1e-10,
                                                                maxit = 100,
                                                                trace = FALSE)) {
  D_k <- wk$D_k
  if (!is.null(weights)) {
    if (length(weights) != length(Y)) stop("weights must have one value per internal subject.")
    if (any(!is.finite(weights)) || any(weights < 0)) stop("weights must be finite and nonnegative.")
  }

  fit <- if (is.null(weights)) {
    stats::glm.fit(
      x = D_k,
      y = Y,
      family = family,
      control = control
    )
  } else {
    stats::glm.fit(
      x = D_k,
      y = Y,
      weights = weights,
      family = family,
      control = control
    )
  }
  gamma_start <- stats::coef(fit)
  if (anyNA(gamma_start)) gamma_start[is.na(gamma_start)] <- 0

  intercept_idx <- retrospective_intercept_index(colnames(D_k))
  if (length(intercept_idx) > 0L) {
    gamma_start <- to_retrospective_from_glm(
      coef_glm = gamma_start,
      rho_sample = rho_for_intercept,
      intercept_index = intercept_idx
    )
  }

  beta_idx <- wk$q_alpha_k + seq_len(wk$p_beta_k)
  if (!is.null(beta_tilde_k)) gamma_start[beta_idx] <- beta_tilde_k
  names(gamma_start) <- wk$gamma_names
  gamma_start
}

build_gamma_from_alpha_beta <- function(wk, alpha, beta) {
  q_alpha <- wk$q_alpha_k
  p_beta <- wk$p_beta_k
  if (length(alpha) != q_alpha) stop("alpha has the wrong length for this working model.")
  if (length(beta) != p_beta) stop("beta has the wrong length for this working model.")

  gamma <- numeric(ncol(wk$D_k))
  gamma[seq_len(q_alpha)] <- as.numeric(alpha)
  gamma[q_alpha + seq_len(p_beta)] <- as.numeric(beta)
  names(gamma) <- wk$gamma_names
  gamma
}

profile_alpha_gmm_start <- function(wk,
                                    beta_fixed,
                                    alpha_start = NULL,
                                    rho_working,
                                    Delta_hat,
                                    p_hat = NULL,
                                    moment_type = c("onepop", "reference"),
                                    optim_control = list(maxit = 200),
                                    bracket_width = 20,
                                    max_expand = 6L,
                                    alpha_tol = 1e-10) {
  moment_type <- match.arg(moment_type)
  if (length(rho_working) != 1L || !is.finite(rho_working) || rho_working <= 0) {
    stop("rho_working must be a positive finite scalar.")
  }
  if (length(Delta_hat) != nrow(wk$D_k)) {
    stop("Delta_hat must have one value per row of the working-model design.")
  }
  if (any(!is.finite(Delta_hat)) || any(Delta_hat <= 0)) {
    stop("Delta_hat must be finite and positive.")
  }
  if (identical(moment_type, "onepop")) {
    if (is.null(p_hat)) stop("p_hat is required for one-population GMM alpha profiling.")
    if (length(p_hat) != nrow(wk$D_k)) stop("p_hat must have one value per internal subject.")
  }

  q_alpha <- wk$q_alpha_k
  p_beta <- wk$p_beta_k
  alpha_idx <- seq_len(q_alpha)
  beta_idx <- q_alpha + seq_len(p_beta)
  beta_fixed <- as.numeric(beta_fixed)
  if (length(beta_fixed) != p_beta) {
    stop("beta_fixed has the wrong length for this working model.")
  }

  if (q_alpha == 0L) {
    gamma_start <- build_gamma_from_alpha_beta(wk = wk, alpha = numeric(0), beta = beta_fixed)
    return(list(
      gamma_start = gamma_start,
      alpha_start_initial = numeric(0),
      alpha_start_profile = numeric(0),
      beta_start = beta_fixed,
      objective_initial = 0,
      objective_profile = 0,
      max_abs_alpha_moment = 0,
      convergence = 0L,
      method = "no_alpha_terms",
      message = NA_character_,
      moment_type = moment_type,
      rho_working = rho_working
    ))
  }

  if (is.null(alpha_start)) alpha_start <- rep(0, q_alpha)
  alpha_start <- as.numeric(alpha_start)
  if (length(alpha_start) != q_alpha) {
    stop("alpha_start has the wrong length for this working model.")
  }
  if (any(!is.finite(alpha_start))) alpha_start[!is.finite(alpha_start)] <- 0

  alpha_moment <- function(alpha) {
    gamma <- build_gamma_from_alpha_beta(wk = wk, alpha = alpha, beta = beta_fixed)
    if (identical(moment_type, "onepop")) {
      m <- colSums(gmm_moment_onepop_cc(
        gamma_k = gamma,
        wk = wk,
        rho_working = rho_working,
        Delta_hat = Delta_hat,
        p_hat = p_hat
      ))
    } else {
      m <- colMeans(psi2_case_control_reference(
        gamma_k = gamma,
        wk = wk,
        rho_k = rho_working,
        Delta_ref = Delta_hat
      ))
    }
    m[alpha_idx]
  }

  objective_alpha <- function(alpha) {
    m_alpha <- alpha_moment(alpha)
    sum(m_alpha * m_alpha)
  }

  obj_initial <- objective_alpha(alpha_start)
  alpha_profile <- alpha_start
  convergence <- NA_integer_
  method_used <- "optim_BFGS"
  message <- NA_character_

  if (q_alpha == 1L) {

    root_fun <- function(a) as.numeric(alpha_moment(a))
    base <- alpha_start[1]
    bw <- bracket_width
    lo <- base - bw
    hi <- base + bw
    f_lo <- root_fun(lo)
    f_hi <- root_fun(hi)
    expand_iter <- 0L
    while (is.finite(f_lo) && is.finite(f_hi) && f_lo * f_hi > 0 &&
           expand_iter < as.integer(max_expand)) {
      bw <- bw * 2
      lo <- base - bw
      hi <- base + bw
      f_lo <- root_fun(lo)
      f_hi <- root_fun(hi)
      expand_iter <- expand_iter + 1L
    }
    if (is.finite(f_lo) && is.finite(f_hi) && f_lo * f_hi <= 0) {
      alpha_profile <- stats::uniroot(root_fun, lower = lo, upper = hi,
                                      tol = alpha_tol)$root
      convergence <- 0L
      method_used <- "uniroot_alpha_score"
      message <- paste0("bracket_expansions=", expand_iter)
    }
  }

  if (!identical(method_used, "uniroot_alpha_score")) {
    opt <- stats::optim(
      par = alpha_start,
      fn = objective_alpha,
      method = "BFGS",
      control = optim_control
    )
    alpha_profile <- opt$par
    convergence <- opt$convergence
    method_used <- "optim_BFGS_alpha_moment"
    message <- if (!is.null(opt$message)) opt$message else NA_character_
  }

  gamma_start <- build_gamma_from_alpha_beta(
    wk = wk,
    alpha = alpha_profile,
    beta = beta_fixed
  )
  obj_profile <- objective_alpha(alpha_profile)
  max_abs_alpha_moment <- max(abs(alpha_moment(alpha_profile)))

  list(
    gamma_start = gamma_start,
    alpha_start_initial = alpha_start,
    alpha_start_profile = alpha_profile,
    beta_start = beta_fixed,
    objective_initial = obj_initial,
    objective_profile = obj_profile,
    max_abs_alpha_moment = max_abs_alpha_moment,
    convergence = convergence,
    method = method_used,
    message = message,
    moment_type = moment_type,
    rho_working = rho_working
  )
}

make_gmm_start_onepop <- function(Y,
                                  wk,
                                  rho_E,
                                  rho_I,
                                  beta_tilde_k,
                                  Delta_hat,
                                  p_hat,
                                  gmm_start_method = c("profile_alpha", "glm_replace_beta", "glm"),
                                  optim_control = list(maxit = 200)) {
  gmm_start_method <- match.arg(gmm_start_method)
  Y <- check_binary_y(Y)
  case_weight <- rho_E / rho_I
  w_gmm_start <- ifelse(Y == 1, case_weight, 1)

  glm_start <- working_start_from_glm(
    Y = Y,
    wk = wk,
    rho_for_intercept = rho_E,
    beta_tilde_k = NULL,
    weights = w_gmm_start,
    family = stats::quasibinomial()
  )

  if (identical(gmm_start_method, "glm")) {
    return(list(
      gamma_start = glm_start,
      diagnostics = list(
        method = "glm",
        description = "weighted internal glm start; beta not replaced by beta_tilde",
        case_weight = case_weight,
        rho_E = rho_E,
        rho_I = rho_I
      )
    ))
  }

  if (identical(gmm_start_method, "glm_replace_beta")) {
    gamma_start <- glm_start
    beta_idx_local <- wk$q_alpha_k + seq_len(wk$p_beta_k)
    gamma_start[beta_idx_local] <- beta_tilde_k
    names(gamma_start) <- wk$gamma_names
    return(list(
      gamma_start = gamma_start,
      diagnostics = list(
        method = "glm_replace_beta",
        description = "weighted internal glm alpha with beta replaced by beta_tilde",
        case_weight = case_weight,
        rho_E = rho_E,
        rho_I = rho_I
      )
    ))
  }

  alpha_idx <- seq_len(wk$q_alpha_k)
  prof <- profile_alpha_gmm_start(
    wk = wk,
    beta_fixed = beta_tilde_k,
    alpha_start = glm_start[alpha_idx],
    rho_working = rho_E,
    Delta_hat = Delta_hat,
    p_hat = p_hat,
    moment_type = "onepop",
    optim_control = optim_control
  )
  prof$diagnostics <- list(
    method = "profile_alpha",
    description = "weighted internal glm alpha used only as initial value; alpha profiled conditional on beta_tilde using one-population GMM alpha moment",
    profile = prof[setdiff(names(prof), c("gamma_start", "diagnostics"))],
    case_weight = case_weight,
    rho_E = rho_E,
    rho_I = rho_I
  )
  list(gamma_start = prof$gamma_start, diagnostics = prof$diagnostics)
}

make_gmm_start_reference <- function(Y,
                                     wk_internal,
                                     wk_reference,
                                     rho_I,
                                     rho_E,
                                     beta_tilde_k,
                                     Delta_ref,
                                     gmm_start_method_reference = c("profile_alpha", "old"),
                                     optim_control = list(maxit = 200)) {
  gmm_start_method_reference <- match.arg(gmm_start_method_reference)
  Y <- check_binary_y(Y)

  glm_start_internal <- working_start_from_glm(
    Y = Y,
    wk = wk_internal,
    rho_for_intercept = rho_I,
    beta_tilde_k = NULL
  )
  glm_start_reference <- align_gamma_start_to_working(
    gamma_start = glm_start_internal,
    wk = wk_reference,
    what = "internal glm start for reference GMM"
  )

  if (identical(gmm_start_method_reference, "old")) {
    gamma_start <- glm_start_reference
    beta_idx_local <- wk_reference$q_alpha_k + seq_len(wk_reference$p_beta_k)
    gamma_start[beta_idx_local] <- beta_tilde_k
    gamma_start <- calibrate_external_working_intercept_start(
      wk_reference = wk_reference,
      gamma_start = gamma_start,
      Delta_ref = Delta_ref,
      rho_E = rho_E
    )
    names(gamma_start) <- wk_reference$gamma_names
    return(list(
      gamma_start = gamma_start,
      diagnostics = list(
        method = "old",
        description = "old reference start: internal glm alpha, beta replaced by beta_tilde, then external intercept recalibrated only",
        rho_E = rho_E,
        rho_I = rho_I
      )
    ))
  }

  alpha_idx <- seq_len(wk_reference$q_alpha_k)
  prof <- profile_alpha_gmm_start(
    wk = wk_reference,
    beta_fixed = beta_tilde_k,
    alpha_start = glm_start_reference[alpha_idx],
    rho_working = rho_E,
    Delta_hat = Delta_ref,
    p_hat = NULL,
    moment_type = "reference",
    optim_control = optim_control
  )
  prof$diagnostics <- list(
    method = "profile_alpha",
    description = "new reference start: internal glm alpha used only as initial value; all reference-population nuisance alpha terms profiled conditional on beta_tilde",
    profile = prof[setdiff(names(prof), c("gamma_start", "diagnostics"))],
    rho_E = rho_E,
    rho_I = rho_I
  )
  list(gamma_start = prof$gamma_start, diagnostics = prof$diagnostics)
}

psi2_case_control <- function(gamma_k, wk, rho_k, rho_1, Delta_hat) {
  D_k <- wk$D_k
  eta_k <- drop(D_k %*% gamma_k)
  delta_k <- safe_exp(eta_k)

  c_k <- rho_k * (Delta_hat - delta_k) / (1 + rho_k * delta_k)
  v_k <- D_k * c_k
  denom <- 1 + rho_1 * Delta_hat
  v_k / denom
}

psi2_case_control_reference <- function(gamma_k, wk, rho_k, Delta_ref) {
  D_k <- wk$D_k
  eta_k <- drop(D_k %*% gamma_k)
  delta_k <- safe_exp(eta_k)
  c_k <- rho_k * (Delta_ref - delta_k) / (1 + rho_k * delta_k)
  D_k * c_k
}

gmm_moment_onepop_cc <- function(gamma_k,
                                wk,
                                rho_working = NULL,
                                Delta_hat = NULL,
                                p_hat = NULL,
                                rho_I = NULL) {

  if (is.null(rho_working)) rho_working <- rho_I
  if (is.null(rho_working) || length(rho_working) != 1L ||
      !is.finite(rho_working) || rho_working <= 0) {
    stop("rho_working must be a positive finite scalar.")
  }
  if (is.null(Delta_hat)) stop("Delta_hat must be supplied.")

  D_k <- wk$D_k
  if (ncol(D_k) != length(gamma_k)) {
    stop("gamma_k length does not match the working-model design matrix.")
  }
  if (length(Delta_hat) != nrow(D_k)) {
    stop("Delta_hat must have one value per internal subject.")
  }

  eta_k <- drop(D_k %*% gamma_k)
  delta_k <- safe_exp(eta_k)

  scalar <- rho_working * (Delta_hat - delta_k) /
    (1 + rho_working * delta_k)
  out <- D_k * as.numeric(scalar)

  if (!is.null(p_hat)) {
    if (length(p_hat) != nrow(D_k)) stop("p_hat must have one value per internal subject.")
    out <- out * as.numeric(p_hat)
  }
  out
}

gmm_moment_sum_onepop <- function(gamma_k, wk, rho_working, Delta_hat, p_hat) {
  colSums(gmm_moment_onepop_cc(
    gamma_k = gamma_k,
    wk = wk,
    rho_working = rho_working,
    Delta_hat = Delta_hat,
    p_hat = p_hat
  ))
}

gmm_jacobian_onepop <- function(gamma_k, wk, rho_working, Delta_hat, p_hat) {
  D_k <- wk$D_k
  if (ncol(D_k) != length(gamma_k)) {
    stop("gamma_k length does not match the working-model design matrix.")
  }
  if (length(Delta_hat) != nrow(D_k) || length(p_hat) != nrow(D_k)) {
    stop("Delta_hat and p_hat must have one value per internal subject.")
  }

  eta_k <- drop(D_k %*% gamma_k)
  delta_k <- safe_exp(eta_k)
  ddelta_k <- safe_exp_derivative(eta_k)
  deriv <- -rho_working * ddelta_k * (1 + rho_working * Delta_hat) /
    (1 + rho_working * delta_k)^2

  crossprod(D_k, D_k * as.numeric(p_hat * deriv))
}

gmm_jacobian_rank_diagnostics <- function(J, rank_tol = 1e-10) {
  J <- as.matrix(J)
  if (nrow(J) != ncol(J)) stop("The GMM Jacobian must be square.")
  if (length(rank_tol) != 1L || !is.finite(rank_tol) || rank_tol <= 0) {
    stop("rank_tol must be a positive finite scalar.")
  }
  if (any(!is.finite(J))) {
    return(list(
      jacobian_rank = 0L,
      jacobian_dimension = ncol(J),
      jacobian_rank_tol = rank_tol,
      jacobian_rank_cutoff = NA_real_,
      jacobian_min_singular_value = NA_real_,
      jacobian_max_singular_value = NA_real_,
      jacobian_condition_number = Inf,
      jacobian_singular_values = rep(NA_real_, ncol(J))
    ))
  }

  singular_values <- svd(J, nu = 0, nv = 0)$d
  max_sv <- if (length(singular_values) > 0L) max(singular_values) else 0
  cutoff <- rank_tol * max(max_sv, .Machine$double.eps)
  jacobian_rank <- sum(singular_values > cutoff)
  min_sv <- if (length(singular_values) > 0L) min(singular_values) else 0
  condition_number <- if (jacobian_rank == ncol(J) && min_sv > 0) {
    max_sv / min_sv
  } else {
    Inf
  }

  list(
    jacobian_rank = jacobian_rank,
    jacobian_dimension = ncol(J),
    jacobian_rank_tol = rank_tol,
    jacobian_rank_cutoff = cutoff,
    jacobian_min_singular_value = min_sv,
    jacobian_max_singular_value = max_sv,
    jacobian_condition_number = condition_number,
    jacobian_singular_values = singular_values
  )
}

validate_gmm_onepop_solution <- function(out,
                                         wk,
                                         rho_working,
                                         Delta_hat,
                                         p_hat,
                                         score_tol = 1e-8,
                                         rank_tol = 1e-10,
                                         stop_on_failure = TRUE) {
  if (is.null(out$gamma_hat)) stop("GMM solver did not return gamma_hat.")
  if (length(score_tol) != 1L || !is.finite(score_tol) || score_tol <= 0) {
    stop("score_tol must be a positive finite scalar.")
  }
  working_eta <- drop(wk$D_k %*% out$gamma_hat)
  exp_diagnostic <- check_exp_range(
    eta = working_eta,
    context = "Final one-population GMM working model",
    action = "none"
  )

  moment <- gmm_moment_sum_onepop(
    gamma_k = out$gamma_hat,
    wk = wk,
    rho_working = rho_working,
    Delta_hat = Delta_hat,
    p_hat = p_hat
  )
  J <- gmm_jacobian_onepop(
    gamma_k = out$gamma_hat,
    wk = wk,
    rho_working = rho_working,
    Delta_hat = Delta_hat,
    p_hat = p_hat
  )
  rank_diag <- gmm_jacobian_rank_diagnostics(J, rank_tol = rank_tol)

  max_abs_moment <- max(abs(moment))
  moment_norm <- sqrt(sum(moment * moment))
  moment_ok <- is.finite(max_abs_moment) && max_abs_moment <= score_tol
  rank_ok <- rank_diag$jacobian_rank == rank_diag$jacobian_dimension
  exp_ok <- !isTRUE(exp_diagnostic$clipping_active)
  valid <- moment_ok && rank_ok && exp_ok && all(is.finite(out$gamma_hat))

  out$moment <- moment
  out$max_abs_moment <- max_abs_moment
  out$moment_norm <- moment_norm
  out$score_tol <- score_tol
  out$jacobian <- J
  out$jacobian_rank <- rank_diag$jacobian_rank
  out$jacobian_dimension <- rank_diag$jacobian_dimension
  out$jacobian_rank_tol <- rank_diag$jacobian_rank_tol
  out$jacobian_rank_cutoff <- rank_diag$jacobian_rank_cutoff
  out$jacobian_min_singular_value <- rank_diag$jacobian_min_singular_value
  out$jacobian_max_singular_value <- rank_diag$jacobian_max_singular_value
  out$jacobian_condition_number <- rank_diag$jacobian_condition_number
  out$jacobian_singular_values <- rank_diag$jacobian_singular_values
  out$exp_diagnostic <- exp_diagnostic
  out$solution_validated <- valid

  if (!valid && isTRUE(stop_on_failure)) {
    problems <- character(0)
    if (!moment_ok) {
      problems <- c(
        problems,
        sprintf("max(abs(moment)) = %.6g exceeds %.6g", max_abs_moment, score_tol)
      )
    }
    if (!rank_ok) {
      problems <- c(
        problems,
        sprintf(
          "Jacobian rank is %d of %d at relative tolerance %.3g",
          rank_diag$jacobian_rank,
          rank_diag$jacobian_dimension,
          rank_tol
        )
      )
    }
    if (!exp_ok) {
      problems <- c(
        problems,
        sprintf(
          "exponential clipping is active (n_clipped = %d, n_nonfinite = %d, max(abs(eta)) = %.6g)",
          exp_diagnostic$n_clipped,
          exp_diagnostic$n_nonfinite,
          exp_diagnostic$max_abs_eta
        )
      )
    }
    if (any(!is.finite(out$gamma_hat))) {
      problems <- c(problems, "gamma_hat contains non-finite values")
    }
    stop("One-population GMM solution rejected: ", paste(problems, collapse = "; "))
  }

  out
}

gmm_moment_mean_reference <- function(gamma_k, wk, rho_k, Delta_ref) {
  colMeans(psi2_case_control_reference(
    gamma_k = gamma_k,
    wk = wk,
    rho_k = rho_k,
    Delta_ref = Delta_ref
  ))
}

gmm_jacobian_reference <- function(gamma_k, wk, rho_k, Delta_ref) {
  D_k <- wk$D_k
  if (ncol(D_k) != length(gamma_k)) {
    stop("gamma_k length does not match the reference working-model design matrix.")
  }
  if (length(Delta_ref) != nrow(D_k)) {
    stop("Delta_ref must have one value per reference control.")
  }

  eta_k <- drop(D_k %*% gamma_k)
  delta_k <- safe_exp(eta_k)
  ddelta_k <- safe_exp_derivative(eta_k)
  deriv <- -rho_k * ddelta_k * (1 + rho_k * Delta_ref) /
    (1 + rho_k * delta_k)^2

  crossprod(D_k, D_k * as.numeric(deriv)) / nrow(D_k)
}

validate_gmm_reference_solution <- function(out,
                                            wk,
                                            rho_k,
                                            Delta_ref,
                                            score_tol = 1e-8,
                                            rank_tol = ZIM_DEFAULT_BREAD_RANK_TOL,
                                            condition_warn = ZIM_DEFAULT_CONDITION_WARN,
                                            stop_on_failure = TRUE) {
  if (is.null(out$gamma_hat)) stop("Reference GMM solver did not return gamma_hat.")
  score_tol <- validate_positive_tolerance(score_tol, "score_tol")
  rank_tol <- validate_positive_tolerance(rank_tol, "rank_tol")
  condition_warn <- validate_positive_tolerance(condition_warn, "condition_warn")

  working_eta <- drop(wk$D_k %*% out$gamma_hat)
  exp_diagnostic <- check_exp_range(
    eta = working_eta,
    context = "Final reference-control GMM working model",
    action = "none"
  )
  moment <- gmm_moment_mean_reference(
    gamma_k = out$gamma_hat,
    wk = wk,
    rho_k = rho_k,
    Delta_ref = Delta_ref
  )
  J <- gmm_jacobian_reference(
    gamma_k = out$gamma_hat,
    wk = wk,
    rho_k = rho_k,
    Delta_ref = Delta_ref
  )
  rank_diag <- matrix_rank_diagnostics(J, rank_tol = rank_tol)
  if (isTRUE(rank_diag$full_rank)) {
    warn_on_condition_number(
      rank_diag,
      threshold = condition_warn,
      context = "Reference-control GMM Jacobian"
    )
  }

  max_abs_moment <- max(abs(moment))
  moment_norm <- sqrt(sum(moment * moment))
  moment_ok <- is.finite(max_abs_moment) && max_abs_moment <= score_tol
  rank_ok <- isTRUE(rank_diag$full_rank)
  exp_ok <- !isTRUE(exp_diagnostic$clipping_active)
  finite_ok <- all(is.finite(out$gamma_hat))
  valid <- moment_ok && rank_ok && exp_ok && finite_ok

  out$moment <- moment
  out$max_abs_moment <- max_abs_moment
  out$moment_norm <- moment_norm
  out$score_tol <- score_tol
  out$jacobian <- J
  out$jacobian_rank <- rank_diag$rank
  out$jacobian_dimension <- rank_diag$dimension
  out$jacobian_rank_tol <- rank_diag$rank_tol
  out$jacobian_rank_cutoff <- rank_diag$rank_cutoff
  out$jacobian_min_singular_value <- rank_diag$min_singular_value
  out$jacobian_max_singular_value <- rank_diag$max_singular_value
  out$jacobian_condition_number <- rank_diag$condition_number
  out$jacobian_singular_values <- rank_diag$singular_values
  out$jacobian_rank_diagnostics <- rank_diag
  out$exp_diagnostic <- exp_diagnostic
  out$solution_validated <- valid

  if (!valid && isTRUE(stop_on_failure)) {
    problems <- character(0)
    if (!moment_ok) {
      problems <- c(problems, sprintf(
        "max(abs(mean moment)) = %.6g exceeds %.6g",
        max_abs_moment, score_tol
      ))
    }
    if (!rank_ok) {
      problems <- c(problems, sprintf(
        "raw Jacobian rank is %d of %d at relative tolerance %.3g",
        rank_diag$rank, rank_diag$dimension, rank_tol
      ))
    }
    if (!exp_ok) {
      problems <- c(problems, sprintf(
        "exponential clipping is active (n_clipped = %d, n_nonfinite = %d, max(abs(eta)) = %.6g)",
        exp_diagnostic$n_clipped,
        exp_diagnostic$n_nonfinite,
        exp_diagnostic$max_abs_eta
      ))
    }
    if (!finite_ok) problems <- c(problems, "gamma_hat contains non-finite values")
    stop("Reference-control GMM solution rejected: ", paste(problems, collapse = "; "))
  }

  out
}

solve_moment_newton <- function(gamma_start,
                                moment_fun,
                                jacobian_fun,
                                ridge,
                                maxit = 100,
                                score_tol = 1e-8,
                                step_tol = 1e-10,
                                max_halving = 30L) {
  ridge <- validate_ridge(ridge)
  gamma <- as.numeric(gamma_start)
  names(gamma) <- names(gamma_start)

  m <- moment_fun(gamma)
  obj <- sum(m * m)
  max_abs <- max(abs(m))
  converged <- is.finite(max_abs) && max_abs <= score_tol
  iter <- 0L
  message <- NA_character_
  step_factor_last <- NA_real_

  while (!converged && iter < as.integer(maxit)) {
    iter <- iter + 1L

    J <- jacobian_fun(gamma)
    J_reg <- regularize_jacobian(J, ridge = ridge)
    step <- tryCatch(
      solve(J_reg, m),
      error = function(e) rep(NA_real_, length(gamma))
    )

    if (any(!is.finite(step))) {
      message <- "Newton linear solve failed."
      break
    }

    step_norm <- sqrt(sum(step * step))
    if (!is.finite(step_norm) || step_norm <= step_tol) {
      message <- "Newton step below step_tol."
      break
    }

    accepted <- FALSE
    step_factor <- 1
    for (h in seq_len(as.integer(max_halving) + 1L)) {
      gamma_new <- gamma - step_factor * step
      m_new <- moment_fun(gamma_new)
      obj_new <- sum(m_new * m_new)

      if (all(is.finite(m_new)) && is.finite(obj_new) && obj_new < obj) {
        gamma <- gamma_new
        m <- m_new
        obj <- obj_new
        max_abs <- max(abs(m_new))
        accepted <- TRUE
        step_factor_last <- step_factor
        break
      }
      step_factor <- step_factor / 2
    }

    if (!accepted) {
      message <- "Newton line search failed to reduce the moment norm."
      break
    }

    converged <- is.finite(max_abs) && max_abs <= score_tol
  }

  if (!converged && is.na(message) && iter >= as.integer(maxit)) {
    message <- "Newton reached maxit before satisfying score_tol."
  }

  list(
    gamma_hat = gamma,
    convergence = if (converged) 0L else 1L,
    objective = obj,
    max_abs_moment = max_abs,
    iterations = iter,
    ridge = ridge,
    score_tol = score_tol,
    step_tol = step_tol,
    newton_step_factor_last = step_factor_last,
    method = "damped_Newton_analytic_Jacobian",
    message = message
  )
}

solve_gamma_gmm_onepop_k_newton <- function(wk,
                                            rho_working = NULL,
                                            Delta_hat = NULL,
                                            p_hat = NULL,
                                            gamma_start = NULL,
                                            optim_control = list(maxit = 1000),
                                            ridge,
                                            score_tol = 1e-8,
                                            step_tol = 1e-10,
                                            max_halving = 30L,
                                            rho_I = NULL) {
  ridge <- validate_ridge(ridge)
  if (is.null(rho_working)) rho_working <- rho_I
  if (is.null(rho_working) || length(rho_working) != 1L ||
      !is.finite(rho_working) || rho_working <= 0) {
    stop("rho_working must be a positive finite scalar.")
  }
  if (is.null(Delta_hat) || is.null(p_hat) || is.null(gamma_start)) {
    stop("Delta_hat, p_hat, and gamma_start must be supplied.")
  }

  out <- solve_moment_newton(
    gamma_start = gamma_start,
    moment_fun = function(gamma) {
      gmm_moment_sum_onepop(
        gamma_k = gamma,
        wk = wk,
        rho_working = rho_working,
        Delta_hat = Delta_hat,
        p_hat = p_hat
      )
    },
    jacobian_fun = function(gamma) {
      gmm_jacobian_onepop(
        gamma_k = gamma,
        wk = wk,
        rho_working = rho_working,
        Delta_hat = Delta_hat,
        p_hat = p_hat
      )
    },
    ridge = ridge,
    maxit = zim_null_default(optim_control$maxit, 100L),
    score_tol = score_tol,
    step_tol = step_tol,
    max_halving = max_halving
  )
  names(out$gamma_hat) <- wk$gamma_names
  out$rho_working <- rho_working
  out
}

solve_gamma_gmm_reference_k_newton <- function(wk,
                                               rho_k,
                                               Delta_ref,
                                               gamma_start,
                                               optim_control = list(maxit = 1000),
                                               ridge,
                                               score_tol = 1e-8,
                                               step_tol = 1e-10,
                                               max_halving = 30L) {
  ridge <- validate_ridge(ridge)
  if (is.null(rho_k) || length(rho_k) != 1L || !is.finite(rho_k) || rho_k <= 0) {
    stop("rho_k must be a positive finite scalar.")
  }
  if (is.null(Delta_ref) || is.null(gamma_start)) {
    stop("Delta_ref and gamma_start must be supplied.")
  }

  out <- solve_moment_newton(
    gamma_start = gamma_start,
    moment_fun = function(gamma) {
      gmm_moment_mean_reference(
        gamma_k = gamma,
        wk = wk,
        rho_k = rho_k,
        Delta_ref = Delta_ref
      )
    },
    jacobian_fun = function(gamma) {
      gmm_jacobian_reference(
        gamma_k = gamma,
        wk = wk,
        rho_k = rho_k,
        Delta_ref = Delta_ref
      )
    },
    ridge = ridge,
    maxit = zim_null_default(optim_control$maxit, 100L),
    score_tol = score_tol,
    step_tol = step_tol,
    max_halving = max_halving
  )
  names(out$gamma_hat) <- wk$gamma_names
  out$rho_working <- rho_k
  out
}

extract_package_gmm_convergence <- function(fit) {
  out <- NA_integer_
  if (!is.null(fit$algoInfo) && !is.null(fit$algoInfo$convergence)) {
    out <- fit$algoInfo$convergence
  } else if (!is.null(fit$convCode)) {
    out <- fit$convCode
  } else if (!is.null(fit$convergence)) {
    out <- fit$convergence
  }
  if (length(out) == 0L) NA_integer_ else as.integer(out[1])
}

solve_gamma_gmm_onepop_k_package_gmm <- function(wk,
                                                 rho_working = NULL,
                                                 Delta_hat = NULL,
                                                 p_hat = NULL,
                                                 gamma_start = NULL,
                                                 optim_control = list(maxit = 1000),
                                                 rho_I = NULL,
                                                 gmm_type = c("twoStep", "cue", "iterative"),
                                                 gmm_wmatrix = c("optimal", "ident")) {
  gmm_type <- match.arg(gmm_type)
  gmm_wmatrix <- match.arg(gmm_wmatrix)
  if (is.null(rho_working)) rho_working <- rho_I
  if (is.null(rho_working) || length(rho_working) != 1L ||
      !is.finite(rho_working) || rho_working <= 0) {
    stop("rho_working must be a positive finite scalar.")
  }
  if (is.null(Delta_hat) || is.null(p_hat) || is.null(gamma_start)) {
    stop("Delta_hat, p_hat, and gamma_start must be supplied.")
  }
  if (!requireNamespace("gmm", quietly = TRUE)) {
    stop("gmm_solver = 'gmm' requires the CRAN package 'gmm'. Install it or use gmm_solver = 'optim' or 'newton'.")
  }

  D_k <- wk$D_k
  p_gamma <- ncol(D_k)
  x <- cbind(D_k, Delta_hat = Delta_hat, p_hat = p_hat)

  gfun <- function(theta, x) {
    Q <- x[, seq_len(p_gamma), drop = FALSE]
    Delta <- x[, p_gamma + 1L]
    p_w <- x[, p_gamma + 2L]
    delta <- safe_exp(drop(Q %*% theta))
    scalar <- rho_working * (Delta - delta) / (1 + rho_working * delta)

    Q * as.numeric(nrow(Q) * p_w * scalar)
  }

  grad_fun <- function(theta, x) {
    Q <- x[, seq_len(p_gamma), drop = FALSE]
    Delta <- x[, p_gamma + 1L]
    p_w <- x[, p_gamma + 2L]
    eta <- drop(Q %*% theta)
    delta <- safe_exp(eta)
    ddelta <- safe_exp_derivative(eta)
    deriv <- -rho_working * ddelta * (1 + rho_working * Delta) /
      (1 + rho_working * delta)^2
    crossprod(Q, Q * as.numeric(p_w * deriv))
  }

  fit <- gmm::gmm(
    g = gfun,
    x = x,
    t0 = gamma_start,
    gradv = grad_fun,
    type = gmm_type,
    wmatrix = gmm_wmatrix,
    vcov = "iid",
    prewhite = 0,
    centeredVcov = TRUE,
    model = FALSE,
    X = FALSE,
    Y = FALSE,
    control = optim_control
  )

  gamma_hat <- as.numeric(stats::coef(fit))
  names(gamma_hat) <- wk$gamma_names
  moment_bar <- gmm_moment_sum_onepop(
    gamma_k = gamma_hat,
    wk = wk,
    rho_working = rho_working,
    Delta_hat = Delta_hat,
    p_hat = p_hat
  )

  list(
    gamma_hat = gamma_hat,
    convergence = extract_package_gmm_convergence(fit),
    objective = sum(moment_bar * moment_bar),
    max_abs_moment = max(abs(moment_bar)),
    iterations = NA_integer_,
    rho_working = rho_working,
    method = paste0("package_gmm_", gmm_type, "_", gmm_wmatrix),
    message = NA_character_,
    gmm_object = fit
  )
}

solve_gamma_gmm_reference_k_package_gmm <- function(wk,
                                                    rho_k,
                                                    Delta_ref,
                                                    gamma_start,
                                                    optim_control = list(maxit = 1000),
                                                    gmm_type = c("twoStep", "cue", "iterative"),
                                                    gmm_wmatrix = c("optimal", "ident")) {
  gmm_type <- match.arg(gmm_type)
  gmm_wmatrix <- match.arg(gmm_wmatrix)
  if (is.null(rho_k) || length(rho_k) != 1L || !is.finite(rho_k) || rho_k <= 0) {
    stop("rho_k must be a positive finite scalar.")
  }
  if (is.null(Delta_ref) || is.null(gamma_start)) {
    stop("Delta_ref and gamma_start must be supplied.")
  }
  if (!requireNamespace("gmm", quietly = TRUE)) {
    stop("gmm_solver = 'gmm' requires the CRAN package 'gmm'. Install it or use gmm_solver = 'optim' or 'newton'.")
  }

  D_k <- wk$D_k
  p_gamma <- ncol(D_k)
  x <- cbind(D_k, Delta_ref = Delta_ref)

  gfun <- function(theta, x) {
    Q <- x[, seq_len(p_gamma), drop = FALSE]
    Delta <- x[, p_gamma + 1L]
    delta <- safe_exp(drop(Q %*% theta))
    scalar <- rho_k * (Delta - delta) / (1 + rho_k * delta)
    Q * as.numeric(scalar)
  }

  grad_fun <- function(theta, x) {
    Q <- x[, seq_len(p_gamma), drop = FALSE]
    Delta <- x[, p_gamma + 1L]
    eta <- drop(Q %*% theta)
    delta <- safe_exp(eta)
    ddelta <- safe_exp_derivative(eta)
    deriv <- -rho_k * ddelta * (1 + rho_k * Delta) /
      (1 + rho_k * delta)^2
    crossprod(Q, Q * as.numeric(deriv)) / nrow(Q)
  }

  fit <- gmm::gmm(
    g = gfun,
    x = x,
    t0 = gamma_start,
    gradv = grad_fun,
    type = gmm_type,
    wmatrix = gmm_wmatrix,
    vcov = "iid",
    prewhite = 0,
    centeredVcov = TRUE,
    model = FALSE,
    X = FALSE,
    Y = FALSE,
    control = optim_control
  )

  gamma_hat <- as.numeric(stats::coef(fit))
  names(gamma_hat) <- wk$gamma_names
  moment_bar <- gmm_moment_mean_reference(
    gamma_k = gamma_hat,
    wk = wk,
    rho_k = rho_k,
    Delta_ref = Delta_ref
  )

  list(
    gamma_hat = gamma_hat,
    convergence = extract_package_gmm_convergence(fit),
    objective = sum(moment_bar * moment_bar),
    max_abs_moment = max(abs(moment_bar)),
    iterations = NA_integer_,
    rho_working = rho_k,
    method = paste0("package_gmm_", gmm_type, "_", gmm_wmatrix),
    message = NA_character_,
    gmm_object = fit
  )
}

solve_gamma_gmm_onepop_k_dispatch <- function(wk,
                                              rho_working = NULL,
                                              Delta_hat = NULL,
                                              p_hat = NULL,
                                              gamma_start = NULL,
                                              optim_control = list(maxit = 1000),
                                              rho_I = NULL,
                                              ridge,
                                              gmm_solver = c("newton", "optim", "gmm"),
                                              score_tol = 1e-8,
                                              rank_tol = 1e-10,
                                              newton_step_tol = 1e-10,
                                              newton_max_halving = 30L,
                                              newton_fallback_to_optim = TRUE,
                                              gmm_type = c("twoStep", "cue", "iterative"),
                                              gmm_wmatrix = c("optimal", "ident")) {
  ridge <- validate_ridge(ridge)
  gmm_solver <- match.arg(gmm_solver)
  gmm_type <- match.arg(gmm_type)
  gmm_wmatrix <- match.arg(gmm_wmatrix)
  if (is.null(rho_working)) rho_working <- rho_I

  if (identical(gmm_solver, "optim")) {
    out <- solve_gamma_gmm_onepop_k(
      wk = wk,
      rho_working = rho_working,
      Delta_hat = Delta_hat,
      p_hat = p_hat,
      gamma_start = gamma_start,
      optim_control = optim_control,
      rho_I = rho_I
    )
    out$method <- paste0("optim_BFGS_", out$method)
  } else if (identical(gmm_solver, "newton")) {
    out <- solve_gamma_gmm_onepop_k_newton(
      wk = wk,
      rho_working = rho_working,
      Delta_hat = Delta_hat,
      p_hat = p_hat,
      gamma_start = gamma_start,
      optim_control = optim_control,
      ridge = ridge,
      score_tol = score_tol,
      step_tol = newton_step_tol,
      max_halving = newton_max_halving,
      rho_I = rho_I
    )
    if (isTRUE(newton_fallback_to_optim) &&
        (!is.finite(out$max_abs_moment) || out$max_abs_moment > score_tol)) {
      fallback_start <- if (all(is.finite(out$gamma_hat))) out$gamma_hat else gamma_start
      out2 <- solve_gamma_gmm_onepop_k(
        wk = wk,
        rho_working = rho_working,
        Delta_hat = Delta_hat,
        p_hat = p_hat,
        gamma_start = fallback_start,
        optim_control = optim_control,
        rho_I = rho_I
      )
      out2$method <- "optim_BFGS_fallback_after_Newton"
      out2$newton_result <- out
      out <- out2
    }
  } else {
    out <- solve_gamma_gmm_onepop_k_package_gmm(
      wk = wk,
      rho_working = rho_working,
      Delta_hat = Delta_hat,
      p_hat = p_hat,
      gamma_start = gamma_start,
      optim_control = optim_control,
      rho_I = rho_I,
      gmm_type = gmm_type,
      gmm_wmatrix = gmm_wmatrix
    )
  }

  validate_gmm_onepop_solution(
    out = out,
    wk = wk,
    rho_working = rho_working,
    Delta_hat = Delta_hat,
    p_hat = p_hat,
    score_tol = score_tol,
    rank_tol = rank_tol,
    stop_on_failure = TRUE
  )
}

solve_gamma_gmm_reference_k_dispatch <- function(wk,
                                                 rho_k,
                                                 Delta_ref,
                                                 gamma_start,
                                                 optim_control = list(maxit = 1000),
                                                 ridge,
                                                 gmm_solver = c("newton", "optim", "gmm"),
                                                 score_tol = 1e-8,
                                                 rank_tol = ZIM_DEFAULT_BREAD_RANK_TOL,
                                                 condition_warn = ZIM_DEFAULT_CONDITION_WARN,
                                                 newton_step_tol = 1e-10,
                                                 newton_max_halving = 30L,
                                                 newton_fallback_to_optim = TRUE,
                                                 gmm_type = c("twoStep", "cue", "iterative"),
                                                 gmm_wmatrix = c("optimal", "ident")) {
  ridge <- validate_ridge(ridge)
  gmm_solver <- match.arg(gmm_solver)
  gmm_type <- match.arg(gmm_type)
  gmm_wmatrix <- match.arg(gmm_wmatrix)

  if (identical(gmm_solver, "optim")) {
    out <- solve_gamma_gmm_reference_k(
      wk = wk,
      rho_k = rho_k,
      Delta_ref = Delta_ref,
      gamma_start = gamma_start,
      optim_control = optim_control
    )
    out$method <- "optim_BFGS_GMM_reference_control_external_ratio_target"
    out$rho_working <- rho_k
  } else if (identical(gmm_solver, "newton")) {
    out <- solve_gamma_gmm_reference_k_newton(
      wk = wk,
      rho_k = rho_k,
      Delta_ref = Delta_ref,
      gamma_start = gamma_start,
      optim_control = optim_control,
      ridge = ridge,
      score_tol = score_tol,
      step_tol = newton_step_tol,
      max_halving = newton_max_halving
    )
    if (isTRUE(newton_fallback_to_optim) &&
        (!is.finite(out$max_abs_moment) || out$max_abs_moment > score_tol)) {
      fallback_start <- if (all(is.finite(out$gamma_hat))) out$gamma_hat else gamma_start
      out2 <- solve_gamma_gmm_reference_k(
        wk = wk,
        rho_k = rho_k,
        Delta_ref = Delta_ref,
        gamma_start = fallback_start,
        optim_control = optim_control
      )
      out2$method <- "optim_BFGS_fallback_after_Newton_reference"
      out2$rho_working <- rho_k
      out2$newton_result <- out
      out <- out2
    }
  } else {
    out <- solve_gamma_gmm_reference_k_package_gmm(
      wk = wk,
      rho_k = rho_k,
      Delta_ref = Delta_ref,
      gamma_start = gamma_start,
      optim_control = optim_control,
      gmm_type = gmm_type,
      gmm_wmatrix = gmm_wmatrix
    )
  }

  validate_gmm_reference_solution(
    out = out,
    wk = wk,
    rho_k = rho_k,
    Delta_ref = Delta_ref,
    score_tol = score_tol,
    rank_tol = rank_tol,
    condition_warn = condition_warn,
    stop_on_failure = TRUE
  )
}

solve_gamma_gmm_onepop_k <- function(wk,
                                     rho_working = NULL,
                                     Delta_hat = NULL,
                                     p_hat = NULL,
                                     gamma_start = NULL,
                                     optim_control = list(maxit = 1000),
                                     rho_I = NULL) {
  if (is.null(rho_working)) rho_working <- rho_I
  if (is.null(rho_working) || length(rho_working) != 1L ||
      !is.finite(rho_working) || rho_working <= 0) {
    stop("rho_working must be a positive finite scalar.")
  }
  if (is.null(Delta_hat) || is.null(p_hat) || is.null(gamma_start)) {
    stop("Delta_hat, p_hat, and gamma_start must be supplied.")
  }

  moment_fun <- function(gamma_k) {
    gmm_moment_sum_onepop(
      gamma_k = gamma_k,
      wk = wk,
      rho_working = rho_working,
      Delta_hat = Delta_hat,
      p_hat = p_hat
    )
  }
  objective <- function(gamma_k) {
    m <- moment_fun(gamma_k)
    sum(m * m)
  }
  gradient <- function(gamma_k) {
    m <- moment_fun(gamma_k)
    J <- gmm_jacobian_onepop(
      gamma_k = gamma_k,
      wk = wk,
      rho_working = rho_working,
      Delta_hat = Delta_hat,
      p_hat = p_hat
    )
    drop(2 * crossprod(J, m))
  }

  opt <- stats::optim(
    par = gamma_start,
    fn = objective,
    gr = gradient,
    method = "BFGS",
    control = optim_control
  )

  gamma_hat <- opt$par
  names(gamma_hat) <- wk$gamma_names
  moment_bar <- colSums(gmm_moment_onepop_cc(
    gamma_k = gamma_hat,
    wk = wk,
    rho_working = rho_working,
    Delta_hat = Delta_hat,
    p_hat = p_hat
  ))

  list(
    gamma_hat = gamma_hat,
    convergence = opt$convergence,
    objective = opt$value,
    max_abs_moment = max(abs(moment_bar)),
    rho_working = rho_working,
    method = "GMM_onepop_external_ratio_target"
  )
}

solve_gamma_gmm_k <- function(wk, rho_k, rho_1, Delta_hat, gamma_start,
                              optim_control = list(maxit = 1000)) {
  objective <- function(gamma_k) {
    m <- colMeans(psi2_case_control(gamma_k, wk, rho_k, rho_1, Delta_hat))
    sum(m * m)
  }

  opt <- stats::optim(
    par = gamma_start,
    fn = objective,
    method = "BFGS",
    control = optim_control
  )

  gamma_hat <- opt$par
  names(gamma_hat) <- wk$gamma_names
  moment_bar <- colMeans(psi2_case_control(gamma_hat, wk, rho_k, rho_1, Delta_hat))

  list(
    gamma_hat = gamma_hat,
    convergence = opt$convergence,
    objective = opt$value,
    max_abs_moment = max(abs(moment_bar))
  )
}

solve_gamma_gmm_reference_k <- function(wk, rho_k, Delta_ref, gamma_start,
                                        optim_control = list(maxit = 1000)) {
  moment_fun <- function(gamma_k) {
    gmm_moment_mean_reference(gamma_k, wk, rho_k, Delta_ref)
  }
  objective <- function(gamma_k) {
    m <- moment_fun(gamma_k)
    sum(m * m)
  }
  gradient <- function(gamma_k) {
    m <- moment_fun(gamma_k)
    J <- gmm_jacobian_reference(gamma_k, wk, rho_k, Delta_ref)
    drop(2 * crossprod(J, m))
  }

  opt <- stats::optim(
    par = gamma_start,
    fn = objective,
    gr = gradient,
    method = "BFGS",
    control = optim_control
  )

  gamma_hat <- opt$par
  names(gamma_hat) <- wk$gamma_names
  moment_bar <- moment_fun(gamma_hat)

  list(
    gamma_hat = gamma_hat,
    convergence = opt$convergence,
    objective = opt$value,
    max_abs_moment = max(abs(moment_bar)),
    iterations = if (!is.null(opt$counts[["function"]])) opt$counts[["function"]] else NA_integer_,
    message = if (!is.null(opt$message)) opt$message else NA_character_
  )
}

make_reference_target_case_control <- function(
    reference_data, M_formula, theta_hat,
    normalization_tol = 1e-10) {
  normalization_tol <- validate_positive_tolerance(
    normalization_tol, "normalization_tol"
  )
  X_ref <- rhs_model_matrix(M_formula, reference_data, keep_intercept = FALSE)
  if (ncol(X_ref) != length(theta_hat)) {
    stop("reference_data target model matrix has the wrong number of columns.")
  }
  theta_names <- names(theta_hat)
  if (!is.null(theta_names) && length(theta_names) == length(theta_hat) &&
      all(!is.na(theta_names)) && all(nzchar(theta_names))) {
    if (all(theta_names %in% colnames(X_ref))) {
      X_ref <- X_ref[, theta_names, drop = FALSE]
    } else if (!identical(colnames(X_ref), theta_names)) {
      stop(
        "reference_data target model matrix columns do not match the internal target columns. ",
        "Use the same factor levels/contrasts and M_formula in both data sets."
      )
    }
  }
  m_ref <- drop(X_ref %*% theta_hat)
  if (any(!is.finite(m_ref))) {
    stop("The external reference target linear predictor is non-finite.")
  }
  log_mean_exp_ref <- log_mean_exp(m_ref)
  tau_ref <- -log_mean_exp_ref
  log_Delta_ref <- m_ref - log_mean_exp_ref
  exp_diagnostic <- check_exp_range(
    eta = log_Delta_ref,
    context = "Final normalized external target tilt",
    action = "stop"
  )
  Delta_ref <- exp(log_Delta_ref)
  normalization_error <- abs(mean(Delta_ref) - 1)
  if (!is.finite(normalization_error) || normalization_error > normalization_tol) {
    stop(
      "External tilt normalization failed: abs(mean(Delta_ref) - 1) = ",
      signif(normalization_error, 8), ", tolerance = ",
      signif(normalization_tol, 8), "."
    )
  }

  list(
    X_theta = X_ref,
    tau_hat = tau_ref,
    theta_hat = theta_hat,
    log_Delta_hat = log_Delta_ref,
    Delta_hat = Delta_ref,
    diagnostics = list(
      log_mean_exp = log_mean_exp_ref,
      normalization_mean = mean(Delta_ref),
      normalization_error = normalization_error,
      normalization_tol = normalization_tol,
      Delta_min = min(Delta_ref),
      Delta_max = max(Delta_ref),
      exp_diagnostic = exp_diagnostic
    )
  )
}

case_control_target_influence <- function(
    Y, target, ridge,
    fixed_case_control = TRUE,
    rank_tol = ZIM_DEFAULT_BREAD_RANK_TOL,
    condition_warn = ZIM_DEFAULT_CONDITION_WARN) {
  ridge <- validate_ridge(ridge)
  rank_tol <- validate_positive_tolerance(rank_tol, "rank_tol")
  condition_warn <- validate_positive_tolerance(condition_warn, "condition_warn")
  n <- length(Y)
  Z <- target$Z
  q_hat <- target$q_hat
  psi1 <- Z * as.numeric(Y - q_hat)
  A_lambda_lambda <- -crossprod(Z, Z * as.numeric(q_hat * (1 - q_hat))) / n
  bread_diagnostics <- matrix_rank_diagnostics(A_lambda_lambda, rank_tol = rank_tol)
  if (!isTRUE(bread_diagnostics$full_rank)) {
    stop(
      "Internal target bread is rank deficient: rank ", bread_diagnostics$rank,
      " of ", bread_diagnostics$dimension, "."
    )
  }
  warn_on_condition_number(
    bread_diagnostics,
    threshold = condition_warn,
    context = "Internal target bread"
  )
  A_lambda_lambda_reg <- regularize_negative_bread(A_lambda_lambda, ridge = ridge)
  bread_regularized_diagnostics <- matrix_rank_diagnostics(
    A_lambda_lambda_reg, rank_tol = rank_tol
  )
  IF_lambda <- -right_solve_transpose(A_lambda_lambda_reg, psi1)
  colnames(IF_lambda) <- colnames(Z)
  IF_theta <- IF_lambda[, -1, drop = FALSE]

  list(
    IF_lambda = IF_lambda,
    IF_theta = IF_theta,
    H_theta_theta = symmetrize(case_control_cov(IF_theta, Y, fixed_case_control)),
    A_lambda_lambda = A_lambda_lambda,
    A_lambda_lambda_reg = A_lambda_lambda_reg,
    bread_rank_diagnostics = bread_diagnostics,
    bread_regularized_rank_diagnostics = bread_regularized_diagnostics
  )
}

influence_internal_gmm <- function(Y, target, working, gamma_hat_list,
                                   rho_1, rho_k, ridge,
                                   fixed_case_control = TRUE,
                                   joint_psd_tol = ZIM_DEFAULT_PSD_TOL) {
  ridge <- validate_ridge(ridge)
  joint_psd_tol <- validate_positive_tolerance(joint_psd_tol, "joint_psd_tol")
  Y <- check_binary_y(Y)
  n <- length(Y)
  Z <- target$Z
  q_hat <- target$q_hat
  Delta_hat <- target$Delta_hat
  if (length(rho_k) == 1L) rho_k <- rep(rho_k, length(working))
  if (length(rho_k) != length(working)) {
    stop("rho_k must have length 1 or length equal to the number of working models.")
  }
  if (any(!is.finite(rho_k)) || any(rho_k <= 0)) {
    stop("rho_k must contain positive finite ratios.")
  }

  psi1 <- Z * as.numeric(Y - q_hat)
  A_lambda_lambda <- -crossprod(Z, Z * as.numeric(q_hat * (1 - q_hat))) / n
  IF_lambda <- -right_solve_transpose(A_lambda_lambda, psi1, ridge = ridge)
  colnames(IF_lambda) <- colnames(Z)

  IF_beta_parts <- vector("list", length(working))
  J_gamma_eta_list <- vector("list", length(working))
  J_gamma_gamma_list <- vector("list", length(working))

  for (k in seq_along(working)) {
    wk <- working[[k]]
    gamma_k <- gamma_hat_list[[k]]
    D_k <- wk$D_k
    eta_k <- drop(D_k %*% gamma_k)
    delta_k <- safe_exp(eta_k)
    denom_k <- 1 + rho_k[k] * delta_k
    denom_target <- 1 + rho_1 * Delta_hat

    c_k <- rho_k[k] * (Delta_hat - delta_k) / denom_k
    v_k <- D_k * c_k
    psi2_k <- v_k / denom_target

    dc_dm_k <- -rho_k[k] * delta_k * (1 + rho_k[k] * Delta_hat) / (denom_k^2)
    A_gamma_gamma <- crossprod(D_k, D_k * as.numeric(dc_dm_k / denom_target)) / n

    dc_dDelta_times_dDelta_dm <- rho_k[k] * Delta_hat / denom_k
    da_target_dm <- -rho_1 * Delta_hat / (denom_target^2)
    dpsi2_dm_target <- dc_dDelta_times_dDelta_dm / denom_target + c_k * da_target_dm
    A_gamma_lambda <- crossprod(D_k, Z * as.numeric(dpsi2_dm_target)) / n

    tmp <- psi2_k + IF_lambda %*% t(A_gamma_lambda)
    IF_gamma <- -right_solve_transpose(A_gamma_gamma, tmp, ridge = ridge)
    colnames(IF_gamma) <- wk$gamma_names

    beta_idx <- wk$q_alpha_k + seq_len(wk$p_beta_k)
    IF_beta_parts[[k]] <- IF_gamma[, beta_idx, drop = FALSE]
    J_gamma_eta_list[[k]] <- A_gamma_lambda
    J_gamma_gamma_list[[k]] <- A_gamma_gamma
  }

  IF_theta <- IF_lambda[, -1, drop = FALSE]
  IF_beta <- do.call(cbind, IF_beta_parts)
  theta_names <- names(target$theta_hat)
  beta_names <- unlist(lapply(working, `[[`, "beta_names"), use.names = FALSE)
  colnames(IF_theta) <- theta_names
  colnames(IF_beta) <- beta_names

  H_theta_theta <- symmetrize(case_control_cov(IF_theta, Y, fixed_case_control))
  H_theta_beta <- case_control_crosscov(IF_theta, IF_beta, Y, fixed_case_control)
  H_beta_beta <- symmetrize(case_control_cov(IF_beta, Y, fixed_case_control))
  dimnames(H_theta_theta) <- list(theta_names, theta_names)
  dimnames(H_theta_beta) <- list(theta_names, beta_names)
  dimnames(H_beta_beta) <- list(beta_names, beta_names)
  joint_psd_diagnostics <- validate_joint_covariance_blocks(
    H_theta_theta, H_theta_beta, H_beta_beta,
    tol = joint_psd_tol,
    context = "GMM sample-IF joint internal covariance"
  )

  list(
    IF_lambda = IF_lambda,
    IF_theta = IF_theta,
    IF_beta = IF_beta,
    H_theta_theta = H_theta_theta,
    H_theta_beta = H_theta_beta,
    H_beta_beta = H_beta_beta,
    H_joint_psd_diagnostics = joint_psd_diagnostics,
    H_joint_min_eigenvalue = joint_psd_diagnostics$minimum_eigenvalue,
    H_joint_max_eigenvalue = joint_psd_diagnostics$maximum_eigenvalue,
    H_joint_psd_cutoff = joint_psd_diagnostics$psd_cutoff,
    H_joint_is_psd = joint_psd_diagnostics$is_psd,
    J_eta_eta_sample = A_lambda_lambda,
    J_gamma_eta_sample = J_gamma_eta_list,
    J_gamma_gamma_sample = J_gamma_gamma_list,
    B_centering = if (fixed_case_control) {
      "observed_case_control_stratum_centering"
    } else {
      "observed_pooled_centering"
    },
    max_abs_center_resid = NA_real_,
    rho_I = rho_1,
    rho_k = rho_k,
    method = "GMM_onepop_external_ratio_target_sample_if",
    H_theta_theta_source = "GMM_sample_IF_case_control_covariance",
    H_theta_beta_source = "GMM_sample_IF_case_control_crosscovariance",
    H_beta_beta_source = "GMM_sample_IF_case_control_covariance"
  )
}

influence_reference_gmm_case_control <- function(
    Y, target_internal, target_reference,
    working_reference, gamma_hat_list, rho_k, ridge,
    fixed_case_control = TRUE,
    rank_tol = ZIM_DEFAULT_BREAD_RANK_TOL,
    condition_warn = ZIM_DEFAULT_CONDITION_WARN,
    joint_psd_tol = ZIM_DEFAULT_PSD_TOL,
    profile_equivalence_tol = 1e-8) {
  ridge <- validate_ridge(ridge)
  rank_tol <- validate_positive_tolerance(rank_tol, "rank_tol")
  condition_warn <- validate_positive_tolerance(condition_warn, "condition_warn")
  joint_psd_tol <- validate_positive_tolerance(joint_psd_tol, "joint_psd_tol")
  profile_equivalence_tol <- validate_positive_tolerance(
    profile_equivalence_tol, "profile_equivalence_tol"
  )
  Y <- check_binary_y(Y)
  n <- length(Y)
  m_ref <- length(target_reference$Delta_hat)
  if (length(rho_k) == 1L) rho_k <- rep(rho_k, length(working_reference))
  if (length(rho_k) != length(working_reference)) {
    stop("rho_k must have length one or the number of reference working models.")
  }

  X_theta_ref <- target_reference$X_theta
  Delta_ref <- target_reference$Delta_hat
  u_tau <- Delta_ref - 1
  b_E <- colMeans(X_theta_ref * as.numeric(Delta_ref))
  h_E <- sweep(X_theta_ref, 2, b_E, "-")
  normalized_derivative_residual <- colMeans(h_E * as.numeric(Delta_ref))

  target_if <- case_control_target_influence(
    Y = Y,
    target = target_internal,
    ridge = ridge,
    fixed_case_control = fixed_case_control,
    rank_tol = rank_tol,
    condition_warn = condition_warn
  )
  IF_theta <- target_if$IF_theta
  H_theta_theta <- target_if$H_theta_theta

  G_parts <- vector("list", length(working_reference))
  IF_beta_reference_parts <- vector("list", length(working_reference))
  profile_objects <- vector("list", length(working_reference))
  D_rank_diagnostics <- vector("list", length(working_reference))
  D_regularized_rank_diagnostics <- vector("list", length(working_reference))
  profile_equivalence_diagnostics <- vector("list", length(working_reference))

  for (k in seq_along(working_reference)) {
    wk <- working_reference[[k]]
    gamma_k <- align_gamma_start_to_working(
      gamma_start = gamma_hat_list[[k]],
      wk = wk,
      what = paste0("gamma_hat_list[[", k, "]]" )
    )
    D_design <- wk$D_k
    eta_k <- drop(D_design %*% gamma_k)
    exp_diag <- check_exp_range(
      eta_k,
      context = paste0("Final reference GMM working model ", k),
      action = "stop"
    )
    delta_k <- exp(eta_k)
    denom_k <- 1 + rho_k[k] * delta_k

    g_k <- D_design * as.numeric(
      rho_k[k] * (Delta_ref - delta_k) / denom_k
    )
    deriv_gamma <- -rho_k[k] * delta_k *
      (1 + rho_k[k] * Delta_ref) / denom_k^2
    D_raw <- crossprod(
      D_design,
      D_design * as.numeric(deriv_gamma)
    ) / m_ref
    D_diag <- matrix_rank_diagnostics(D_raw, rank_tol = rank_tol)
    D_rank_diagnostics[[k]] <- D_diag
    if (!isTRUE(D_diag$full_rank)) {
      stop(
        "Reference working Jacobian D_", k, " is rank deficient: rank ",
        D_diag$rank, " of ", D_diag$dimension, "."
      )
    }
    warn_on_condition_number(
      D_diag,
      threshold = condition_warn,
      context = paste0("Reference working Jacobian D_", k)
    )
    D_reg <- regularize_negative_bread(D_raw, ridge = ridge)
    D_regularized_rank_diagnostics[[k]] <- matrix_rank_diagnostics(
      D_reg, rank_tol = rank_tol
    )

    deriv_tau <- rho_k[k] * Delta_ref / denom_k
    a_tau <- colMeans(D_design * as.numeric(deriv_tau))
    C0 <- crossprod(
      D_design,
      X_theta_ref * as.numeric(deriv_tau)
    ) / m_ref
    C_profile_direct <- crossprod(
      D_design,
      h_E * as.numeric(deriv_tau)
    ) / m_ref
    C_profile_identity <- C0 - tcrossprod(a_tau, b_E)
    C_profile_diff <- max(abs(C_profile_direct - C_profile_identity))
    C_profile <- C_profile_direct
    r_k <- g_k - tcrossprod(u_tau, a_tau)

    p_gamma <- ncol(D_design)
    p_beta_k <- wk$p_beta_k
    beta_local <- wk$q_alpha_k + seq_len(p_beta_k)
    S_k <- matrix(0, nrow = p_beta_k, ncol = p_gamma)
    S_k[cbind(seq_len(p_beta_k), beta_local)] <- 1
    dimnames(S_k) <- list(wk$beta_names, wk$gamma_names)
    L_beta <- right_solve(D_reg, S_k)
    G_k <- -L_beta %*% C_profile
    IF_beta_ref_k <- -r_k %*% t(L_beta)
    colnames(IF_beta_ref_k) <- wk$beta_names
    rownames(G_k) <- wk$beta_names
    colnames(G_k) <- names(target_internal$theta_hat)

    A_aug_reg <- rbind(
      c(1, rep(0, p_gamma)),
      cbind(a_tau, D_reg)
    )
    A_theta_aug <- rbind(b_E, C0)
    S_aug <- cbind(0, S_k)
    psi_aug <- cbind(u_tau, g_k)
    G_joint <- -S_aug %*% solve(A_aug_reg, A_theta_aug)
    IF_beta_ref_joint <- -psi_aug %*% solve(t(A_aug_reg), t(S_aug))
    G_joint_diff <- max(abs(G_joint - G_k))
    IF_joint_diff <- max(abs(IF_beta_ref_joint - IF_beta_ref_k))
    equivalence_ok <- C_profile_diff <= profile_equivalence_tol &&
      G_joint_diff <= profile_equivalence_tol &&
      IF_joint_diff <= profile_equivalence_tol
    if (!isTRUE(equivalence_ok)) {
      stop(
        "Joint/profile reference implementation mismatch for model ", k,
        ": C diff = ", signif(C_profile_diff, 6),
        ", G diff = ", signif(G_joint_diff, 6),
        ", IF diff = ", signif(IF_joint_diff, 6),
        ", tolerance = ", signif(profile_equivalence_tol, 6), "."
      )
    }

    G_parts[[k]] <- G_k
    IF_beta_reference_parts[[k]] <- IF_beta_ref_k
    profile_equivalence_diagnostics[[k]] <- list(
      model = k,
      C_profile_max_abs_diff = C_profile_diff,
      G_joint_profile_max_abs_diff = G_joint_diff,
      IF_joint_profile_max_abs_diff = IF_joint_diff,
      tolerance = profile_equivalence_tol,
      equivalent = equivalence_ok
    )
    profile_objects[[k]] <- list(
      model = k,
      D_raw = D_raw,
      D_reg = D_reg,
      C0 = C0,
      C = C_profile,
      a_tau = a_tau,
      b_E = b_E,
      g = g_k,
      r = r_k,
      u_tau = u_tau,
      G = G_k,
      IF_beta_ref = IF_beta_ref_k,
      S_beta = S_k,
      L_beta = L_beta,
      exp_diagnostic = exp_diag,
      mean_g = colMeans(g_k),
      mean_r = colMeans(r_k)
    )
  }

  G <- do.call(rbind, G_parts)
  IF_beta_reference <- do.call(cbind, IF_beta_reference_parts)
  beta_names <- unlist(lapply(working_reference, `[[`, "beta_names"), use.names = FALSE)
  theta_names <- names(target_internal$theta_hat)
  rownames(G) <- beta_names
  colnames(G) <- theta_names
  colnames(IF_beta_reference) <- beta_names
  IF_beta_internal <- IF_theta %*% t(G)
  colnames(IF_beta_internal) <- beta_names

  H_theta_beta <- H_theta_theta %*% t(G)
  H_beta_beta_internal <- G %*% H_theta_theta %*% t(G)
  IF_beta_reference_centered <- sweep(
    IF_beta_reference, 2, colMeans(IF_beta_reference), "-"
  )
  H_beta_beta_reference <- (n / m_ref) *
    crossprod(IF_beta_reference_centered) / m_ref
  H_beta_beta <- symmetrize(H_beta_beta_internal + H_beta_beta_reference)

  dimnames(H_theta_theta) <- list(theta_names, theta_names)
  dimnames(H_theta_beta) <- list(theta_names, beta_names)
  dimnames(H_beta_beta_internal) <- list(beta_names, beta_names)
  dimnames(H_beta_beta_reference) <- list(beta_names, beta_names)
  dimnames(H_beta_beta) <- list(beta_names, beta_names)

  H_theta_beta_row_check <- case_control_crosscov(
    IF_theta, IF_beta_internal, Y, fixed_case_control
  )
  H_beta_beta_internal_row_check <- case_control_cov(
    IF_beta_internal, Y, fixed_case_control
  )
  identity_diagnostics <- list(
    H_theta_beta_max_abs_diff = max(abs(H_theta_beta - H_theta_beta_row_check)),
    H_beta_beta_internal_max_abs_diff = max(
      abs(H_beta_beta_internal - H_beta_beta_internal_row_check)
    ),
    normalized_tilt_derivative_max_abs = max(abs(normalized_derivative_residual))
  )

  joint_psd_diagnostics <- validate_joint_covariance_blocks(
    H_theta_theta = H_theta_theta,
    H_theta_beta = H_theta_beta,
    H_beta_beta = H_beta_beta,
    tol = joint_psd_tol,
    context = "Reference-control GMM joint internal/reference covariance"
  )

  list(
    IF_theta = IF_theta,
    IF_beta_internal = IF_beta_internal,
    IF_beta_reference = IF_beta_reference,
    G = G,
    G_parts = G_parts,
    H_theta_theta = H_theta_theta,
    H_theta_beta = H_theta_beta,
    H_beta_beta = H_beta_beta,
    H_beta_beta_internal = symmetrize(H_beta_beta_internal),
    H_beta_beta_reference = symmetrize(H_beta_beta_reference),
    H_theta_beta_row_check = H_theta_beta_row_check,
    H_beta_beta_internal_row_check = H_beta_beta_internal_row_check,
    identity_diagnostics = identity_diagnostics,
    IF_beta_reference_col_means = colMeans(IF_beta_reference),
    tau_reference = target_reference$tau_hat,
    b_E = b_E,
    h_E = h_E,
    u_tau = u_tau,
    n_reference = m_ref,
    D_rank_diagnostics = D_rank_diagnostics,
    D_regularized_rank_diagnostics = D_regularized_rank_diagnostics,
    target_bread_rank_diagnostics = target_if$bread_rank_diagnostics,
    target_bread_regularized_rank_diagnostics = target_if$bread_regularized_rank_diagnostics,
    profile_objects = profile_objects,
    profile_equivalence_diagnostics = profile_equivalence_diagnostics,
    profile_equivalence_max_abs = max(vapply(
      profile_equivalence_diagnostics,
      function(x) max(
        x$C_profile_max_abs_diff,
        x$G_joint_profile_max_abs_diff,
        x$IF_joint_profile_max_abs_diff
      ), numeric(1)
    )),
    H_joint_psd_diagnostics = joint_psd_diagnostics,
    H_joint_min_eigenvalue = joint_psd_diagnostics$minimum_eigenvalue,
    H_joint_max_eigenvalue = joint_psd_diagnostics$maximum_eigenvalue,
    H_joint_psd_cutoff = joint_psd_diagnostics$psd_cutoff,
    H_joint_is_psd = joint_psd_diagnostics$is_psd,
    B_centering = "ordinary_reference_control_centering",
    max_abs_center_resid = max(abs(colMeans(IF_beta_reference))),
    method = "two_population_reference_control_gmm_profiled_beta_block_IF"
  )
}

fit_internal_wkm_case_control <- function(Y,
                                          target = NULL,
                                          working,
                                          rho_I = NULL,
                                          rho_1 = NULL,
                                          rho_E = NULL,
                                          ridge,
                                          fixed_case_control = TRUE,
                                          glm_control = stats::glm.control(
                                            epsilon = 1e-10,
                                            maxit = 100,
                                            trace = FALSE
                                          ),
                                          score_tol = 1e-6,
                                          require_intercept = TRUE,
                                          rank_tol = ZIM_DEFAULT_BREAD_RANK_TOL,
                                          condition_warn = ZIM_DEFAULT_CONDITION_WARN,
                                          joint_psd_tol = ZIM_DEFAULT_PSD_TOL) {
  ridge <- validate_ridge(ridge)
  rank_tol <- validate_positive_tolerance(rank_tol, "rank_tol")
  condition_warn <- validate_positive_tolerance(condition_warn, "condition_warn")
  joint_psd_tol <- validate_positive_tolerance(joint_psd_tol, "joint_psd_tol")
  Y <- check_binary_y(Y)
  if (is.null(rho_I)) rho_I <- rho_1
  if (is.null(rho_I) || length(rho_I) != 1L || !is.finite(rho_I) || rho_I <= 0) {
    stop("rho_I must be supplied as a positive finite scalar.")
  }
  if (is.null(rho_E) || length(rho_E) != 1L || !is.finite(rho_E) || rho_E <= 0) {
    stop("rho_E must be supplied as a positive finite scalar for weighted WKM.")
  }
  if (is.null(target)) {
    stop("target must be supplied so weighted WKM can return joint sample IF/H matrices.")
  }

  n <- length(Y)
  n0 <- sum(Y == 0)
  n1 <- sum(Y == 1)
  if (n0 <= 0 || n1 <= 0) stop("Internal data must contain both cases and controls.")

  case_weight <- rho_E / rho_I
  w_wkm <- ifelse(Y == 1, case_weight, 1)

  Z <- target$Z
  q_hat <- target$q_hat
  psi_lambda <- Z * as.numeric(Y - q_hat)
  A_lambda_lambda <- -crossprod(Z, Z * as.numeric(q_hat * (1 - q_hat))) / n
  A_lambda_lambda_reg <- regularize_negative_bread(A_lambda_lambda, ridge = ridge)
  A_lambda_inv <- solve(A_lambda_lambda_reg)
  IF_lambda <- -right_solve_transpose(A_lambda_lambda, psi_lambda, ridge = ridge)
  colnames(IF_lambda) <- colnames(Z)
  IF_theta <- IF_lambda[, -1, drop = FALSE]
  H_theta_theta_shortcut <- -A_lambda_inv[-1, -1, drop = FALSE]

  gamma_hat_list <- vector("list", length(working))
  gamma_pros_list <- vector("list", length(working))
  beta_hat_parts <- vector("list", length(working))
  IF_beta_parts <- vector("list", length(working))
  H_theta_beta_shortcut_parts <- vector("list", length(working))
  J_gamma_eta_sample_list <- vector("list", length(working))
  J_gamma_gamma_sample_list <- vector("list", length(working))
  working_information_diagnostics <- vector("list", length(working))
  convergence <- vector("list", length(working))
  problem_models <- integer(0)

  for (k in seq_along(working)) {
    wk <- working[[k]]
    D_k <- wk$D_k
    intercept_idx <- retrospective_intercept_index(colnames(D_k))
    if (require_intercept && length(intercept_idx) == 0L) {
      stop(
        "Working model ", k, " has no free intercept in its alpha part. ",
        "Weighted WKM for case-control data requires an intercept so that ",
        "the prospective external-ratio intercept can be converted by subtracting log(rho_E)."
      )
    }

    fit_k <- stats::glm.fit(
      x = D_k,
      y = Y,
      weights = w_wkm,
      family = stats::quasibinomial(),
      control = glm_control
    )
    gamma_sample <- stats::coef(fit_k)
    if (!isTRUE(fit_k$converged)) {
      stop("Weighted internal working model ", k, " did not converge.")
    }
    if (anyNA(gamma_sample) || any(!is.finite(gamma_sample))) {
      stop("Weighted internal working model ", k,
           " produced non-finite or aliased coefficients.")
    }
    if (isTRUE(fit_k$boundary)) {
      warning("Weighted internal working model ", k,
              " terminated on the GLM boundary; inspect separation and fitted probabilities.")
    }
    names(gamma_sample) <- wk$gamma_names
    gamma_pros_list[[k]] <- gamma_sample

    gamma_hat <- gamma_sample
    if (length(intercept_idx) > 0L) {
      gamma_hat <- to_retrospective_from_glm(
        coef_glm = gamma_sample,
        rho_sample = rho_E,
        intercept_index = intercept_idx
      )
    }
    names(gamma_hat) <- wk$gamma_names
    working_exp_diagnostic <- check_exp_range(
      eta = drop(D_k %*% gamma_hat),
      context = paste0("Final WKM working model ", k),
      action = "stop"
    )
    gamma_hat_list[[k]] <- gamma_hat

    beta_idx <- wk$q_alpha_k + seq_len(wk$p_beta_k)
    beta_hat_parts[[k]] <- gamma_hat[beta_idx]

    pieces <- working_score_pieces(D_k, gamma_hat, rho = rho_E)
    score_mat <- pieces$phi0 * as.numeric(Y == 0) +
      pieces$phi1 * as.numeric((Y == 1) * case_weight)
    score_sum <- colSums(score_mat)
    score_scaled <- score_sum / n0
    max_abs_scaled <- max(abs(score_scaled))

    eta_sample <- drop(D_k %*% gamma_sample)
    pi_k <- stats::quasibinomial()$linkinv(eta_sample)
    psi_gamma <- D_k * as.numeric(w_wkm * (Y - pi_k))
    A_gamma_gamma <- -crossprod(
      D_k,
      D_k * as.numeric(w_wkm * pi_k * (1 - pi_k))
    ) / n
    working_info_diag <- matrix_rank_diagnostics(A_gamma_gamma, rank_tol = rank_tol)
    working_information_diagnostics[[k]] <- working_info_diag
    if (!isTRUE(working_info_diag$full_rank)) {
      stop("Weighted WKM information for model ", k,
           " is rank deficient: rank ", working_info_diag$rank,
           " of ", working_info_diag$dimension, ".")
    }
    warn_on_condition_number(working_info_diag, threshold = condition_warn,
                             context = paste0("Weighted WKM information for model ", k))

    J_gamma_eta_k <- crossprod(psi_gamma, psi_lambda) / n
    A_gamma_gamma_reg <- regularize_negative_bread(A_gamma_gamma, ridge = ridge)
    P_beta_t <- matrix(0, nrow = ncol(D_k), ncol = length(beta_idx))
    P_beta_t[cbind(beta_idx, seq_along(beta_idx))] <- 1
    H_eta_beta_k <- solve(t(A_lambda_lambda_reg), t(J_gamma_eta_k)) %*%
      solve(t(A_gamma_gamma_reg), P_beta_t)
    H_theta_beta_shortcut_parts[[k]] <- H_eta_beta_k[-1, , drop = FALSE]
    J_gamma_eta_sample_list[[k]] <- J_gamma_eta_k
    J_gamma_gamma_sample_list[[k]] <- A_gamma_gamma

    IF_gamma <- -right_solve_transpose(A_gamma_gamma, psi_gamma, ridge = ridge)
    colnames(IF_gamma) <- wk$gamma_names
    IF_beta_parts[[k]] <- IF_gamma[, beta_idx, drop = FALSE]

    converged <- isTRUE(fit_k$converged)
    if (!is.finite(max_abs_scaled) || max_abs_scaled > score_tol) {
      problem_models <- c(problem_models, k)
    }

    convergence[[k]] <- list(
      converged = converged,
      max_abs_retrospective_score = max(abs(score_sum)),
      mean_abs_retrospective_score = mean(abs(score_sum)),
      max_abs_scaled_retrospective_score = max_abs_scaled,
      mean_abs_scaled_retrospective_score = mean(abs(score_scaled)),
      score_scaled_by = "n_control_internal",
      score_tol = score_tol,
      glm_iterations = if (!is.null(fit_k$iter)) fit_k$iter else NA_integer_,
      glm_boundary = if (!is.null(fit_k$boundary)) fit_k$boundary else NA,
      information_rank = working_info_diag$rank,
      information_dimension = working_info_diag$dimension,
      information_rank_tol = working_info_diag$rank_tol,
      information_rank_cutoff = working_info_diag$rank_cutoff,
      information_min_singular_value = working_info_diag$min_singular_value,
      information_max_singular_value = working_info_diag$max_singular_value,
      information_condition_number = working_info_diag$condition_number,
      information_singular_values = working_info_diag$singular_values,
      case_weight = case_weight,
      rho_I = rho_I,
      rho_target = rho_E,
      exp_diagnostic = working_exp_diagnostic,
      intercept_scale = "weighted glm prospective intercept converted by subtracting log(rho_E)"
    )
  }

  if (length(problem_models) > 0L) {
    shown <- paste(utils::head(problem_models, 20L), collapse = ", ")
    warning(
      "Weighted WKM may not have solved the external-ratio retrospective score for model(s): ",
      shown,
      if (length(problem_models) > 20L) ", ..." else "",
      ". Inspect fit$convergence[[k]]$max_abs_scaled_retrospective_score."
    )
  }

  IF_beta <- do.call(cbind, IF_beta_parts)
  beta_names <- unlist(lapply(working, `[[`, "beta_names"), use.names = FALSE)
  colnames(IF_beta) <- beta_names

  H_theta_theta_sample_if <- symmetrize(case_control_cov(IF_theta, Y, fixed_case_control))
  H_theta_beta_sample_if <- case_control_crosscov(IF_theta, IF_beta, Y, fixed_case_control)
  H_beta_beta <- symmetrize(case_control_cov(IF_beta, Y, fixed_case_control))

  H_theta_theta <- H_theta_theta_sample_if
  H_theta_beta <- H_theta_beta_sample_if
  H_theta_beta_shortcut <- do.call(cbind, H_theta_beta_shortcut_parts)
  joint_psd_diagnostics <- validate_joint_covariance_blocks(
    H_theta_theta, H_theta_beta, H_beta_beta,
    tol = joint_psd_tol,
    context = "WKM sample-IF joint internal covariance"
  )

  dimnames(H_theta_theta) <- list(names(target$theta_hat), names(target$theta_hat))
  dimnames(H_theta_beta) <- list(names(target$theta_hat), beta_names)
  dimnames(H_beta_beta) <- list(beta_names, beta_names)
  dimnames(H_theta_theta_sample_if) <- list(names(target$theta_hat), names(target$theta_hat))
  dimnames(H_theta_beta_sample_if) <- list(names(target$theta_hat), beta_names)

  list(
    beta_hat_WKM = unlist(beta_hat_parts, use.names = TRUE),
    gamma_hat_WKM = gamma_hat_list,
    gamma_pros_WKM = gamma_pros_list,
    convergence = convergence,
    IF_lambda = IF_lambda,
    IF_theta = IF_theta,
    IF_beta = IF_beta,
    H_theta_theta = H_theta_theta,
    H_theta_beta = H_theta_beta,
    H_beta_beta = H_beta_beta,
    H_theta_theta_sample_if = H_theta_theta_sample_if,
    H_theta_beta_sample_if = H_theta_beta_sample_if,
    H_theta_theta_shortcut = symmetrize(H_theta_theta_shortcut),
    H_theta_beta_shortcut = H_theta_beta_shortcut,
    H_theta_theta_source = "WKM_sample_weighted_glm_IF_case_control_covariance",
    H_theta_beta_source = "WKM_sample_weighted_glm_IF_case_control_crosscovariance",
    H_beta_beta_source = "WKM_sample_weighted_glm_IF_case_control_covariance",
    H_theta_theta_max_abs_diff_from_shortcut = max(abs(H_theta_theta - H_theta_theta_shortcut)),
    H_theta_beta_max_abs_diff_from_shortcut = max(abs(H_theta_beta - H_theta_beta_shortcut)),

    H_theta_theta_max_abs_diff_from_sample_if = 0,
    H_theta_beta_max_abs_diff_from_sample_if = 0,
    J_eta_eta_sample = A_lambda_lambda,
    J_gamma_eta_sample = J_gamma_eta_sample_list,
    J_gamma_gamma_sample = J_gamma_gamma_sample_list,
    working_information_diagnostics = working_information_diagnostics,
    H_joint_psd_diagnostics = joint_psd_diagnostics,
    H_joint_min_eigenvalue = joint_psd_diagnostics$minimum_eigenvalue,
    H_joint_max_eigenvalue = joint_psd_diagnostics$maximum_eigenvalue,
    H_joint_psd_cutoff = joint_psd_diagnostics$psd_cutoff,
    H_joint_is_psd = joint_psd_diagnostics$is_psd,
    B_centering = if (fixed_case_control) {
      "observed_case_control_stratum_centering"
    } else {
      "observed_pooled_centering"
    },
    max_abs_center_resid = NA_real_,
    case_weight = case_weight,
    rho_I = rho_I,
    rho_E = rho_E,
    method = "WKM_weighted_glm_external_ratio_target_sample_if"
  )
}

# Fit case-control working models to external individual-level data.
#
# Arguments:
#   external_data: data frame containing the outcome and working-model variables.
#   Y: name of the binary outcome column (0/1 or a two-level factor).
#   working_models: nonempty list of beta formulas or lists with alpha and beta
#     formulas; alpha defaults to ~ 1.
#   n_internal: optional internal sample size for the Sigma_beta_beta output.
#   rho_k: optional case/control ratio, either a scalar or one value per model.
#   covariance_method: 'prospective_centered' or 'formula6'; the latter requires
#     target_delta. See args(fit_external_summary_case_control) for all controls.
#
# Returns a list containing beta_tilde, gamma_tilde, external sample counts,
# V_beta_tilde, Sigma_beta_beta, and fit diagnostics.
fit_external_summary_case_control <- function(external_data, Y,
                                              working_models,
                                              n_internal = NULL,
                                              rho_k = NULL,
                                              ridge = ZIM_DEFAULT_RIDGE,
                                              fixed_case_control = TRUE,
                                              covariance_method = c("prospective_centered", "formula6"),
                                              target_delta = NULL,
                                              rank_tol = ZIM_DEFAULT_BREAD_RANK_TOL,
                                              condition_warn = ZIM_DEFAULT_CONDITION_WARN,
                                              psd_tol = ZIM_DEFAULT_PSD_TOL) {
  ridge <- validate_ridge(ridge)
  rank_tol <- validate_positive_tolerance(rank_tol, "rank_tol")
  condition_warn <- validate_positive_tolerance(condition_warn, "condition_warn")
  psd_tol <- validate_positive_tolerance(psd_tol, "psd_tol")
  covariance_method <- match.arg(covariance_method)
  Y_ext <- check_binary_y(external_data[[Y]])
  N <- length(Y_ext)
  rho_ext <- sum(Y_ext == 1) / sum(Y_ext == 0)

  working_ext <- pack_working_models(working_models, external_data)
  K <- length(working_ext)
  if (is.null(rho_k)) rho_k <- rep(rho_ext, K)
  if (length(rho_k) == 1) rho_k <- rep(rho_k, K)
  if (length(rho_k) != K) stop("rho_k must have length 1 or length K.")

  beta_tilde_parts <- vector("list", K)
  IF_beta_parts <- vector("list", K)
  gamma_tilde_list <- vector("list", K)
  working_information_diagnostics <- vector("list", K)

  for (k in seq_along(working_ext)) {
    wk <- working_ext[[k]]
    D_k <- wk$D_k

    fit_k <- stats::glm.fit(x = D_k, y = Y_ext, family = stats::binomial())
    gamma_sample <- stats::coef(fit_k)
    if (!isTRUE(fit_k$converged)) stop("External working model ", k, " did not converge.")
    if (anyNA(gamma_sample) || any(!is.finite(gamma_sample))) {
      stop("External working model ", k, " produced non-finite or aliased coefficients.")
    }
    if (isTRUE(fit_k$boundary)) {
      warning("External working model ", k,
              " terminated on the GLM boundary; inspect separation and fitted probabilities.")
    }

    gamma_tilde <- gamma_sample
    intercept_col <- grepl(":\\(Intercept\\)$", colnames(D_k))
    if (any(intercept_col)) {
      gamma_tilde[which(intercept_col)[1]] <- gamma_tilde[which(intercept_col)[1]] -
        log(rho_k[k])
    }
    names(gamma_tilde) <- wk$gamma_names
    gamma_tilde_list[[k]] <- gamma_tilde

    beta_idx <- wk$q_alpha_k + seq_len(wk$p_beta_k)
    beta_tilde_parts[[k]] <- gamma_tilde[beta_idx]

    eta_k_sample <- drop(D_k %*% gamma_sample)
    q_k <- stats::binomial()$linkinv(eta_k_sample)
    psi3_k <- D_k * as.numeric(Y_ext - q_k)
    A_gamma_gamma <- -crossprod(D_k, D_k * as.numeric(q_k * (1 - q_k))) / N
    info_diag <- matrix_rank_diagnostics(A_gamma_gamma, rank_tol = rank_tol)
    working_information_diagnostics[[k]] <- info_diag
    if (!isTRUE(info_diag$full_rank)) {
      stop("External working information for model ", k,
           " is rank deficient: rank ", info_diag$rank,
           " of ", info_diag$dimension, ".")
    }
    warn_on_condition_number(info_diag, threshold = condition_warn,
                             context = paste0("External working information for model ", k))
    IF_gamma <- -right_solve_transpose(A_gamma_gamma, psi3_k, ridge = ridge)
    colnames(IF_gamma) <- wk$gamma_names
    IF_beta_parts[[k]] <- IF_gamma[, beta_idx, drop = FALSE]
  }

  beta_tilde <- unlist(beta_tilde_parts, use.names = TRUE)
  IF_beta <- do.call(cbind, IF_beta_parts)

  H_ext_beta_beta <- if (identical(covariance_method, "formula6")) {
    if (is.null(target_delta)) {
      stop("covariance_method='formula6' requires target_delta.")
    }
    formula6_beta_cov_case_control(
      Y = Y_ext,
      working = working_ext,
      gamma_tilde_list = gamma_tilde_list,
      rho_k = rho_k,
      target_delta = target_delta,
      ridge = ridge
    )
  } else {
    symmetrize(case_control_cov(IF_beta, Y_ext, fixed_case_control))
  }
  H_ext_psd_diagnostics <- matrix_psd_diagnostics(H_ext_beta_beta, tol = psd_tol)
  if (!isTRUE(H_ext_psd_diagnostics$is_psd)) {
    stop("External fitted-summary beta covariance is not positive semidefinite: minimum eigenvalue = ",
         signif(H_ext_psd_diagnostics$minimum_eigenvalue, 6), ".")
  }
  V_beta_tilde <- H_ext_beta_beta / N
  Sigma_beta_beta <- NULL
  if (!is.null(n_internal)) Sigma_beta_beta <- (n_internal / N) * H_ext_beta_beta

  list(
    beta_tilde = beta_tilde,
    gamma_tilde = gamma_tilde_list,
    ridge = ridge,
    rho_external = rho_ext,
    external_n_case = sum(Y_ext == 1),
    external_n_control = sum(Y_ext == 0),
    n_external = N,
    rho_k = rho_k,
    covariance_method = covariance_method,
    H_ext_beta_beta = H_ext_beta_beta,
    V_beta_tilde = V_beta_tilde,
    Sigma_beta_beta = Sigma_beta_beta,
    IF_beta_external = IF_beta,
    working_information_diagnostics = working_information_diagnostics,
    H_ext_psd_diagnostics = H_ext_psd_diagnostics
  )
}

# Fit ZIM for a one- or two-population case-control study.
#
# Arguments:
#   internal_data: data frame containing the outcome and target/working variables.
#   Y: name of the binary outcome column (0/1 or a two-level factor).
#   M_formula: one-sided formula for target covariates, for example ~ prs + age.
#   working_models: nonempty list of beta formulas or lists with alpha and beta
#     formulas, for example list(list(alpha = ~ 1, beta = ~ snp1)).
#   beta_tilde: external working-model coefficient estimates, ordered or named
#     to match the beta terms in working_models.
#   Sigma_beta_beta: optional covariance input for one population; use NULL for
#     program-estimated covariance. Must be NULL when reference_data is supplied.
#   external_n_case, external_n_control: external case and control sample sizes;
#     required when Sigma_beta_beta is NULL.
#   rho_k: optional external case/control ratio, scalar or one value per model.
#   reference_data: optional independent external-population control data frame
#     containing target and working covariates; requires beta_method = 'GMM'.
#   W: optional user-supplied weight matrix.
#   ridge: nonnegative numeric value; default ZIM_DEFAULT_RIDGE.
#   fixed_case_control: logical; defaults to TRUE.
#   beta_method: 'GMM' or 'WKM'.
#   inverse_method: 'solve', 'pca', 'hbb_pca', 'sigma_pca', 'm_factor_diag',
#     'hbb_factor_diag', 'shrink', or 'shrinkage'. 'shrink' is unavailable with
#     reference_data; 'shrinkage' requires GMM and shrink_lambda in [0, 1].
#   pca_rule: 'variance', 'condition', 'condition_or_variance', or 'fixed'.
#   pc_var: variance-fraction threshold for applicable PCA rules.
#   pc_max_condition: condition-number threshold for applicable PCA rules.
#   n_pc: requested component count; implies 'fixed' when pca_rule is omitted.
#   pca_scale: 'covariance' or 'correlation'.
#   wkm_mode: 'sample', 'el', or 'hybrid'.
#   gmm_H_mode: 'el' or 'sample_if'.
#   gmm_solver: 'newton', 'optim', or 'gmm'; 'gmm' requires the gmm package.
#   gmm_start_method_onepop: 'profile_alpha', 'glm_replace_beta', or 'glm'.
#   gmm_start_method_reference: 'profile_alpha' or 'old'.
#   return_gmm_el_pseudo_if, return_reference_if: logical output controls.
#   Additional solver and diagnostic controls are listed by args(zim_case_control).
#
# Returns a list containing theta_ZIM, se_ZIM, theta_MLE, se_MLE, W_hat,
# inverse_diagnostics, convergence, and additional fit components.
#
# Example using existing internal data and external summaries:
#   fit <- zim_case_control(
#     internal_data = dat_internal, Y = 'Y', M_formula = ~ prs + age,
#     working_models = working_models, beta_tilde = beta_external,
#     external_n_case = n_external_case, external_n_control = n_external_control,
#     Sigma_beta_beta = NULL, beta_method = 'GMM')
#   fit$theta_ZIM
#   fit$se_ZIM
zim_case_control <- function(internal_data, Y, M_formula, working_models,
                             beta_tilde,
                             Sigma_beta_beta = NULL,
                             external_n_case = NULL,
                             external_n_control = NULL,
                             rho_k = NULL,
                             reference_data = NULL,
                             W = NULL,
                             ridge = ZIM_DEFAULT_RIDGE,
                             fixed_case_control = TRUE,
                             optim_control = list(maxit = 1000),
                             beta_method = c("GMM", "WKM"),
                             inverse_method = c("solve", "pca", "hbb_pca", "sigma_pca", "m_factor_diag", "hbb_factor_diag", "shrink", "shrinkage"),
                             pc_var = 0.99,
                             n_pc = NULL,
                             eig_tol = 1e-10,
                             pca_rule = c("variance", "condition", "condition_or_variance", "fixed"),
                             pc_max_condition = 1e4,
                             pca_scale = c("covariance", "correlation"),
                             shrink_structured = TRUE,
                             H_B_centering = c("separate", "compact"),
                             sigma_B_centering = c("separate", "compact"),
                             el_mass_tol = ZIM_DEFAULT_EL_MASS_TOL,
                             el_mass_action = c("auto", "stop", "warn"),
                             compact_centering_tol = ZIM_DEFAULT_COMPACT_CENTERING_TOL,
                             joint_psd_tol = ZIM_DEFAULT_PSD_TOL,
                             wkm_rank_tol = ZIM_DEFAULT_BREAD_RANK_TOL,
                             bread_condition_warn = ZIM_DEFAULT_CONDITION_WARN,
                             reference_normalization_tol = 1e-10,
                             profile_equivalence_tol = 1e-8,
                             gmm_start_method_onepop = c("profile_alpha", "glm_replace_beta", "glm"),
                             gmm_start_method_reference = c("profile_alpha", "old"),
                             wkm_mode = c("sample", "el", "hybrid"),
                             gmm_H_mode = c("el", "sample_if"),
                             gmm_solver = c("newton", "optim", "gmm"),
                             gmm_score_tol = 1e-8,
                             gmm_rank_tol = 1e-10,
                             gmm_newton_step_tol = 1e-10,
                             gmm_newton_max_halving = 30L,
                             gmm_newton_fallback_to_optim = TRUE,
                             gmm_package_type = c("twoStep", "cue", "iterative"),
                             gmm_package_wmatrix = c("optimal", "ident"),
                             return_gmm_el_pseudo_if = FALSE,
                             return_reference_if = FALSE,
                             shrink_lambda = NULL) {
  ridge <- validate_ridge(ridge)
  pca_rule_was_missing <- missing(pca_rule)
  beta_method <- match.arg(beta_method)
  inverse_method <- match.arg(inverse_method)
  shrink_lambda <- validate_shrinkage_controls(
    inverse_method, shrink_lambda, beta_method, shrink_structured
  )
  if (identical(inverse_method, "shrinkage") && !is.null(W)) {
    stop("Choose either inverse_method = 'shrinkage' or a user-supplied W, not both.")
  }
  pca_controls <- validate_pca_controls(
    pca_rule = pca_rule,
    pc_var = pc_var,
    n_pc = n_pc,
    pc_max_condition = pc_max_condition,
    pca_scale = pca_scale,
    pca_rule_was_missing = pca_rule_was_missing,
    active = is.null(W) && is_pca_family_method(inverse_method)
  )
  pca_rule <- pca_controls$pca_rule
  pc_var <- pca_controls$pc_var
  n_pc <- pca_controls$n_pc
  pc_max_condition <- pca_controls$pc_max_condition
  pca_scale <- pca_controls$pca_scale
  H_B_centering <- match.arg(H_B_centering)
  sigma_B_centering <- match.arg(sigma_B_centering)
  el_mass_action <- match.arg(el_mass_action)
  el_mass_tol <- validate_positive_tolerance(el_mass_tol, "el_mass_tol")
  compact_centering_tol <- validate_positive_tolerance(
    compact_centering_tol, "compact_centering_tol"
  )
  joint_psd_tol <- validate_positive_tolerance(joint_psd_tol, "joint_psd_tol")
  wkm_rank_tol <- validate_positive_tolerance(wkm_rank_tol, "wkm_rank_tol")
  bread_condition_warn <- validate_positive_tolerance(
    bread_condition_warn, "bread_condition_warn"
  )
  reference_normalization_tol <- validate_positive_tolerance(
    reference_normalization_tol, "reference_normalization_tol"
  )
  profile_equivalence_tol <- validate_positive_tolerance(
    profile_equivalence_tol, "profile_equivalence_tol"
  )
  if (!is.logical(return_reference_if) || length(return_reference_if) != 1L ||
      is.na(return_reference_if)) {
    stop("return_reference_if must be TRUE or FALSE.")
  }
  gmm_start_method_onepop <- match.arg(gmm_start_method_onepop)
  gmm_start_method_reference <- match.arg(gmm_start_method_reference)
  wkm_mode <- match.arg(wkm_mode)
  gmm_H_mode <- match.arg(gmm_H_mode)
  gmm_solver <- match.arg(gmm_solver)
  gmm_package_type <- match.arg(gmm_package_type)
  gmm_package_wmatrix <- match.arg(gmm_package_wmatrix)
  if (length(gmm_score_tol) != 1L || !is.finite(gmm_score_tol) || gmm_score_tol <= 0) {
    stop("gmm_score_tol must be a positive finite scalar.")
  }
  if (length(gmm_rank_tol) != 1L || !is.finite(gmm_rank_tol) || gmm_rank_tol <= 0) {
    stop("gmm_rank_tol must be a positive finite scalar.")
  }
  if (length(gmm_newton_step_tol) != 1L || !is.finite(gmm_newton_step_tol) || gmm_newton_step_tol <= 0) {
    stop("gmm_newton_step_tol must be a positive finite scalar.")
  }
  use_reference <- !is.null(reference_data)
  if (use_reference && !is.null(Sigma_beta_beta)) {
    stop(
      "For the two-population reference-control path, Sigma_beta_beta is ",
      "reconstructed from reference controls and external case/control counts; ",
      "do not supply Sigma_beta_beta."
    )
  }
  if (use_reference && identical(inverse_method, "shrink")) {
    stop(
      "inverse_method = 'shrink' is not supported for the two-population ",
      "reference-control path. Use solve or a PCA-family method."
    )
  }
  if (use_reference && beta_method != "GMM") {
    stop("reference_data is X-only, so beta_method = 'WKM' cannot be fit; use beta_method = 'GMM'.")
  }

  strict_el_mass_required <- !use_reference && (
    is.null(Sigma_beta_beta) ||
      (identical(beta_method, "GMM") && identical(gmm_H_mode, "el")) ||
      (identical(beta_method, "WKM") && wkm_mode %in% c("el", "hybrid"))
  )
  el_mass_action_used <- if (identical(el_mass_action, "auto")) {
    if (strict_el_mass_required) "stop" else "warn"
  } else el_mass_action

  Y_int <- check_binary_y(internal_data[[Y]])
  n <- length(Y_int)
  n1_I <- sum(Y_int == 1)
  n0_I <- sum(Y_int == 0)
  if (n0_I <= 0 || n1_I <= 0) stop("Internal data must contain both cases and controls.")
  rho_I <- n1_I / n0_I

  X_theta <- rhs_model_matrix(M_formula, internal_data, keep_intercept = FALSE)
  target <- fit_profile_logistic(
    Y_int,
    X_theta,
    rho = rho_I,
    rank_tol = if (identical(beta_method, "GMM")) gmm_rank_tol else wkm_rank_tol,
    condition_warn = bread_condition_warn
  )
  el <- compute_internal_el_weights(
    Y = Y_int,
    Delta_hat = target$Delta_hat,
    rho_I = rho_I,
    mass_tol = el_mass_tol,
    mass_action = el_mass_action_used
  )
  p_hat <- el$p_hat

  working <- pack_working_models(working_models, internal_data)
  working_reference <- NULL
  target_reference <- NULL
  if (use_reference) {
    working_reference <- pack_working_models(working_models, reference_data)
    if (length(working_reference) != length(working)) {
      stop("reference_data produced the wrong number of working-model design objects.")
    }
    for (k in seq_along(working)) {
      if (!identical(working_reference[[k]]$gamma_names, working[[k]]$gamma_names)) {
        stop(
          "reference_data working model ", k,
          " has design columns that do not match internal_data. ",
          "Use the same factor levels/contrasts and formulas in both data sets."
        )
      }
    }
    target_reference <- make_reference_target_case_control(
      reference_data = reference_data,
      M_formula = M_formula,
      theta_hat = target$theta_hat,
      normalization_tol = reference_normalization_tol
    )
  }
  K <- length(working)

  has_one_external_count <- !is.null(external_n_case) || !is.null(external_n_control)
  external_info <- NULL
  rho_E <- NA_real_
  n_external <- NA_real_
  if (has_one_external_count) {
    external_info <- validate_external_case_control_counts(external_n_case, external_n_control)
    rho_E <- external_info$rho_E
    n_external <- external_info$n_external
  }
  if (!use_reference && !has_one_external_count) {
    stop("For the one-population one-external-study path, provide external_n_case and external_n_control. These counts define rho_E for the internal beta target and for Sigma_beta_beta.")
  }
  if (use_reference && !has_one_external_count) {
    stop(
      "The two-population reference-control path requires external_n_case and ",
      "external_n_control so rho_E and the reconstructed Sigma_beta_beta can be ",
      "defined on the internal sqrt(n) scale."
    )
  }

  if (use_reference) {
    if (is.null(rho_k)) {
      if (has_one_external_count) {
        rho_k <- rep(rho_E, K)
      } else {
        stop("reference_data path requires rho_k or external_n_case/external_n_control.")
      }
    }
    if (length(rho_k) == 1L) rho_k <- rep(rho_k, K)
    if (length(rho_k) != K) stop("rho_k must have length 1 or length K.")
    if (any(!is.finite(rho_k)) || any(rho_k <= 0)) {
      stop("rho_k must contain positive finite ratios.")
    }
    if (has_one_external_count && max(abs(rho_k - rho_E)) > 1e-8) {
      stop("Supplied rho_k is inconsistent with external_n_case/external_n_control.")
    }
    if (max(abs(rho_k - rho_k[1])) > 1e-8) {
      stop("The reference-control implementation assumes all beta summaries come from one external case-control study with a common rho_E.")
    }
  } else {
    if (!is.null(rho_k)) {
      if (length(rho_k) != 1L && length(rho_k) != K) {
        stop("rho_k is deprecated in the one-external-study path; if supplied, it must have length 1 or length K.")
      }
      rho_k_tmp <- if (length(rho_k) == 1L) rep(rho_k, K) else rho_k
      if (max(abs(rho_k_tmp - rho_k_tmp[1])) > 1e-8) {
        stop("For one external study, all entries of rho_k must be equal. Prefer external_n_case/external_n_control.")
      }
      if (has_one_external_count && abs(rho_k_tmp[1] - rho_E) > 1e-8) {
        stop("Supplied rho_k is inconsistent with external_n_case/external_n_control.")
      }
      if (!has_one_external_count) {
        stop("rho_k alone is no longer sufficient in the one-external-study path; provide external_n_case and external_n_control.")
      }
    }
  }

  beta_names <- unlist(lapply(working, `[[`, "beta_names"), use.names = FALSE)
  p_beta <- length(beta_names)
  beta_tilde <- align_beta_tilde_to_working(beta_tilde, beta_names)

  if (is.null(Sigma_beta_beta) && !has_one_external_count) {
    if (use_reference) {
      stop("When Sigma_beta_beta is not supplied with reference_data, provide external_n_case and external_n_control so the reference-control external sandwich can estimate Sigma_beta_beta on the sqrt(n_internal) scale.")
    } else {
      stop("When Sigma_beta_beta is not supplied, provide external_n_case and external_n_control for the one external study.")
    }
  }

  gamma_hat_GMM <- NULL
  gamma_start_GMM <- NULL
  gmm_start_diagnostics <- NULL
  gamma_hat_WKM <- NULL
  gamma_pros_WKM <- NULL
  beta_hat_GMM <- NULL
  beta_hat_WKM <- NULL
  sigma_external <- NULL
  H <- NULL
  H_weight <- NULL
  H_sample <- NULL
  H_el <- NULL
  weight_source <- NA_character_
  variance_source <- NA_character_
  V_W_internal <- NULL
  V_W_external <- NULL
  V_W_method <- NA_character_

  if (beta_method == "GMM") {
    gamma_hat_list <- vector("list", K)
    gamma_start_list <- vector("list", K)
    gmm_start_diagnostics <- vector("list", K)
    beta_hat_parts <- vector("list", K)
    convergence <- vector("list", K)
    offset <- 0L

    for (k in seq_along(working)) {
      wk <- working[[k]]
      wk_gmm <- if (use_reference) working_reference[[k]] else wk
      beta_idx <- offset + seq_len(wk$p_beta_k)
      if (use_reference) {
        start_k <- make_gmm_start_reference(
          Y = Y_int,
          wk_internal = wk,
          wk_reference = wk_gmm,
          rho_I = rho_I,
          rho_E = rho_k[k],
          beta_tilde_k = beta_tilde[beta_idx],
          Delta_ref = target_reference$Delta_hat,
          gmm_start_method_reference = gmm_start_method_reference,
          optim_control = optim_control
        )
      } else {
        start_k <- make_gmm_start_onepop(
          Y = Y_int,
          wk = wk,
          rho_E = rho_E,
          rho_I = rho_I,
          beta_tilde_k = beta_tilde[beta_idx],
          Delta_hat = target$Delta_hat,
          p_hat = p_hat,
          gmm_start_method = gmm_start_method_onepop,
          optim_control = optim_control
        )
      }
      gamma_start <- start_k$gamma_start
      gamma_start_list[[k]] <- gamma_start
      gmm_start_diagnostics[[k]] <- start_k$diagnostics

      if (use_reference) {
        sol_k <- solve_gamma_gmm_reference_k_dispatch(
          wk = wk_gmm,
          rho_k = rho_k[k],
          Delta_ref = target_reference$Delta_hat,
          gamma_start = gamma_start,
          optim_control = optim_control,
          ridge = ridge,
          gmm_solver = gmm_solver,
          score_tol = gmm_score_tol,
          rank_tol = gmm_rank_tol,
          condition_warn = bread_condition_warn,
          newton_step_tol = gmm_newton_step_tol,
          newton_max_halving = gmm_newton_max_halving,
          newton_fallback_to_optim = gmm_newton_fallback_to_optim,
          gmm_type = gmm_package_type,
          gmm_wmatrix = gmm_package_wmatrix
        )
      } else {
        sol_k <- solve_gamma_gmm_onepop_k_dispatch(
          wk = wk,
          rho_working = rho_E,
          Delta_hat = target$Delta_hat,
          p_hat = p_hat,
          gamma_start = gamma_start,
          optim_control = optim_control,
          ridge = ridge,
          gmm_solver = gmm_solver,
          score_tol = gmm_score_tol,
          rank_tol = gmm_rank_tol,
          newton_step_tol = gmm_newton_step_tol,
          newton_max_halving = gmm_newton_max_halving,
          newton_fallback_to_optim = gmm_newton_fallback_to_optim,
          gmm_type = gmm_package_type,
          gmm_wmatrix = gmm_package_wmatrix
        )
      }

      gamma_hat_list[[k]] <- sol_k$gamma_hat
      beta_local_idx <- wk_gmm$q_alpha_k + seq_len(wk_gmm$p_beta_k)
      beta_hat_parts[[k]] <- sol_k$gamma_hat[beta_local_idx]
      convergence[[k]] <- list(
        convergence = sol_k$convergence,
        objective = sol_k$objective,
        max_abs_moment = sol_k$max_abs_moment,
        moment_norm = if (!is.null(sol_k$moment_norm)) sol_k$moment_norm else NA_real_,
        solution_validated = if (!is.null(sol_k$solution_validated)) sol_k$solution_validated else NA,
        jacobian_rank = if (!is.null(sol_k$jacobian_rank)) sol_k$jacobian_rank else NA_integer_,
        jacobian_dimension = if (!is.null(sol_k$jacobian_dimension)) sol_k$jacobian_dimension else NA_integer_,
        jacobian_rank_tol = if (!is.null(sol_k$jacobian_rank_tol)) sol_k$jacobian_rank_tol else NA_real_,
        jacobian_rank_cutoff = if (!is.null(sol_k$jacobian_rank_cutoff)) sol_k$jacobian_rank_cutoff else NA_real_,
        jacobian_min_singular_value = if (!is.null(sol_k$jacobian_min_singular_value)) sol_k$jacobian_min_singular_value else NA_real_,
        jacobian_max_singular_value = if (!is.null(sol_k$jacobian_max_singular_value)) sol_k$jacobian_max_singular_value else NA_real_,
        jacobian_condition_number = if (!is.null(sol_k$jacobian_condition_number)) sol_k$jacobian_condition_number else NA_real_,
        exp_diagnostic = if (!is.null(sol_k$exp_diagnostic)) sol_k$exp_diagnostic else NULL,
        iterations = if (!is.null(sol_k$iterations)) sol_k$iterations else NA_integer_,
        method = if (!is.null(sol_k$method)) sol_k$method else NA_character_,
        message = if (!is.null(sol_k$message)) sol_k$message else NA_character_,
        rho_working = if (!is.null(sol_k$rho_working)) sol_k$rho_working else NA_real_,
        score_tol = gmm_score_tol,
        gmm_start = start_k$diagnostics,
        newton_result = if (!is.null(sol_k$newton_result)) sol_k$newton_result else NULL,
        package_gmm_object = if (!is.null(sol_k$gmm_object)) sol_k$gmm_object else NULL
      )
      offset <- offset + wk$p_beta_k
    }

    beta_hat_internal <- unlist(beta_hat_parts, use.names = TRUE)
    names(beta_hat_internal) <- beta_names
    beta_hat_GMM <- beta_hat_internal
    gamma_hat_GMM <- gamma_hat_list
    gamma_start_GMM <- gamma_start_list

    if (use_reference) {
      H <- influence_reference_gmm_case_control(
        Y = Y_int,
        target_internal = target,
        target_reference = target_reference,
        working_reference = working_reference,
        gamma_hat_list = gamma_hat_list,
        rho_k = rho_k,
        ridge = ridge,
        fixed_case_control = fixed_case_control,
        rank_tol = gmm_rank_tol,
        condition_warn = bread_condition_warn,
        joint_psd_tol = joint_psd_tol,
        profile_equivalence_tol = profile_equivalence_tol
      )
    } else if (identical(gmm_H_mode, "sample_if")) {

      H <- influence_internal_gmm(
        Y = Y_int,
        target = target,
        working = working,
        gamma_hat_list = gamma_hat_list,
        rho_1 = rho_I,
        rho_k = rep(rho_E, K),
        ridge = ridge,
        fixed_case_control = fixed_case_control,
        joint_psd_tol = joint_psd_tol
      )
    } else {
      H <- estimate_H_onepop_cc(
        method = "GMM",
        target = target,
        working = working,
        gamma_hat_list = gamma_hat_list,
        rho_I = rho_I,
        p_hat = p_hat,
        Delta_hat = target$Delta_hat,
        rho_working = rho_E,
        ridge = ridge,
        B_centering = H_B_centering,
        return_el_pseudo_if = return_gmm_el_pseudo_if,
        mass_tol = el_mass_tol,
        mass_action = el_mass_action_used,
        compact_tol = compact_centering_tol,
        rank_tol = gmm_rank_tol,
        condition_warn = bread_condition_warn,
        joint_psd_tol = joint_psd_tol
      )
    }
  } else {
    if (use_reference) stop("beta_method='WKM' is not available with reference_data.")
    wkm_fit <- fit_internal_wkm_case_control(
      Y = Y_int,
      target = target,
      working = working,
      rho_I = rho_I,
      rho_E = rho_E,
      ridge = ridge,
      fixed_case_control = fixed_case_control,
      rank_tol = wkm_rank_tol,
      condition_warn = bread_condition_warn,
      joint_psd_tol = joint_psd_tol
    )
    beta_hat_internal <- wkm_fit$beta_hat_WKM
    names(beta_hat_internal) <- beta_names
    beta_hat_WKM <- beta_hat_internal
    gamma_hat_WKM <- wkm_fit$gamma_hat_WKM
    gamma_pros_WKM <- wkm_fit$gamma_pros_WKM
    convergence <- wkm_fit$convergence

    H_sample <- wkm_fit
    if (wkm_mode %in% c("el", "hybrid")) {
      H_el <- estimate_H_onepop_cc(
        method = "WKM",
        target = target,
        working = working,
        gamma_hat_list = gamma_hat_WKM,
        rho_I = rho_I,
        p_hat = p_hat,
        Delta_hat = target$Delta_hat,
        rho_working = rho_E,
        ridge = ridge,
        B_centering = H_B_centering,
        mass_tol = el_mass_tol,
        mass_action = el_mass_action_used,
        compact_tol = compact_centering_tol,
        rank_tol = wkm_rank_tol,
        condition_warn = bread_condition_warn,
        joint_psd_tol = joint_psd_tol
      )
    }

    if (identical(wkm_mode, "sample")) {

      H <- H_sample
      H_weight <- H_sample
      weight_source <- "sample_if"
      variance_source <- "sample_if"
    } else if (identical(wkm_mode, "el")) {

      H <- H_el
      H_weight <- H_el
      weight_source <- "el"
      variance_source <- "el"
    } else {

      H <- H_sample
      H_weight <- H_el
      weight_source <- "el"
      variance_source <- "sample_if"
    }
  }
  if (is.null(H_weight)) {
    H_weight <- H
    if (beta_method == "GMM") {
      if (use_reference) {
        weight_source <- "reference_gmm_if"
        variance_source <- "reference_gmm_if"
      } else if (identical(gmm_H_mode, "sample_if")) {
        weight_source <- "gmm_sample_if"
        variance_source <- "sample_if"
      } else {
        weight_source <- "gmm_el"
        variance_source <- "gmm_el"
      }
    }
  }

  if (is.null(Sigma_beta_beta)) {
    gamma_plugin_list <- if (beta_method == "GMM") gamma_hat_GMM else gamma_hat_WKM
    if (use_reference) {
      sigma_external <- estimate_sigma_betabeta_reference_cc(
        working_reference = working_reference,
        gamma_plugin_list = gamma_plugin_list,
        Delta_ref = target_reference$Delta_hat,
        n_internal = n,
        external_n_case = external_n_case,
        external_n_control = external_n_control,
        ridge = ridge,
        B_centering = sigma_B_centering,
        compact_tol = compact_centering_tol,
        rank_tol = gmm_rank_tol,
        condition_warn = bread_condition_warn,
        psd_tol = joint_psd_tol
      )
      sigma_external$gamma_plugin_source <- "GMM_reference_control_external_ratio_target"
    } else {
      sigma_external <- estimate_sigma_betabeta_one_external_cc(
        working = working,
        gamma_plugin_list = gamma_plugin_list,
        Delta_hat = target$Delta_hat,
        p_hat = p_hat,
        n_internal = n,
        external_n_case = external_n_case,
        external_n_control = external_n_control,
        ridge = ridge,
        B_centering = sigma_B_centering,
        mass_tol = el_mass_tol,
        mass_action = el_mass_action_used,
        compact_tol = compact_centering_tol,
        rank_tol = if (identical(beta_method, "GMM")) gmm_rank_tol else wkm_rank_tol,
        condition_warn = bread_condition_warn,
        psd_tol = joint_psd_tol
      )
      sigma_external$gamma_plugin_source <- if (beta_method == "GMM") {
        "GMM_onepop_external_ratio_target"
      } else {
        "WKM_weighted_glm_external_ratio_target"
      }
    }
    Sigma_beta_beta <- sigma_external$Sigma_beta_beta
    rho_E <- sigma_external$rho_E
    n_external <- sigma_external$n_external
    external_info <- list(
      external_n_case = sigma_external$external_n_case,
      external_n_control = sigma_external$external_n_control,
      n_external = sigma_external$n_external,
      rho_E = sigma_external$rho_E
    )
    sigma_method <- sigma_external$method
  } else {
    Sigma_beta_beta <- align_sigma_betabeta_to_working(Sigma_beta_beta, beta_names)
    sigma_method <- "user_supplied"
  }

  H_theta_theta <- H$H_theta_theta
  H_theta_beta <- H$H_theta_beta
  H_beta_beta <- H$H_beta_beta
  dimnames(H_theta_theta) <- list(names(target$theta_hat), names(target$theta_hat))
  dimnames(H_theta_beta) <- list(names(target$theta_hat), beta_names)
  dimnames(H_beta_beta) <- list(beta_names, beta_names)

  H_weight_theta_theta <- H_weight$H_theta_theta
  H_weight_theta_beta <- H_weight$H_theta_beta
  H_weight_beta_beta <- H_weight$H_beta_beta
  dimnames(H_weight_theta_theta) <- list(names(target$theta_hat), names(target$theta_hat))
  dimnames(H_weight_theta_beta) <- list(names(target$theta_hat), beta_names)
  dimnames(H_weight_beta_beta) <- list(beta_names, beta_names)

  IF_beta_for_shrink <- NULL
  if (!is.null(H_weight$IF_beta)) IF_beta_for_shrink <- H_weight$IF_beta
  if (is.null(IF_beta_for_shrink) && !is.null(H_weight$IF_beta_internal)) {
    IF_beta_for_shrink <- H_weight$IF_beta_internal
  }
  IF_theta_for_shrink <- H_weight$IF_theta
  if (inverse_method == "shrink" &&
      (is.null(IF_beta_for_shrink) || is.null(IF_theta_for_shrink))) {
    stop("inverse_method='shrink' requires per-observation IF rows for the W source. For WKM use wkm_mode='sample'; for one-population GMM use gmm_H_mode='sample_if'; otherwise use inverse_method='solve', 'pca', 'hbb_pca', 'sigma_pca', 'm_factor_diag', or 'hbb_factor_diag'.")
  }
  M_beta_raw <- symmetrize(H_beta_beta + Sigma_beta_beta)
  M_beta <- regularize_matrix(M_beta_raw, ridge = ridge)
  M_beta_for_weight_raw <- symmetrize(H_weight_beta_beta + Sigma_beta_beta)
  M_beta_for_weight_regularized <- regularize_matrix(
    M_beta_for_weight_raw, ridge = ridge
  )
  M_beta_for_weight <- M_beta_for_weight_regularized
  M_beta_weight_model <- NULL
  M_beta_weight_model_raw <- NULL
  H_weight_beta_beta_model <- NULL
  H_weight_theta_beta_model <- NULL
  M_beta_weight_reduced <- NULL
  M_beta_weight_reduced_raw <- NULL
  beta_reduction_matrix <- NULL
  H_weight_theta_beta_reduced <- NULL
  if (is.null(W)) {
    inv <- compute_zim_weight(
      H_theta_beta = H_weight_theta_beta,
      H_beta_beta = H_weight_beta_beta,
      Sigma_beta_beta = Sigma_beta_beta,
      ridge = ridge,
      inverse_method = inverse_method,
      pc_var = pc_var,
      n_pc = n_pc,
      eig_tol = eig_tol,
      pca_rule = pca_rule,
      pc_max_condition = pc_max_condition,
      pca_scale = pca_scale,
      H_theta_theta = H_weight_theta_theta,
      IF_theta = IF_theta_for_shrink,
      IF_beta = IF_beta_for_shrink,
      Y = Y_int,
      fixed_case_control = fixed_case_control,
      shrink_structured = shrink_structured,
      shrink_lambda = shrink_lambda
    )
    M_beta_for_weight <- if (inverse_method %in% c("m_factor_diag", "hbb_factor_diag", "shrinkage") &&
                             !is.null(inv$M_beta_weight_model)) {
      inv$M_beta_weight_model
    } else {
      inv$M_beta
    }
    M_beta_weight_model <- inv$M_beta_weight_model
    M_beta_weight_model_raw <- inv$M_beta_weight_model_raw
    H_weight_beta_beta_model <- inv$H_beta_beta_weight_model
    H_weight_theta_beta_model <- inv$H_theta_beta_weight_model
    M_beta_weight_reduced <- inv$M_beta_weight_reduced
    M_beta_weight_reduced_raw <- inv$M_beta_weight_reduced_raw
    beta_reduction_matrix <- inv$beta_reduction_matrix
    H_weight_theta_beta_reduced <- inv$H_theta_beta_weight_reduced
    W_hat <- inv$W_hat
    weight_method <- paste0("estimated_", inverse_method)
    if (beta_method == "WKM") {
      weight_method <- paste0(weight_method, "_wkm_mode_", wkm_mode)
    } else if (beta_method == "GMM" && !use_reference) {
      weight_method <- paste0(weight_method, "_gmm_H_mode_", gmm_H_mode)
    }
    inverse_diagnostics <- inv$diagnostics
    if (is_pca_family_method(inverse_method)) {
      inverse_diagnostics$pca_rule_inferred_from_n_pc <-
        pca_controls$inferred_from_n_pc
    }
  } else {
    W_hat <- as.matrix(W)
    if (!all(dim(W_hat) == c(length(target$theta_hat), p_beta))) {
      stop("User-supplied W must be a d by p_beta matrix, where d = length(theta).")
    }
    weight_method <- "user_supplied"
    inverse_diagnostics <- list(
      inverse_method = "none_user_supplied_W",
      ridge = ridge,
      p_beta = p_beta,
      pca_rule = NA_character_,
      pca_scale = NA_character_,
      pca_rule_inferred_from_n_pc = FALSE,
      n_pc_requested = NA_integer_,
      n_pc_used = NA_integer_,
      pc_var_requested = NA_real_,
      pc_var_used = NA_real_,
      pc_max_condition_requested = NA_real_,
      eig_tol_requested = NA_real_,
      n_positive_eigenvalues = NA_integer_,
      numerical_eigen_cutoff = NA_real_,
      rule_eigen_cutoff = NA_real_,
      effective_eigen_cutoff = NA_real_,
      retained_eigen_min = NA_real_,
      first_discarded_eigenvalue = NA_real_,
      eigen_min = NA_real_,
      eigen_max = NA_real_,
      condition_number = NA_real_,
      retained_condition_number = NA_real_,
      covariance_condition_number = NA_real_,
      beta_scale_sd_min = NA_real_,
      beta_scale_sd_max = NA_real_,
      factor_diag_basis = NA_character_,
      factor_diag_diagonal_target = NA_character_,
      factor_diag_diagonal_preserved = NA,
      factor_diag_max_abs_diagonal_error = NA_real_,
      factor_diag_total_M_max_abs_diagonal_error = NA_real_,
      factor_diag_residual_variance_min = NA_real_,
      factor_diag_residual_variance_max = NA_real_,
      factor_diag_offdiag_removed_frobenius = NA_real_,
      factor_diag_basis_raw_min_eigenvalue = NA_real_,
      factor_diag_basis_material_negative_eigenvalues = NA_integer_,
      factor_diag_matrix_min_eigenvalue = NA_real_,
      factor_diag_matrix_condition_number = NA_real_,
      factor_diag_H_theta_beta_projected = NA,
      factor_diag_H_theta_beta_projection_max_abs_change = NA_real_,
      shrink_lambda = NA_real_,
      shrink_structured = NA
    )
  }
  dimnames(W_hat) <- list(names(target$theta_hat), beta_names)

  theta_ZIM <- as.numeric(target$theta_hat - W_hat %*% (beta_hat_internal - beta_tilde))
  names(theta_ZIM) <- names(target$theta_hat)

  if (identical(variance_source, "sample_if") &&
      !is.null(H$IF_theta) && !is.null(H$IF_beta)) {
    g_if <- H$IF_theta - H$IF_beta %*% t(W_hat)
    V_W_internal <- symmetrize(case_control_cov(g_if, Y_int, fixed_case_control))
    V_W_external <- symmetrize(W_hat %*% Sigma_beta_beta %*% t(W_hat))
    V_W <- symmetrize(V_W_internal + V_W_external)
    V_W_method <- "sample_if_residual_plus_external_summary"
  } else {
    V_W_internal <- symmetrize(
      H_theta_theta -
        H_theta_beta %*% t(W_hat) -
        W_hat %*% t(H_theta_beta) +
        W_hat %*% H_beta_beta %*% t(W_hat)
    )
    V_W_external <- symmetrize(W_hat %*% Sigma_beta_beta %*% t(W_hat))
    V_W <- symmetrize(V_W_internal + V_W_external)
    V_W_method <- "expanded_zim_sandwich"
  }
  dimnames(V_W) <- dimnames(H_theta_theta)
  dimnames(V_W_internal) <- dimnames(H_theta_theta)
  dimnames(V_W_external) <- dimnames(H_theta_theta)

  V_W_source_decomposed <- NULL
  V_W_source_internal_case <- NULL
  V_W_source_reference <- NULL
  V_W_source_external <- NULL
  V_W_source_max_abs_diff <- NA_real_
  V_W_source_relative_diff <- NA_real_
  if (use_reference) {
    g_internal_source <- H$IF_theta - H$IF_beta_internal %*% t(W_hat)
    g_reference_source <- -H$IF_beta_reference %*% t(W_hat)
    V_W_source_internal_case <- symmetrize(
      case_control_cov(g_internal_source, Y_int, fixed_case_control)
    )
    V_W_source_reference <- symmetrize(
      reference_centered_cov_scaled(
        g_reference_source,
        scale_factor = n / H$n_reference
      )
    )
    V_W_source_external <- symmetrize(
      W_hat %*% Sigma_beta_beta %*% t(W_hat)
    )
    V_W_source_decomposed <- symmetrize(
      V_W_source_internal_case +
        V_W_source_reference +
        V_W_source_external
    )
    V_W_source_max_abs_diff <- max(abs(V_W - V_W_source_decomposed))
    source_scale <- max(1, max(abs(V_W)), max(abs(V_W_source_decomposed)))
    V_W_source_relative_diff <- V_W_source_max_abs_diff / source_scale
    source_tol <- 100 * sqrt(.Machine$double.eps)
    if (!is.finite(V_W_source_relative_diff) ||
        V_W_source_relative_diff > source_tol) {
      warning(
        "Expanded and source-decomposed reference-control variances differ: ",
        "relative max difference = ", signif(V_W_source_relative_diff, 6),
        ", tolerance = ", signif(source_tol, 6), "."
      )
    }
    V_W_method <- paste0(
      V_W_method,
      "_with_reference_source_decomposition_check"
    )
    dimnames(V_W_source_decomposed) <- dimnames(H_theta_theta)
    dimnames(V_W_source_internal_case) <- dimnames(H_theta_theta)
    dimnames(V_W_source_reference) <- dimnames(H_theta_theta)
    dimnames(V_W_source_external) <- dimnames(H_theta_theta)
  }

  se_ZIM <- sqrt(pmax(diag(V_W) / n, 0))
  z_ZIM <- theta_ZIM / se_ZIM
  p_ZIM <- 2 * stats::pnorm(abs(z_ZIM), lower.tail = FALSE)

  se_MLE <- sqrt(pmax(diag(H_theta_theta) / n, 0))
  H_theta_theta_glm <- target_glm_H_theta_theta(target, ridge = ridge)
  se_MLE_glm <- sqrt(pmax(diag(H_theta_theta_glm) / n, 0))
  H_theta_theta_glm_max_abs_diff <- max(abs(H_theta_theta - H_theta_theta_glm))
  H_theta_theta_glm_max_rel_diff <- max(
    abs(H_theta_theta - H_theta_theta_glm) /
      pmax(abs(H_theta_theta_glm), .Machine$double.eps)
  )
  rho_k_return <- if (use_reference) {
    rho_k
  } else if (is.finite(rho_E)) {
    rep(rho_E, K)
  } else {
    rep(NA_real_, K)
  }

  list(
    theta_ZIM = theta_ZIM,
    se_ZIM = se_ZIM,
    ridge = ridge,
    el_mass_tol = el_mass_tol,
    el_mass_action_requested = el_mass_action,
    el_mass_action_used = el_mass_action_used,
    compact_centering_tol = compact_centering_tol,
    joint_psd_tol = joint_psd_tol,
    wkm_rank_tol = wkm_rank_tol,
    bread_condition_warn = bread_condition_warn,
    reference_normalization_tol = reference_normalization_tol,
    profile_equivalence_tol = profile_equivalence_tol,
    z_ZIM = z_ZIM,
    p_ZIM = p_ZIM,
    theta_MLE = target$theta_hat,
    se_MLE = se_MLE,
    se_MLE_glm = se_MLE_glm,
    tau_hat = target$tau_hat,
    beta_method = beta_method,
    weight_method = weight_method,
    reference_data_used = use_reference,
    n_reference = if (use_reference) nrow(reference_data) else NA_integer_,
    tau_reference = if (use_reference) target_reference$tau_hat else NA_real_,
    Delta_reference = if (use_reference) target_reference$Delta_hat else NULL,
    reference_target = if (use_reference) target_reference else NULL,
    beta_hat = beta_hat_internal,
    beta_hat_GMM = beta_hat_GMM,
    beta_hat_WKM = beta_hat_WKM,
    beta_tilde = beta_tilde,
    gamma_hat_GMM = gamma_hat_GMM,
    gamma_start_GMM = gamma_start_GMM,
    gmm_start_method = if (beta_method == "GMM") {
      if (use_reference) gmm_start_method_reference else gmm_start_method_onepop
    } else NA_character_,
    gmm_active_start_method = if (beta_method == "GMM") {
      if (use_reference) paste0("reference_", gmm_start_method_reference) else paste0("onepop_", gmm_start_method_onepop)
    } else NA_character_,
    gmm_start_method_onepop = if (beta_method == "GMM" && !use_reference) gmm_start_method_onepop else NA_character_,
    gmm_start_method_reference = if (beta_method == "GMM" && use_reference) gmm_start_method_reference else NA_character_,
    gmm_solver = if (beta_method == "GMM") gmm_solver else NA_character_,
    gmm_score_tol = if (beta_method == "GMM") gmm_score_tol else NA_real_,
    gmm_rank_tol = if (beta_method == "GMM") gmm_rank_tol else NA_real_,
    gmm_package_type = if (beta_method == "GMM" && identical(gmm_solver, "gmm")) gmm_package_type else NA_character_,
    gmm_package_wmatrix = if (beta_method == "GMM" && identical(gmm_solver, "gmm")) gmm_package_wmatrix else NA_character_,
    gmm_start_diagnostics = if (beta_method == "GMM") gmm_start_diagnostics else NULL,
    reference_gmm_diagnostics = if (beta_method == "GMM" && use_reference) convergence else NULL,
    gamma_hat_WKM = gamma_hat_WKM,
    gamma_pros_WKM = gamma_pros_WKM,
    W_hat = W_hat,
    V_W = V_W,
    H_theta_theta = H_theta_theta,
    H_theta_theta_glm = H_theta_theta_glm,
    H_theta_theta_glm_max_abs_diff = H_theta_theta_glm_max_abs_diff,
    H_theta_theta_glm_max_rel_diff = H_theta_theta_glm_max_rel_diff,
    H_theta_beta = H_theta_beta,
    H_beta_beta = H_beta_beta,
    G = if (use_reference) H$G else NULL,
    IF_beta_ref = if (use_reference && isTRUE(return_reference_if)) {
      H$IF_beta_reference
    } else NULL,
    IF_beta_internal_reference_path = if (use_reference && isTRUE(return_reference_if)) {
      H$IF_beta_internal
    } else NULL,
    reference_profile_equivalence_diagnostics = if (use_reference) {
      H$profile_equivalence_diagnostics
    } else NULL,
    reference_profile_equivalence_max_abs = if (use_reference) {
      H$profile_equivalence_max_abs
    } else NA_real_,
    reference_identity_diagnostics = if (use_reference) H$identity_diagnostics else NULL,
    reference_bread_rank_diagnostics = if (use_reference) {
      list(
        target = H$target_bread_rank_diagnostics,
        working = H$D_rank_diagnostics,
        external = if (!is.null(sigma_external)) sigma_external$A_E_rank_diagnostics else NULL,
        external_regularized = if (!is.null(sigma_external)) {
          sigma_external$A_E_regularized_rank_diagnostics
        } else NULL
      )
    } else NULL,
    H_method = H$method,
    gmm_H_mode = if (beta_method == "GMM" && !use_reference) gmm_H_mode else NA_character_,
    gmm_W_source = if (beta_method == "GMM") weight_source else NA_character_,
    gmm_variance_source = if (beta_method == "GMM") variance_source else NA_character_,
    wkm_mode = if (beta_method == "WKM") wkm_mode else NA_character_,
    wkm_W_source = if (beta_method == "WKM") weight_source else NA_character_,
    wkm_variance_source = if (beta_method == "WKM") variance_source else NA_character_,
    H_weight_method = if (!is.null(H_weight$method)) H_weight$method else NA_character_,
    H_B_centering = if (!is.null(H$B_centering)) H$B_centering else NA_character_,
    H_weight_B_centering = if (!is.null(H_weight$B_centering)) H_weight$B_centering else NA_character_,
    H_max_abs_center_resid = if (!is.null(H$max_abs_center_resid)) H$max_abs_center_resid else NA_real_,
    H_weight_max_abs_center_resid = if (!is.null(H_weight$max_abs_center_resid)) H_weight$max_abs_center_resid else NA_real_,
    H_joint_is_psd = if (!is.null(H$H_joint_is_psd)) H$H_joint_is_psd else NA,
    H_joint_min_eigenvalue = if (!is.null(H$H_joint_min_eigenvalue)) H$H_joint_min_eigenvalue else NA_real_,
    H_joint_psd_cutoff = if (!is.null(H$H_joint_psd_cutoff)) H$H_joint_psd_cutoff else NA_real_,
    H_weight_joint_is_psd = if (!is.null(H_weight$H_joint_is_psd)) H_weight$H_joint_is_psd else NA,
    H_weight_joint_min_eigenvalue = if (!is.null(H_weight$H_joint_min_eigenvalue)) H_weight$H_joint_min_eigenvalue else NA_real_,
    H_full = H,
    H_weight_full = H_weight,
    H_sample_full = H_sample,
    H_el_full = H_el,
    Sigma_beta_beta = Sigma_beta_beta,
    Sigma_beta_beta_method = sigma_method,
    reference_sigma_reconstructed = use_reference,
    Sigma_B_centering = if (!is.null(sigma_external) && !is.null(sigma_external$B_centering)) sigma_external$B_centering else if (identical(sigma_method, "user_supplied")) "user_supplied" else NA_character_,
    Sigma_max_abs_center_resid = if (!is.null(sigma_external) && !is.null(sigma_external$max_abs_center_resid_E)) sigma_external$max_abs_center_resid_E else NA_real_,
    sigma_external = sigma_external,
    external_psd_diagnostics = if (!is.null(sigma_external)) {
      list(
        H_ext = sigma_external$H_ext_psd_diagnostics,
        Sigma = sigma_external$Sigma_psd_diagnostics
      )
    } else NULL,
    external_blockwise_beta_IF = if (use_reference && !is.null(sigma_external)) {
      sigma_external$blockwise_beta_IF
    } else NA,
    external_full_gamma_sandwich_formed = if (
      use_reference && !is.null(sigma_external)
    ) {
      sigma_external$full_gamma_sandwich_formed
    } else NA,
    M_beta_raw = M_beta_raw,
    M_beta = M_beta,
    M_beta_regularized = M_beta,
    M_beta_for_weight_raw = M_beta_for_weight_raw,
    M_beta_for_weight_regularized = M_beta_for_weight_regularized,
    M_beta_for_weight = M_beta_for_weight,
    M_beta_weight_model = M_beta_weight_model,
    M_beta_weight_model_raw = M_beta_weight_model_raw,
    H_weight_beta_beta_model = H_weight_beta_beta_model,
    H_weight_theta_beta_model = H_weight_theta_beta_model,
    M_beta_weight_reduced = M_beta_weight_reduced,
    M_beta_weight_reduced_raw = M_beta_weight_reduced_raw,
    beta_reduction_matrix = beta_reduction_matrix,
    H_weight_theta_beta_reduced = H_weight_theta_beta_reduced,
    V_W_internal = V_W_internal,
    V_W_external = V_W_external,
    V_W_source_decomposed = V_W_source_decomposed,
    V_W_source_internal_case = V_W_source_internal_case,
    V_W_source_reference = V_W_source_reference,
    V_W_source_external = V_W_source_external,
    V_W_source_max_abs_diff = V_W_source_max_abs_diff,
    V_W_source_relative_diff = V_W_source_relative_diff,
    V_W_method = V_W_method,
    inverse_method = if (is.null(W)) inverse_method else "user_supplied_W",
    pca_rule = if (is.null(W) && is_pca_family_method(inverse_method)) pca_rule else NA_character_,
    pca_scale = if (is.null(W) && is_pca_family_method(inverse_method)) pca_scale else NA_character_,
    pca_rule_inferred_from_n_pc = if (is.null(W) && is_pca_family_method(inverse_method)) {
      pca_controls$inferred_from_n_pc
    } else FALSE,
    pc_var = if (is.null(W) && is_pca_family_method(inverse_method) &&
                 pca_rule %in% c("variance", "condition_or_variance")) pc_var else NA_real_,
    pc_max_condition = if (is.null(W) && is_pca_family_method(inverse_method) &&
                           pca_rule %in% c("condition", "condition_or_variance")) pc_max_condition else NA_real_,
    n_pc = if (is.null(W) && is_pca_family_method(inverse_method) &&
               identical(pca_rule, "fixed")) n_pc else NA_integer_,
    inverse_diagnostics = inverse_diagnostics,
    shrink_lambda = inverse_diagnostics$shrink_lambda,
    IF_theta = if (!is.null(H$IF_theta)) H$IF_theta else IF_theta_for_shrink,
    IF_beta = if (!is.null(H$IF_beta)) H$IF_beta else IF_beta_for_shrink,
    IF_theta_for_weight = IF_theta_for_shrink,
    IF_beta_for_weight = IF_beta_for_shrink,
    rho_I = rho_I,
    rho_1 = rho_I,
    rho_E = rho_E,
    rho_k = rho_k_return,
    n_internal = n,
    n_case_internal = n1_I,
    n_control_internal = n0_I,
    external_n_case = if (!is.null(external_info)) external_info$external_n_case else external_n_case,
    external_n_control = if (!is.null(external_info)) external_info$external_n_control else external_n_control,
    n_external = n_external,
    p_hat = p_hat,
    target_diagnostics = target$diagnostics,
    reference_target_diagnostics = if (use_reference) target_reference$diagnostics else NULL,
    el_diagnostics = el$diagnostics,
    convergence = convergence,
    call = match.call()
  )
}

# Convenience wrapper for zim_case_control(..., beta_method = 'WKM').
# Arguments: ... supplies zim_case_control arguments other than beta_method;
# ridge is a nonnegative numeric value with default ZIM_DEFAULT_RIDGE.
# Returns the same result structure as zim_case_control().
zim_case_control_wkm <- function(..., ridge = ZIM_DEFAULT_RIDGE) {
  ridge <- validate_ridge(ridge)
  zim_case_control(..., beta_method = "WKM", ridge = ridge)
}

make_stratified_folds <- function(Y, K, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  fold <- integer(length(Y))
  for (y in c(0, 1)) {
    idx <- which(Y == y)
    if (length(idx) > 0) fold[idx] <- sample(rep_len(seq_len(K), length(idx)))
  }
  fold
}

make_plain_folds <- function(n, K, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  sample(rep_len(seq_len(K), n))
}

crossfit_weight_on_train <- function(IF_theta, IF_beta, Y, train_idx,
                                     Sigma_beta_beta, ridge, inverse_method,
                                     pc_var, n_pc, eig_tol,
                                     pca_rule, pc_max_condition, pca_scale,
                                     fixed_case_control, shrink_structured,
                                     shrink_lambda = NULL) {
  ridge <- validate_ridge(ridge)
  IFt <- IF_theta[train_idx, , drop = FALSE]
  IFb <- IF_beta[train_idx, , drop = FALSE]
  Ytr <- Y[train_idx]
  Htt <- symmetrize(case_control_cov(IFt, Ytr, fixed_case_control))
  Htb <- case_control_crosscov(IFt, IFb, Ytr, fixed_case_control)
  Hbb <- symmetrize(case_control_cov(IFb, Ytr, fixed_case_control))
  inv <- compute_zim_weight(
    H_theta_beta = Htb,
    H_beta_beta = Hbb,
    Sigma_beta_beta = Sigma_beta_beta,
    ridge = ridge,
    inverse_method = inverse_method,
    pc_var = pc_var,
    n_pc = n_pc,
    eig_tol = eig_tol,
    pca_rule = pca_rule,
    pc_max_condition = pc_max_condition,
    pca_scale = pca_scale,
    H_theta_theta = Htt,
    IF_theta = IFt,
    IF_beta = IFb,
    Y = Ytr,
    fixed_case_control = fixed_case_control,
    shrink_structured = shrink_structured,
    shrink_lambda = shrink_lambda
  )
  inv$W_hat
}

el_pseudo_crosscov <- function(IF0_left, IF1_left,
                               IF0_right, IF1_right,
                               p_hat, Delta_hat, rho_I,
                               idx = NULL) {
  if (is.null(idx)) idx <- seq_along(p_hat)
  p0 <- p_hat[idx]
  p1 <- p_hat[idx] * Delta_hat[idx]
  mass0 <- sum(p0)
  mass1 <- sum(p1)
  if (!is.finite(mass0) || mass0 <= 0 || !is.finite(mass1) || mass1 <= 0) {
    stop("EL pseudo-row masses must be positive in every training split.")
  }

  L0 <- IF0_left[idx, , drop = FALSE]
  L1 <- IF1_left[idx, , drop = FALSE]
  R0 <- IF0_right[idx, , drop = FALSE]
  R1 <- IF1_right[idx, , drop = FALSE]

  muL0 <- colSums(L0 * as.numeric(p0)) / mass0
  muR0 <- colSums(R0 * as.numeric(p0)) / mass0
  muL1 <- colSums(L1 * as.numeric(p1)) / mass1
  muR1 <- colSums(R1 * as.numeric(p1)) / mass1

  L0c <- sweep(L0, 2, muL0, "-")
  R0c <- sweep(R0, 2, muR0, "-")
  L1c <- sweep(L1, 2, muL1, "-")
  R1c <- sweep(R1, 2, muR1, "-")

  w0 <- (p0 / mass0) / (1 + rho_I)
  w1 <- rho_I * (p1 / mass1) / (1 + rho_I)

  crossprod(L0c, R0c * as.numeric(w0)) +
    crossprod(L1c, R1c * as.numeric(w1))
}

el_pseudo_cov <- function(IF0, IF1, p_hat, Delta_hat, rho_I, idx = NULL) {
  symmetrize(el_pseudo_crosscov(
    IF0_left = IF0,
    IF1_left = IF1,
    IF0_right = IF0,
    IF1_right = IF1,
    p_hat = p_hat,
    Delta_hat = Delta_hat,
    rho_I = rho_I,
    idx = idx
  ))
}

el_pseudo_fold_mass <- function(p_hat, Delta_hat, rho_I, idx = NULL) {
  if (is.null(idx)) idx <- seq_along(p_hat)
  sum(p_hat[idx] / (1 + rho_I) + rho_I * p_hat[idx] * Delta_hat[idx] / (1 + rho_I))
}

crossfit_weight_on_train_el_gmm <- function(EL_IF0_theta, EL_IF1_theta,
                                            EL_IF0_beta, EL_IF1_beta,
                                            p_hat, Delta_hat, rho_I,
                                            train_idx,
                                            Sigma_beta_beta, ridge,
                                            inverse_method, pc_var, n_pc,
                                            eig_tol, pca_rule,
                                            pc_max_condition, pca_scale,
                                            shrink_lambda = NULL) {
  ridge <- validate_ridge(ridge)
  if (identical(inverse_method, "shrink")) {
    stop("GMM EL cross-fitting does not support automatic 'shrink'; use 'shrinkage' with shrink_lambda for target-preserving residual shrinkage.")
  }
  Htt <- el_pseudo_cov(EL_IF0_theta, EL_IF1_theta, p_hat, Delta_hat, rho_I, train_idx)
  Htb <- el_pseudo_crosscov(
    IF0_left = EL_IF0_theta,
    IF1_left = EL_IF1_theta,
    IF0_right = EL_IF0_beta,
    IF1_right = EL_IF1_beta,
    p_hat = p_hat,
    Delta_hat = Delta_hat,
    rho_I = rho_I,
    idx = train_idx
  )
  Hbb <- el_pseudo_cov(EL_IF0_beta, EL_IF1_beta, p_hat, Delta_hat, rho_I, train_idx)
  inv <- compute_zim_weight(
    H_theta_beta = Htb,
    H_beta_beta = Hbb,
    Sigma_beta_beta = Sigma_beta_beta,
    ridge = ridge,
    inverse_method = inverse_method,
    pc_var = pc_var,
    n_pc = n_pc,
    eig_tol = eig_tol,
    pca_rule = pca_rule,
    pc_max_condition = pc_max_condition,
    pca_scale = pca_scale,
    H_theta_theta = Htt,
    shrink_lambda = shrink_lambda
  )
  inv$W_hat
}

reference_centered_cov_scaled <- function(X, scale_factor) {
  Xc <- sweep(X, 2, colMeans(X), "-")
  scale_factor * crossprod(Xc) / nrow(Xc)
}

crossfit_weight_on_train_reference_gmm <- function(IF_theta,
                                                   IF_beta_internal,
                                                   IF_beta_reference,
                                                   Y,
                                                   train_internal_idx,
                                                   train_reference_idx,
                                                   n_internal,
                                                   n_reference,
                                                   Sigma_beta_beta,
                                                   ridge,
                                                   inverse_method,
                                                   pc_var,
                                                   n_pc,
                                                   eig_tol,
                                                   pca_rule,
                                                   pc_max_condition,
                                                   pca_scale,
                                                   fixed_case_control,
                                                   shrink_lambda = NULL) {
  ridge <- validate_ridge(ridge)
  if (identical(inverse_method, "shrink")) {
    stop("Reference-control GMM cross-fitting does not support automatic 'shrink'; use 'shrinkage' with shrink_lambda for target-preserving residual shrinkage.")
  }
  IFt <- IF_theta[train_internal_idx, , drop = FALSE]
  IFb_int <- IF_beta_internal[train_internal_idx, , drop = FALSE]
  Ytr <- Y[train_internal_idx]
  IFb_ref <- IF_beta_reference[train_reference_idx, , drop = FALSE]

  Htt <- symmetrize(case_control_cov(IFt, Ytr, fixed_case_control))
  Htb <- case_control_crosscov(IFt, IFb_int, Ytr, fixed_case_control)
  Hbb_internal <- symmetrize(case_control_cov(IFb_int, Ytr, fixed_case_control))
  Hbb_reference <- symmetrize(reference_centered_cov_scaled(
    IFb_ref,
    scale_factor = n_internal / n_reference
  ))
  Hbb <- symmetrize(Hbb_internal + Hbb_reference)

  inv <- compute_zim_weight(
    H_theta_beta = Htb,
    H_beta_beta = Hbb,
    Sigma_beta_beta = Sigma_beta_beta,
    ridge = ridge,
    inverse_method = inverse_method,
    pc_var = pc_var,
    n_pc = n_pc,
    eig_tol = eig_tol,
    pca_rule = pca_rule,
    pc_max_condition = pc_max_condition,
    pca_scale = pca_scale,
    H_theta_theta = Htt,
    shrink_lambda = shrink_lambda
  )
  inv$W_hat
}

# Fit cross-fitted ZIM for a case-control study.
#
# Data, summary, reference, and tuning arguments follow zim_case_control().
# K is the number of folds; n_reps is the number of repeated splits; fold_seed
# is an optional random seed. beta_method defaults to 'WKM'; WKM requires
# wkm_mode = 'sample' and no reference_data. GMM supports either population
# setting. See args(zim_case_control_crossfit) for the full argument list.
#
# Returns theta_ZIM, se_ZIM, theta_MLE, se_MLE, fold assignments, per-repeat
# results, and fit_full. theta_ZIM_plugin and se_ZIM_plugin contain the
# corresponding non-cross-fitted results.
zim_case_control_crossfit <- function(internal_data, Y, M_formula, working_models,
                                      beta_tilde,
                                      Sigma_beta_beta = NULL,
                                      external_n_case = NULL,
                                      external_n_control = NULL,
                                      rho_k = NULL,
                                      K = 5,
                                      reference_data = NULL,
                                      ridge = ZIM_DEFAULT_RIDGE,
                                      fixed_case_control = TRUE,
                                      optim_control = list(maxit = 1000),
                                      beta_method = c("WKM", "GMM"),
                                      wkm_mode = c("sample", "el", "hybrid"),
                                      inverse_method = c("solve", "pca", "hbb_pca", "sigma_pca", "m_factor_diag", "hbb_factor_diag", "shrink", "shrinkage"),
                                      pc_var = 0.99,
                                      n_pc = NULL,
                                      eig_tol = 1e-10,
                                      pca_rule = c("variance", "condition", "condition_or_variance", "fixed"),
                                      pc_max_condition = 1e4,
                                      pca_scale = c("covariance", "correlation"),
                                      shrink_structured = TRUE,
                                      H_B_centering = c("separate", "compact"),
                                      sigma_B_centering = c("separate", "compact"),
                                      el_mass_tol = ZIM_DEFAULT_EL_MASS_TOL,
                                      el_mass_action = c("auto", "stop", "warn"),
                                      compact_centering_tol = ZIM_DEFAULT_COMPACT_CENTERING_TOL,
                                      joint_psd_tol = ZIM_DEFAULT_PSD_TOL,
                                      wkm_rank_tol = ZIM_DEFAULT_BREAD_RANK_TOL,
                                      bread_condition_warn = ZIM_DEFAULT_CONDITION_WARN,
                                      gmm_start_method_onepop = c("profile_alpha", "glm_replace_beta", "glm"),
                                      gmm_start_method_reference = c("profile_alpha", "old"),
                                      gmm_H_mode = c("el", "sample_if"),
                                      gmm_solver = c("newton", "optim", "gmm"),
                                      gmm_score_tol = 1e-8,
                                      gmm_rank_tol = 1e-10,
                                      gmm_newton_step_tol = 1e-10,
                                      gmm_newton_max_halving = 30L,
                                      gmm_newton_fallback_to_optim = TRUE,
                                      gmm_package_type = c("twoStep", "cue", "iterative"),
                                      gmm_package_wmatrix = c("optimal", "ident"),
                                      n_reps = 1L,
                                      fold_seed = NULL,
                                      shrink_lambda = NULL) {
  ridge <- validate_ridge(ridge)
  gmm_solver_was_missing <- missing(gmm_solver)
  pca_rule_was_missing <- missing(pca_rule)
  beta_method <- match.arg(beta_method)
  wkm_mode <- match.arg(wkm_mode)
  inverse_method <- match.arg(inverse_method)
  shrink_lambda <- validate_shrinkage_controls(
    inverse_method, shrink_lambda, beta_method, shrink_structured
  )
  pca_controls <- validate_pca_controls(
    pca_rule = pca_rule,
    pc_var = pc_var,
    n_pc = n_pc,
    pc_max_condition = pc_max_condition,
    pca_scale = pca_scale,
    pca_rule_was_missing = pca_rule_was_missing,
    active = is_pca_family_method(inverse_method)
  )
  pca_rule <- pca_controls$pca_rule
  pc_var <- pca_controls$pc_var
  n_pc <- pca_controls$n_pc
  pc_max_condition <- pca_controls$pc_max_condition
  pca_scale <- pca_controls$pca_scale
  H_B_centering <- match.arg(H_B_centering)
  sigma_B_centering <- match.arg(sigma_B_centering)
  el_mass_action <- match.arg(el_mass_action)
  el_mass_tol <- validate_positive_tolerance(el_mass_tol, "el_mass_tol")
  compact_centering_tol <- validate_positive_tolerance(
    compact_centering_tol, "compact_centering_tol"
  )
  joint_psd_tol <- validate_positive_tolerance(joint_psd_tol, "joint_psd_tol")
  wkm_rank_tol <- validate_positive_tolerance(wkm_rank_tol, "wkm_rank_tol")
  bread_condition_warn <- validate_positive_tolerance(
    bread_condition_warn, "bread_condition_warn"
  )
  gmm_start_method_onepop <- match.arg(gmm_start_method_onepop)
  gmm_start_method_reference <- match.arg(gmm_start_method_reference)
  gmm_H_mode <- match.arg(gmm_H_mode)
  gmm_solver <- match.arg(gmm_solver)
  gmm_package_type <- match.arg(gmm_package_type)
  gmm_package_wmatrix <- match.arg(gmm_package_wmatrix)

  use_reference <- !is.null(reference_data)
  if (gmm_solver_was_missing && use_reference) {
    gmm_solver <- "optim"
  }
  if (identical(beta_method, "WKM")) {
    if (use_reference) stop("WKM cross-fitting is not available with X-only reference_data.")
    if (!identical(wkm_mode, "sample")) {
      stop("WKM cross-fitting is implemented only for wkm_mode = 'sample'.")
    }
  }
  if (identical(beta_method, "GMM") && identical(inverse_method, "shrink") &&
      (use_reference || !identical(gmm_H_mode, "sample_if"))) {
    stop("GMM cross-fitting with shrinkage is available only for the one-population gmm_H_mode = 'sample_if' path.")
  }
  if (identical(beta_method, "GMM") && !use_reference && identical(gmm_H_mode, "el") &&
      !identical(H_B_centering, "separate")) {
    stop("One-population GMM EL cross-fitting requires H_B_centering = 'separate'.")
  }

  K <- as.integer(K)
  if (is.na(K) || K < 2L) stop("K must be an integer at least 2.")
  n_reps <- as.integer(n_reps)
  if (is.na(n_reps) || n_reps < 1L) stop("n_reps must be a positive integer.")

  Y_int <- check_binary_y(internal_data[[Y]])
  n <- length(Y_int)
  n0 <- sum(Y_int == 0)
  n1 <- sum(Y_int == 1)
  if (n0 <= 1L || n1 <= 1L) {
    stop("Cross-fitting requires at least two cases and two controls.")
  }
  if (K > min(n0, n1)) {
    stop("K must be no larger than the smaller internal case/control stratum size.")
  }
  if (use_reference && K > nrow(reference_data)) {
    stop("K must be no larger than the number of reference controls.")
  }

  fit_full <- zim_case_control(
    internal_data = internal_data,
    Y = Y,
    M_formula = M_formula,
    working_models = working_models,
    beta_tilde = beta_tilde,
    Sigma_beta_beta = Sigma_beta_beta,
    external_n_case = external_n_case,
    external_n_control = external_n_control,
    rho_k = rho_k,
    reference_data = reference_data,
    ridge = ridge,
    fixed_case_control = fixed_case_control,
    optim_control = optim_control,
    beta_method = beta_method,
    wkm_mode = wkm_mode,
    inverse_method = inverse_method,
    pc_var = pc_var,
    n_pc = n_pc,
    eig_tol = eig_tol,
    pca_rule = pca_rule,
    pc_max_condition = pc_max_condition,
    pca_scale = pca_scale,
    shrink_structured = shrink_structured,
    shrink_lambda = shrink_lambda,
    H_B_centering = H_B_centering,
    sigma_B_centering = sigma_B_centering,
    el_mass_tol = el_mass_tol,
    el_mass_action = el_mass_action,
    compact_centering_tol = compact_centering_tol,
    joint_psd_tol = joint_psd_tol,
    wkm_rank_tol = wkm_rank_tol,
    bread_condition_warn = bread_condition_warn,
    gmm_start_method_onepop = gmm_start_method_onepop,
    gmm_start_method_reference = gmm_start_method_reference,
    gmm_H_mode = gmm_H_mode,
    gmm_solver = gmm_solver,
    gmm_score_tol = gmm_score_tol,
    gmm_rank_tol = gmm_rank_tol,
    gmm_newton_step_tol = gmm_newton_step_tol,
    gmm_newton_max_halving = gmm_newton_max_halving,
    gmm_newton_fallback_to_optim = gmm_newton_fallback_to_optim,
    gmm_package_type = gmm_package_type,
    gmm_package_wmatrix = gmm_package_wmatrix,
    return_gmm_el_pseudo_if = identical(beta_method, "GMM") && !use_reference &&
      identical(gmm_H_mode, "el")
  )

  if (is_pca_family_method(inverse_method)) {
    fit_full$inverse_diagnostics$pca_rule_inferred_from_n_pc <-
      pca_controls$inferred_from_n_pc
    fit_full$pca_rule_inferred_from_n_pc <- pca_controls$inferred_from_n_pc
  }

  d <- length(fit_full$theta_MLE)
  p_beta <- length(fit_full$beta_hat)
  Sigma_bb <- fit_full$Sigma_beta_beta

  theta_reps <- matrix(NA_real_, nrow = n_reps, ncol = d)
  colnames(theta_reps) <- names(fit_full$theta_MLE)
  V_reps <- vector("list", n_reps)
  V_internal_reps <- vector("list", n_reps)
  V_reference_reps <- vector("list", n_reps)
  V_external_reps <- vector("list", n_reps)
  W_bar_reps <- vector("list", n_reps)
  W_folds_reps <- vector("list", n_reps)
  fold_reps <- vector("list", n_reps)
  reference_fold_reps <- vector("list", n_reps)

  for (r in seq_len(n_reps)) {
    seed_r <- NULL
    if (!is.null(fold_seed)) seed_r <- as.integer(fold_seed) + r - 1L
    fold <- make_stratified_folds(Y_int, K, seed = seed_r)
    ref_fold <- NULL
    if (use_reference) {
      seed_ref <- if (is.null(seed_r)) NULL else seed_r + 100000L
      ref_fold <- make_plain_folds(nrow(reference_data), K, seed = seed_ref)
    }

    W_bar <- matrix(0, nrow = d, ncol = p_beta)
    dimnames(W_bar) <- dimnames(fit_full$W_hat)
    W_folds <- vector("list", K)

    V_reference <- matrix(0, nrow = d, ncol = d)
    dimnames(V_reference) <- dimnames(fit_full$V_W)

    if (identical(beta_method, "WKM") ||
        (identical(beta_method, "GMM") && !use_reference &&
           identical(gmm_H_mode, "sample_if"))) {
      IF_theta <- fit_full$IF_theta
      IF_beta <- fit_full$IF_beta
      if (is.null(IF_theta) || is.null(IF_beta)) {
        stop("Full-data sample-IF fit did not return per-observation IF rows; cannot cross-fit.")
      }
      g <- matrix(0, nrow = n, ncol = d)
      colnames(g) <- names(fit_full$theta_MLE)

      for (k in seq_len(K)) {
        train_idx <- which(fold != k)
        test_idx <- which(fold == k)
        if (length(test_idx) == 0L) next
        if (length(unique(Y_int[train_idx])) < 2L) {
          stop("Each training split must contain both cases and controls; reduce K.")
        }
        W_k <- crossfit_weight_on_train(
          IF_theta = IF_theta,
          IF_beta = IF_beta,
          Y = Y_int,
          train_idx = train_idx,
          Sigma_beta_beta = Sigma_bb,
          ridge = ridge,
          inverse_method = inverse_method,
          pc_var = pc_var,
          n_pc = n_pc,
          eig_tol = eig_tol,
          pca_rule = pca_rule,
          pc_max_condition = pc_max_condition,
          pca_scale = pca_scale,
          fixed_case_control = fixed_case_control,
          shrink_structured = shrink_structured,
          shrink_lambda = shrink_lambda
        )
        dimnames(W_k) <- dimnames(fit_full$W_hat)
        W_folds[[k]] <- W_k
        W_bar <- W_bar + (length(test_idx) / n) * W_k
        g[test_idx, ] <- IF_theta[test_idx, , drop = FALSE] -
          IF_beta[test_idx, , drop = FALSE] %*% t(W_k)
      }
      V_internal <- symmetrize(case_control_cov(g, Y_int, fixed_case_control))
      crossfit_method <- if (identical(beta_method, "WKM")) {
        "wkm_sample_if_out_of_fold_weight_residual_variance"
      } else {
        "gmm_onepop_sample_if_out_of_fold_weight_residual_variance"
      }

    } else if (!use_reference) {
      H_el <- fit_full$H_full
      required <- c("EL_IF0_theta", "EL_IF1_theta", "EL_IF0_beta", "EL_IF1_beta",
                    "EL_p_hat", "EL_Delta_hat", "rho_I")
      if (!all(required %in% names(H_el)) ||
          is.null(H_el$EL_IF0_theta) || is.null(H_el$EL_IF1_theta) ||
          is.null(H_el$EL_IF0_beta) || is.null(H_el$EL_IF1_beta) ||
          is.null(H_el$EL_p_hat) || is.null(H_el$EL_Delta_hat)) {
        stop("Full-data GMM fit did not return EL pseudo-influence rows; cannot cross-fit.")
      }
      EL_IF0_theta <- H_el$EL_IF0_theta
      EL_IF1_theta <- H_el$EL_IF1_theta
      EL_IF0_beta <- H_el$EL_IF0_beta
      EL_IF1_beta <- H_el$EL_IF1_beta
      p_hat <- H_el$EL_p_hat
      Delta_hat <- H_el$EL_Delta_hat
      rho_I <- H_el$rho_I
      g0 <- matrix(0, nrow = n, ncol = d)
      g1 <- matrix(0, nrow = n, ncol = d)
      colnames(g0) <- colnames(g1) <- names(fit_full$theta_MLE)
      total_mass <- el_pseudo_fold_mass(p_hat, Delta_hat, rho_I)

      for (k in seq_len(K)) {
        train_idx <- which(fold != k)
        test_idx <- which(fold == k)
        if (length(test_idx) == 0L) next
        W_k <- crossfit_weight_on_train_el_gmm(
          EL_IF0_theta = EL_IF0_theta,
          EL_IF1_theta = EL_IF1_theta,
          EL_IF0_beta = EL_IF0_beta,
          EL_IF1_beta = EL_IF1_beta,
          p_hat = p_hat,
          Delta_hat = Delta_hat,
          rho_I = rho_I,
          train_idx = train_idx,
          Sigma_beta_beta = Sigma_bb,
          ridge = ridge,
          inverse_method = inverse_method,
          pc_var = pc_var,
          n_pc = n_pc,
          eig_tol = eig_tol,
          pca_rule = pca_rule,
          pc_max_condition = pc_max_condition,
          pca_scale = pca_scale,
          shrink_lambda = shrink_lambda
        )
        dimnames(W_k) <- dimnames(fit_full$W_hat)
        W_folds[[k]] <- W_k
        fold_mass <- el_pseudo_fold_mass(p_hat, Delta_hat, rho_I, test_idx)
        W_bar <- W_bar + (fold_mass / total_mass) * W_k
        g0[test_idx, ] <- EL_IF0_theta[test_idx, , drop = FALSE] -
          EL_IF0_beta[test_idx, , drop = FALSE] %*% t(W_k)
        g1[test_idx, ] <- EL_IF1_theta[test_idx, , drop = FALSE] -
          EL_IF1_beta[test_idx, , drop = FALSE] %*% t(W_k)
      }
      V_internal <- el_pseudo_cov(g0, g1, p_hat, Delta_hat, rho_I)
      crossfit_method <- "gmm_onepop_el_pseudo_row_out_of_fold_weight_residual_variance"

    } else {
      H_ref <- fit_full$H_full
      if (is.null(H_ref$IF_theta) || is.null(H_ref$IF_beta_internal) ||
          is.null(H_ref$IF_beta_reference)) {
        stop("Full-data reference-control GMM fit did not return IF rows; cannot cross-fit.")
      }
      IF_theta <- H_ref$IF_theta
      IF_beta_internal <- H_ref$IF_beta_internal
      IF_beta_reference <- H_ref$IF_beta_reference
      m_ref <- nrow(IF_beta_reference)
      g_internal <- matrix(0, nrow = n, ncol = d)
      g_reference <- matrix(0, nrow = m_ref, ncol = d)
      colnames(g_internal) <- colnames(g_reference) <- names(fit_full$theta_MLE)

      for (k in seq_len(K)) {
        train_idx <- which(fold != k)
        test_idx <- which(fold == k)
        train_ref_idx <- which(ref_fold != k)
        test_ref_idx <- which(ref_fold == k)
        if (length(test_idx) == 0L || length(test_ref_idx) == 0L) next
        if (length(unique(Y_int[train_idx])) < 2L) {
          stop("Each internal training split must contain both cases and controls; reduce K.")
        }
        if (length(train_ref_idx) < 2L) {
          stop("Each reference training split must contain at least two controls; reduce K.")
        }
        W_k <- crossfit_weight_on_train_reference_gmm(
          IF_theta = IF_theta,
          IF_beta_internal = IF_beta_internal,
          IF_beta_reference = IF_beta_reference,
          Y = Y_int,
          train_internal_idx = train_idx,
          train_reference_idx = train_ref_idx,
          n_internal = n,
          n_reference = m_ref,
          Sigma_beta_beta = Sigma_bb,
          ridge = ridge,
          inverse_method = inverse_method,
          pc_var = pc_var,
          n_pc = n_pc,
          eig_tol = eig_tol,
          pca_rule = pca_rule,
          pc_max_condition = pc_max_condition,
          pca_scale = pca_scale,
          fixed_case_control = fixed_case_control,
          shrink_lambda = shrink_lambda
        )
        dimnames(W_k) <- dimnames(fit_full$W_hat)
        W_folds[[k]] <- W_k
        W_bar <- W_bar + (length(test_idx) / n) * W_k
        g_internal[test_idx, ] <- IF_theta[test_idx, , drop = FALSE] -
          IF_beta_internal[test_idx, , drop = FALSE] %*% t(W_k)
        g_reference[test_ref_idx, ] <- -IF_beta_reference[test_ref_idx, , drop = FALSE] %*% t(W_k)
      }
      V_internal_case <- symmetrize(case_control_cov(g_internal, Y_int, fixed_case_control))
      V_reference <- symmetrize(reference_centered_cov_scaled(
        g_reference,
        scale_factor = n / m_ref
      ))
      V_internal <- symmetrize(V_internal_case + V_reference)
      crossfit_method <- "gmm_reference_control_out_of_fold_weight_internal_and_reference_residual_variance"
    }

    V_external <- symmetrize(W_bar %*% Sigma_bb %*% t(W_bar))
    V_cf <- symmetrize(V_internal + V_external)
    dimnames(V_internal) <- dimnames(fit_full$V_W)
    dimnames(V_external) <- dimnames(fit_full$V_W)
    dimnames(V_cf) <- dimnames(fit_full$V_W)

    theta_r <- as.numeric(fit_full$theta_MLE -
                            W_bar %*% (fit_full$beta_hat - fit_full$beta_tilde))
    names(theta_r) <- names(fit_full$theta_MLE)

    theta_reps[r, ] <- theta_r
    V_reps[[r]] <- V_cf
    V_internal_reps[[r]] <- V_internal
    V_reference_reps[[r]] <- V_reference
    V_external_reps[[r]] <- V_external
    W_bar_reps[[r]] <- W_bar
    W_folds_reps[[r]] <- W_folds
    fold_reps[[r]] <- fold
    reference_fold_reps[r] <- list(ref_fold)
  }

  theta_ZIM <- colMeans(theta_reps)
  names(theta_ZIM) <- names(fit_full$theta_MLE)
  V_internal <- Reduce(`+`, V_internal_reps) / n_reps
  V_reference <- Reduce(`+`, V_reference_reps) / n_reps
  V_external <- Reduce(`+`, V_external_reps) / n_reps
  V_mean <- Reduce(`+`, V_reps) / n_reps
  if (n_reps > 1L) {
    theta_centered <- sweep(theta_reps, 2, theta_ZIM, "-")
    V_between_split <- symmetrize(n * crossprod(theta_centered) / (n_reps - 1L))
  } else {
    V_between_split <- matrix(0, nrow = d, ncol = d)
  }
  dimnames(V_between_split) <- dimnames(fit_full$V_W)
  V_cf <- symmetrize(V_mean + V_between_split)
  dimnames(V_cf) <- dimnames(fit_full$V_W)

  W_bar <- Reduce(`+`, W_bar_reps) / n_reps
  dimnames(W_bar) <- dimnames(fit_full$W_hat)

  se_ZIM <- sqrt(pmax(diag(V_cf) / n, 0))
  z_ZIM <- theta_ZIM / se_ZIM
  p_ZIM <- 2 * stats::pnorm(abs(z_ZIM), lower.tail = FALSE)

  list(
    theta_ZIM = theta_ZIM,
    se_ZIM = se_ZIM,
    ridge = ridge,
    el_mass_tol = el_mass_tol,
    el_mass_action = el_mass_action,
    compact_centering_tol = compact_centering_tol,
    joint_psd_tol = joint_psd_tol,
    wkm_rank_tol = wkm_rank_tol,
    bread_condition_warn = bread_condition_warn,
    z_ZIM = z_ZIM,
    p_ZIM = p_ZIM,
    theta_MLE = fit_full$theta_MLE,
    se_MLE = fit_full$se_MLE,
    V_W = V_cf,
    V_W_internal = V_internal,
    V_W_reference = V_reference,
    V_W_external = V_external,
    V_W_between_split = V_between_split,
    V_W_reps = V_reps,
    theta_ZIM_reps = theta_reps,
    W_bar = W_bar,
    W_bar_reps = W_bar_reps,
    W_folds_reps = W_folds_reps,
    fold_reps = fold_reps,
    reference_fold_reps = reference_fold_reps,
    fold = fold_reps[[1L]],
    reference_fold = reference_fold_reps[[1L]],
    n_reps = n_reps,
    K = K,
    crossfit = TRUE,
    crossfit_method = crossfit_method,
    beta_method = beta_method,
    gmm_H_mode = if (identical(beta_method, "GMM") && !use_reference) gmm_H_mode else NA_character_,
    gmm_solver = if (identical(beta_method, "GMM")) gmm_solver else NA_character_,
    gmm_score_tol = if (identical(beta_method, "GMM")) gmm_score_tol else NA_real_,
    gmm_rank_tol = if (identical(beta_method, "GMM") && !use_reference) gmm_rank_tol else NA_real_,
    gmm_package_type = if (identical(beta_method, "GMM") && identical(gmm_solver, "gmm")) gmm_package_type else NA_character_,
    gmm_package_wmatrix = if (identical(beta_method, "GMM") && identical(gmm_solver, "gmm")) gmm_package_wmatrix else NA_character_,
    wkm_mode = if (identical(beta_method, "WKM")) wkm_mode else NA_character_,
    inverse_method = inverse_method,
    shrink_lambda = if (identical(inverse_method, "shrinkage")) shrink_lambda else NA_real_,
    pca_rule = if (is_pca_family_method(inverse_method)) pca_rule else NA_character_,
    pca_scale = if (is_pca_family_method(inverse_method)) pca_scale else NA_character_,
    pca_rule_inferred_from_n_pc = if (is_pca_family_method(inverse_method)) {
      pca_controls$inferred_from_n_pc
    } else FALSE,
    pc_var = if (is_pca_family_method(inverse_method) &&
                 pca_rule %in% c("variance", "condition_or_variance")) pc_var else NA_real_,
    pc_max_condition = if (is_pca_family_method(inverse_method) &&
                           pca_rule %in% c("condition", "condition_or_variance")) pc_max_condition else NA_real_,
    n_pc = if (is_pca_family_method(inverse_method) &&
               identical(pca_rule, "fixed")) n_pc else NA_integer_,
    reference_data_used = use_reference,

    theta_ZIM_plugin = fit_full$theta_ZIM,
    se_ZIM_plugin = fit_full$se_ZIM,
    W_hat_plugin = fit_full$W_hat,
    V_W_plugin = fit_full$V_W,
    fit_full = fit_full,
    call = match.call()
  )
}
