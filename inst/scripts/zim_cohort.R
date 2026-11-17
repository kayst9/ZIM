`%||%` <- function(x, y) if (is.null(x)) y else x

cohort_symmetrize <- function(A) {
  A <- as.matrix(A)
  if (nrow(A) != ncol(A)) stop("A must be square.")
  if (length(A) == 0L) return(A)
  (A + t(A)) / 2
}

cohort_matrix_scale <- function(A) {
  A <- as.matrix(A)
  if (nrow(A) == 0L) return(1)
  max(mean(abs(diag(A))), 1)
}

cohort_regularize_covariance <- function(A, ridge = 0) {
  A <- cohort_symmetrize(A)
  if (!is.numeric(ridge) || length(ridge) != 1L || !is.finite(ridge) || ridge < 0) {
    stop("ridge must be a finite nonnegative scalar.")
  }
  if (nrow(A) == 0L || ridge == 0) return(A)
  A + diag(ridge * cohort_matrix_scale(A), nrow(A))
}

cohort_regularize_negative_bread <- function(A, ridge = 0) {
  A <- cohort_symmetrize(A)
  if (!is.numeric(ridge) || length(ridge) != 1L || !is.finite(ridge) || ridge < 0) {
    stop("ridge must be a finite nonnegative scalar.")
  }
  if (nrow(A) == 0L || ridge == 0) return(A)
  A - diag(ridge * cohort_matrix_scale(A), nrow(A))
}

cohort_right_solve <- function(A, X) {
  A <- as.matrix(A)
  X <- as.matrix(X)
  if (nrow(A) != ncol(A) || ncol(X) != nrow(A)) {
    stop("For X %*% solve(A), A must be square and ncol(X) must equal nrow(A).")
  }
  if (nrow(A) == 0L) return(X)
  t(solve(t(A), t(X)))
}

cohort_right_solve_transpose <- function(A, X) {
  A <- as.matrix(A)
  X <- as.matrix(X)
  if (nrow(A) != ncol(A) || ncol(X) != nrow(A)) {
    stop("For X %*% solve(t(A)), A must be square and ncol(X) must equal nrow(A).")
  }
  if (nrow(A) == 0L) return(X)
  t(solve(A, t(X)))
}

cohort_center_rows <- function(X) {
  X <- as.matrix(X)
  if (ncol(X) == 0L) return(X)
  sweep(X, 2L, colMeans(X), FUN = "-")
}

cohort_crosscov <- function(A, B) {
  A <- as.matrix(A)
  B <- as.matrix(B)
  if (nrow(A) != nrow(B)) stop("A and B must have the same number of rows.")
  if (nrow(A) == 0L) stop("Covariance requires at least one row.")
  crossprod(cohort_center_rows(A), cohort_center_rows(B)) / nrow(A)
}

cohort_cov <- function(A) {
  cohort_symmetrize(cohort_crosscov(A, A))
}

cohort_matrix_diagnostics <- function(A, tol = 1e-10) {
  A <- as.matrix(A)
  if (nrow(A) != ncol(A)) stop("Matrix diagnostics require a square matrix.")
  p <- nrow(A)
  if (p == 0L) {
    return(list(
      dimension = 0L, rank = 0L, full_rank = TRUE,
      singular_values = numeric(0), cutoff = 0,
      condition_number = 1, min_singular_value = Inf,
      max_singular_value = 0
    ))
  }
  s <- svd(A, nu = 0L, nv = 0L)$d
  smax <- max(s)
  cutoff <- tol * max(smax, .Machine$double.eps)
  rank <- sum(s > cutoff)
  smin <- min(s)
  condition <- if (smin > 0) smax / smin else Inf
  list(
    dimension = p,
    rank = rank,
    full_rank = rank == p,
    singular_values = s,
    cutoff = cutoff,
    condition_number = condition,
    min_singular_value = smin,
    max_singular_value = smax
  )
}

cohort_assert_full_rank <- function(A, label, rank_tol = 1e-10,
                                    condition_warn = 1e12) {
  d <- cohort_matrix_diagnostics(A, rank_tol)
  if (!d$full_rank) {
    stop(
      label, " is rank deficient: rank ", d$rank,
      " of ", d$dimension, "."
    )
  }
  if (is.finite(d$condition_number) && d$condition_number > condition_warn) {
    warning(
      label, " is poorly conditioned (condition number ",
      format(d$condition_number, scientific = TRUE), ").",
      call. = FALSE
    )
  }
  d
}

cohort_psd_diagnostics <- function(A, tol = 1e-8, label = "matrix",
                                   action = c("stop", "warn", "none")) {
  action <- match.arg(action)
  A <- cohort_symmetrize(A)
  if (nrow(A) == 0L) {
    return(list(
      valid = TRUE, min_eigenvalue = Inf, max_eigenvalue = 0,
      cutoff = 0, eigenvalues = numeric(0), dimension = 0L
    ))
  }
  values <- eigen(A, symmetric = TRUE, only.values = TRUE)$values
  cutoff <- tol * max(1, max(abs(values)))
  valid <- min(values) >= -cutoff
  if (!valid) {
    msg <- paste0(
      label, " is not positive semidefinite within tolerance: minimum eigenvalue = ",
      format(min(values), scientific = TRUE), ", allowed cutoff = -",
      format(cutoff, scientific = TRUE), "."
    )
    if (action == "stop") stop(msg)
    if (action == "warn") warning(msg, call. = FALSE)
  }
  list(
    valid = valid,
    min_eigenvalue = min(values),
    max_eigenvalue = max(values),
    cutoff = cutoff,
    eigenvalues = values,
    dimension = nrow(A)
  )
}

cohort_check_finite_matrix <- function(X, label) {
  X <- as.matrix(X)
  if (any(!is.finite(X))) stop(label, " contains NA, NaN, Inf, or -Inf.")
  storage.mode(X) <- "double"
  X
}

if (!exists("symmetrize", mode = "function")) symmetrize <- cohort_symmetrize
if (!exists("right_solve", mode = "function")) {
  right_solve <- function(A, X, ridge = 0) {
    if (ridge > 0) A <- cohort_regularize_negative_bread(A, ridge)
    cohort_right_solve(A, X)
  }
}
if (!exists("right_solve_transpose", mode = "function")) {
  right_solve_transpose <- function(A, X, ridge = 0) {
    if (ridge > 0) A <- cohort_regularize_negative_bread(A, ridge)
    cohort_right_solve_transpose(A, X)
  }
}
if (!exists("regularize_matrix", mode = "function")) {
  regularize_matrix <- cohort_regularize_covariance
}

cohort_normalize_working_models <- function(working_models) {
  if (inherits(working_models, "formula")) working_models <- list(working_models)
  if (!is.list(working_models) || length(working_models) == 0L) {
    stop("working_models must be a nonempty list of formulas or alpha/beta lists.")
  }
  working_models
}

cohort_working_formulas <- function(working_models) {
  working_models <- cohort_normalize_working_models(working_models)
  out <- list()
  for (k in seq_along(working_models)) {
    spec <- working_models[[k]]
    if (inherits(spec, "formula")) {
      out[[length(out) + 1L]] <- spec
    } else if (is.list(spec) && inherits(spec$beta, "formula")) {
      out[[length(out) + 1L]] <- spec$alpha %||% stats::as.formula("~ 1")
      out[[length(out) + 1L]] <- spec$beta
    } else {
      stop(
        "working_models[[", k,
        "]] must be a formula or a list with beta = <formula> and optional alpha = <formula>."
      )
    }
  }
  out
}

cohort_rhs_variables <- function(formula, data) {
  if (!inherits(formula, "formula")) stop("All model specifications must be formulas.")
  trm <- stats::delete.response(stats::terms(formula, data = data))
  vars <- all.vars(stats::formula(trm))
  if (grepl("(^|[^[:alnum:]_])\\.([^[:alnum:]_]|$)", paste(deparse(formula), collapse = ""))) {
    vars <- union(vars, names(data))
  }
  vars
}

cohort_collect_required_variables <- function(data, Y = NULL, M_formula = NULL,
                                              working_models = NULL) {
  vars <- character(0)
  if (!is.null(Y)) vars <- c(vars, Y)
  formulas <- list()
  if (!is.null(M_formula)) formulas[[length(formulas) + 1L]] <- M_formula
  if (!is.null(working_models)) formulas <- c(formulas, cohort_working_formulas(working_models))
  for (f in formulas) vars <- c(vars, cohort_rhs_variables(f, data))
  vars <- unique(vars)
  missing <- setdiff(vars, names(data))
  if (length(missing)) {
    stop("Missing required variables in data: ", paste(missing, collapse = ", "), ".")
  }
  vars
}

cohort_prepare_analysis_data <- function(data, Y = NULL, M_formula = NULL,
                                         working_models = NULL,
                                         data_label = "data") {
  if (!is.data.frame(data)) data <- as.data.frame(data)
  vars <- cohort_collect_required_variables(
    data = data, Y = Y, M_formula = M_formula, working_models = working_models
  )
  mask <- if (length(vars)) {
    stats::complete.cases(data[, vars, drop = FALSE])
  } else {
    rep(TRUE, nrow(data))
  }
  n_removed <- sum(!mask)
  if (n_removed > 0L) {
    warning(
      data_label, ": removed ", n_removed,
      " row(s) using one common complete-case mask across all requested models.",
      call. = FALSE
    )
  }
  out <- data[mask, , drop = FALSE]
  if (nrow(out) == 0L) stop(data_label, " has no complete observations.")
  list(
    data = out,
    complete_case = mask,
    n_original = nrow(data),
    n_used = nrow(out),
    n_removed = n_removed,
    variables = vars
  )
}

cohort_rhs_model_matrix <- function(formula, data, keep_intercept = FALSE,
                                    label = "model matrix") {
  if (!inherits(formula, "formula")) stop(label, " formula is invalid.")
  trm <- stats::delete.response(stats::terms(formula, data = data))
  mm <- stats::model.matrix(trm, data = data, na.action = stats::na.fail)
  if (!keep_intercept && "(Intercept)" %in% colnames(mm)) {
    mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
  }
  cohort_check_finite_matrix(mm, label)
}

if (!exists("rhs_model_matrix", mode = "function")) rhs_model_matrix <- cohort_rhs_model_matrix

cohort_pack_working_models <- function(working_models, data) {
  working_models <- cohort_normalize_working_models(working_models)
  labels <- names(working_models)
  if (is.null(labels)) labels <- rep("", length(working_models))
  empty <- is.na(labels) | labels == ""
  labels[empty] <- paste0("wm", which(empty))

  packed <- vector("list", length(working_models))
  all_beta_raw <- character(0)
  beta_lengths <- integer(length(working_models))

  for (k in seq_along(working_models)) {
    spec <- working_models[[k]]
    if (inherits(spec, "formula")) {
      alpha_formula <- stats::as.formula("~ 1")
      beta_formula <- spec
    } else {
      alpha_formula <- spec$alpha %||% stats::as.formula("~ 1")
      beta_formula <- spec$beta
      if (!inherits(alpha_formula, "formula") || !inherits(beta_formula, "formula")) {
        stop("working_models[[", k, "]] has invalid alpha or beta formula.")
      }
    }

    A_k <- cohort_rhs_model_matrix(
      alpha_formula, data, keep_intercept = TRUE,
      label = paste0("working model ", k, " alpha design")
    )
    B_k <- cohort_rhs_model_matrix(
      beta_formula, data, keep_intercept = FALSE,
      label = paste0("working model ", k, " beta design")
    )
    if (ncol(B_k) == 0L) stop("working model ", k, " contributes no reported beta columns.")

    D_k <- cohort_check_finite_matrix(
      cbind(A_k, B_k), paste0("working model ", k, " design")
    )
    alpha_names <- colnames(A_k)
    beta_raw <- colnames(B_k)
    gamma_names <- make.unique(c(alpha_names, beta_raw), sep = "__")
    colnames(D_k) <- gamma_names

    packed[[k]] <- list(
      D_k = D_k,
      A_k = A_k,
      B_k = B_k,
      q_alpha_k = ncol(A_k),
      p_beta_k = ncol(B_k),
      alpha_names = alpha_names,
      beta_names_raw = beta_raw,
      beta_names = beta_raw,
      gamma_names = gamma_names,
      alpha_formula = alpha_formula,
      beta_formula = beta_formula,
      model_label = labels[k],
      model_index = k
    )
    all_beta_raw <- c(all_beta_raw, beta_raw)
    beta_lengths[k] <- length(beta_raw)
  }

  duplicated_anywhere <- duplicated(all_beta_raw) | duplicated(all_beta_raw, fromLast = TRUE)
  all_beta_names <- all_beta_raw
  if (any(duplicated_anywhere)) {
    offset <- 0L
    for (k in seq_along(packed)) {
      idx <- offset + seq_len(beta_lengths[k])
      use_prefix <- duplicated_anywhere[idx]
      all_beta_names[idx[use_prefix]] <- paste0(
        packed[[k]]$model_label, "::", all_beta_raw[idx[use_prefix]]
      )
      offset <- offset + beta_lengths[k]
    }
  }
  all_beta_names <- make.unique(all_beta_names, sep = "__")

  offset <- 0L
  for (k in seq_along(packed)) {
    idx <- offset + seq_len(beta_lengths[k])
    packed[[k]]$beta_names <- all_beta_names[idx]
    offset <- offset + beta_lengths[k]
  }
  attr(packed, "beta_names") <- all_beta_names
  packed
}

if (!exists("pack_working_models", mode = "function")) pack_working_models <- cohort_pack_working_models

cohort_align_beta_vector <- function(x, beta_names, object_name = "beta_tilde") {
  original_names <- names(x)
  x_num <- as.numeric(x)
  if (is.null(original_names)) {
    if (length(x_num) != length(beta_names)) {
      stop(object_name, " has length ", length(x_num), "; expected ", length(beta_names), ".")
    }
    warning(object_name, " is unnamed; using positional alignment.", call. = FALSE)
    return(stats::setNames(x_num, beta_names))
  }
  if (anyDuplicated(original_names)) stop(object_name, " contains duplicated names.")
  missing_names <- setdiff(beta_names, original_names)
  extra_names <- setdiff(original_names, beta_names)
  if (length(missing_names) || length(extra_names)) {
    stop(
      object_name, " names do not match the internal beta names. Missing: ",
      paste(missing_names, collapse = ", "), "; extra: ",
      paste(extra_names, collapse = ", "), "."
    )
  }
  stats::setNames(x_num[match(beta_names, original_names)], beta_names)
}

cohort_align_beta_covariance <- function(Sigma, beta_names,
                                         object_name = "Sigma_beta_beta",
                                         psd_tol = 1e-8) {
  Sigma <- as.matrix(Sigma)
  p <- length(beta_names)
  if (!all(dim(Sigma) == c(p, p))) {
    stop(object_name, " must be a ", p, " by ", p, " matrix.")
  }
  rn <- rownames(Sigma)
  cn <- colnames(Sigma)
  if (is.null(rn) && is.null(cn)) {
    warning(object_name, " is unnamed; using positional alignment.", call. = FALSE)
  } else {
    if (is.null(rn) || is.null(cn)) {
      stop(object_name, " must have both row and column names, or neither.")
    }
    if (anyDuplicated(rn) || anyDuplicated(cn)) stop(object_name, " has duplicated dimnames.")
    if (!setequal(rn, beta_names) || !setequal(cn, beta_names)) {
      stop(object_name, " dimnames do not match the internal beta names.")
    }
    Sigma <- Sigma[beta_names, beta_names, drop = FALSE]
  }
  if (any(!is.finite(Sigma))) stop(object_name, " contains nonfinite values.")
  asym <- max(abs(Sigma - t(Sigma)))
  scale <- max(1, max(abs(Sigma)))
  if (asym > psd_tol * scale) {
    stop(object_name, " is materially asymmetric (max absolute asymmetry = ", asym, ").")
  }
  Sigma <- cohort_symmetrize(Sigma)
  dimnames(Sigma) <- list(beta_names, beta_names)
  psd <- cohort_psd_diagnostics(Sigma, psd_tol, object_name, action = "stop")
  list(matrix = Sigma, diagnostics = c(psd, list(max_abs_asymmetry = asym)))
}

cohort_align_weight_matrix <- function(W, theta_names, beta_names) {
  W <- as.matrix(W)
  if (!all(dim(W) == c(length(theta_names), length(beta_names)))) {
    stop(
      "User-supplied W must be a ", length(theta_names), " by ",
      length(beta_names), " matrix."
    )
  }
  rn <- rownames(W)
  cn <- colnames(W)
  if (!is.null(rn)) {
    if (anyDuplicated(rn) || !setequal(rn, theta_names)) {
      stop("User-supplied W row names do not match target coefficient names.")
    }
    W <- W[theta_names, , drop = FALSE]
  }
  if (!is.null(cn)) {
    if (anyDuplicated(cn) || !setequal(cn, beta_names)) {
      stop("User-supplied W column names do not match beta names.")
    }
    W <- W[, beta_names, drop = FALSE]
  }
  if (any(!is.finite(W))) stop("User-supplied W contains nonfinite values.")
  dimnames(W) <- list(theta_names, beta_names)
  W
}

normalize_outcome_type <- function(outcome_type) {
  match.arg(outcome_type, c("binary", "continuous", "poisson"))
}

family_from_outcome_type <- function(outcome_type) {
  outcome_type <- normalize_outcome_type(outcome_type)
  switch(
    outcome_type,
    binary = stats::binomial(link = "logit"),
    continuous = stats::gaussian(link = "identity"),
    poisson = stats::poisson(link = "log")
  )
}

as_glm_family <- function(family) {
  if (is.function(family)) family <- family()
  required <- c("family", "link", "linkinv", "mu.eta", "variance")
  if (!is.list(family) || !all(required %in% names(family))) {
    stop("family must be a GLM family object, such as binomial(), gaussian(), or poisson().")
  }
  family
}

validate_canonical_cohort_family <- function(family) {
  family <- as_glm_family(family)
  ok <- (family$family == "binomial" && family$link == "logit") ||
    (family$family == "gaussian" && family$link == "identity") ||
    (family$family == "poisson" && family$link == "log")
  if (!ok) {
    stop(
      "Only canonical-link cohort GLMs are supported: ",
      "binomial(logit), gaussian(identity), and poisson(log)."
    )
  }
  family
}

outcome_type_from_family <- function(family) {
  if (family$family == "binomial") return("binary")
  if (family$family == "gaussian") return("continuous")
  if (family$family == "poisson") return("poisson")
  stop("Unsupported family.")
}

resolve_cohort_family <- function(family = NULL,
                                  outcome_type = c("binary", "continuous", "poisson")) {
  if (is.null(family)) {
    outcome_type <- normalize_outcome_type(outcome_type)
    family <- family_from_outcome_type(outcome_type)
  } else {
    family <- validate_canonical_cohort_family(family)
    outcome_type <- outcome_type_from_family(family)
  }
  list(family = family, outcome_type = outcome_type)
}

cohort_check_binary_y <- function(y) {
  if (is.logical(y)) y <- as.integer(y)
  if (is.factor(y)) {
    lev <- levels(y)
    if (length(lev) != 2L) stop("Binary outcome factor must have exactly two levels.")
    if (identical(lev, c("0", "1"))) {
      y <- as.integer(as.character(y))
    } else {
      warning(
        "Mapping the first binary factor level to 0 and the second to 1: ",
        paste(lev, collapse = " -> "), ".",
        call. = FALSE
      )
      y <- as.integer(y) - 1L
    }
  }
  y <- as.numeric(y)
  if (any(!is.finite(y)) || any(!(y %in% c(0, 1)))) {
    stop("Binary outcomes must be finite and coded 0/1.")
  }
  if (length(unique(y)) != 2L) stop("Both binary outcome classes must be present.")
  y
}

if (!exists("check_binary_y", mode = "function")) check_binary_y <- cohort_check_binary_y

check_outcome <- function(y, family) {
  if (family$family == "binomial") return(cohort_check_binary_y(y))
  y <- as.numeric(y)
  if (any(!is.finite(y))) stop(family$family, " outcomes must be finite.")
  if (family$family == "poisson") {
    if (any(y < 0)) stop("Poisson outcomes must be nonnegative.")
    if (any(abs(y - round(y)) > 1e-8)) {
      warning(
        "Poisson outcomes are usually integer counts; proceeding with the supplied non-integer values.",
        call. = FALSE
      )
    }
  }
  y
}

cohort_linkinv <- function(eta, family, poisson_clip = 35) {
  eta <- as.numeric(eta)
  if (family$family == "binomial") return(stats::plogis(eta))
  if (family$family == "gaussian") return(eta)
  if (family$family == "poisson") {
    return(exp(pmin(pmax(eta, -poisson_clip), poisson_clip)))
  }
  as.numeric(family$linkinv(eta))
}

cohort_mu_eta <- function(eta, family, poisson_clip = 35) {
  eta <- as.numeric(eta)
  if (family$family == "binomial") {
    mu <- stats::plogis(eta)
    return(mu * (1 - mu))
  }
  if (family$family == "gaussian") return(rep(1, length(eta)))
  if (family$family == "poisson") {
    inside <- eta > -poisson_clip & eta < poisson_clip
    out <- numeric(length(eta))
    out[inside] <- exp(eta[inside])
    return(out)
  }
  as.numeric(family$mu.eta(eta))
}

cohort_clipping_diagnostic <- function(eta, family, poisson_clip = 35) {
  eta <- as.numeric(eta)
  if (family$family != "poisson") {
    return(list(active = FALSE, n_clipped = 0L, max_abs_eta = max(abs(eta))))
  }
  active <- abs(eta) >= poisson_clip
  list(
    active = any(active),
    n_clipped = sum(active),
    max_abs_eta = max(abs(eta)),
    clip_limit = poisson_clip
  )
}

cohort_fit_target <- function(Y, X_no_intercept, family,
                              ridge = 0,
                              rank_tol = 1e-10,
                              bread_condition_warn = 1e12,
                              poisson_clip = 35) {
  X_no_intercept <- cohort_check_finite_matrix(X_no_intercept, "target slope design")
  Z <- cbind("(Intercept)" = 1, X_no_intercept)
  fit <- stats::glm.fit(x = Z, y = Y, family = family)

  converged <- if (family$family %in% c("binomial", "poisson")) {
    isTRUE(fit$converged)
  } else {
    TRUE
  }
  if (!converged) stop("Target GLM did not converge.")

  lambda_hat <- as.numeric(stats::coef(fit))
  if (length(lambda_hat) != ncol(Z) || any(!is.finite(lambda_hat))) {
    stop("Target model produced nonfinite or aliased coefficients.")
  }
  names(lambda_hat) <- colnames(Z)

  eta_hat <- drop(Z %*% lambda_hat)
  clipping <- cohort_clipping_diagnostic(eta_hat, family, poisson_clip)
  if (isTRUE(clipping$active)) stop("Target Poisson fit requires final clipping.")
  mu_hat <- cohort_linkinv(eta_hat, family, poisson_clip)
  mu_eta_hat <- cohort_mu_eta(eta_hat, family, poisson_clip)

  score_rows <- Z * as.numeric(Y - mu_hat)
  bread_raw <- -crossprod(Z, Z * as.numeric(mu_eta_hat)) / length(Y)
  bread_diag <- cohort_assert_full_rank(
    bread_raw, "Target raw bread", rank_tol, bread_condition_warn
  )
  bread_used <- cohort_regularize_negative_bread(bread_raw, ridge)
  IF_lambda <- -cohort_right_solve_transpose(bread_used, score_rows)
  colnames(IF_lambda) <- colnames(Z)
  IF_theta <- IF_lambda[, -1, drop = FALSE]
  theta_hat <- lambda_hat[-1]
  names(theta_hat) <- colnames(X_no_intercept)
  colnames(IF_theta) <- names(theta_hat)

  boundary <- isTRUE(fit$boundary)
  if (boundary) {
    warning("Target GLM reported a boundary fit.", call. = FALSE)
  }

  list(
    fit = fit,
    Z = Z,
    lambda_hat = lambda_hat,
    theta_hat = theta_hat,
    eta_hat = eta_hat,
    mu_hat = mu_hat,
    mu_eta_hat = mu_eta_hat,
    score_rows = score_rows,
    bread_raw = bread_raw,
    bread_used = bread_used,
    IF_lambda = IF_lambda,
    IF_theta = IF_theta,
    family = family,
    diagnostics = list(
      converged = converged,
      boundary = boundary,
      max_abs_score = max(abs(colMeans(score_rows))),
      score_norm = sqrt(sum(colMeans(score_rows)^2)),
      fitted_mean_range = range(mu_hat),
      clipping = clipping,
      bread = bread_diag,
      residual_scale = if (family$family == "gaussian") {
        residual <- Y - mu_hat
        sqrt(sum(residual^2) / max(1, length(Y) - ncol(Z)))
      } else {
        NA_real_
      }
    )
  )
}

fit_target_cohort <- function(Y, X_no_intercept, family,
                              ridge = 0,
                              rank_tol = 1e-10,
                              bread_condition_warn = 1e12,
                              poisson_clip = 35) {
  cohort_fit_target(
    Y, X_no_intercept, family, ridge, rank_tol,
    bread_condition_warn, poisson_clip
  )
}

make_reference_target_cohort <- function(reference_data, M_formula,
                                         lambda_hat, family,
                                         poisson_clip = 35) {
  X_ref <- cohort_rhs_model_matrix(
    M_formula, reference_data, keep_intercept = FALSE,
    label = "reference target slope design"
  )
  if (ncol(X_ref) != length(lambda_hat) - 1L) {
    stop("reference_data target model matrix has the wrong number of columns.")
  }
  Z_ref <- cbind("(Intercept)" = 1, X_ref)
  if (!identical(colnames(Z_ref), names(lambda_hat))) {

    if (!identical(colnames(Z_ref)[-1], names(lambda_hat)[-1])) {
      stop("reference_data target columns do not match the internal target design.")
    }
  }
  eta_ref <- drop(Z_ref %*% as.numeric(lambda_hat))
  clipping <- cohort_clipping_diagnostic(eta_ref, family, poisson_clip)
  if (isTRUE(clipping$active)) stop("Reference target evaluation requires Poisson clipping.")
  mu_ref <- cohort_linkinv(eta_ref, family, poisson_clip)
  mu_eta_ref <- cohort_mu_eta(eta_ref, family, poisson_clip)

  list(
    Z = Z_ref,
    lambda_hat = lambda_hat,
    theta_hat = lambda_hat[-1],
    eta_hat = eta_ref,
    mu_hat = mu_ref,
    mu_eta_hat = mu_eta_ref,
    family = family,
    clipping = clipping
  )
}

score_glm_cohort <- function(Y, D, gamma, family, poisson_clip = 35) {
  eta <- drop(D %*% gamma)
  mu <- cohort_linkinv(eta, family, poisson_clip)
  D * as.numeric(Y - mu)
}

sensitivity_glm_cohort <- function(D, gamma, family, poisson_clip = 35) {
  eta <- drop(D %*% gamma)
  v <- cohort_mu_eta(eta, family, poisson_clip)
  -crossprod(D, D * as.numeric(v)) / nrow(D)
}

target_influence_cohort <- function(Y, target, family,
                                    ridge = 0,
                                    rank_tol = 1e-10,
                                    bread_condition_warn = 1e12,
                                    poisson_clip = 35) {
  if (!is.null(target$IF_lambda) && !is.null(target$bread_raw)) {
    return(list(
      IF_lambda = target$IF_lambda,
      IF_theta = target$IF_theta,
      A_lambda_lambda = target$bread_raw,
      A_lambda_lambda_used = target$bread_used,
      diagnostics = target$diagnostics$bread
    ))
  }
  Z <- target$Z
  score_rows <- Z * as.numeric(Y - target$mu_hat)
  v <- cohort_mu_eta(target$eta_hat, family, poisson_clip)
  bread_raw <- -crossprod(Z, Z * as.numeric(v)) / length(Y)
  diagnostics <- cohort_assert_full_rank(
    bread_raw, "Target raw bread", rank_tol, bread_condition_warn
  )
  bread_used <- cohort_regularize_negative_bread(bread_raw, ridge)
  IF_lambda <- -cohort_right_solve_transpose(bread_used, score_rows)
  colnames(IF_lambda) <- colnames(Z)
  list(
    IF_lambda = IF_lambda,
    IF_theta = IF_lambda[, -1, drop = FALSE],
    A_lambda_lambda = bread_raw,
    A_lambda_lambda_used = bread_used,
    diagnostics = diagnostics
  )
}

psi2_cohort <- function(gamma_k, wk, target, family, poisson_clip = 35) {
  D_k <- wk$D_k
  eta_k <- drop(D_k %*% gamma_k)
  mu_k <- cohort_linkinv(eta_k, family, poisson_clip)
  D_k * as.numeric(target$mu_hat - mu_k)
}

cohort_gmm_eval <- function(gamma, D, target_mu, family, poisson_clip = 35) {
  eta <- drop(D %*% gamma)
  mu <- cohort_linkinv(eta, family, poisson_clip)
  mu_eta <- cohort_mu_eta(eta, family, poisson_clip)
  residual <- as.numeric(target_mu - mu)
  rows <- D * residual
  moment <- colMeans(rows)
  jacobian <- -crossprod(D, D * as.numeric(mu_eta)) / nrow(D)
  list(
    eta = eta,
    mu = mu,
    mu_eta = mu_eta,
    residual = residual,
    rows = rows,
    moment = moment,
    jacobian = jacobian,
    objective = sum(moment^2),
    clipping = cohort_clipping_diagnostic(eta, family, poisson_clip)
  )
}

fit_glm_start_cohort <- function(Y, wk, family, beta_tilde_k = NULL,
                                 poisson_clip = 35) {
  gamma_start <- tryCatch({
    if (family$family == "gaussian") {
      as.numeric(stats::coef(stats::lm.fit(x = wk$D_k, y = Y)))
    } else {
      as.numeric(stats::coef(stats::glm.fit(x = wk$D_k, y = Y, family = family)))
    }
  }, error = function(e) rep(NA_real_, ncol(wk$D_k)))
  if (length(gamma_start) != ncol(wk$D_k) || any(!is.finite(gamma_start))) {
    gamma_start <- rep(0, ncol(wk$D_k))
  }
  if (!is.null(beta_tilde_k)) {
    beta_idx <- wk$q_alpha_k + seq_len(wk$p_beta_k)
    gamma_start[beta_idx] <- as.numeric(beta_tilde_k)
  }
  names(gamma_start) <- wk$gamma_names
  gamma_start
}

cohort_solve_gmm_newton <- function(D, target_mu, family, gamma_start,
                                    ridge = 0,
                                    score_tol = 1e-8,
                                    step_tol = 1e-10,
                                    maxit = 100L,
                                    max_halving = 30L,
                                    poisson_clip = 35) {
  gamma <- as.numeric(gamma_start)
  last <- cohort_gmm_eval(gamma, D, target_mu, family, poisson_clip)
  converged <- max(abs(last$moment)) <= score_tol
  accepted_steps <- 0L
  total_halving <- 0L
  stop_reason <- if (converged) "initial_value_satisfies_tolerance" else "maxit"

  if (!converged) {
    for (iter in seq_len(as.integer(maxit))) {
      J_used <- cohort_regularize_negative_bread(last$jacobian, ridge)
      step <- tryCatch(
        as.numeric(solve(J_used, last$moment)),
        error = function(e) rep(NA_real_, length(gamma))
      )
      if (any(!is.finite(step))) {
        stop_reason <- "nonfinite_newton_step"
        break
      }

      step_scale <- 1
      accepted <- FALSE
      candidate <- last
      for (h in 0:as.integer(max_halving)) {
        gamma_try <- gamma - step_scale * step
        candidate <- cohort_gmm_eval(gamma_try, D, target_mu, family, poisson_clip)
        if (is.finite(candidate$objective) &&
            candidate$objective <= last$objective * (1 + 1e-12)) {
          accepted <- TRUE
          gamma <- gamma_try
          total_halving <- total_halving + h
          accepted_steps <- accepted_steps + 1L
          break
        }
        step_scale <- step_scale / 2
      }
      if (!accepted) {
        stop_reason <- "line_search_failed"
        break
      }

      last <- candidate
      max_moment <- max(abs(last$moment))
      if (max_moment <= score_tol) {
        converged <- TRUE
        stop_reason <- "score_tolerance"
        break
      }
      if (max(abs(step_scale * step)) <= step_tol) {
        stop_reason <- "step_tolerance_before_score_tolerance"
        break
      }
    }
  }

  list(
    gamma_hat = gamma,
    eval = last,
    convergence = if (converged) 0L else 1L,
    converged = converged,
    stop_reason = stop_reason,
    accepted_steps = accepted_steps,
    total_halving = total_halving,
    solver = "newton"
  )
}

cohort_solve_gmm_optim <- function(D, target_mu, family, gamma_start,
                                   optim_control = list(maxit = 1000),
                                   poisson_clip = 35) {
  objective <- function(gamma) {
    cohort_gmm_eval(gamma, D, target_mu, family, poisson_clip)$objective
  }
  gradient <- function(gamma) {
    ev <- cohort_gmm_eval(gamma, D, target_mu, family, poisson_clip)
    as.numeric(2 * crossprod(ev$jacobian, ev$moment))
  }
  opt <- stats::optim(
    par = as.numeric(gamma_start),
    fn = objective,
    gr = gradient,
    method = "BFGS",
    control = optim_control
  )
  ev <- cohort_gmm_eval(opt$par, D, target_mu, family, poisson_clip)
  list(
    gamma_hat = opt$par,
    eval = ev,
    convergence = opt$convergence,
    converged = opt$convergence == 0L,
    stop_reason = opt$message %||% "optim_completed",
    accepted_steps = NA_integer_,
    total_halving = NA_integer_,
    solver = "optim",
    optim = opt
  )
}

cohort_solve_gmm <- function(D, target_mu, family, gamma_start,
                             solver = c("newton", "optim"),
                             ridge = 0,
                             score_tol = 1e-8,
                             step_tol = 1e-10,
                             maxit = 100L,
                             max_halving = 30L,
                             fallback_to_optim = TRUE,
                             optim_control = list(maxit = 1000),
                             rank_tol = 1e-10,
                             bread_condition_warn = 1e12,
                             poisson_clip = 35,
                             label = "GMM model") {
  solver <- match.arg(solver)

  if (family$family == "gaussian") {
    fit <- stats::lm.fit(x = D, y = target_mu)
    gamma <- as.numeric(stats::coef(fit))
    if (length(gamma) != ncol(D) || any(!is.finite(gamma))) {
      stop(label, " Gaussian closed-form fit produced nonfinite or aliased coefficients.")
    }
    result <- list(
      gamma_hat = gamma,
      eval = cohort_gmm_eval(gamma, D, target_mu, family, poisson_clip),
      convergence = 0L,
      converged = TRUE,
      stop_reason = "gaussian_closed_form",
      accepted_steps = 1L,
      total_halving = 0L,
      solver = "gaussian_closed_form"
    )
  } else if (solver == "newton") {
    result <- cohort_solve_gmm_newton(
      D, target_mu, family, gamma_start,
      ridge = ridge,
      score_tol = score_tol,
      step_tol = step_tol,
      maxit = maxit,
      max_halving = max_halving,
      poisson_clip = poisson_clip
    )
    if (!result$converged && isTRUE(fallback_to_optim)) {
      optim_result <- cohort_solve_gmm_optim(
        D, target_mu, family,
        gamma_start = result$gamma_hat,
        optim_control = optim_control,
        poisson_clip = poisson_clip
      )
      newton_max <- max(abs(result$eval$moment))
      optim_max <- max(abs(optim_result$eval$moment))

      if ((optim_max <= score_tol && newton_max > score_tol) ||
          (optim_max < newton_max && !(newton_max <= score_tol))) {
        result <- optim_result
      }
    }
  } else {
    result <- cohort_solve_gmm_optim(
      D, target_mu, family, gamma_start,
      optim_control = optim_control,
      poisson_clip = poisson_clip
    )
  }

  gamma <- as.numeric(result$gamma_hat)
  if (length(gamma) != ncol(D) || any(!is.finite(gamma))) {
    stop(label, " returned nonfinite coefficients.")
  }
  ev <- cohort_gmm_eval(gamma, D, target_mu, family, poisson_clip)
  max_abs_moment <- max(abs(ev$moment))
  if (!is.finite(max_abs_moment) || max_abs_moment > score_tol) {
    stop(
      label, " did not satisfy the unregularized moment tolerance: max |moment| = ",
      format(max_abs_moment, scientific = TRUE), ", tolerance = ",
      format(score_tol, scientific = TRUE), "."
    )
  }
  jac_diag <- cohort_assert_full_rank(
    ev$jacobian, paste0(label, " raw Jacobian"),
    rank_tol, bread_condition_warn
  )
  if (isTRUE(ev$clipping$active)) {
    stop(label, " requires final Poisson clipping.")
  }

  bread_used <- cohort_regularize_negative_bread(ev$jacobian, ridge)
  names(gamma) <- colnames(D)
  list(
    gamma_hat = gamma,
    eta_hat = ev$eta,
    mu_hat = ev$mu,
    mu_eta_hat = ev$mu_eta,
    residual = ev$residual,
    moment = ev$moment,
    objective = ev$objective,
    max_abs_moment = max_abs_moment,
    moment_norm = sqrt(sum(ev$moment^2)),
    bread_raw = ev$jacobian,
    bread_used = bread_used,
    jacobian_diagnostics = jac_diag,
    clipping = ev$clipping,
    solution_validated = TRUE,
    solver = result$solver,
    solver_convergence = result$convergence,
    solver_stop_reason = result$stop_reason,
    accepted_steps = result$accepted_steps,
    total_halving = result$total_halving
  )
}

solve_gamma_gmm_cohort_k <- function(wk, target, family, gamma_start,
                                     optim_control = list(maxit = 1000),
                                     gmm_solver = c("newton", "optim"),
                                     ridge = 0,
                                     gmm_score_tol = 1e-8,
                                     gmm_step_tol = 1e-10,
                                     gmm_maxit = 100L,
                                     gmm_max_halving = 30L,
                                     gmm_newton_fallback_to_optim = TRUE,
                                     rank_tol = 1e-10,
                                     bread_condition_warn = 1e12,
                                     poisson_clip = 35) {
  sol <- cohort_solve_gmm(
    D = wk$D_k,
    target_mu = target$mu_hat,
    family = family,
    gamma_start = gamma_start,
    solver = match.arg(gmm_solver),
    ridge = ridge,
    score_tol = gmm_score_tol,
    step_tol = gmm_step_tol,
    maxit = gmm_maxit,
    max_halving = gmm_max_halving,
    fallback_to_optim = gmm_newton_fallback_to_optim,
    optim_control = optim_control,
    rank_tol = rank_tol,
    bread_condition_warn = bread_condition_warn,
    poisson_clip = poisson_clip,
    label = paste0("Internal GMM model ", wk$model_index %||% "")
  )
  list(
    gamma_hat = sol$gamma_hat,
    convergence = sol$solver_convergence,
    objective = sol$objective,
    max_abs_moment = sol$max_abs_moment,
    moment_norm = sol$moment_norm,
    solution_validated = sol$solution_validated,
    jacobian_rank = sol$jacobian_diagnostics$rank,
    jacobian_dimension = sol$jacobian_diagnostics$dimension,
    jacobian_condition_number = sol$jacobian_diagnostics$condition_number,
    solver = sol$solver,
    details = sol
  )
}

profile_working_system <- function(D_full_reg, q_alpha, C_full = NULL,
                                   D_full_raw = NULL,
                                   rank_tol = 1e-10,
                                   bread_condition_warn = 1e12,
                                   label = "working bread") {
  D_full_reg <- as.matrix(D_full_reg)
  r <- nrow(D_full_reg)
  q_alpha <- as.integer(q_alpha)
  if (q_alpha < 0L || q_alpha >= r) {
    stop(label, ": q_alpha must be positive and smaller than the full dimension.")
  }
  alpha_idx <- seq_len(q_alpha)
  beta_idx <- q_alpha + seq_len(r - q_alpha)

  make_profile <- function(D, prefix) {
    D_aa <- D[alpha_idx, alpha_idx, drop = FALSE]
    D_ab <- D[alpha_idx, beta_idx, drop = FALSE]
    D_ba <- D[beta_idx, alpha_idx, drop = FALSE]
    D_bb <- D[beta_idx, beta_idx, drop = FALSE]
    aa_diag <- cohort_assert_full_rank(
      D_aa, paste0(label, " ", prefix, " alpha block"),
      rank_tol, bread_condition_warn
    )
    T <- cohort_right_solve(D_aa, D_ba)
    S <- D_bb - T %*% D_ab
    S_diag <- cohort_assert_full_rank(
      S, paste0(label, " ", prefix, " Schur complement"),
      rank_tol, bread_condition_warn
    )
    list(T = T, S = S, alpha = aa_diag, schur = S_diag)
  }

  used <- make_profile(D_full_reg, "used")
  raw <- NULL
  if (!is.null(D_full_raw)) raw <- make_profile(as.matrix(D_full_raw), "raw")

  C_prof <- NULL
  if (!is.null(C_full)) {
    C_full <- as.matrix(C_full)
    if (nrow(C_full) != r) stop(label, ": C_full has incompatible dimensions.")
    C_alpha <- C_full[alpha_idx, , drop = FALSE]
    C_beta <- C_full[beta_idx, , drop = FALSE]
    C_prof <- C_beta - used$T %*% C_alpha
  }

  list(
    alpha_idx = alpha_idx,
    beta_idx = beta_idx,
    T_beta_alpha = used$T,
    S_beta = used$S,
    C_prof = C_prof,
    C_full = C_full,
    raw_profile_diagnostics = if (is.null(raw)) NULL else list(
      alpha = raw$alpha, schur = raw$schur
    ),
    used_profile_diagnostics = list(
      alpha = used$alpha, schur = used$schur
    )
  )
}

cohort_fit_working_wkm_k <- function(Y, wk, family,
                                     ridge = 0,
                                     rank_tol = 1e-10,
                                     bread_condition_warn = 1e12,
                                     poisson_clip = 35,
                                     model_index = wk$model_index) {
  D <- wk$D_k
  fit <- stats::glm.fit(x = D, y = Y, family = family)
  converged <- if (family$family %in% c("binomial", "poisson")) {
    isTRUE(fit$converged)
  } else {
    TRUE
  }
  if (!converged) stop("Internal WKM model ", model_index, " did not converge.")

  gamma_hat <- as.numeric(stats::coef(fit))
  if (length(gamma_hat) != ncol(D) || any(!is.finite(gamma_hat))) {
    stop("Internal WKM model ", model_index, " produced nonfinite or aliased coefficients.")
  }
  names(gamma_hat) <- wk$gamma_names

  eta_hat <- drop(D %*% gamma_hat)
  clipping <- cohort_clipping_diagnostic(eta_hat, family, poisson_clip)
  if (isTRUE(clipping$active)) {
    stop("Internal WKM Poisson model ", model_index, " requires final clipping.")
  }
  mu_hat <- cohort_linkinv(eta_hat, family, poisson_clip)
  mu_eta_hat <- cohort_mu_eta(eta_hat, family, poisson_clip)
  residual <- as.numeric(Y - mu_hat)
  score_rows <- D * residual
  score_bar <- colMeans(score_rows)

  bread_raw <- -crossprod(D, D * as.numeric(mu_eta_hat)) / length(Y)
  bread_diag <- cohort_assert_full_rank(
    bread_raw, paste0("Internal WKM model ", model_index, " raw bread"),
    rank_tol, bread_condition_warn
  )
  bread_used <- cohort_regularize_negative_bread(bread_raw, ridge)
  profile <- profile_working_system(
    D_full_reg = bread_used,
    D_full_raw = bread_raw,
    q_alpha = wk$q_alpha_k,
    C_full = NULL,
    rank_tol = rank_tol,
    bread_condition_warn = bread_condition_warn,
    label = paste0("Internal WKM model ", model_index, " bread")
  )

  beta_hat <- gamma_hat[profile$beta_idx]
  names(beta_hat) <- wk$beta_names
  boundary <- isTRUE(fit$boundary)
  if (boundary) {
    warning("Internal WKM model ", model_index, " reported a boundary fit.", call. = FALSE)
  }

  list(
    method = "WKM",
    wk = wk,
    fit = fit,
    gamma_hat = gamma_hat,
    beta_hat = beta_hat,
    eta_hat = eta_hat,
    mu_hat = mu_hat,
    mu_eta_hat = mu_eta_hat,
    residual = residual,
    score_rows = score_rows,
    bread_raw = bread_raw,
    bread_used = bread_used,
    bread_diagnostics = bread_diag,
    C_full = NULL,
    profile = profile,
    convergence = list(
      converged = converged,
      boundary = boundary,
      max_abs_score = max(abs(score_bar)),
      score_norm = sqrt(sum(score_bar^2)),
      fitted_mean_range = range(mu_hat),
      clipping = clipping,
      bread = bread_diag,
      residual_scale = if (family$family == "gaussian") {
        sqrt(sum(residual^2) / max(1, length(Y) - ncol(D)))
      } else {
        NA_real_
      }
    )
  )
}

cohort_fit_working_gmm_k <- function(Y, wk, target, family,
                                     beta_tilde_k = NULL,
                                     ridge = 0,
                                     rank_tol = 1e-10,
                                     bread_condition_warn = 1e12,
                                     gmm_solver = c("newton", "optim"),
                                     gmm_score_tol = 1e-8,
                                     gmm_step_tol = 1e-10,
                                     gmm_maxit = 100L,
                                     gmm_max_halving = 30L,
                                     gmm_newton_fallback_to_optim = TRUE,
                                     optim_control = list(maxit = 1000),
                                     poisson_clip = 35,
                                     model_index = wk$model_index) {
  gmm_solver <- match.arg(gmm_solver)
  gamma_start <- fit_glm_start_cohort(
    Y = Y,
    wk = wk,
    family = family,
    beta_tilde_k = beta_tilde_k,
    poisson_clip = poisson_clip
  )

  sol <- cohort_solve_gmm(
    D = wk$D_k,
    target_mu = target$mu_hat,
    family = family,
    gamma_start = gamma_start,
    solver = gmm_solver,
    ridge = ridge,
    score_tol = gmm_score_tol,
    step_tol = gmm_step_tol,
    maxit = gmm_maxit,
    max_halving = gmm_max_halving,
    fallback_to_optim = gmm_newton_fallback_to_optim,
    optim_control = optim_control,
    rank_tol = rank_tol,
    bread_condition_warn = bread_condition_warn,
    poisson_clip = poisson_clip,
    label = paste0("Internal GMM model ", model_index)
  )

  C_full <- crossprod(
    wk$D_k,
    target$Z * as.numeric(target$mu_eta_hat)
  ) / length(Y)
  profile <- profile_working_system(
    D_full_reg = sol$bread_used,
    D_full_raw = sol$bread_raw,
    q_alpha = wk$q_alpha_k,
    C_full = C_full,
    rank_tol = rank_tol,
    bread_condition_warn = bread_condition_warn,
    label = paste0("Internal GMM model ", model_index, " bread")
  )

  beta_hat <- sol$gamma_hat[profile$beta_idx]
  names(beta_hat) <- wk$beta_names

  list(
    method = "GMM",
    wk = wk,
    gamma_start = gamma_start,
    gamma_hat = sol$gamma_hat,
    beta_hat = beta_hat,
    eta_hat = sol$eta_hat,
    mu_hat = sol$mu_hat,
    mu_eta_hat = sol$mu_eta_hat,
    residual = sol$residual,
    moment_rows = wk$D_k * as.numeric(sol$residual),
    bread_raw = sol$bread_raw,
    bread_used = sol$bread_used,
    bread_diagnostics = sol$jacobian_diagnostics,
    C_full = C_full,
    profile = profile,
    convergence = list(
      solution_validated = TRUE,
      solver = sol$solver,
      solver_convergence = sol$solver_convergence,
      solver_stop_reason = sol$solver_stop_reason,
      max_abs_moment = sol$max_abs_moment,
      moment_norm = sol$moment_norm,
      jacobian_rank = sol$jacobian_diagnostics$rank,
      jacobian_dimension = sol$jacobian_diagnostics$dimension,
      jacobian_condition_number = sol$jacobian_diagnostics$condition_number,
      accepted_steps = sol$accepted_steps,
      total_halving = sol$total_halving,
      clipping = sol$clipping
    )
  )
}

cohort_prepare_working_states <- function(Y, target, working, family,
                                          beta_method = c("GMM", "WKM"),
                                          beta_tilde = NULL,
                                          ridge = 0,
                                          rank_tol = 1e-10,
                                          bread_condition_warn = 1e12,
                                          gmm_solver = c("newton", "optim"),
                                          gmm_score_tol = 1e-8,
                                          gmm_step_tol = 1e-10,
                                          gmm_maxit = 100L,
                                          gmm_max_halving = 30L,
                                          gmm_newton_fallback_to_optim = TRUE,
                                          optim_control = list(maxit = 1000),
                                          poisson_clip = 35) {
  beta_method <- match.arg(beta_method)
  gmm_solver <- match.arg(gmm_solver)
  states <- vector("list", length(working))
  offset <- 0L
  for (k in seq_along(working)) {
    wk <- working[[k]]
    if (beta_method == "WKM") {
      states[[k]] <- cohort_fit_working_wkm_k(
        Y = Y,
        wk = wk,
        family = family,
        ridge = ridge,
        rank_tol = rank_tol,
        bread_condition_warn = bread_condition_warn,
        poisson_clip = poisson_clip,
        model_index = k
      )
    } else {
      idx <- offset + seq_len(wk$p_beta_k)
      states[[k]] <- cohort_fit_working_gmm_k(
        Y = Y,
        wk = wk,
        target = target,
        family = family,
        beta_tilde_k = if (is.null(beta_tilde)) NULL else beta_tilde[idx],
        ridge = ridge,
        rank_tol = rank_tol,
        bread_condition_warn = bread_condition_warn,
        gmm_solver = gmm_solver,
        gmm_score_tol = gmm_score_tol,
        gmm_step_tol = gmm_step_tol,
        gmm_maxit = gmm_maxit,
        gmm_max_halving = gmm_max_halving,
        gmm_newton_fallback_to_optim = gmm_newton_fallback_to_optim,
        optim_control = optim_control,
        poisson_clip = poisson_clip,
        model_index = k
      )
    }
    offset <- offset + wk$p_beta_k
  }
  states
}

cohort_make_if_beta_columns <- function(state, local_beta_indices,
                                        target,
                                        influence_method = c("profiled", "full")) {
  influence_method <- match.arg(influence_method)
  wk <- state$wk
  local_beta_indices <- as.integer(local_beta_indices)
  if (length(local_beta_indices) == 0L ||
      any(local_beta_indices < 1L | local_beta_indices > wk$p_beta_k) ||
      anyDuplicated(local_beta_indices)) {
    stop("Invalid local beta-coordinate request for working model ", wk$model_index, ".")
  }
  m <- length(local_beta_indices)

  if (influence_method == "profiled") {
    p <- wk$p_beta_k
    E <- matrix(0, nrow = p, ncol = m)
    E[cbind(local_beta_indices, seq_len(m))] <- 1
    X <- solve(t(state$profile$S_beta), E)

    design_times_X <- wk$B_k %*% X
    if (wk$q_alpha_k > 0L) {
      design_times_X <- design_times_X -
        wk$A_k %*% (t(state$profile$T_beta_alpha) %*% X)
    }
    effective_times_X <- sweep(
      design_times_X, 1L, as.numeric(state$residual), FUN = "*"
    )
    if (state$method == "GMM") {
      target_map <- t(state$profile$C_prof) %*% X
      effective_times_X <- effective_times_X + target$IF_lambda %*% target_map
    }
    IF_part <- -effective_times_X
  } else {
    r <- ncol(wk$D_k)
    E <- matrix(0, nrow = r, ncol = m)
    full_beta_idx <- wk$q_alpha_k + local_beta_indices
    E[cbind(full_beta_idx, seq_len(m))] <- 1
    X <- solve(t(state$bread_used), E)

    effective_times_X <- sweep(
      wk$D_k %*% X, 1L, as.numeric(state$residual), FUN = "*"
    )
    if (state$method == "GMM") {
      target_map <- t(state$C_full) %*% X
      effective_times_X <- effective_times_X + target$IF_lambda %*% target_map
    }
    IF_part <- -effective_times_X
  }

  colnames(IF_part) <- wk$beta_names[local_beta_indices]
  IF_part
}

cohort_beta_index_map <- function(states) {
  rows <- vector("list", length(states))
  offset <- 0L
  for (k in seq_along(states)) {
    p <- states[[k]]$wk$p_beta_k
    rows[[k]] <- data.frame(
      global_index = offset + seq_len(p),
      model_index = k,
      local_beta_index = seq_len(p),
      beta_name = states[[k]]$wk$beta_names,
      stringsAsFactors = FALSE
    )
    offset <- offset + p
  }
  do.call(rbind, rows)
}

cohort_make_beta_blocks <- function(p_beta, block_size = 256L) {
  p_beta <- as.integer(p_beta)
  block_size <- as.integer(block_size)
  if (length(block_size) != 1L || is.na(block_size) || block_size <= 0L) {
    stop("H_block_size must be a positive integer.")
  }
  if (p_beta <= 0L) return(list())
  split(seq_len(p_beta), ceiling(seq_len(p_beta) / block_size))
}

cohort_make_if_beta_block <- function(global_indices, states, beta_map,
                                      target,
                                      influence_method = c("profiled", "full")) {
  influence_method <- match.arg(influence_method)
  global_indices <- as.integer(global_indices)
  subset_map <- beta_map[match(global_indices, beta_map$global_index), , drop = FALSE]
  if (anyNA(subset_map$model_index)) stop("Invalid global beta block indices.")
  n <- nrow(states[[1]]$wk$D_k)
  out <- matrix(NA_real_, nrow = n, ncol = length(global_indices))
  for (k in unique(subset_map$model_index)) {
    positions <- which(subset_map$model_index == k)
    local <- subset_map$local_beta_index[positions]
    out[, positions] <- cohort_make_if_beta_columns(
      state = states[[k]],
      local_beta_indices = local,
      target = target,
      influence_method = influence_method
    )
  }
  colnames(out) <- subset_map$beta_name
  out
}

cohort_make_joint_H <- function(H_theta_theta, H_theta_beta, H_beta_beta,
                                theta_names, beta_names,
                                joint_psd_tol = 1e-8,
                                label = "Joint internal H covariance") {
  d <- nrow(H_theta_theta)
  p <- nrow(H_beta_beta)
  if (d == 0L) {
    H_joint <- H_beta_beta
    dimnames(H_joint) <- list(beta_names, beta_names)
  } else if (p == 0L) {
    H_joint <- H_theta_theta
    dimnames(H_joint) <- list(theta_names, theta_names)
  } else {
    H_joint <- rbind(
      cbind(H_theta_theta, H_theta_beta),
      cbind(t(H_theta_beta), H_beta_beta)
    )
    all_names <- c(theta_names, beta_names)
    dimnames(H_joint) <- list(all_names, all_names)
  }
  psd <- cohort_psd_diagnostics(
    H_joint, tol = joint_psd_tol, label = label, action = "stop"
  )
  list(H_joint = H_joint, diagnostics = psd)
}

cohort_build_H_dense <- function(IF_theta, states, target,
                                 influence_method = c("profiled", "full"),
                                 beta_names = NULL,
                                 joint_psd_tol = 1e-8,
                                 return_if = FALSE) {
  influence_method <- match.arg(influence_method)
  beta_map <- cohort_beta_index_map(states)
  p <- nrow(beta_map)
  if (is.null(beta_names)) beta_names <- beta_map$beta_name
  IF_beta <- cohort_make_if_beta_block(
    seq_len(p), states, beta_map, target, influence_method
  )
  colnames(IF_beta) <- beta_names

  IF_theta <- as.matrix(IF_theta)
  theta_names <- colnames(IF_theta) %||% character(ncol(IF_theta))
  IF_joint_c <- cohort_center_rows(cbind(IF_theta, IF_beta))
  H_joint_raw <- cohort_symmetrize(crossprod(IF_joint_c) / nrow(IF_joint_c))
  d <- ncol(IF_theta)

  Htt <- if (d) H_joint_raw[seq_len(d), seq_len(d), drop = FALSE] else matrix(0, 0, 0)
  Htb <- if (d) {
    H_joint_raw[seq_len(d), d + seq_len(p), drop = FALSE]
  } else {
    matrix(0, 0, p)
  }
  Hbb <- H_joint_raw[d + seq_len(p), d + seq_len(p), drop = FALSE]
  dimnames(Htt) <- list(theta_names, theta_names)
  dimnames(Htb) <- list(theta_names, beta_names)
  dimnames(Hbb) <- list(beta_names, beta_names)
  joint <- cohort_make_joint_H(
    Htt, Htb, Hbb, theta_names, beta_names,
    joint_psd_tol = joint_psd_tol
  )

  list(
    H_theta_theta = Htt,
    H_theta_beta = Htb,
    H_beta_beta = Hbb,
    H_joint = joint$H_joint,
    H_joint_psd_diagnostics = joint$diagnostics,
    beta_index_map = beta_map,
    beta_blocks = list(seq_len(p)),
    H_compute_diagnostics = list(
      method = "dense",
      influence_method = influence_method,
      n_blocks = 1L,
      block_sizes = p,
      cache = "none"
    ),
    IF_beta = if (return_if) IF_beta else NULL
  )
}

cohort_build_H_blocked <- function(IF_theta, states, target,
                                   influence_method = c("profiled", "full"),
                                   H_block_size = 256L,
                                   H_block_cache = c("disk", "recompute", "memory"),
                                   H_cache_dir = tempdir(),
                                   beta_names = NULL,
                                   joint_psd_tol = 1e-8,
                                   return_if = FALSE) {
  influence_method <- match.arg(influence_method)
  H_block_cache <- match.arg(H_block_cache)
  beta_map <- cohort_beta_index_map(states)
  p <- nrow(beta_map)
  if (is.null(beta_names)) beta_names <- beta_map$beta_name
  beta_blocks <- cohort_make_beta_blocks(p, H_block_size)
  B <- length(beta_blocks)

  IF_theta <- as.matrix(IF_theta)
  n <- nrow(IF_theta)
  d <- ncol(IF_theta)
  theta_names <- colnames(IF_theta) %||% character(d)
  IF_theta_c <- cohort_center_rows(IF_theta)

  Htt <- if (d) cohort_symmetrize(crossprod(IF_theta_c) / n) else matrix(0, 0, 0)
  Htb <- matrix(0, nrow = d, ncol = p)
  Hbb <- matrix(0, nrow = p, ncol = p)
  dimnames(Htt) <- list(theta_names, theta_names)
  dimnames(Htb) <- list(theta_names, beta_names)
  dimnames(Hbb) <- list(beta_names, beta_names)
  IF_beta_return <- if (return_if) matrix(NA_real_, nrow = n, ncol = p) else NULL
  if (return_if) colnames(IF_beta_return) <- beta_names

  cache_memory <- vector("list", if (H_block_cache == "memory") B else 0L)
  cache_paths <- character(if (H_block_cache == "disk") B else 0L)
  run_dir <- NULL
  if (H_block_cache == "disk") {
    if (!dir.exists(H_cache_dir)) {
      dir.create(H_cache_dir, recursive = TRUE, showWarnings = FALSE)
    }
    run_dir <- tempfile("zim_cohort_H_blocks_", tmpdir = H_cache_dir)
    if (!dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)) {
      stop("Could not create H block cache directory: ", run_dir)
    }
    on.exit(unlink(run_dir, recursive = TRUE, force = TRUE), add = TRUE)
  }

  make_uncentered <- function(b) {
    cohort_make_if_beta_block(
      beta_blocks[[b]], states, beta_map, target, influence_method
    )
  }
  make_centered <- function(b) cohort_center_rows(make_uncentered(b))

  for (b in seq_len(B)) {
    idx_b <- beta_blocks[[b]]
    IF_b <- make_uncentered(b)
    IF_b_c <- cohort_center_rows(IF_b)
    if (d) Htb[, idx_b] <- crossprod(IF_theta_c, IF_b_c) / n
    Hbb[idx_b, idx_b] <- cohort_symmetrize(crossprod(IF_b_c) / n)
    if (return_if) IF_beta_return[, idx_b] <- IF_b

    if (H_block_cache == "memory") {
      cache_memory[[b]] <- IF_b_c
    } else if (H_block_cache == "disk") {
      path <- file.path(run_dir, sprintf("if_beta_%05d.rds", b))
      saveRDS(IF_b_c, path, compress = FALSE)
      cache_paths[b] <- path
    }
    rm(IF_b, IF_b_c)
  }

  get_centered <- function(b) {
    if (H_block_cache == "memory") return(cache_memory[[b]])
    if (H_block_cache == "disk") return(readRDS(cache_paths[b]))
    make_centered(b)
  }

  if (B > 1L) {
    for (b in 2:B) {
      idx_b <- beta_blocks[[b]]
      IF_b_c <- get_centered(b)
      for (c_block in seq_len(b - 1L)) {
        idx_c <- beta_blocks[[c_block]]
        IF_c_c <- get_centered(c_block)
        H_cb <- crossprod(IF_c_c, IF_b_c) / n
        Hbb[idx_c, idx_b] <- H_cb
        Hbb[idx_b, idx_c] <- t(H_cb)
        rm(IF_c_c)
      }
      rm(IF_b_c)
    }
  }
  Hbb <- cohort_symmetrize(Hbb)
  dimnames(Hbb) <- list(beta_names, beta_names)

  joint <- cohort_make_joint_H(
    Htt, Htb, Hbb, theta_names, beta_names,
    joint_psd_tol = joint_psd_tol
  )

  list(
    H_theta_theta = Htt,
    H_theta_beta = Htb,
    H_beta_beta = Hbb,
    H_joint = joint$H_joint,
    H_joint_psd_diagnostics = joint$diagnostics,
    beta_index_map = beta_map,
    beta_blocks = beta_blocks,
    H_compute_diagnostics = list(
      method = "blocked",
      influence_method = influence_method,
      n_blocks = B,
      block_sizes = vapply(beta_blocks, length, integer(1)),
      cache = H_block_cache,
      peak_beta_if_columns = max(vapply(beta_blocks, length, integer(1))),
      preserves_off_diagonal_blocks = TRUE
    ),
    IF_beta = IF_beta_return
  )
}

cohort_build_H <- function(IF_theta, states, target,
                           influence_method = c("profiled", "full"),
                           H_compute = c("blocked", "dense"),
                           H_block_size = 256L,
                           H_block_cache = c("disk", "recompute", "memory"),
                           H_cache_dir = tempdir(),
                           beta_names = NULL,
                           joint_psd_tol = 1e-8,
                           return_if = FALSE) {
  influence_method <- match.arg(influence_method)
  H_compute <- match.arg(H_compute)
  H_block_cache <- match.arg(H_block_cache)
  if (H_compute == "dense") {
    return(cohort_build_H_dense(
      IF_theta = IF_theta,
      states = states,
      target = target,
      influence_method = influence_method,
      beta_names = beta_names,
      joint_psd_tol = joint_psd_tol,
      return_if = return_if
    ))
  }
  cohort_build_H_blocked(
    IF_theta = IF_theta,
    states = states,
    target = target,
    influence_method = influence_method,
    H_block_size = H_block_size,
    H_block_cache = H_block_cache,
    H_cache_dir = H_cache_dir,
    beta_names = beta_names,
    joint_psd_tol = joint_psd_tol,
    return_if = return_if
  )
}

influence_internal_gmm_cohort <- function(Y, target, working, gamma_hat_list,
                                          family, ridge = 0,
                                          influence_method = c("profiled", "full"),
                                          H_compute = c("blocked", "dense"),
                                          H_block_size = 256L,
                                          H_block_cache = c("disk", "recompute", "memory"),
                                          H_cache_dir = tempdir(),
                                          rank_tol = 1e-10,
                                          bread_condition_warn = 1e12,
                                          joint_psd_tol = 1e-8,
                                          return_if = FALSE,
                                          poisson_clip = 35) {
  influence_method <- match.arg(influence_method)
  H_compute <- match.arg(H_compute)
  states <- vector("list", length(working))
  for (k in seq_along(working)) {
    wk <- working[[k]]
    gamma <- as.numeric(gamma_hat_list[[k]])
    names(gamma) <- wk$gamma_names
    ev <- cohort_gmm_eval(gamma, wk$D_k, target$mu_hat, family, poisson_clip)
    bread_diag <- cohort_assert_full_rank(
      ev$jacobian, paste0("GMM model ", k, " raw bread"),
      rank_tol, bread_condition_warn
    )
    bread_used <- cohort_regularize_negative_bread(ev$jacobian, ridge)
    C_full <- crossprod(
      wk$D_k,
      target$Z * as.numeric(target$mu_eta_hat %||%
        cohort_mu_eta(target$eta_hat, family, poisson_clip))
    ) / length(Y)
    profile <- profile_working_system(
      D_full_reg = bread_used,
      D_full_raw = ev$jacobian,
      q_alpha = wk$q_alpha_k,
      C_full = C_full,
      rank_tol = rank_tol,
      bread_condition_warn = bread_condition_warn,
      label = paste0("GMM model ", k, " bread")
    )
    beta_hat <- gamma[profile$beta_idx]
    names(beta_hat) <- wk$beta_names
    states[[k]] <- list(
      method = "GMM", wk = wk, gamma_hat = gamma, beta_hat = beta_hat,
      residual = ev$residual, bread_raw = ev$jacobian,
      bread_used = bread_used, bread_diagnostics = bread_diag,
      C_full = C_full, profile = profile,
      convergence = list(max_abs_moment = max(abs(ev$moment)))
    )
  }
  H <- cohort_build_H(
    IF_theta = target$IF_theta,
    states = states,
    target = target,
    influence_method = influence_method,
    H_compute = H_compute,
    H_block_size = H_block_size,
    H_block_cache = H_block_cache,
    H_cache_dir = H_cache_dir,
    beta_names = unlist(lapply(working, `[[`, "beta_names"), use.names = FALSE),
    joint_psd_tol = joint_psd_tol,
    return_if = return_if
  )
  c(H, list(IF_lambda = target$IF_lambda, IF_theta = target$IF_theta))
}

fit_internal_wkm_cohort <- function(Y, target, working, family,
                                    ridge = 0,
                                    influence_method = c("profiled", "full"),
                                    H_compute = c("blocked", "dense"),
                                    H_block_size = 256L,
                                    H_block_cache = c("disk", "recompute", "memory"),
                                    H_cache_dir = tempdir(),
                                    rank_tol = 1e-10,
                                    bread_condition_warn = 1e12,
                                    joint_psd_tol = 1e-8,
                                    return_if = FALSE,
                                    poisson_clip = 35) {
  influence_method <- match.arg(influence_method)
  H_compute <- match.arg(H_compute)
  states <- lapply(seq_along(working), function(k) {
    cohort_fit_working_wkm_k(
      Y = Y, wk = working[[k]], family = family,
      ridge = ridge, rank_tol = rank_tol,
      bread_condition_warn = bread_condition_warn,
      poisson_clip = poisson_clip, model_index = k
    )
  })
  beta_hat <- unlist(lapply(states, `[[`, "beta_hat"), use.names = TRUE)
  H <- cohort_build_H(
    IF_theta = target$IF_theta,
    states = states,
    target = target,
    influence_method = influence_method,
    H_compute = H_compute,
    H_block_size = H_block_size,
    H_block_cache = H_block_cache,
    H_cache_dir = H_cache_dir,
    beta_names = unlist(lapply(working, `[[`, "beta_names"), use.names = FALSE),
    joint_psd_tol = joint_psd_tol,
    return_if = return_if
  )
  c(H, list(
    beta_hat_WKM = beta_hat,
    gamma_hat_WKM = lapply(states, `[[`, "gamma_hat"),
    convergence = lapply(states, `[[`, "convergence"),
    IF_lambda = target$IF_lambda,
    IF_theta = target$IF_theta,
    states = states
  ))
}

cohort_reference_if_beta_columns <- function(state, local_beta_indices,
                                             target_internal,
                                             target_reference,
                                             influence_method = c("profiled", "full")) {
  influence_method <- match.arg(influence_method)
  wk <- state$wk
  local_beta_indices <- as.integer(local_beta_indices)
  if (length(local_beta_indices) == 0L ||
      any(local_beta_indices < 1L | local_beta_indices > wk$p_beta_k) ||
      anyDuplicated(local_beta_indices)) {
    stop("Invalid reference-path local beta-coordinate request for working model ",
         wk$model_index, ".")
  }
  m <- length(local_beta_indices)

  if (influence_method == "profiled") {
    E <- matrix(0, nrow = wk$p_beta_k, ncol = m)
    E[cbind(local_beta_indices, seq_len(m))] <- 1
    X <- solve(t(state$profile$S_beta), E)

    internal_effective <- target_internal$IF_lambda %*%
      (t(state$profile$C_prof) %*% X)

    reference_design <- wk$B_k %*% X
    if (wk$q_alpha_k > 0L) {
      reference_design <- reference_design -
        wk$A_k %*% (t(state$profile$T_beta_alpha) %*% X)
    }
    reference_effective <- sweep(
      reference_design, 1L, state$residual, FUN = "*"
    )
  } else {
    r <- ncol(wk$D_k)
    E <- matrix(0, nrow = r, ncol = m)
    full_idx <- wk$q_alpha_k + local_beta_indices
    E[cbind(full_idx, seq_len(m))] <- 1
    X <- solve(t(state$bread_used), E)
    internal_effective <- target_internal$IF_lambda %*%
      (t(state$C_full) %*% X)
    reference_effective <- sweep(
      wk$D_k %*% X, 1L, state$residual, FUN = "*"
    )
  }

  internal <- -internal_effective
  reference <- -reference_effective
  colnames(internal) <- wk$beta_names[local_beta_indices]
  colnames(reference) <- wk$beta_names[local_beta_indices]
  list(internal = internal, reference = reference)
}

cohort_make_reference_if_beta_block <- function(
    global_indices, states, beta_map,
    target_internal, target_reference,
    influence_method = c("profiled", "full")) {
  influence_method <- match.arg(influence_method)
  global_indices <- as.integer(global_indices)
  subset_map <- beta_map[
    match(global_indices, beta_map$global_index), , drop = FALSE
  ]
  if (anyNA(subset_map$model_index)) {
    stop("Invalid reference-path global beta block indices.")
  }

  n <- nrow(target_internal$IF_lambda)
  m_ref <- nrow(target_reference$Z)
  internal <- matrix(NA_real_, nrow = n, ncol = length(global_indices))
  reference <- matrix(NA_real_, nrow = m_ref, ncol = length(global_indices))
  for (k in unique(subset_map$model_index)) {
    positions <- which(subset_map$model_index == k)
    local <- subset_map$local_beta_index[positions]
    cols <- cohort_reference_if_beta_columns(
      state = states[[k]],
      local_beta_indices = local,
      target_internal = target_internal,
      target_reference = target_reference,
      influence_method = influence_method
    )
    internal[, positions] <- cols$internal
    reference[, positions] <- cols$reference
  }
  colnames(internal) <- subset_map$beta_name
  colnames(reference) <- subset_map$beta_name
  list(internal = internal, reference = reference)
}

cohort_build_H_reference_dense <- function(
    states, target_internal, target_reference,
    influence_method = c("profiled", "full"),
    beta_names = NULL,
    joint_psd_tol = 1e-8,
    return_if = FALSE) {
  influence_method <- match.arg(influence_method)
  beta_map <- cohort_beta_index_map(states)
  p <- nrow(beta_map)
  if (is.null(beta_names)) beta_names <- beta_map$beta_name
  IF <- cohort_make_reference_if_beta_block(
    seq_len(p), states, beta_map,
    target_internal, target_reference, influence_method
  )
  colnames(IF$internal) <- beta_names
  colnames(IF$reference) <- beta_names

  IF_theta <- as.matrix(target_internal$IF_theta)
  n <- nrow(IF_theta)
  m_ref <- nrow(IF$reference)
  theta_names <- colnames(IF_theta) %||% character(ncol(IF_theta))
  Htt <- cohort_cov(IF_theta)
  Htb <- cohort_crosscov(IF_theta, IF$internal)
  Hbb_internal <- cohort_cov(IF$internal)
  Hbb_reference <- (n / m_ref) * cohort_cov(IF$reference)
  Hbb <- cohort_symmetrize(Hbb_internal + Hbb_reference)

  dimnames(Htt) <- list(theta_names, theta_names)
  dimnames(Htb) <- list(theta_names, beta_names)
  dimnames(Hbb_internal) <- list(beta_names, beta_names)
  dimnames(Hbb_reference) <- list(beta_names, beta_names)
  dimnames(Hbb) <- list(beta_names, beta_names)
  joint <- cohort_make_joint_H(
    Htt, Htb, Hbb, theta_names, beta_names,
    joint_psd_tol = joint_psd_tol,
    label = "Two-population reference-path joint H covariance"
  )

  list(
    H_theta_theta = Htt,
    H_theta_beta = Htb,
    H_beta_beta = Hbb,
    H_beta_beta_internal = Hbb_internal,
    H_beta_beta_reference = Hbb_reference,
    H_joint = joint$H_joint,
    H_joint_psd_diagnostics = joint$diagnostics,
    beta_index_map = beta_map,
    beta_blocks = list(seq_len(p)),
    H_compute_diagnostics = list(
      method = "dense_reference",
      influence_method = influence_method,
      n_blocks = 1L,
      block_sizes = p,
      cache = "none",
      n_reference = m_ref,
      reference_scale = n / m_ref,
      preserves_off_diagonal_blocks = TRUE
    ),
    IF_beta_internal = if (return_if) IF$internal else NULL,
    IF_beta_reference = if (return_if) IF$reference else NULL
  )
}

cohort_build_H_reference_blocked <- function(
    states, target_internal, target_reference,
    influence_method = c("profiled", "full"),
    H_block_size = 256L,
    H_block_cache = c("disk", "recompute", "memory"),
    H_cache_dir = tempdir(),
    beta_names = NULL,
    joint_psd_tol = 1e-8,
    return_if = FALSE) {
  influence_method <- match.arg(influence_method)
  H_block_cache <- match.arg(H_block_cache)
  beta_map <- cohort_beta_index_map(states)
  p <- nrow(beta_map)
  if (is.null(beta_names)) beta_names <- beta_map$beta_name
  beta_blocks <- cohort_make_beta_blocks(p, H_block_size)
  B <- length(beta_blocks)

  IF_theta <- as.matrix(target_internal$IF_theta)
  n <- nrow(IF_theta)
  m_ref <- nrow(target_reference$Z)
  d <- ncol(IF_theta)
  theta_names <- colnames(IF_theta) %||% character(d)
  IF_theta_c <- cohort_center_rows(IF_theta)
  reference_scale <- n / m_ref

  Htt <- if (d) {
    cohort_symmetrize(crossprod(IF_theta_c) / n)
  } else {
    matrix(0, 0, 0)
  }
  Htb <- matrix(0, nrow = d, ncol = p)
  Hbb_internal <- matrix(0, nrow = p, ncol = p)
  Hbb_reference <- matrix(0, nrow = p, ncol = p)
  IF_internal_return <- if (return_if) {
    matrix(NA_real_, nrow = n, ncol = p)
  } else NULL
  IF_reference_return <- if (return_if) {
    matrix(NA_real_, nrow = m_ref, ncol = p)
  } else NULL
  if (return_if) {
    colnames(IF_internal_return) <- beta_names
    colnames(IF_reference_return) <- beta_names
  }

  cache_memory <- vector("list", if (H_block_cache == "memory") B else 0L)
  cache_paths <- character(if (H_block_cache == "disk") B else 0L)
  run_dir <- NULL
  if (H_block_cache == "disk") {
    if (!dir.exists(H_cache_dir)) {
      dir.create(H_cache_dir, recursive = TRUE, showWarnings = FALSE)
    }
    run_dir <- tempfile("zim_cohort_reference_H_blocks_", tmpdir = H_cache_dir)
    if (!dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)) {
      stop("Could not create reference-path H block cache directory: ", run_dir)
    }
    on.exit(unlink(run_dir, recursive = TRUE, force = TRUE), add = TRUE)
  }

  make_uncentered <- function(b) {
    cohort_make_reference_if_beta_block(
      beta_blocks[[b]], states, beta_map,
      target_internal, target_reference, influence_method
    )
  }
  make_centered <- function(b) {
    block <- make_uncentered(b)
    list(
      internal = cohort_center_rows(block$internal),
      reference = cohort_center_rows(block$reference)
    )
  }

  for (b in seq_len(B)) {
    idx_b <- beta_blocks[[b]]
    IF_b <- make_uncentered(b)
    IF_b_c <- list(
      internal = cohort_center_rows(IF_b$internal),
      reference = cohort_center_rows(IF_b$reference)
    )
    if (d) {
      Htb[, idx_b] <- crossprod(IF_theta_c, IF_b_c$internal) / n
    }
    Hbb_internal[idx_b, idx_b] <- cohort_symmetrize(
      crossprod(IF_b_c$internal) / n
    )
    Hbb_reference[idx_b, idx_b] <- reference_scale * cohort_symmetrize(
      crossprod(IF_b_c$reference) / m_ref
    )
    if (return_if) {
      IF_internal_return[, idx_b] <- IF_b$internal
      IF_reference_return[, idx_b] <- IF_b$reference
    }

    if (H_block_cache == "memory") {
      cache_memory[[b]] <- IF_b_c
    } else if (H_block_cache == "disk") {
      path <- file.path(run_dir, sprintf("if_beta_reference_%05d.rds", b))
      saveRDS(IF_b_c, path, compress = FALSE)
      cache_paths[b] <- path
    }
    rm(IF_b, IF_b_c)
  }

  get_centered <- function(b) {
    if (H_block_cache == "memory") return(cache_memory[[b]])
    if (H_block_cache == "disk") return(readRDS(cache_paths[b]))
    make_centered(b)
  }

  if (B > 1L) {
    for (b in 2:B) {
      idx_b <- beta_blocks[[b]]
      IF_b_c <- get_centered(b)
      for (c_block in seq_len(b - 1L)) {
        idx_c <- beta_blocks[[c_block]]
        IF_c_c <- get_centered(c_block)
        H_i_cb <- crossprod(IF_c_c$internal, IF_b_c$internal) / n
        H_r_cb <- reference_scale *
          crossprod(IF_c_c$reference, IF_b_c$reference) / m_ref
        Hbb_internal[idx_c, idx_b] <- H_i_cb
        Hbb_internal[idx_b, idx_c] <- t(H_i_cb)
        Hbb_reference[idx_c, idx_b] <- H_r_cb
        Hbb_reference[idx_b, idx_c] <- t(H_r_cb)
        rm(IF_c_c)
      }
      rm(IF_b_c)
    }
  }

  Hbb_internal <- cohort_symmetrize(Hbb_internal)
  Hbb_reference <- cohort_symmetrize(Hbb_reference)
  Hbb <- cohort_symmetrize(Hbb_internal + Hbb_reference)
  dimnames(Htt) <- list(theta_names, theta_names)
  dimnames(Htb) <- list(theta_names, beta_names)
  dimnames(Hbb_internal) <- list(beta_names, beta_names)
  dimnames(Hbb_reference) <- list(beta_names, beta_names)
  dimnames(Hbb) <- list(beta_names, beta_names)
  joint <- cohort_make_joint_H(
    Htt, Htb, Hbb, theta_names, beta_names,
    joint_psd_tol = joint_psd_tol,
    label = "Two-population blockwise reference-path joint H covariance"
  )

  list(
    H_theta_theta = Htt,
    H_theta_beta = Htb,
    H_beta_beta = Hbb,
    H_beta_beta_internal = Hbb_internal,
    H_beta_beta_reference = Hbb_reference,
    H_joint = joint$H_joint,
    H_joint_psd_diagnostics = joint$diagnostics,
    beta_index_map = beta_map,
    beta_blocks = beta_blocks,
    H_compute_diagnostics = list(
      method = "blocked_reference",
      influence_method = influence_method,
      n_blocks = B,
      block_sizes = vapply(beta_blocks, length, integer(1)),
      cache = H_block_cache,
      peak_beta_if_columns_per_source = max(
        vapply(beta_blocks, length, integer(1))
      ),
      n_reference = m_ref,
      reference_scale = reference_scale,
      preserves_off_diagonal_blocks = TRUE,
      covariance_sources = c("internal", "reference")
    ),
    IF_beta_internal = IF_internal_return,
    IF_beta_reference = IF_reference_return
  )
}

cohort_build_H_reference <- function(
    states, target_internal, target_reference,
    influence_method = c("profiled", "full"),
    H_compute = c("blocked", "dense"),
    H_block_size = 256L,
    H_block_cache = c("disk", "recompute", "memory"),
    H_cache_dir = tempdir(),
    beta_names = NULL,
    joint_psd_tol = 1e-8,
    return_if = FALSE) {
  influence_method <- match.arg(influence_method)
  H_compute <- match.arg(H_compute)
  H_block_cache <- match.arg(H_block_cache)
  if (H_compute == "dense") {
    return(cohort_build_H_reference_dense(
      states = states,
      target_internal = target_internal,
      target_reference = target_reference,
      influence_method = influence_method,
      beta_names = beta_names,
      joint_psd_tol = joint_psd_tol,
      return_if = return_if
    ))
  }
  cohort_build_H_reference_blocked(
    states = states,
    target_internal = target_internal,
    target_reference = target_reference,
    influence_method = influence_method,
    H_block_size = H_block_size,
    H_block_cache = H_block_cache,
    H_cache_dir = H_cache_dir,
    beta_names = beta_names,
    joint_psd_tol = joint_psd_tol,
    return_if = return_if
  )
}

cohort_fit_reference_gmm <- function(Y, target_internal, target_reference,
                                     working_internal, working_reference,
                                     beta_tilde, family,
                                     influence_method = c("profiled", "full"),
                                     H_compute = c("blocked", "dense"),
                                     H_block_size = 256L,
                                     H_block_cache = c("disk", "recompute", "memory"),
                                     H_cache_dir = tempdir(),
                                     ridge = 0,
                                     rank_tol = 1e-10,
                                     bread_condition_warn = 1e12,
                                     joint_psd_tol = 1e-8,
                                     gmm_solver = c("newton", "optim"),
                                     gmm_score_tol = 1e-8,
                                     gmm_step_tol = 1e-10,
                                     gmm_maxit = 100L,
                                     gmm_max_halving = 30L,
                                     gmm_newton_fallback_to_optim = TRUE,
                                     optim_control = list(maxit = 1000),
                                     poisson_clip = 35,
                                     return_if = FALSE) {
  influence_method <- match.arg(influence_method)
  H_compute <- match.arg(H_compute)
  H_block_cache <- match.arg(H_block_cache)
  gmm_solver <- match.arg(gmm_solver)
  states <- vector("list", length(working_reference))
  beta_parts <- vector("list", length(working_reference))
  offset <- 0L

  for (k in seq_along(working_reference)) {
    wk_ref <- working_reference[[k]]
    wk_int <- working_internal[[k]]
    idx <- offset + seq_len(wk_ref$p_beta_k)
    if (family$family == "gaussian") {
      gamma_start <- as.numeric(stats::coef(
        stats::lm.fit(x = wk_ref$D_k, y = target_reference$mu_hat)
      ))
      if (length(gamma_start) != ncol(wk_ref$D_k) || any(!is.finite(gamma_start))) {
        stop("Reference GMM model ", k, " has an invalid Gaussian starting value.")
      }
      names(gamma_start) <- wk_ref$gamma_names
    } else {
      gamma_start <- fit_glm_start_cohort(
        Y = Y,
        wk = wk_int,
        family = family,
        beta_tilde_k = beta_tilde[idx],
        poisson_clip = poisson_clip
      )
      names(gamma_start) <- wk_ref$gamma_names
    }

    sol <- cohort_solve_gmm(
      D = wk_ref$D_k,
      target_mu = target_reference$mu_hat,
      family = family,
      gamma_start = gamma_start,
      solver = gmm_solver,
      ridge = ridge,
      score_tol = gmm_score_tol,
      step_tol = gmm_step_tol,
      maxit = gmm_maxit,
      max_halving = gmm_max_halving,
      fallback_to_optim = gmm_newton_fallback_to_optim,
      optim_control = optim_control,
      rank_tol = rank_tol,
      bread_condition_warn = bread_condition_warn,
      poisson_clip = poisson_clip,
      label = paste0("Reference GMM model ", k)
    )
    C_full <- crossprod(
      wk_ref$D_k,
      target_reference$Z * as.numeric(target_reference$mu_eta_hat)
    ) / nrow(wk_ref$D_k)
    profile <- profile_working_system(
      D_full_reg = sol$bread_used,
      D_full_raw = sol$bread_raw,
      q_alpha = wk_ref$q_alpha_k,
      C_full = C_full,
      rank_tol = rank_tol,
      bread_condition_warn = bread_condition_warn,
      label = paste0("Reference GMM model ", k, " bread")
    )
    beta <- sol$gamma_hat[profile$beta_idx]
    names(beta) <- wk_ref$beta_names
    beta_parts[[k]] <- beta
    states[[k]] <- list(
      method = "GMM",
      wk = wk_ref,
      gamma_hat = sol$gamma_hat,
      beta_hat = beta,
      residual = sol$residual,
      bread_raw = sol$bread_raw,
      bread_used = sol$bread_used,
      bread_diagnostics = sol$jacobian_diagnostics,
      C_full = C_full,
      profile = profile,
      convergence = list(
        solution_validated = TRUE,
        solver = sol$solver,
        max_abs_moment = sol$max_abs_moment,
        jacobian_condition_number = sol$jacobian_diagnostics$condition_number
      )
    )
    offset <- offset + wk_ref$p_beta_k
  }

  beta_hat <- unlist(beta_parts, use.names = TRUE)
  beta_names <- unlist(lapply(states, function(s) s$wk$beta_names), use.names = FALSE)
  H <- cohort_build_H_reference(
    states = states,
    target_internal = target_internal,
    target_reference = target_reference,
    influence_method = influence_method,
    H_compute = H_compute,
    H_block_size = H_block_size,
    H_block_cache = H_block_cache,
    H_cache_dir = H_cache_dir,
    beta_names = beta_names,
    joint_psd_tol = joint_psd_tol,
    return_if = return_if
  )

  c(H, list(
    beta_hat = beta_hat,
    gamma_hat = lapply(states, `[[`, "gamma_hat"),
    states = states,
    convergence = lapply(states, `[[`, "convergence")
  ))
}

influence_reference_gmm_cohort <- function(Y, target_internal, target_reference,
                                           working_reference, gamma_hat_list,
                                           family, ridge = 0,
                                           influence_method = c("profiled", "full"),
                                           H_compute = c("blocked", "dense"),
                                           H_block_size = 256L,
                                           H_block_cache = c("disk", "recompute", "memory"),
                                           H_cache_dir = tempdir(),
                                           rank_tol = 1e-10,
                                           bread_condition_warn = 1e12,
                                           joint_psd_tol = 1e-8,
                                           return_if = FALSE,
                                           poisson_clip = 35) {
  influence_method <- match.arg(influence_method)
  H_compute <- match.arg(H_compute)
  H_block_cache <- match.arg(H_block_cache)
  states <- vector("list", length(working_reference))
  for (k in seq_along(working_reference)) {
    wk <- working_reference[[k]]
    gamma <- as.numeric(gamma_hat_list[[k]])
    names(gamma) <- wk$gamma_names
    ev <- cohort_gmm_eval(gamma, wk$D_k, target_reference$mu_hat, family, poisson_clip)
    bread_diag <- cohort_assert_full_rank(
      ev$jacobian, paste0("Reference GMM model ", k, " raw bread"),
      rank_tol, bread_condition_warn
    )
    bread_used <- cohort_regularize_negative_bread(ev$jacobian, ridge)
    C_full <- crossprod(
      wk$D_k,
      target_reference$Z * as.numeric(target_reference$mu_eta_hat)
    ) / nrow(wk$D_k)
    profile <- profile_working_system(
      D_full_reg = bread_used, D_full_raw = ev$jacobian,
      q_alpha = wk$q_alpha_k, C_full = C_full,
      rank_tol = rank_tol, bread_condition_warn = bread_condition_warn,
      label = paste0("Reference GMM model ", k, " bread")
    )
    beta <- gamma[profile$beta_idx]
    names(beta) <- wk$beta_names
    states[[k]] <- list(
      method = "GMM", wk = wk, gamma_hat = gamma, beta_hat = beta,
      residual = ev$residual, bread_raw = ev$jacobian,
      bread_used = bread_used, bread_diagnostics = bread_diag,
      C_full = C_full, profile = profile,
      convergence = list(max_abs_moment = max(abs(ev$moment)))
    )
  }
  beta_names <- unlist(lapply(states, function(s) s$wk$beta_names), use.names = FALSE)
  H <- cohort_build_H_reference(
    states = states,
    target_internal = target_internal,
    target_reference = target_reference,
    influence_method = influence_method,
    H_compute = H_compute,
    H_block_size = H_block_size,
    H_block_cache = H_block_cache,
    H_cache_dir = H_cache_dir,
    beta_names = beta_names,
    joint_psd_tol = joint_psd_tol,
    return_if = return_if
  )
  c(H, list(
    IF_lambda = target_internal$IF_lambda,
    IF_theta = target_internal$IF_theta,
    n_reference = nrow(target_reference$Z)
  ))
}

# fit_external_summary_cohort(): fit working models to a cohort dataset.
# external_data: data frame; Y: outcome column name.
# working_models: list of formulas or lists with alpha and beta formulas.
# outcome_type: "binary", "continuous", or "poisson"; family overrides it.
# n_internal: optional internal sample size for the returned Sigma_beta_beta.
# return_if: whether to include IF_beta_external.
# Returns beta_tilde, gamma_tilde, covariance matrices, sample size, and diagnostics.

fit_external_summary_cohort <- function(external_data, Y, working_models,
                                        outcome_type = c("binary", "continuous", "poisson"),
                                        family = NULL,
                                        n_internal = NULL,
                                        ridge = 0,
                                        influence_method = c("profiled", "full"),
                                        H_compute = c("blocked", "dense"),
                                        H_block_size = 256L,
                                        H_block_cache = c("disk", "recompute", "memory"),
                                        H_cache_dir = tempdir(),
                                        rank_tol = 1e-10,
                                        bread_condition_warn = 1e12,
                                        joint_psd_tol = 1e-8,
                                        return_if = FALSE,
                                        poisson_clip = 35) {
  fam <- resolve_cohort_family(family = family, outcome_type = outcome_type)
  family <- fam$family
  outcome_type <- fam$outcome_type
  influence_method <- match.arg(influence_method)
  H_compute <- match.arg(H_compute)
  H_block_cache <- match.arg(H_block_cache)

  prep <- cohort_prepare_analysis_data(
    data = external_data,
    Y = Y,
    working_models = working_models,
    data_label = "external_data"
  )
  external_data <- prep$data
  Y_ext <- check_outcome(external_data[[Y]], family)
  N <- length(Y_ext)
  working_ext <- cohort_pack_working_models(working_models, external_data)
  beta_names <- attr(working_ext, "beta_names")

  states <- lapply(seq_along(working_ext), function(k) {
    cohort_fit_working_wkm_k(
      Y = Y_ext,
      wk = working_ext[[k]],
      family = family,
      ridge = ridge,
      rank_tol = rank_tol,
      bread_condition_warn = bread_condition_warn,
      poisson_clip = poisson_clip,
      model_index = k
    )
  })
  beta_tilde <- unlist(lapply(states, `[[`, "beta_hat"), use.names = TRUE)
  beta_tilde <- beta_tilde[beta_names]

  empty_target <- list(IF_lambda = matrix(0, nrow = N, ncol = 0L))
  H <- cohort_build_H(
    IF_theta = matrix(0, nrow = N, ncol = 0L),
    states = states,
    target = empty_target,
    influence_method = influence_method,
    H_compute = H_compute,
    H_block_size = H_block_size,
    H_block_cache = H_block_cache,
    H_cache_dir = H_cache_dir,
    beta_names = beta_names,
    joint_psd_tol = joint_psd_tol,
    return_if = return_if
  )
  H_ext_beta_beta <- H$H_beta_beta
  V_beta_tilde <- H_ext_beta_beta / N
  Sigma_beta_beta <- NULL
  if (!is.null(n_internal)) {
    if (!is.numeric(n_internal) || length(n_internal) != 1L ||
        !is.finite(n_internal) || n_internal <= 0) {
      stop("n_internal must be a positive finite scalar.")
    }
    Sigma_beta_beta <- (n_internal / N) * H_ext_beta_beta
    dimnames(Sigma_beta_beta) <- list(beta_names, beta_names)
  }

  list(
    beta_tilde = beta_tilde,
    gamma_tilde = lapply(states, `[[`, "gamma_hat"),
    H_ext_beta_beta = H_ext_beta_beta,
    V_beta_tilde = V_beta_tilde,
    Sigma_beta_beta = Sigma_beta_beta,
    IF_beta_external = if (return_if) H$IF_beta else NULL,
    family = family,
    outcome_type = outcome_type,
    n_external = N,
    data_diagnostics = prep,
    convergence = lapply(states, `[[`, "convergence"),
    bread_diagnostics = lapply(states, `[[`, "bread_diagnostics"),
    profile_diagnostics = lapply(states, function(s) {
      list(
        raw = s$profile$raw_profile_diagnostics,
        used = s$profile$used_profile_diagnostics
      )
    }),
    H_joint_psd_diagnostics = H$H_joint_psd_diagnostics,
    H_compute_diagnostics = H$H_compute_diagnostics,
    beta_index_map = H$beta_index_map,
    beta_block_map = H$beta_blocks
  )
}

cohort_validate_external_sample_size <- function(n_external) {
  if (!is.numeric(n_external) || length(n_external) != 1L ||
      !is.finite(n_external) || n_external <= 0) {
    stop("n_external must be a positive finite scalar when estimating Sigma_beta_beta.")
  }
  as.numeric(n_external)
}

estimate_sigma_betabeta_internal_cohort <- function(
    internal_data, Y, working_models, n_external,
    outcome_type = c("binary", "continuous", "poisson"), family = NULL,
    influence_method = c("profiled", "full"),
    H_compute = c("blocked", "dense"), H_block_size = 256L,
    H_block_cache = c("disk", "recompute", "memory"),
    H_cache_dir = tempdir(), rank_tol = 1e-10, bread_condition_warn = 1e12,
    joint_psd_tol = 1e-8, return_if = FALSE, poisson_clip = 35) {
  n_external <- cohort_validate_external_sample_size(n_external)
  prep <- cohort_prepare_analysis_data(
    internal_data, Y = Y, working_models = working_models,
    data_label = "internal_data"
  )

  marginal <- fit_external_summary_cohort(
    external_data = prep$data, Y = Y, working_models = working_models,
    outcome_type = outcome_type, family = family, n_internal = NULL,
    ridge = 0, influence_method = match.arg(influence_method),
    H_compute = match.arg(H_compute), H_block_size = H_block_size,
    H_block_cache = match.arg(H_block_cache), H_cache_dir = H_cache_dir,
    rank_tol = rank_tol, bread_condition_warn = bread_condition_warn,
    joint_psd_tol = joint_psd_tol, return_if = return_if,
    poisson_clip = poisson_clip
  )
  n <- marginal$n_external
  H <- marginal$H_ext_beta_beta
  list(
    Sigma_beta_beta = (n / n_external) * H,
    V_beta_tilde = H / n_external,
    H_internal_beta_beta = H,
    beta_hat_internal_marginal = marginal$beta_tilde,
    gamma_hat_internal_marginal = marginal$gamma_tilde,
    n_internal = n, n_external = n_external,
    method = "internal_marginal_sandwich", source = "internal_data",
    scaling_factor = n / n_external,
    data_diagnostics = prep,
    convergence = marginal$convergence,
    bread_diagnostics = marginal$bread_diagnostics,
    profile_diagnostics = marginal$profile_diagnostics,
    H_compute_diagnostics = marginal$H_compute_diagnostics,
    IF_beta_internal_marginal = marginal$IF_beta_external
  )
}

estimate_sigma_betabeta_reference_cohort <- function(
    states, target_internal, target_reference, n_external,
    H_compute = c("blocked", "dense"), H_block_size = 256L,
    H_block_cache = c("disk", "recompute", "memory"),
    H_cache_dir = tempdir(), joint_psd_tol = 1e-8, return_if = FALSE) {
  n_external <- cohort_validate_external_sample_size(n_external)
  H_compute <- match.arg(H_compute)
  H_block_cache <- match.arg(H_block_cache)
  if (!is.logical(return_if) || length(return_if) != 1L || is.na(return_if)) {
    stop("return_if must be TRUE or FALSE.")
  }
  n <- nrow(target_internal$IF_lambda)
  m_ref <- nrow(target_reference$Z)
  mu <- as.numeric(target_reference$mu_hat)
  family_name <- target_internal$family$family
  conditional_variance <- switch(
    family_name,
    binomial = mu * (1 - mu),
    gaussian = rep(target_internal$diagnostics$residual_scale^2, m_ref),
    poisson = mu,
    stop("Unsupported family for reference covariance estimation.")
  )
  if (length(conditional_variance) != m_ref ||
      any(!is.finite(conditional_variance)) || any(conditional_variance < 0)) {
    stop("Invalid target conditional variance for reference covariance estimation.")
  }
  variance_assumption <- switch(
    family_name,
    binomial = "Transported Bernoulli conditional outcome law",
    gaussian = paste("Transported homoskedastic Gaussian conditional outcome law;",
                     "variance estimated from internal target RSS / (n - q)"),
    poisson = "Transported equidispersed Poisson conditional outcome law"
  )
  beta_map <- cohort_beta_index_map(states)
  p <- nrow(beta_map)
  beta_names <- beta_map$beta_name
  beta_blocks <- if (H_compute == "dense") list(seq_len(p)) else {
    cohort_make_beta_blocks(p, H_block_size)
  }
  B <- length(beta_blocks)

  make_block <- function(b) {
    idx <- beta_blocks[[b]]
    subset_map <- beta_map[match(idx, beta_map$global_index), , drop = FALSE]
    noise <- mean_if <- matrix(NA_real_, nrow = m_ref, ncol = length(idx))
    for (k in unique(subset_map$model_index)) {
      state <- states[[k]]
      wk <- state$wk
      positions <- which(subset_map$model_index == k)
      local <- subset_map$local_beta_index[positions]
      E <- matrix(0, nrow = ncol(wk$D_k), ncol = length(local))
      E[cbind(wk$q_alpha_k + local, seq_along(local))] <- 1

      loadings <- -wk$D_k %*% solve(t(state$bread_raw), E)
      noise[, positions] <- loadings * sqrt(conditional_variance)
      mean_if[, positions] <- loadings * as.numeric(state$residual)
    }
    colnames(noise) <- colnames(mean_if) <- subset_map$beta_name

    rbind(noise, cohort_center_rows(mean_if))
  }

  cache <- if (H_compute == "dense") "memory" else H_block_cache
  cache_memory <- vector("list", if (cache == "memory") B else 0L)
  cache_paths <- character(if (cache == "disk") B else 0L)
  if (cache == "disk") {
    if (!dir.exists(H_cache_dir)) {
      dir.create(H_cache_dir, recursive = TRUE, showWarnings = FALSE)
    }
    run_dir <- tempfile("zim_cohort_reference_sigma_blocks_", tmpdir = H_cache_dir)
    if (!dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)) {
      stop("Could not create reference covariance block cache directory: ", run_dir)
    }
    on.exit(unlink(run_dir, recursive = TRUE, force = TRUE), add = TRUE)
  }
  H <- matrix(0, nrow = p, ncol = p, dimnames = list(beta_names, beta_names))
  noise_return <- mean_return <- if (return_if) {
    matrix(NA_real_, nrow = m_ref, ncol = p, dimnames = list(NULL, beta_names))
  } else NULL
  for (b in seq_len(B)) {
    idx <- beta_blocks[[b]]
    block <- make_block(b)
    H[idx, idx] <- crossprod(block) / m_ref
    if (return_if) {
      noise_return[, idx] <- block[seq_len(m_ref), , drop = FALSE]
      mean_return[, idx] <- block[m_ref + seq_len(m_ref), , drop = FALSE]
    }
    if (cache == "memory") {
      cache_memory[[b]] <- block
    } else if (cache == "disk") {
      cache_paths[b] <- file.path(run_dir, sprintf("sigma_beta_%05d.rds", b))
      saveRDS(block, cache_paths[b], compress = FALSE)
    }
    rm(block)
  }
  get_block <- function(b) {
    if (cache == "memory") return(cache_memory[[b]])
    if (cache == "disk") return(readRDS(cache_paths[b]))
    make_block(b)
  }
  if (B > 1L) {
    for (b in 2:B) {
      idx_b <- beta_blocks[[b]]
      block_b <- get_block(b)
      for (c_block in seq_len(b - 1L)) {
        idx_c <- beta_blocks[[c_block]]
        block_c <- get_block(c_block)
        H_cb <- crossprod(block_c, block_b) / m_ref
        H[idx_c, idx_b] <- H_cb
        H[idx_b, idx_c] <- t(H_cb)
        rm(block_c)
      }
      rm(block_b)
    }
  }
  H <- cohort_symmetrize(H)
  diagnostics <- cohort_psd_diagnostics(
    H, tol = joint_psd_tol,
    label = "Reference-model external marginal covariance", action = "stop"
  )
  list(
    Sigma_beta_beta = (n / n_external) * H,
    V_beta_tilde = H / n_external,
    H_external_beta_beta = H,
    n_internal = n, n_external = n_external, n_reference = m_ref,
    method = "reference_model_based_sandwich",
    source = "internal_and_reference_data", scaling_factor = n / n_external,
    conditional_variance_assumption = variance_assumption,
    conditional_variance_range = range(conditional_variance),
    convergence = lapply(states, `[[`, "convergence"),
    bread_diagnostics = lapply(states, `[[`, "bread_diagnostics"),
    H_compute_diagnostics = list(
      method = paste0(H_compute, "_reference_external_covariance"),
      n_blocks = B, block_sizes = vapply(beta_blocks, length, integer(1)),
      cache = if (H_compute == "dense") "none" else cache,
      preserves_off_diagonal_blocks = TRUE,
      noise_loadings_centered = FALSE, mean_influence_centered = TRUE,
      covariance_denominator = m_ref, raw_bread_used = TRUE
    ),
    covariance_diagnostics = diagnostics,
    noise_loadings = noise_return,
    mean_influence = mean_return
  )
}

cohort_select_pc_components <- function(values,
                                        pca_rule = c(
                                          "variance", "condition",
                                          "condition_or_variance", "fixed"
                                        ),
                                        pc_var = 0.99,
                                        n_pc = NULL,
                                        eig_tol = 1e-10,
                                        pc_max_condition = 1e4,
                                        label = "PCA basis") {
  pca_rule <- match.arg(pca_rule)
  values <- as.numeric(values)
  if (!length(values) || any(!is.finite(values))) {
    stop(label, " has no finite eigenvalues.")
  }
  lambda <- pmax(values, 0)
  lambda1 <- max(lambda)
  positive_cutoff <- eig_tol * max(lambda1, .Machine$double.eps)
  q <- sum(lambda > positive_cutoff)
  if (q < 1L) stop(label, " has no numerically positive eigenvalues.")
  lambda_pos <- lambda[seq_len(q)]

  if (!is.numeric(pc_var) || length(pc_var) != 1L || !is.finite(pc_var) ||
      pc_var <= 0 || pc_var > 1) {
    stop("pc_var must lie in (0, 1].")
  }
  if (!is.numeric(pc_max_condition) || length(pc_max_condition) != 1L ||
      !is.finite(pc_max_condition) || pc_max_condition < 1) {
    stop("pc_max_condition must be a finite number at least 1.")
  }

  cumulative <- cumsum(lambda_pos) / sum(lambda_pos)
  r_var <- which(cumulative >= pc_var)[1]
  admissible_condition <- which(lambda_pos[1] / lambda_pos <= pc_max_condition)
  r_cond <- if (length(admissible_condition)) max(admissible_condition) else 1L

  if (pca_rule == "fixed") {
    if (is.null(n_pc) || length(n_pc) != 1L || !is.finite(n_pc)) {
      stop("n_pc must be supplied for pca_rule = 'fixed'.")
    }
    r <- as.integer(n_pc)
    if (r < 1L || r > q) {
      stop("n_pc must be between 1 and the numerical positive rank (", q, ").")
    }
  } else if (pca_rule == "variance") {
    r <- r_var
  } else if (pca_rule == "condition") {
    r <- r_cond
  } else {
    r <- min(r_var, r_cond)
  }

  list(
    rank = as.integer(r),
    positive_rank = as.integer(q),
    positive_cutoff = positive_cutoff,
    pc_var_used = cumulative[r],
    retained_condition_number = lambda_pos[1] / lambda_pos[r],
    variance_rank = as.integer(r_var),
    condition_rank = as.integer(r_cond),
    working_eigenvalues = lambda,
    selected_values = lambda_pos[seq_len(r)]
  )
}

cohort_build_pc_basis <- function(A_raw,
                                  pca_scale = c("covariance", "correlation"),
                                  regularize_spectral = FALSE,
                                  ridge = 0,
                                  pca_rule = c(
                                    "variance", "condition",
                                    "condition_or_variance", "fixed"
                                  ),
                                  pc_var = 0.99,
                                  n_pc = NULL,
                                  eig_tol = 1e-10,
                                  pc_max_condition = 1e4,
                                  label = "PCA basis") {
  pca_scale <- match.arg(pca_scale)
  pca_rule <- match.arg(pca_rule)
  A_raw <- cohort_symmetrize(A_raw)
  p <- nrow(A_raw)
  if (p < 1L) stop(label, " is empty.")

  if (pca_scale == "correlation") {
    diagonal <- diag(A_raw)
    if (any(!is.finite(diagonal)) || any(diagonal <= 0)) {
      stop(label, " requires strictly positive finite diagonal entries on correlation scale.")
    }
    scale_sd <- sqrt(diagonal)
    A_spec_raw <- sweep(sweep(A_raw, 1L, scale_sd, FUN = "/"),
                        2L, scale_sd, FUN = "/")
  } else {
    scale_sd <- rep(1, p)
    A_spec_raw <- A_raw
  }
  A_spec_raw <- cohort_symmetrize(A_spec_raw)
  A_spec <- if (regularize_spectral) {
    cohort_regularize_covariance(A_spec_raw, ridge)
  } else {
    A_spec_raw
  }

  eigen_raw <- eigen(A_spec_raw, symmetric = TRUE, only.values = TRUE)$values
  eig <- eigen(A_spec, symmetric = TRUE)
  selection <- cohort_select_pc_components(
    eig$values,
    pca_rule = pca_rule,
    pc_var = pc_var,
    n_pc = n_pc,
    eig_tol = eig_tol,
    pc_max_condition = pc_max_condition,
    label = label
  )
  r <- selection$rank
  U <- eig$vectors[, seq_len(r), drop = FALSE]
  B <- if (pca_scale == "correlation") {
    sweep(U, 1L, scale_sd, FUN = "/")
  } else {
    U
  }
  material_cutoff <- eig_tol * max(1, max(abs(eigen_raw)))

  list(
    B = B,
    U = U,
    selected_values = selection$selected_values,
    scale_sd = scale_sd,
    pca_scale = pca_scale,
    A_raw = A_raw,
    A_spectral_raw = A_spec_raw,
    A_spectral = A_spec,
    eigenvalues = eig$values,
    selection = c(selection, list(
      raw_min_eigenvalue = min(eigen_raw),
      raw_max_eigenvalue = max(eigen_raw),
      material_negative_eigenvalues = sum(eigen_raw < -material_cutoff)
    ))
  )
}

cohort_build_factor_diag_covariance <- function(A_raw, basis,
                                                tolerance = 1e-10) {
  A_raw <- cohort_symmetrize(A_raw)
  U <- basis$U
  lambda <- basis$selected_values
  low_spectral <- U %*% diag(lambda, nrow = length(lambda)) %*% t(U)
  residual_spectral <- diag(basis$A_spectral_raw) - diag(low_spectral)
  scale <- max(1, max(abs(diag(basis$A_spectral_raw))))
  if (min(residual_spectral) < -tolerance * scale) {
    warning(
      "Factor-plus-diagonal residual variance was slightly negative; truncating at zero.",
      call. = FALSE
    )
  }
  residual_spectral <- pmax(residual_spectral, 0)
  model_spectral <- cohort_symmetrize(
    low_spectral + diag(residual_spectral, length(residual_spectral))
  )
  if (basis$pca_scale == "correlation") {
    D <- diag(basis$scale_sd, length(basis$scale_sd))
    model_raw <- D %*% model_spectral %*% D
  } else {
    model_raw <- model_spectral
  }
  model_raw <- cohort_symmetrize(model_raw)
  diagonal_error <- max(abs(diag(model_raw) - diag(A_raw)))
  offdiag <- function(A) A - diag(diag(A), nrow(A))

  list(
    covariance = model_raw,
    covariance_spectral = model_spectral,
    low_rank_spectral = low_spectral,
    residual_diag_spectral = residual_spectral,
    diagonal_error = diagonal_error,
    diagonal_preserved = diagonal_error <= tolerance * max(1, max(abs(diag(A_raw)))),
    offdiag_removed_frobenius = sqrt(sum((offdiag(A_raw) - offdiag(model_raw))^2))
  )
}

cohort_project_crosscov_to_factor_space <- function(C, basis) {
  C <- as.matrix(C)
  U <- basis$U
  if (basis$pca_scale == "covariance") {
    return(C %*% U %*% t(U))
  }
  sd <- basis$scale_sd
  C_std <- sweep(C, 2L, sd, FUN = "/")
  C_proj_std <- C_std %*% U %*% t(U)
  sweep(C_proj_std, 2L, sd, FUN = "*")
}

cohort_compute_zim_weight <- function(H_theta_beta,
                                      H_beta_beta,
                                      Sigma_beta_beta,
                                      ridge = 0,
                                      inverse_method = c(
                                        "solve", "pca", "hbb_pca", "sigma_pca",
                                        "m_factor_diag", "hbb_factor_diag", "shrink"
                                      ),
                                      pc_var = 0.99,
                                      n_pc = NULL,
                                      eig_tol = 1e-10,
                                      pca_rule = c(
                                        "variance", "condition",
                                        "condition_or_variance", "fixed"
                                      ),
                                      pc_max_condition = 1e4,
                                      pca_scale = c("covariance", "correlation")) {
  inverse_method <- match.arg(inverse_method)
  pca_rule <- match.arg(pca_rule)
  pca_scale <- match.arg(pca_scale)
  H_theta_beta <- as.matrix(H_theta_beta)
  H_beta_beta <- cohort_symmetrize(H_beta_beta)
  Sigma_beta_beta <- cohort_symmetrize(Sigma_beta_beta)
  M_raw <- cohort_symmetrize(H_beta_beta + Sigma_beta_beta)
  M_for_weight <- cohort_regularize_covariance(M_raw, ridge)
  p <- nrow(M_raw)

  M_model <- NULL
  M_model_raw <- NULL
  H_model <- NULL
  C_model <- H_theta_beta
  B <- NULL
  M_reduced <- NULL
  M_reduced_raw <- NULL
  C_reduced <- NULL
  factor_details <- NULL
  pca_details <- NULL

  if (inverse_method == "shrink") {
    stop(
      "inverse_method = 'shrink' is not implemented in this cohort revision. ",
      "Use solve/PCA/factor-plus-diagonal. Cross-fitting is also intentionally out of scope."
    )
  }

  if (inverse_method == "solve") {
    W <- cohort_right_solve(M_for_weight, H_theta_beta)
  } else if (inverse_method == "pca") {
    basis <- cohort_build_pc_basis(
      A_raw = M_raw,
      pca_scale = pca_scale,
      regularize_spectral = TRUE,
      ridge = ridge,
      pca_rule = pca_rule,
      pc_var = pc_var,
      n_pc = n_pc,
      eig_tol = eig_tol,
      pc_max_condition = pc_max_condition,
      label = "Total-covariance PCA basis"
    )
    B <- basis$B
    lambda <- basis$selected_values
    C_reduced <- H_theta_beta %*% B
    M_reduced <- diag(lambda, nrow = length(lambda))
    M_reduced_raw <- t(B) %*% M_raw %*% B
    W <- C_reduced %*% diag(1 / lambda, nrow = length(lambda)) %*% t(B)
    pca_details <- basis
  } else if (inverse_method %in% c("hbb_pca", "sigma_pca")) {
    A_basis <- if (inverse_method == "hbb_pca") H_beta_beta else Sigma_beta_beta
    basis <- cohort_build_pc_basis(
      A_raw = A_basis,
      pca_scale = pca_scale,
      regularize_spectral = FALSE,
      ridge = ridge,
      pca_rule = pca_rule,
      pc_var = pc_var,
      n_pc = n_pc,
      eig_tol = eig_tol,
      pc_max_condition = pc_max_condition,
      label = if (inverse_method == "hbb_pca") {
        "Internal-covariance PCA basis"
      } else {
        "External-covariance PCA basis"
      }
    )
    B <- basis$B
    C_reduced <- H_theta_beta %*% B
    M_reduced_raw <- cohort_symmetrize(t(B) %*% M_raw %*% B)
    M_reduced <- cohort_symmetrize(t(B) %*% M_for_weight %*% B)
    W <- cohort_right_solve(M_reduced, C_reduced) %*% t(B)
    pca_details <- basis
  } else if (inverse_method == "m_factor_diag") {
    basis <- cohort_build_pc_basis(
      A_raw = M_raw,
      pca_scale = pca_scale,
      regularize_spectral = FALSE,
      ridge = ridge,
      pca_rule = pca_rule,
      pc_var = pc_var,
      n_pc = n_pc,
      eig_tol = eig_tol,
      pc_max_condition = pc_max_condition,
      label = "Total-covariance factor basis"
    )
    factor_details <- cohort_build_factor_diag_covariance(M_raw, basis)
    M_model_raw <- factor_details$covariance
    M_model <- cohort_regularize_covariance(M_model_raw, ridge)
    W <- cohort_right_solve(M_model, H_theta_beta)
    pca_details <- basis
  } else if (inverse_method == "hbb_factor_diag") {
    basis <- cohort_build_pc_basis(
      A_raw = H_beta_beta,
      pca_scale = pca_scale,
      regularize_spectral = FALSE,
      ridge = ridge,
      pca_rule = pca_rule,
      pc_var = pc_var,
      n_pc = n_pc,
      eig_tol = eig_tol,
      pc_max_condition = pc_max_condition,
      label = "Internal-covariance factor basis"
    )
    factor_details <- cohort_build_factor_diag_covariance(H_beta_beta, basis)
    H_model <- factor_details$covariance
    M_model_raw <- cohort_symmetrize(H_model + Sigma_beta_beta)
    M_model <- cohort_regularize_covariance(M_model_raw, ridge)
    C_model <- cohort_project_crosscov_to_factor_space(H_theta_beta, basis)
    W <- cohort_right_solve(M_model, C_model)
    pca_details <- basis
  }

  dimnames(W) <- list(rownames(H_theta_beta), colnames(H_theta_beta))
  if (!is.null(B)) rownames(B) <- colnames(H_theta_beta)

  M_diag <- cohort_matrix_diagnostics(M_for_weight, eig_tol)
  diagnostics <- list(
    inverse_method = inverse_method,
    p_beta = p,
    ridge = ridge,
    pca_rule = if (inverse_method == "solve") NA_character_ else pca_rule,
    pca_scale = if (inverse_method == "solve") NA_character_ else pca_scale,
    pca_basis = if (inverse_method == "pca") {
      "total_covariance"
    } else if (inverse_method %in% c("hbb_pca", "hbb_factor_diag")) {
      "internal_covariance"
    } else if (inverse_method == "sigma_pca") {
      "external_covariance"
    } else if (inverse_method == "m_factor_diag") {
      "total_covariance"
    } else {
      NA_character_
    },
    n_pc_requested = if (is.null(n_pc)) NA_integer_ else as.integer(n_pc),
    n_pc_used = if (is.null(pca_details)) NA_integer_ else pca_details$selection$rank,
    n_positive_eigenvalues = if (is.null(pca_details)) NA_integer_ else
      pca_details$selection$positive_rank,
    pc_var_requested = if (is.null(pca_details)) NA_real_ else pc_var,
    pc_var_used = if (is.null(pca_details)) NA_real_ else
      pca_details$selection$pc_var_used,
    retained_condition_number = if (is.null(pca_details)) NA_real_ else
      pca_details$selection$retained_condition_number,
    M_beta_condition_number = M_diag$condition_number,
    pca_basis_raw_min_eigenvalue = if (is.null(pca_details)) NA_real_ else
      pca_details$selection$raw_min_eigenvalue,
    pca_basis_material_negative_eigenvalues = if (is.null(pca_details)) NA_integer_ else
      pca_details$selection$material_negative_eigenvalues,
    pca_reduced_total_condition_number = if (is.null(M_reduced)) NA_real_ else
      cohort_matrix_diagnostics(M_reduced, eig_tol)$condition_number,
    pca_reduction_matrix_rank = if (is.null(B)) NA_integer_ else
      cohort_matrix_diagnostics(crossprod(B), eig_tol)$rank,
    factor_diag_diagonal_preserved = if (is.null(factor_details)) NA else
      factor_details$diagonal_preserved,
    factor_diag_max_abs_diagonal_error = if (is.null(factor_details)) NA_real_ else
      factor_details$diagonal_error,
    factor_diag_residual_variance_min = if (is.null(factor_details)) NA_real_ else
      min(factor_details$residual_diag_spectral),
    factor_diag_offdiag_removed_frobenius = if (is.null(factor_details)) NA_real_ else
      factor_details$offdiag_removed_frobenius,
    factor_diag_H_theta_beta_projected = identical(inverse_method, "hbb_factor_diag")
  )

  list(
    W_hat = W,
    M_beta_raw = M_raw,
    M_beta_for_weight = M_for_weight,
    M_beta = M_for_weight,
    M_beta_weight_model = M_model,
    M_beta_weight_model_raw = M_model_raw,
    H_weight_beta_beta_model = H_model,
    H_weight_theta_beta_model = C_model,
    beta_reduction_matrix = B,
    M_beta_weight_reduced = M_reduced,
    M_beta_weight_reduced_raw = M_reduced_raw,
    H_weight_theta_beta_reduced = C_reduced,
    diagnostics = diagnostics
  )
}

# compute_zim_weight(): alias of cohort_compute_zim_weight with the same arguments
# and returned list.

compute_zim_weight <- cohort_compute_zim_weight

cohort_accumulate_point_adjustment <- function(theta_hat, W_hat,
                                               beta_hat, beta_tilde,
                                               beta_blocks) {
  theta <- as.numeric(theta_hat)
  names(theta) <- names(theta_hat)
  for (idx in beta_blocks) {
    theta <- theta - as.numeric(
      W_hat[, idx, drop = FALSE] %*%
        (beta_hat[idx] - beta_tilde[idx])
    )
  }
  names(theta) <- names(theta_hat)
  theta
}

# zim_cohort(): fit ZIM for one- or two-population cohort studies.
# internal_data: data frame containing the outcome and model covariates.
# Y: name of the outcome column.
# M_formula: formula specifying target-model covariates.
# working_models: list of formulas or lists with alpha and beta formulas.
# beta_tilde: external coefficient vector matching working-model beta terms.
# Sigma_beta_beta: optional covariance input on the n_internal * Cov(beta_tilde)
# scale; leave NULL to use the default estimator.
# n_external: sample size underlying beta_tilde; required when Sigma_beta_beta
# is NULL.
# outcome_type: "binary", "continuous", or "poisson".
# family: optional binomial(logit), gaussian(identity), or poisson(log) family.
# reference_data: external-population covariate data frame, or NULL for one
# population; reference data require beta_method = "GMM".
# W: optional weight matrix; ridge: nonnegative regularization parameter.
# beta_method: "GMM" or "WKM"; influence_method: "profiled" or "full".
# H_compute: "blocked" or "dense"; H_block_size: number of columns per block.
# H_block_cache: "disk", "recompute", or "memory"; H_cache_dir: cache directory.
# return_if: whether to include influence-function matrices in the result.
# gmm_solver: "newton" or "optim"; optim_control: optimizer control list.
# gmm_score_tol, gmm_step_tol, gmm_maxit, gmm_max_halving, and
# gmm_newton_fallback_to_optim: GMM solver controls.
# rank_tol, bread_condition_warn, joint_psd_tol, and final_psd_tol:
# diagnostic tolerances and warning threshold.
# inverse_method: "solve", "pca", "hbb_pca", "sigma_pca", "m_factor_diag",
# or "hbb_factor_diag"; "shrink" is reserved and currently unavailable.
# pc_var: retained-variance threshold; n_pc: optional fixed component count.
# eig_tol: eigenvalue tolerance; pc_max_condition: condition-number threshold.
# pca_rule: "variance", "condition", "condition_or_variance", or "fixed".
# pca_scale: "covariance" or "correlation"; poisson_clip: Poisson predictor limit.
# Returns a list containing theta_ZIM, se_ZIM, z_ZIM, p_ZIM, theta_MLE, se_MLE,
# sample sizes, fitted coefficients, covariance matrices, and diagnostics.

zim_cohort <- function(internal_data, Y, M_formula, working_models,
                       beta_tilde,
                       Sigma_beta_beta = NULL,
                       outcome_type = c("binary", "continuous", "poisson"),
                       family = NULL,
                       reference_data = NULL,
                       W = NULL,
                       ridge = 0,
                       optim_control = list(maxit = 1000),
                       beta_method = c("GMM", "WKM"),
                       influence_method = c("profiled", "full"),
                       H_compute = c("blocked", "dense"),
                       H_block_size = 256L,
                       H_block_cache = c("disk", "recompute", "memory"),
                       H_cache_dir = tempdir(),
                       return_if = FALSE,
                       gmm_solver = c("newton", "optim"),
                       gmm_score_tol = 1e-8,
                       gmm_step_tol = 1e-10,
                       gmm_maxit = 100L,
                       gmm_max_halving = 30L,
                       gmm_newton_fallback_to_optim = TRUE,
                       rank_tol = 1e-10,
                       bread_condition_warn = 1e12,
                       joint_psd_tol = 1e-8,
                       final_psd_tol = 1e-8,
                       inverse_method = c(
                         "solve", "pca", "hbb_pca", "sigma_pca",
                         "m_factor_diag", "hbb_factor_diag", "shrink"
                       ),
                       pc_var = 0.99,
                       n_pc = NULL,
                       eig_tol = 1e-10,
                       pca_rule = c(
                         "variance", "condition",
                         "condition_or_variance", "fixed"
                       ),
                       pc_max_condition = 1e4,
                       pca_scale = c("covariance", "correlation"),
                       poisson_clip = 35,
                       n_external = NULL) {
  pca_rule_was_missing <- missing(pca_rule)
  fam <- resolve_cohort_family(family = family, outcome_type = outcome_type)
  family <- fam$family
  outcome_type <- fam$outcome_type
  beta_method <- match.arg(beta_method)
  influence_method <- match.arg(influence_method)
  H_compute <- match.arg(H_compute)
  H_block_cache <- match.arg(H_block_cache)
  gmm_solver <- match.arg(gmm_solver)
  inverse_method <- match.arg(inverse_method)
  pca_scale <- match.arg(pca_scale)
  if (!is.null(n_pc) && pca_rule_was_missing) pca_rule <- "fixed"
  pca_rule <- match.arg(pca_rule)

  use_reference <- !is.null(reference_data)
  if (!is.null(n_external)) {
    n_external <- cohort_validate_external_sample_size(n_external)
  }
  if (is.null(Sigma_beta_beta) && is.null(n_external)) {
    stop("Supply n_external when Sigma_beta_beta is omitted; it is the sample size underlying beta_tilde.")
  }

  if (!is.character(Y) || length(Y) != 1L || !nzchar(Y)) {
    stop("Y must be the name of one outcome column.")
  }
  if (!is.logical(return_if) || length(return_if) != 1L || is.na(return_if)) {
    stop("return_if must be TRUE or FALSE.")
  }
  if (!is.numeric(poisson_clip) || length(poisson_clip) != 1L ||
      !is.finite(poisson_clip) || poisson_clip <= 0) {
    stop("poisson_clip must be a positive finite scalar.")
  }

  prep <- cohort_prepare_analysis_data(
    data = internal_data,
    Y = Y,
    M_formula = M_formula,
    working_models = working_models,
    data_label = "internal_data"
  )
  internal_used <- prep$data
  Y_int <- check_outcome(internal_used[[Y]], family)
  n <- length(Y_int)

  X_theta <- cohort_rhs_model_matrix(
    M_formula, internal_used, keep_intercept = FALSE,
    label = "target slope design"
  )
  target <- cohort_fit_target(
    Y = Y_int,
    X_no_intercept = X_theta,
    family = family,
    ridge = ridge,
    rank_tol = rank_tol,
    bread_condition_warn = bread_condition_warn,
    poisson_clip = poisson_clip
  )
  working <- cohort_pack_working_models(working_models, internal_used)
  beta_names <- attr(working, "beta_names")
  p_beta <- length(beta_names)

  beta_tilde <- cohort_align_beta_vector(
    beta_tilde, beta_names, object_name = "beta_tilde"
  )
  sigma_estimate <- NULL
  sigma_method <- "user_supplied"
  sigma_source <- "user_supplied"
  gamma_hat_GMM <- NULL
  gamma_hat_WKM <- NULL
  beta_hat_GMM <- NULL
  beta_hat_WKM <- NULL
  profile_diagnostics <- NULL
  bread_diagnostics <- NULL
  reference_prep <- NULL

  if (!use_reference) {
    states <- cohort_prepare_working_states(
      Y = Y_int,
      target = target,
      working = working,
      family = family,
      beta_method = beta_method,
      beta_tilde = beta_tilde,
      ridge = ridge,
      rank_tol = rank_tol,
      bread_condition_warn = bread_condition_warn,
      gmm_solver = gmm_solver,
      gmm_score_tol = gmm_score_tol,
      gmm_step_tol = gmm_step_tol,
      gmm_maxit = gmm_maxit,
      gmm_max_halving = gmm_max_halving,
      gmm_newton_fallback_to_optim = gmm_newton_fallback_to_optim,
      optim_control = optim_control,
      poisson_clip = poisson_clip
    )
    beta_hat_internal <- unlist(lapply(states, `[[`, "beta_hat"), use.names = TRUE)
    beta_hat_internal <- beta_hat_internal[beta_names]
    convergence <- lapply(states, `[[`, "convergence")
    if (beta_method == "GMM") {
      beta_hat_GMM <- beta_hat_internal
      gamma_hat_GMM <- lapply(states, `[[`, "gamma_hat")
    } else {
      beta_hat_WKM <- beta_hat_internal
      gamma_hat_WKM <- lapply(states, `[[`, "gamma_hat")
    }

    H <- cohort_build_H(
      IF_theta = target$IF_theta,
      states = states,
      target = target,
      influence_method = influence_method,
      H_compute = H_compute,
      H_block_size = H_block_size,
      H_block_cache = H_block_cache,
      H_cache_dir = H_cache_dir,
      beta_names = beta_names,
      joint_psd_tol = joint_psd_tol,
      return_if = return_if
    )
    profile_diagnostics <- lapply(states, function(s) {
      list(
        raw = s$profile$raw_profile_diagnostics,
        used = s$profile$used_profile_diagnostics
      )
    })
    bread_diagnostics <- lapply(states, `[[`, "bread_diagnostics")
  } else {
    if (beta_method != "GMM") {
      stop("reference_data is X-only, so beta_method = 'WKM' is unavailable.")
    }
    reference_prep <- cohort_prepare_analysis_data(
      data = reference_data,
      M_formula = M_formula,
      working_models = working_models,
      data_label = "reference_data"
    )
    reference_used <- reference_prep$data
    target_reference <- make_reference_target_cohort(
      reference_data = reference_used,
      M_formula = M_formula,
      lambda_hat = target$lambda_hat,
      family = family,
      poisson_clip = poisson_clip
    )
    working_reference <- cohort_pack_working_models(working_models, reference_used)
    beta_names_reference <- attr(working_reference, "beta_names")
    if (!identical(beta_names_reference, beta_names)) {
      stop("reference_data working-model beta columns do not match internal_data.")
    }
    if (length(working_reference) != length(working)) {
      stop("reference_data produced the wrong number of working-model designs.")
    }
    for (k in seq_along(working)) {
      if (!identical(working_reference[[k]]$gamma_names,
                     working[[k]]$gamma_names)) {
        stop(
          "reference_data working model ", k,
          " has design columns that do not match internal_data. ",
          "Use the same factor levels, contrasts, and formulas in both data sets."
        )
      }
    }
    ref_fit <- cohort_fit_reference_gmm(
      Y = Y_int,
      target_internal = target,
      target_reference = target_reference,
      working_internal = working,
      working_reference = working_reference,
      beta_tilde = beta_tilde,
      family = family,
      influence_method = influence_method,
      H_compute = H_compute,
      H_block_size = H_block_size,
      H_block_cache = H_block_cache,
      H_cache_dir = H_cache_dir,
      ridge = ridge,
      rank_tol = rank_tol,
      bread_condition_warn = bread_condition_warn,
      joint_psd_tol = joint_psd_tol,
      gmm_solver = gmm_solver,
      gmm_score_tol = gmm_score_tol,
      gmm_step_tol = gmm_step_tol,
      gmm_maxit = gmm_maxit,
      gmm_max_halving = gmm_max_halving,
      gmm_newton_fallback_to_optim = gmm_newton_fallback_to_optim,
      optim_control = optim_control,
      poisson_clip = poisson_clip,
      return_if = return_if
    )
    beta_hat_internal <- ref_fit$beta_hat[beta_names]
    beta_hat_GMM <- beta_hat_internal
    gamma_hat_GMM <- ref_fit$gamma_hat
    convergence <- ref_fit$convergence
    states <- ref_fit$states
    H <- ref_fit
    profile_diagnostics <- lapply(states, function(s) {
      list(
        raw = s$profile$raw_profile_diagnostics,
        used = s$profile$used_profile_diagnostics
      )
    })
    bread_diagnostics <- lapply(states, `[[`, "bread_diagnostics")
  }

  if (is.null(Sigma_beta_beta)) {
    if (use_reference) {
      sigma_estimate <- estimate_sigma_betabeta_reference_cohort(
        states = states, target_internal = target,
        target_reference = target_reference, n_external = n_external,
        H_compute = H_compute, H_block_size = H_block_size,
        H_block_cache = H_block_cache, H_cache_dir = H_cache_dir,
        joint_psd_tol = joint_psd_tol
      )
    } else {
      sigma_estimate <- estimate_sigma_betabeta_internal_cohort(
        internal_data = internal_used, Y = Y, working_models = working_models,
        n_external = n_external, family = family,
        influence_method = influence_method, H_compute = H_compute,
        H_block_size = H_block_size, H_block_cache = H_block_cache,
        H_cache_dir = H_cache_dir, rank_tol = rank_tol,
        bread_condition_warn = bread_condition_warn,
        joint_psd_tol = joint_psd_tol, poisson_clip = poisson_clip
      )
    }
    stopifnot(sigma_estimate$n_internal == n)
    Sigma_beta_beta <- sigma_estimate$Sigma_beta_beta
    sigma_method <- sigma_estimate$method
    sigma_source <- sigma_estimate$source
  }
  Sigma_aligned <- cohort_align_beta_covariance(
    Sigma = Sigma_beta_beta, beta_names = beta_names,
    object_name = "Sigma_beta_beta", psd_tol = joint_psd_tol
  )
  Sigma_beta_beta <- Sigma_aligned$matrix

  H_theta_theta <- H$H_theta_theta
  H_theta_beta <- H$H_theta_beta
  H_beta_beta <- H$H_beta_beta
  theta_names <- names(target$theta_hat)
  dimnames(H_theta_theta) <- list(theta_names, theta_names)
  dimnames(H_theta_beta) <- list(theta_names, beta_names)
  dimnames(H_beta_beta) <- list(beta_names, beta_names)

  if (is.null(W)) {
    inv <- cohort_compute_zim_weight(
      H_theta_beta = H_theta_beta,
      H_beta_beta = H_beta_beta,
      Sigma_beta_beta = Sigma_beta_beta,
      ridge = ridge,
      inverse_method = inverse_method,
      pc_var = pc_var,
      n_pc = n_pc,
      eig_tol = eig_tol,
      pca_rule = pca_rule,
      pc_max_condition = pc_max_condition,
      pca_scale = pca_scale
    )
    W_hat <- inv$W_hat
    weight_method <- paste0("estimated_", inverse_method)
    inverse_diagnostics <- inv$diagnostics
    M_beta_raw <- inv$M_beta_raw
    M_beta_for_weight <- inv$M_beta_for_weight
    M_beta_weight_model <- inv$M_beta_weight_model
    M_beta_weight_model_raw <- inv$M_beta_weight_model_raw
    H_weight_beta_beta_model <- inv$H_weight_beta_beta_model
    H_weight_theta_beta_model <- inv$H_weight_theta_beta_model
    beta_reduction_matrix <- inv$beta_reduction_matrix
    M_beta_weight_reduced <- inv$M_beta_weight_reduced
    M_beta_weight_reduced_raw <- inv$M_beta_weight_reduced_raw
    H_weight_theta_beta_reduced <- inv$H_weight_theta_beta_reduced
  } else {
    W_hat <- cohort_align_weight_matrix(W, theta_names, beta_names)
    weight_method <- "user_supplied"
    M_beta_raw <- cohort_symmetrize(H_beta_beta + Sigma_beta_beta)
    M_beta_for_weight <- cohort_regularize_covariance(M_beta_raw, ridge)
    M_beta_weight_model <- NULL
    M_beta_weight_model_raw <- NULL
    H_weight_beta_beta_model <- NULL
    H_weight_theta_beta_model <- NULL
    beta_reduction_matrix <- NULL
    M_beta_weight_reduced <- NULL
    M_beta_weight_reduced_raw <- NULL
    H_weight_theta_beta_reduced <- NULL
    inverse_diagnostics <- list(
      inverse_method = "none_user_supplied_W",
      p_beta = p_beta,
      ridge = ridge
    )
  }
  dimnames(W_hat) <- list(theta_names, beta_names)
  dimnames(M_beta_raw) <- list(beta_names, beta_names)
  dimnames(M_beta_for_weight) <- list(beta_names, beta_names)

  beta_blocks <- H$beta_blocks %||% cohort_make_beta_blocks(p_beta, H_block_size)
  theta_ZIM <- cohort_accumulate_point_adjustment(
    target$theta_hat, W_hat, beta_hat_internal, beta_tilde, beta_blocks
  )

  V_W_internal <- cohort_symmetrize(
    H_theta_theta -
      H_theta_beta %*% t(W_hat) -
      W_hat %*% t(H_theta_beta) +
      W_hat %*% H_beta_beta %*% t(W_hat)
  )
  V_W_external <- cohort_symmetrize(
    W_hat %*% Sigma_beta_beta %*% t(W_hat)
  )
  V_W <- cohort_symmetrize(V_W_internal + V_W_external)
  dimnames(V_W_internal) <- dimnames(H_theta_theta)
  dimnames(V_W_external) <- dimnames(H_theta_theta)
  dimnames(V_W) <- dimnames(H_theta_theta)
  V_W_psd_diagnostics <- cohort_psd_diagnostics(
    V_W,
    tol = final_psd_tol,
    label = "Final ZIM asymptotic covariance",
    action = "stop"
  )

  se_ZIM <- sqrt(pmax(diag(V_W) / n, 0))
  z_ZIM <- theta_ZIM / se_ZIM
  p_ZIM <- 2 * stats::pnorm(abs(z_ZIM), lower.tail = FALSE)
  se_MLE <- sqrt(pmax(diag(H_theta_theta) / n, 0))
  names(se_ZIM) <- names(z_ZIM) <- names(p_ZIM) <- names(theta_ZIM)
  names(se_MLE) <- names(target$theta_hat)

  list(
    theta_ZIM = theta_ZIM,
    se_ZIM = se_ZIM,
    z_ZIM = z_ZIM,
    p_ZIM = p_ZIM,
    theta_MLE = target$theta_hat,
    se_MLE = se_MLE,
    lambda_hat = target$lambda_hat,
    beta_method = beta_method,
    influence_method = influence_method,
    H_compute = H_compute,
    H_block_size = as.integer(H_block_size),
    H_block_cache = if (H_compute == "blocked") H_block_cache else NA_character_,
    weight_method = weight_method,
    family = family,
    outcome_type = outcome_type,
    reference_data_used = use_reference,
    n_internal = n,
    n_external = n_external %||% NA_real_,
    n_reference = if (use_reference) reference_prep$n_used else NA_integer_,
    data_diagnostics = prep,
    reference_data_diagnostics = reference_prep,
    beta_hat = beta_hat_internal,
    beta_hat_GMM = beta_hat_GMM,
    beta_hat_WKM = beta_hat_WKM,
    beta_tilde = beta_tilde,
    gamma_hat_GMM = gamma_hat_GMM,
    gamma_hat_WKM = gamma_hat_WKM,
    W_hat = W_hat,
    V_W = V_W,
    V_W_internal = V_W_internal,
    V_W_external = V_W_external,
    V_W_psd_diagnostics = V_W_psd_diagnostics,
    H_theta_theta = H_theta_theta,
    H_theta_beta = H_theta_beta,
    H_beta_beta = H_beta_beta,
    H_beta_beta_internal = if (use_reference) H$H_beta_beta_internal else NULL,
    H_beta_beta_reference = if (use_reference) H$H_beta_beta_reference else NULL,
    H_joint = H$H_joint,
    H_joint_psd_diagnostics = H$H_joint_psd_diagnostics,
    Sigma_beta_beta = Sigma_beta_beta,
    Sigma_beta_beta_method = sigma_method,
    Sigma_beta_beta_source = sigma_source,
    Sigma_beta_beta_estimation = if (is.null(sigma_estimate)) NULL else {
      metadata_names <- c(
        "method", "source", "n_internal", "n_external", "n_reference",
        "scaling_factor", "conditional_variance_assumption",
        "conditional_variance_range", "convergence", "bread_diagnostics",
        "profile_diagnostics", "H_compute_diagnostics", "covariance_diagnostics"
      )
      sigma_estimate[intersect(metadata_names, names(sigma_estimate))]
    },
    Sigma_beta_beta_diagnostics = Sigma_aligned$diagnostics,
    M_beta_raw = M_beta_raw,
    M_beta_for_weight = M_beta_for_weight,
    M_beta = M_beta_for_weight,
    M_beta_weight_model = M_beta_weight_model,
    M_beta_weight_model_raw = M_beta_weight_model_raw,
    H_weight_beta_beta_model = H_weight_beta_beta_model,
    H_weight_theta_beta_model = H_weight_theta_beta_model,
    beta_reduction_matrix = beta_reduction_matrix,
    M_beta_weight_reduced = M_beta_weight_reduced,
    M_beta_weight_reduced_raw = M_beta_weight_reduced_raw,
    H_weight_theta_beta_reduced = H_weight_theta_beta_reduced,
    inverse_diagnostics = inverse_diagnostics,
    target_diagnostics = target$diagnostics,
    bread_diagnostics = bread_diagnostics,
    profile_diagnostics = profile_diagnostics,
    convergence = convergence,
    H_compute_diagnostics = H$H_compute_diagnostics,
    beta_index_map = H$beta_index_map,
    beta_block_map = H$beta_blocks,
    IF_lambda = if (return_if) target$IF_lambda else NULL,
    IF_theta = if (return_if) target$IF_theta else NULL,
    IF_beta = if (return_if && !use_reference) H$IF_beta else NULL,
    IF_beta_internal = if (return_if && use_reference) H$IF_beta_internal else NULL,
    IF_beta_reference = if (return_if && use_reference) H$IF_beta_reference else NULL,
    call = match.call()
  )
}

# zim_cohort_wkm(): call zim_cohort() with beta_method = "WKM".
# ...: arguments passed to zim_cohort(), excluding beta_method.
# Returns the same list as zim_cohort().

zim_cohort_wkm <- function(...) {
  zim_cohort(..., beta_method = "WKM")
}
