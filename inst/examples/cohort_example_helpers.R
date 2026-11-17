example_package_file <- function(example_dir, ...) {
  installed <- system.file(..., package = "ZIM")
  if (nzchar(installed)) return(installed)
  root <- normalizePath(file.path(example_dir, "..", ".."), mustWork = TRUE)
  file.path(root, "inst", ...)
}

example_setup <- function(example_dir, example_name) {
  args <- commandArgs(trailingOnly = TRUE)
  allowed <- args == "--generate" | grepl("^--(zim-source|output-dir|data-file)=.+$", args)
  if (any(!allowed)) {
    stop("Use --generate, --data-file=PATH, --zim-source=PATH, and/or --output-dir=PATH.")
  }
  if (sum(args == "--generate") > 1L) stop("Repeated option: generate")
  option <- function(name, default) {
    prefix <- paste0("--", name, "=")
    values <- args[startsWith(args, prefix)]
    if (length(values) > 1L) stop("Repeated option: ", name)
    if (length(values)) substring(values, nchar(prefix) + 1L) else default
  }
  default_source <- example_package_file(example_dir, "scripts", "zim_cohort.R")
  source_path <- normalizePath(option("zim-source", default_source), mustWork = TRUE)
  source_md5 <- unname(tools::md5sum(source_path))
  core <- new.env(parent = globalenv())
  sys.source(source_path, envir = core)
  if (!identical(source_md5, unname(tools::md5sum(source_path)))) {
    stop("The program source changed while it was loading.")
  }
  output_dir <- option("output-dir", file.path(getwd(), "zim-example-output", example_name))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_dir <- normalizePath(output_dir, mustWork = TRUE)
  default_data <- if ("--generate" %in% args) {
    file.path(output_dir, paste0(example_name, ".rds"))
  } else {
    example_package_file(example_dir, "extdata", paste0(example_name, ".rds"))
  }
  data_file <- option("data-file", default_data)
  generate <- "--generate" %in% args
  if (generate) {
    dir.create(dirname(data_file), recursive = TRUE, showWarnings = FALSE)
  } else if (!file.exists(data_file)) {
    stop("Example dataset not found: ", data_file, ". Use --generate to create it.")
  }
  data_file <- file.path(normalizePath(dirname(data_file), mustWork = TRUE),
                         basename(data_file))
  script_path <- file.path(example_dir, paste0(example_name, ".R"))
  helper_path <- file.path(example_dir, "cohort_example_helpers.R")
  provenance <- list(
    program_source = source_path,
    program_md5 = source_md5,
    example_source = normalizePath(script_path, mustWork = TRUE),
    example_md5 = unname(tools::md5sum(script_path)),
    helper_md5 = unname(tools::md5sum(helper_path)),
    R_version = R.version.string,
    RNG_kind = paste(RNGkind(), collapse = "; "),
    started_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
  )
  message("Program: ", source_path, " (MD5 ", source_md5, ")")
  list(core = core, output_dir = output_dir, data_file = data_file,
       generate = generate, provenance = provenance)
}

example_design <- function(p = 100L, example_dir) {
  design_file <- example_package_file(example_dir, "extdata", "fixed_snp_design.csv")
  values <- read.csv(design_file, check.names = FALSE)
  if (!all(c("omega", "fre") %in% names(values)) || nrow(values) < p ||
      any(!is.finite(values$omega)) || any(!is.finite(values$fre)) ||
      any(values$fre <= 0 | values$fre >= 1)) {
    stop("Invalid fixed SNP design: ", design_file)
  }
  values <- values[seq_len(p), , drop = FALSE]
  list(omega = values$omega, fre = values$fre)
}

example_draw_cohort <- function(n, rho, design, outcome_model, include_y = TRUE) {
  p <- length(design$fre)
  latent <- matrix(rnorm(n * p), nrow = n, ncol = p)
  innovation_sd <- sqrt(1 - rho^2)
  if (p > 1L) for (j in 2:p) {
    latent[, j] <- rho * latent[, j - 1L] + innovation_sd * latent[, j]
  }
  u <- pnorm(latent)
  rm(latent)
  X <- matrix(0L, nrow = n, ncol = p)
  for (j in seq_len(p)) {
    X[, j] <- as.integer(u[, j] < design$fre[j] * (2 - design$fre[j])) +
      as.integer(u[, j] < design$fre[j]^2)
  }
  rm(u)
  colnames(X) <- paste0("x", seq_len(p))
  S <- drop(X %*% design$omega)
  prob_y <- plogis(outcome_model[["(Intercept)"]] + outcome_model[["S"]] * S)
  out <- data.frame(S = S, X, check.names = FALSE)
  if (include_y) out <- data.frame(Y = rbinom(n, 1L, prob_y), out, check.names = FALSE)
  attr(out, "generation") <- list(
    n = n, rho = rho, outcome_prevalence = if (include_y) mean(out$Y) else NA_real_,
    mean_outcome_probability = mean(prob_y), mean_score = mean(S), sd_score = sd(S)
  )
  out
}

example_external_unadjusted <- function(dat, m) {
  beta <- setNames(numeric(m), paste0("x", seq_len(m)))
  for (j in seq_len(m)) {
    model <- glm.fit(cbind(1, dat[[names(beta)[j]]]), dat$Y, family = binomial())
    if (!isTRUE(model$converged) || any(!is.finite(model$coefficients))) {
      stop("External regression failed for ", names(beta)[j])
    }
    if (isTRUE(model$boundary)) warning("External regression reached a boundary for ", names(beta)[j])
    beta[j] <- unname(tail(model$coefficients, 1L))
    if (j %% 250L == 0L) message("External coefficients: ", j, "/", m)
  }
  beta
}

example_save <- function(fit, setup, settings, generation, beta_tilde, elapsed_seconds) {
  make_rows <- function(method, estimates, standard_errors) {
    standard_errors <- standard_errors[names(estimates)]
    data.frame(
      method = method, term = names(estimates), estimate = unname(estimates),
      standard_error = unname(standard_errors),
      ci_lower = unname(estimates - qnorm(0.975) * standard_errors),
      ci_upper = unname(estimates + qnorm(0.975) * standard_errors),
      p_value = unname(2 * pnorm(abs(estimates / standard_errors), lower.tail = FALSE))
    )
  }
  estimates <- rbind(make_rows("Internal MLE", fit$theta_MLE, fit$se_MLE),
                     make_rows("ZIM-GMM", fit$theta_ZIM, fit$se_ZIM))
  provenance <- c(setup$provenance, list(
    input_dataset = normalizePath(setup$data_file, mustWork = TRUE),
    input_dataset_md5 = unname(tools::md5sum(setup$data_file)),
    input_mode = if (setup$generate) "generated" else "saved_dataset",
    finished_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    elapsed_seconds = elapsed_seconds,
    covariance_source = fit$Sigma_beta_beta_source,
    covariance_method = fit$Sigma_beta_beta_method
  ))
  result <- list(settings = settings, generation = generation, estimates = estimates,
                 beta_tilde = beta_tilde, provenance = provenance,
                 diagnostics = list(inverse = fit$inverse_diagnostics,
                   target = fit$target_diagnostics, convergence = fit$convergence,
                   covariance = fit$Sigma_beta_beta_estimation))
  write.csv(estimates, file.path(setup$output_dir, "estimates.csv"), row.names = FALSE)
  saveRDS(result, file.path(setup$output_dir, "result.rds"))
  write.csv(data.frame(field = names(provenance), value = vapply(provenance,
    function(x) paste(x, collapse = "; "), character(1))),
    file.path(setup$output_dir, "provenance.csv"), row.names = FALSE)
  print(estimates, row.names = FALSE, digits = 6)
  message("Saved results to ", setup$output_dir)
  invisible(estimates)
}
