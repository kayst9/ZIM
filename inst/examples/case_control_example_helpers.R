cc_example_package_file <- function(example_dir, ...) {
  installed <- system.file(..., package = "ZIM")
  if (nzchar(installed)) return(installed)
  root <- normalizePath(file.path(example_dir, "..", ".."), mustWork = TRUE)
  file.path(root, "inst", ...)
}

cc_example_setup <- function(example_dir, example_name) {
  args <- commandArgs(trailingOnly = TRUE)
  allowed <- grepl("^--(zim-source|output-dir|data-file)=.+$", args) | args == "--generate"
  if (any(!allowed)) {
    stop("Use --zim-source=PATH, --output-dir=PATH, --data-file=PATH, and/or --generate.")
  }
  if (sum(args == "--generate") > 1L) stop("Repeated option: generate")
  generate <- "--generate" %in% args
  option <- function(name, default) {
    prefix <- paste0("--", name, "=")
    values <- args[startsWith(args, prefix)]
    if (length(values) > 1L) stop("Repeated option: ", name)
    if (length(values)) substring(values, nchar(prefix) + 1L) else default
  }
  default_source <- cc_example_package_file(example_dir, "scripts", "zim_case_control.R")
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
  default_data <- if (generate) {
    file.path(output_dir, paste0(example_name, ".rds"))
  } else {
    cc_example_package_file(example_dir, "extdata", paste0(example_name, ".rds"))
  }
  data_file <- option("data-file", default_data)
  if (generate) dir.create(dirname(data_file), recursive = TRUE, showWarnings = FALSE)
  data_file <- normalizePath(data_file, mustWork = FALSE)
  if (!generate && !file.exists(data_file)) {
    stop("Input dataset does not exist: ", data_file, ". Use --generate to create it.")
  }
  script_path <- file.path(example_dir, paste0(example_name, ".R"))
  helper_path <- file.path(example_dir, "case_control_example_helpers.R")
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
  list(core = core, output_dir = output_dir, provenance = provenance,
       data_file = data_file, generate = generate)
}

cc_example_design <- function(example_dir, p = 100L) {
  path <- cc_example_package_file(example_dir, "extdata", "fixed_snp_design.csv")
  values <- read.csv(path, check.names = FALSE)
  stopifnot(all(c("omega", "fre") %in% names(values)), nrow(values) >= p,
            all(is.finite(values$omega)), all(is.finite(values$fre)),
            all(values$fre > 0 & values$fre < 1))
  values <- values[seq_len(p), , drop = FALSE]
  list(omega = values$omega, fre = values$fre)
}

cc_example_draw_study <- function(n_case, n_control, rho, design, outcome_model,
                                  batch_size) {
  counts <- c(n_case, n_control)
  stopifnot(length(counts) == 2L, all(is.finite(counts)), all(counts >= 0),
            all(counts == as.integer(counts)), sum(counts) > 0,
            length(batch_size) == 1L, batch_size > 0, batch_size == as.integer(batch_size))
  p <- length(design$fre)
  cases <- controls <- list()
  have_case <- have_control <- 0L
  n_batches <- n_source <- n_source_cases <- 0
  while (have_case < n_case || have_control < n_control) {
    latent <- matrix(rnorm(batch_size * p), nrow = batch_size, ncol = p)
    if (p > 1L) {
      innovation_sd <- sqrt(1 - rho^2)
      for (j in 2:p) latent[, j] <- rho * latent[, j - 1L] + innovation_sd * latent[, j]
    }
    u <- pnorm(latent)
    rm(latent)
    X <- matrix(0L, nrow = batch_size, ncol = p)
    fre1 <- design$fre * design$fre
    fre2 <- design$fre * (2 - design$fre)
    for (j in seq_len(p)) {
      X[, j] <- as.integer(u[, j] < fre2[j]) + as.integer(u[, j] < fre1[j])
    }
    rm(u)
    colnames(X) <- paste0("x", seq_len(p))
    S <- drop(X %*% design$omega)
    prob_y <- plogis(outcome_model[["(Intercept)"]] + outcome_model[["S"]] * S)
    Y <- rbinom(batch_size, size = 1, prob = prob_y)
    n_batches <- n_batches + 1L
    n_source <- n_source + batch_size
    n_source_cases <- n_source_cases + sum(Y)
    if (have_case < n_case) {
      take <- head(which(Y == 1), n_case - have_case)
      if (length(take)) {
        cases[[length(cases) + 1L]] <- data.frame(
          Y = Y[take], S = S[take],
          X[take, , drop = FALSE], check.names = FALSE)
        have_case <- have_case + length(take)
      }
    }
    if (have_control < n_control) {
      take <- head(which(Y == 0), n_control - have_control)
      if (length(take)) {
        controls[[length(controls) + 1L]] <- data.frame(
          Y = if (n_case == 0L) 0 else Y[take], S = S[take],
          X[take, , drop = FALSE], check.names = FALSE)
        have_control <- have_control + length(take)
      }
    }
    if (n_batches %% 5L == 0L) {
      message("Collected ", have_case, "/", n_case, " cases and ", have_control,
              "/", n_control, " controls from ", n_source, " source individuals.")
    }
  }
  out <- do.call(rbind, c(cases, controls))
  rownames(out) <- NULL
  out <- out[sample(seq_len(nrow(out))), , drop = FALSE]
  attr(out, "generation") <- list(n_case = n_case, n_control = n_control,
    rho = rho, batch_size = batch_size, n_batches = n_batches, n_source = n_source,
    n_source_cases = n_source_cases, source_disease_prevalence = n_source_cases / n_source)
  out
}

cc_example_external_unadjusted <- function(dat, m, control = glm.control()) {
  beta <- setNames(numeric(m), paste0("x", seq_len(m)))
  for (j in seq_len(m)) {
    D <- cbind(1, dat[[names(beta)[j]]])
    fit <- glm.fit(D, dat$Y, family = binomial(), control = control)
    if (!isTRUE(fit$converged) || any(!is.finite(fit$coefficients))) {
      stop("External regression failed for ", names(beta)[j])
    }
    if (isTRUE(fit$boundary)) warning("External regression reached a boundary for ", names(beta)[j])
    beta[j] <- unname(fit$coefficients[2L])
    if (j %% 250L == 0L) message("External coefficients: ", j, "/", m)
  }
  beta
}

cc_example_save <- function(fit, setup, settings, generation, beta_tilde, elapsed_seconds) {
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
    covariance_source = if (isTRUE(fit$reference_data_used)) "internal_and_reference_controls" else "internal_data",
    covariance_method = fit$Sigma_beta_beta_method
  ))
  result <- list(settings = settings, generation = generation, estimates = estimates,
                 beta_tilde = beta_tilde, provenance = provenance,
                 diagnostics = list(inverse = fit$inverse_diagnostics,
                   target = fit$target_diagnostics, convergence = fit$convergence,
                   covariance = list(method = fit$Sigma_beta_beta_method,
                     reference_reconstructed = fit$reference_sigma_reconstructed,
                     external_n_case = fit$external_n_case,
                     external_n_control = fit$external_n_control,
                     n_internal = fit$n_internal, n_reference = fit$n_reference)))
  write.csv(estimates, file.path(setup$output_dir, "estimates.csv"), row.names = FALSE)
  saveRDS(result, file.path(setup$output_dir, "result.rds"))
  write.csv(data.frame(field = names(provenance), value = vapply(provenance,
    function(x) paste(x, collapse = "; "), character(1))),
    file.path(setup$output_dir, "provenance.csv"), row.names = FALSE)
  print(estimates, row.names = FALSE, digits = 6)
  message("Saved results to ", setup$output_dir)
  invisible(estimates)
}
