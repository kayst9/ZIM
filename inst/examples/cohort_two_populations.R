example_file <- local({
  installed <- system.file("examples", "cohort_two_populations.R", package = "ZIM")
  if (nzchar(installed)) {
    normalizePath(installed, mustWork = TRUE)
  } else {
    paths <- as.character(unlist(lapply(sys.frames(), function(frame) frame$ofile),
                                 use.names = FALSE))
    paths <- paths[basename(paths) == "cohort_two_populations.R"]
    if (!length(paths)) {
      paths <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
    }
    if (!length(paths)) stop("Cannot locate cohort_two_populations.R.")
    normalizePath(tail(paths, 1L), mustWork = TRUE)
  }
})

source(file.path(dirname(example_file), "cohort_example_helpers.R"), local = TRUE)

run_cohort_two_populations <- function(example_dir = dirname(example_file)) {
  started <- proc.time()[["elapsed"]]
  setup <- example_setup(example_dir, "cohort_two_populations")
  settings <- list(
    source_study = "package-example/cohort-two-populations",
    p = 100L,
    m = 100L,
    n_internal = 2000L,
    n_external = 4000L,
    n_reference = 500L,
    rho_internal = 0.90,
    rho_external = 0.95,
    rho_reference = 0.95,
    reference_population = "external",
    summary_adjustment = "none",
    outcome_model = c("(Intercept)" = -0.014, S = 0.015),
    beta_method = "GMM",
    inverse_method = "hbb_pca",
    pc_var = 0.999,
    pc_max_condition = 500,
    pca_rule = "condition_or_variance",
    pca_scale = "covariance",
    eig_tol = 1e-10,
    ridge = 0,
    influence_method = "profiled",
    H_compute = "blocked",
    H_block_size = 100L,
    H_block_cache = "disk",
    gmm_solver = "newton",
    gmm_score_tol = 1e-8,
    gmm_step_tol = 1e-10,
    gmm_maxit = 100L,
    gmm_max_halving = 30L,
    gmm_newton_fallback_to_optim = TRUE
  )
  pca_settings <- settings[c("pc_var", "pc_max_condition", "pca_rule")]
  if (setup$generate) {
    design <- example_design(settings$p, example_dir = example_dir)
    message("Generating internal cohort (n = ", settings$n_internal, ").")
    internal <- example_draw_cohort(
      settings$n_internal, settings$rho_internal, design,
      settings$outcome_model
    )
    message("Generating external cohort (n = ", settings$n_external, ").")
    external <- example_draw_cohort(
      settings$n_external, settings$rho_external, design,
      settings$outcome_model
    )
    message("Generating external-population reference (n = ",
            settings$n_reference, ").")
    reference <- example_draw_cohort(
      settings$n_reference, settings$rho_reference, design,
      settings$outcome_model, include_y = FALSE
    )
    generation <- list(
      internal = attr(internal, "generation"),
      external = attr(external, "generation"),
      reference = attr(reference, "generation")
    )
    message("Estimating ", settings$m, " unadjusted external SNP coefficients.")
    beta_tilde <- example_external_unadjusted(external, settings$m)
    rm(external)
    saveRDS(list(internal_data = internal, reference_data = reference,
                 beta_tilde = beta_tilde, settings = settings,
                 generation = generation, design = design),
            setup$data_file, compress = "xz")
    message("Saved analysis inputs to ", setup$data_file)
  } else {
    inputs <- readRDS(setup$data_file)
    internal <- inputs$internal_data
    reference <- inputs$reference_data
    beta_tilde <- inputs$beta_tilde
    settings <- inputs$settings
    settings[names(pca_settings)] <- pca_settings
    generation <- inputs$generation
    design <- inputs$design
    rm(inputs)
    message("Loaded analysis inputs from ", setup$data_file)
  }
  invisible(gc())
  working_models <- lapply(paste0("x", seq_len(settings$m)), function(snp) {
    list(alpha = ~ 1, beta = stats::reformulate(snp))
  })
  cache_dir <- tempfile("cohort-two-populations-h-cache-")
  dir.create(cache_dir, recursive = TRUE)
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)
  message("Fitting ZIM-GMM with external-population reference data.")
  fit <- setup$core$zim_cohort(
    internal_data = internal,
    Y = "Y",
    M_formula = ~ S,
    working_models = working_models,
    beta_tilde = beta_tilde,
    Sigma_beta_beta = NULL,
    n_external = settings$n_external,
    outcome_type = "binary",
    reference_data = reference,
    beta_method = settings$beta_method,
    inverse_method = settings$inverse_method,
    pc_var = settings$pc_var,
    pc_max_condition = settings$pc_max_condition,
    pca_rule = settings$pca_rule,
    pca_scale = settings$pca_scale,
    eig_tol = settings$eig_tol,
    ridge = settings$ridge,
    influence_method = settings$influence_method,
    H_compute = settings$H_compute,
    H_block_size = settings$H_block_size,
    H_block_cache = settings$H_block_cache,
    H_cache_dir = cache_dir,
    return_if = FALSE,
    gmm_solver = settings$gmm_solver,
    gmm_score_tol = settings$gmm_score_tol,
    gmm_step_tol = settings$gmm_step_tol,
    gmm_maxit = settings$gmm_maxit,
    gmm_max_halving = settings$gmm_max_halving,
    gmm_newton_fallback_to_optim = settings$gmm_newton_fallback_to_optim
  )
  example_save(fit, setup, settings, generation, beta_tilde,
               elapsed_seconds = proc.time()[["elapsed"]] - started)
  invisible(fit)
}

if (sys.nframe() == 0L) run_cohort_two_populations()
