example_file <- local({
  installed <- system.file("examples", "cohort_one_population.R", package = "ZIM")
  if (nzchar(installed)) {
    normalizePath(installed, mustWork = TRUE)
  } else {
    paths <- as.character(unlist(lapply(sys.frames(), function(frame) frame$ofile),
                                 use.names = FALSE))
    paths <- paths[basename(paths) == "cohort_one_population.R"]
    if (!length(paths)) {
      paths <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
    }
    if (!length(paths)) stop("Cannot locate cohort_one_population.R.")
    normalizePath(tail(paths, 1L), mustWork = TRUE)
  }
})

source(file.path(dirname(example_file), "cohort_example_helpers.R"), local = TRUE)

run_cohort_one_population <- function(example_dir = dirname(example_file)) {
  setup <- example_setup(example_dir, "cohort_one_population")

  settings <- list(
    source_study = "package-example/cohort-one-population",
    p = 100L,
    m = 100L,
    n_internal = 2000L,
    n_external = 4000L,
    rho = 0.90,
    external_summary = "unadjusted",
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

  started <- proc.time()[["elapsed"]]
  if (setup$generate) {
    design <- example_design(p = settings$p, example_dir = example_dir)
    message("Generating internal and external cohorts.")
    internal_data <- example_draw_cohort(
      n = settings$n_internal, rho = settings$rho, design = design,
      outcome_model = settings$outcome_model
    )
    external_data <- example_draw_cohort(
      n = settings$n_external, rho = settings$rho, design = design,
      outcome_model = settings$outcome_model
    )
    generation <- list(internal = attr(internal_data, "generation"),
                       external = attr(external_data, "generation"))
    message("Estimating ", settings$m, " unadjusted external SNP coefficients.")
    beta_tilde <- example_external_unadjusted(external_data, settings$m)
    rm(external_data)
    saveRDS(list(internal_data = internal_data, reference_data = NULL,
                 beta_tilde = beta_tilde, settings = settings,
                 generation = generation, design = design),
            setup$data_file, compress = "xz")
    message("Saved analysis inputs to ", setup$data_file)
  } else {
    inputs <- readRDS(setup$data_file)
    internal_data <- inputs$internal_data
    beta_tilde <- inputs$beta_tilde
    settings <- inputs$settings
    settings[names(pca_settings)] <- pca_settings
    generation <- inputs$generation
    design <- inputs$design
    rm(inputs)
    message("Loaded analysis inputs from ", setup$data_file)
  }
  invisible(gc())
  working_models <- lapply(paste0("x", seq_len(settings$m)), function(xj) {
    list(alpha = ~1, beta = stats::reformulate(xj))
  })
  cache_dir <- tempfile("cohort-onepop-h-cache-")
  if (!dir.create(cache_dir)) stop("Cannot create the temporary H cache directory.")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)

  message("Fitting ZIM-GMM with Hbb PCA variance target ", settings$pc_var,
          " and condition cap ", settings$pc_max_condition, ".")
  fit <- setup$core$zim_cohort(
    internal_data = internal_data,
    Y = "Y",
    M_formula = ~S,
    working_models = working_models,
    beta_tilde = beta_tilde,
    n_external = settings$n_external,
    outcome_type = "binary",
    reference_data = NULL,
    Sigma_beta_beta = NULL,
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
  elapsed_seconds <- proc.time()[["elapsed"]] - started
  example_save(fit, setup, settings, generation, beta_tilde, elapsed_seconds)
}

if (sys.nframe() == 0L) run_cohort_one_population()
