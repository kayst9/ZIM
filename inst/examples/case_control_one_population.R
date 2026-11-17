example_file <- local({
  installed <- system.file("examples", "case_control_one_population.R", package = "ZIM")
  if (nzchar(installed)) {
    normalizePath(installed, mustWork = TRUE)
  } else {
    paths <- as.character(unlist(lapply(sys.frames(), function(frame) frame$ofile),
                                 use.names = FALSE))
    paths <- paths[basename(paths) == "case_control_one_population.R"]
    if (!length(paths)) {
      paths <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
    }
    if (!length(paths)) stop("Cannot locate case_control_one_population.R.")
    normalizePath(tail(paths, 1L), mustWork = TRUE)
  }
})

source(file.path(dirname(example_file), "case_control_example_helpers.R"), local = TRUE)

run_case_control_one_population <- function(example_dir = dirname(example_file)) {
  setup <- cc_example_setup(example_dir, "case_control_one_population")

  settings <- list(
    source_study = "package-example/case-control-one-population",
    p = 100L,
    m = 100L,
    internal_n_case = 1000L,
    internal_n_control = 1000L,
    external_n_case = 2000L,
    external_n_control = 2000L,
    rho = 0.90,
    internal_source_batch_size = 5000L,
    external_source_batch_size = 5000L,
    external_summary = "unadjusted",
    outcome_model = c("(Intercept)" = -0.014, S = 0.015),
    beta_method = "GMM",
    inverse_method = "hbb_pca",
    pc_var = 0.99,
    pc_max_condition = 200,
    pca_rule = "condition_or_variance",
    pca_scale = "covariance",
    eig_tol = 1e-10,
    ridge = 0,
    gmm_H_mode = "el",
    gmm_solver = "newton",
    gmm_score_tol = 1e-8,
    gmm_start_method_onepop = "glm",
    fixed_case_control = TRUE,
    H_B_centering = "separate",
    sigma_B_centering = "separate",
    shrink_structured = TRUE
  )

  started <- proc.time()[["elapsed"]]
  if (setup$generate) {
    message("One-population case-control: internal=", settings$internal_n_case, " cases + ",
            settings$internal_n_control, " controls, external=", settings$external_n_case,
            " cases + ", settings$external_n_control, " controls, p=", settings$p,
            ", m=", settings$m, "; unadjusted summaries.")
    design <- cc_example_design(example_dir, settings$p)
    message("Generating internal and external case-control samples.")
    internal_data <- cc_example_draw_study(
      n_case = settings$internal_n_case, n_control = settings$internal_n_control,
      rho = settings$rho, design = design,
      outcome_model = settings$outcome_model, batch_size = settings$internal_source_batch_size
    )
    external_data <- cc_example_draw_study(
      n_case = settings$external_n_case, n_control = settings$external_n_control,
      rho = settings$rho, design = design,
      outcome_model = settings$outcome_model, batch_size = settings$external_source_batch_size
    )
    generation <- list(internal = attr(internal_data, "generation"),
                       external = attr(external_data, "generation"))
    message("Estimating ", settings$m, " unadjusted external SNP coefficients.")
    beta_tilde <- cc_example_external_unadjusted(external_data, settings$m)
    names(beta_tilde) <- paste0("beta", seq_len(settings$m), ":x", seq_len(settings$m))
    rm(external_data)
    invisible(gc())
    inputs <- list(internal_data = internal_data, reference_data = NULL,
                   beta_tilde = beta_tilde, settings = settings,
                   generation = generation, design = design)
    saveRDS(inputs, setup$data_file, compress = "xz")
    message("Saved input dataset to ", setup$data_file)
    rm(inputs)
  } else {
    inputs <- readRDS(setup$data_file)
    internal_data <- inputs$internal_data
    beta_tilde <- inputs$beta_tilde
    settings <- inputs$settings
    generation <- inputs$generation
    rm(inputs)
    message("Loaded input dataset from ", setup$data_file)
  }
  working_models <- lapply(paste0("x", seq_len(settings$m)), function(xj) {
    list(alpha = ~1, beta = stats::reformulate(xj))
  })

  message("Fitting ZIM-GMM with Hbb PCA variance target 0.99 and condition cap 200.")
  fit <- setup$core$zim_case_control(
    internal_data = internal_data,
    Y = "Y",
    M_formula = ~S,
    working_models = working_models,
    beta_tilde = beta_tilde,
    Sigma_beta_beta = NULL,
    reference_data = NULL,
    external_n_case = settings$external_n_case,
    external_n_control = settings$external_n_control,
    beta_method = settings$beta_method,
    inverse_method = settings$inverse_method,
    pca_scale = settings$pca_scale,
    pca_rule = settings$pca_rule,
    pc_var = settings$pc_var,
    pc_max_condition = settings$pc_max_condition,
    eig_tol = settings$eig_tol,
    ridge = settings$ridge,
    gmm_H_mode = settings$gmm_H_mode,
    gmm_solver = settings$gmm_solver,
    gmm_score_tol = settings$gmm_score_tol,
    gmm_start_method_onepop = settings$gmm_start_method_onepop,
    fixed_case_control = settings$fixed_case_control,
    H_B_centering = settings$H_B_centering,
    sigma_B_centering = settings$sigma_B_centering,
    shrink_structured = settings$shrink_structured
  )
  elapsed_seconds <- proc.time()[["elapsed"]] - started
  cc_example_save(fit, setup, settings, generation, beta_tilde, elapsed_seconds)
}

if (sys.nframe() == 0L) run_case_control_one_population()
