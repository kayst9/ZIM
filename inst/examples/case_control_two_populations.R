example_file <- local({
  installed <- system.file("examples", "case_control_two_populations.R", package = "ZIM")
  if (nzchar(installed)) {
    normalizePath(installed, mustWork = TRUE)
  } else {
    paths <- as.character(unlist(lapply(sys.frames(), function(frame) frame$ofile),
                                 use.names = FALSE))
    paths <- paths[basename(paths) == "case_control_two_populations.R"]
    if (!length(paths)) {
      paths <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
    }
    if (!length(paths)) stop("Cannot locate case_control_two_populations.R.")
    normalizePath(tail(paths, 1L), mustWork = TRUE)
  }
})

source(file.path(dirname(example_file), "case_control_example_helpers.R"), local = TRUE)

run_case_control_two_populations <- function(example_dir = dirname(example_file)) {
  started <- proc.time()[["elapsed"]]
  setup <- cc_example_setup(example_dir, "case_control_two_populations")
  settings <- list(
    source_study = "package-example/case-control-two-populations",
    p = 100L,
    m = 100L,
    n_blocks = 1L,
    block_size = 100L,
    weight_range = c(-0.03, 0.07),
    maf_range = c(0.1, 0.4),
    internal_cases = 1000L,
    internal_controls = 1000L,
    n_internal = 2000L,
    external_cases = 2000L,
    external_controls = 2000L,
    n_external = 4000L,
    n_reference = 500L,
    reference_cases = 0L,
    reference_controls = 500L,
    rho_internal = 0.90,
    rho_external = 0.95,
    rho_reference = 0.95,
    batch_internal = 5000L,
    batch_external = 5000L,
    batch_reference = 2000L,
    reference_population = "external_controls",
    summary_adjustment = "none",
    summary_selection = "All 100 SNPs",
    summary_glm_epsilon = 1e-10,
    summary_glm_maxit = 100L,
    outcome_model = c("(Intercept)" = -0.014, S = 0.015),
    beta_method = "GMM",
    inverse_method = "hbb_pca",
    pc_var = 0.999,
    pc_max_condition = 500,
    pca_rule = "condition_or_variance",
    pca_scale = "covariance",
    eig_tol = 1e-10,
    n_pc = NULL,
    ridge = 0,
    fixed_case_control = TRUE,
    H_B_centering = "separate",
    sigma_B_centering = "separate",
    gmm_start_method_reference = "profile_alpha",
    gmm_H_mode = "el",
    gmm_solver = "newton",
    gmm_score_tol = 1e-8
  )
  pca_settings <- settings[c("pc_var", "pc_max_condition", "pca_rule")]
  if (setup$generate) {
    design <- cc_example_design(example_dir, settings$p)
    message("Generating internal study (", settings$internal_cases,
            " cases and ", settings$internal_controls, " controls).")
    internal <- cc_example_draw_study(
      settings$internal_cases, settings$internal_controls, settings$rho_internal,
      design, settings$outcome_model, settings$batch_internal
    )
    message("Generating external study (", settings$external_cases,
            " cases and ", settings$external_controls, " controls).")
    external <- cc_example_draw_study(
      settings$external_cases, settings$external_controls, settings$rho_external,
      design, settings$outcome_model, settings$batch_external
    )
    message("Generating independent external-population reference (",
            settings$n_reference, " controls).")
    reference <- cc_example_draw_study(
      settings$reference_cases, settings$reference_controls, settings$rho_reference,
      design, settings$outcome_model, settings$batch_reference
    )
    generation <- list(
      internal = attr(internal, "generation"),
      external = attr(external, "generation"),
      reference = attr(reference, "generation")
    )
    reference$Y <- NULL
    message("Estimating ", settings$m, " unadjusted external SNP coefficients.")
    beta_tilde <- cc_example_external_unadjusted(
      external, settings$m,
      control = stats::glm.control(epsilon = settings$summary_glm_epsilon,
                                   maxit = settings$summary_glm_maxit)
    )
    rm(external)
    invisible(gc())
    inputs <- list(internal_data = internal, reference_data = reference,
                   beta_tilde = beta_tilde, settings = settings,
                   generation = generation, design = design)
    saveRDS(inputs, setup$data_file, compress = "xz")
    message("Saved input dataset to ", setup$data_file)
    rm(inputs)
  } else {
    inputs <- readRDS(setup$data_file)
    internal <- inputs$internal_data
    reference <- inputs$reference_data
    beta_tilde <- inputs$beta_tilde
    settings <- inputs$settings
    settings[names(pca_settings)] <- pca_settings
    generation <- inputs$generation
    rm(inputs)
    message("Loaded input dataset from ", setup$data_file)
  }
  working_models <- lapply(paste0("x", seq_len(settings$m)), function(snp) {
    list(alpha = ~ 1, beta = stats::reformulate(snp))
  })
  message("Fitting ZIM-GMM with external-population reference controls.")
  fit <- setup$core$zim_case_control(
    internal_data = internal,
    Y = "Y",
    M_formula = ~ S,
    working_models = working_models,
    beta_tilde = beta_tilde,
    Sigma_beta_beta = NULL,
    external_n_case = settings$external_cases,
    external_n_control = settings$external_controls,
    reference_data = reference,
    beta_method = settings$beta_method,
    inverse_method = settings$inverse_method,
    pca_scale = settings$pca_scale,
    pca_rule = settings$pca_rule,
    pc_var = settings$pc_var,
    pc_max_condition = settings$pc_max_condition,
    n_pc = settings$n_pc,
    ridge = settings$ridge,
    eig_tol = settings$eig_tol,
    fixed_case_control = settings$fixed_case_control,
    H_B_centering = settings$H_B_centering,
    sigma_B_centering = settings$sigma_B_centering,
    gmm_start_method_reference = settings$gmm_start_method_reference,
    gmm_H_mode = settings$gmm_H_mode,
    gmm_solver = settings$gmm_solver,
    gmm_score_tol = settings$gmm_score_tol
  )
  cc_example_save(fit, setup, settings, generation, beta_tilde,
                  elapsed_seconds = proc.time()[["elapsed"]] - started)
  invisible(fit)
}

if (sys.nframe() == 0L) run_case_control_two_populations()
