test_that("the public API and bundled resources are available", {
  expect_true(is.function(zim_cohort))
  expect_true(is.function(zim_cohort_wkm))
  expect_true(is.function(zim_case_control))
  expect_true(is.function(zim_case_control_wkm))
  expect_true(is.function(zim_case_control_crossfit))
  expect_true(file.exists(zim_example_path("cohort_one_population.R")))
  expect_true(nzchar(system.file(
    "extdata", "cohort_one_population.rds", package = "ZIM"
  )))
})

test_that("cohort and case-control implementations load independently", {
  cohort <- ZIM:::.zim_impl_env("cohort")
  case_control <- ZIM:::.zim_impl_env("case_control")
  expect_true(is.function(cohort$zim_cohort))
  expect_true(is.function(case_control$zim_case_control))
  expect_false(identical(cohort, case_control))
})

test_that("bundled examples resolve installed source and data files", {
  cohort_path <- zim_example_path("cohort_one_population.R")
  cohort_environment <- new.env(parent = globalenv())
  sys.source(cohort_path, envir = cohort_environment)
  cohort_setup <- cohort_environment$example_setup(
    dirname(cohort_path), "cohort_one_population"
  )
  expect_true(file.exists(cohort_setup$data_file))
  expect_true(is.function(cohort_setup$core$zim_cohort))

  case_control_path <- zim_example_path("case_control_one_population.R")
  case_control_environment <- new.env(parent = globalenv())
  sys.source(case_control_path, envir = case_control_environment)
  case_control_setup <- case_control_environment$cc_example_setup(
    dirname(case_control_path), "case_control_one_population"
  )
  expect_true(file.exists(case_control_setup$data_file))
  expect_true(is.function(case_control_setup$core$zim_case_control))
})

test_that("bundled examples use the small no-covariate PRS configuration", {
  examples <- c(
    cohort_one_population = 2000L,
    cohort_two_populations = 2000L,
    case_control_one_population = 2000L,
    case_control_two_populations = 2000L
  )
  for (name in names(examples)) {
    inputs <- readRDS(system.file("extdata", paste0(name, ".rds"), package = "ZIM"))
    expect_identical(inputs$settings$p, 100L)
    expect_identical(inputs$settings$m, 100L)
    expect_identical(nrow(inputs$internal_data), unname(examples[[name]]))
    expect_false("Zcov" %in% names(inputs$internal_data))
    expect_false("Zcov" %in% names(inputs$settings$outcome_model))
    expect_length(inputs$beta_tilde, 100L)
    if (!is.null(inputs$reference_data)) {
      expect_false("Zcov" %in% names(inputs$reference_data))
    }
  }
})
