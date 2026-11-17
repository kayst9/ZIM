.zim_impl_cache <- new.env(parent = emptyenv())
.zim_namespace <- environment()

.zim_impl_env <- function(kind) {
  if (!kind %in% c("cohort", "case_control")) {
    stop("Unknown ZIM implementation: ", kind, call. = FALSE)
  }
  if (!exists(kind, envir = .zim_impl_cache, inherits = FALSE)) {
    script <- switch(
      kind,
      cohort = "zim_cohort.R",
      case_control = "zim_case_control.R"
    )
    source_path <- system.file("scripts", script, package = "ZIM")
    if (!nzchar(source_path)) {
      stop("The bundled ZIM implementation is unavailable.", call. = FALSE)
    }
    implementation <- new.env(parent = .zim_namespace)
    sys.source(source_path, envir = implementation)
    assign(kind, implementation, envir = .zim_impl_cache)
  }
  get(kind, envir = .zim_impl_cache, inherits = FALSE)
}
