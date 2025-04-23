#' morphErr: Morphometric Analysis with Measurement Error
#'
#' Provides tools for analysing morphometric data from wildlife populations
#' when measurements are subject to error. Particularly suited for measurements
#' obtained from drone footage, this package includes functions for model fitting,
#' simulation, visualisation, and ratio analysis.
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom stats pnorm qnorm predict vcov
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom Matrix bdiag
#' @importFrom mvtnorm dmvnorm rmvnorm
#' @importFrom nlme lme lmeControl getData groupedData corSymm varIdent
#' @importFrom lmeInfo Fisher_info extract_varcomp
#' @importFrom graphics plot.new plot.window box axis title points lines par abline
#' @importFrom grDevices hcl.colors
#' @importFrom parallel makeCluster stopCluster clusterExport
#' @importFrom pbapply pblapply
#' @importFrom stats nlminb reformulate
## usethis namespace: end
NULL
