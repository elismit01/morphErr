#' morphErr: Morphometric Analysis with Measurement Error
#'
#' Provides tools for analysing morphometric data from wildlife populations
#' when measurements are subject to error. Particularly suited for measurements
#' obtained from drone footage, this package includes functions for model fitting,
#' simulation, visualisation, and ratio analysis.
#'
#' @keywords internal
"_PACKAGE"

#' @importFrom graphics plot.new plot.window box axis title points lines par abline
#' @importFrom grDevices hcl.colors
#' @importFrom lmeInfo Fisher_info extract_varcomp
#' @importFrom Matrix bdiag
#' @importFrom msm deltamethod
#' @importFrom mvtnorm dmvnorm rmvnorm
#' @importFrom nlme lme lmeControl getData groupedData corSymm varIdent
#' @importFrom parallel makeCluster stopCluster clusterExport
#' @importFrom pbapply pblapply
#' @importFrom stats nlminb pnorm qnorm reformulate predict vcov
#' @importFrom utils combn setTxtProgressBar txtProgressBar
NULL
