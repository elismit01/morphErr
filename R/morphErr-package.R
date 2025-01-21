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
#' @importFrom stats model.matrix pnorm dnorm qnorm cov2cor predict vcov
#' @importFrom utils head tail setTxtProgressBar txtProgressBar
#' @importFrom Matrix bdiag
#' @importFrom mvtnorm dmvnorm rmvnorm
#' @importFrom RTMB MakeADFun sdreport
#' @importFrom nlme lme lmeControl fixef getData groupedData corSymm varIdent
#' @importFrom lmeInfo Fisher_info extract_varcomp
#' @importFrom graphics plot.new plot.window box axis title points lines par abline
#' @importFrom grDevices hcl.colors
#' @importFrom parallel makeCluster stopCluster clusterExport clusterEvalQ
#' @importFrom pbapply pblapply
#' @importFrom numDeriv hessian
#' @importFrom stats nlminb reformulate
## usethis namespace: end
NULL
