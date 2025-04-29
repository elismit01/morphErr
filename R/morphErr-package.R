#' morphErr: Measurement error models for morphometric data
#'
#' The `morphErr` package provides functions to analyse morphometric
#' data using the method proposed by Stevenson, Smit, and Setyawan (in
#' submission). View the manuscript \href{../doc/model.pdf}{here}.
#'
#' Two common ways of analysing morphometric data are to fit a linear
#' regression models or find the reduced major axis (i.e., the
#' principal component axis). However, these methods are known to
#' perform poorly when observations are subject to non-negligible
#' measurement error. The model implemented in this package provides
#' the same inference you receive from these alternatives, but
#' explicitly accounts for measurement error.
#'
#' @section Data requirements:
#' 
#' The `morphErr` package might be right for you if your data were
#' collected on photogrammetry surveys (e.g., using drones) and you
#' have multiple photographs of some individuals.
#'
#' The package might be particularly useful relative to alternative
#' options if
#'
#' * you take a lot of different measurements from each photograph
#' (i.e., you observe a lot of "dimensions"),
#'
#' * observations are subject to non-negligible measurement error,
#'
#' * the measurement errors are correlated (e.g., a photograph with
#' positive measurement error for one dimension tends to have positive
#' measurement errors for other dimensions),
#'
#' * you don't necessarily measure every dimension in every photograph,
#'
#' * you'd like to estimate relationships amongst different subsets of
#' the dimensions, or
#'
#' * you'd like to predict different dimensions for different
#' individuals using different subsets of the other dimensions.
#'
#' @section Package overview:
#'
#' The key functions in `morphErr` are as follows:
#'
#' * [`plotmorph()`] to plot morphometric data from a photogrammetry
#' survey.
#'
#' * [`fit.morph()`] to fit the model described by
#' \href{../doc/model.pdf}{Stevenson, Smit, and Setyawan (in
#' submission)}.
#'
#' * [`summary.lme.morph()`] with `type = "pars"` for parameter
#' estimates and standard errors.
#'
#' * [`summary.lme.morph()`] with `type = "betas"` for estimated
#' coefficients of linear relationships to predict one dimension from
#' any subset of the other dimensions.
#'
#' * [`summary.lme.morph()`] with `type = "betas-pca"` for estimated
#' coefficients of the reduced major axis (or pricipal component axis)
#' summarising the relationship between two dimensions.
#'
#' * [`summary.lme.morph()`] with `type = "isometric-pca"` or `type =
#' "isometric-pca-boot"` to test for isometric growth between all
#' pairs of dimensions.
#'
#' * [`plot.lme.morph()`] to plot estimated relationships between
#' dimensions.
#'
#' * [`sim.measurements()`] to simulate morphometric data.
#'
#' * [`sim.morph()`] to simulate multiple data sets and fit a model to
#' each one.
#'
#' * [`extract.sim.morph()`] to extract estimates from the models
#' fitted using [`sim.morph()`].
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

#' Manta Ray Morphometric Data
#'
#' Data collected on a drone photogrammetry survey of the reef manta
#' ray \emph{Mobula alfredi} in Raja Ampat, Indonesia.
#'
#' @format ## `manta`
#' A data frame with four columns:
#' \describe{
#' 
#'   \item{\code{animal.id}}{An individual identification number. Rows
#'                    with the same \code{animal.id} correspond to
#'                    measurements of the same individual manta ray.}
#' 
#'   \item{\code{photo.id}}{A photo identification number. Rows with
#'                   the same \code{photo.id} correspond to
#'                   measurements taken from the same image.}
#' 
#'   \item{\code{dim}}{An integer indicating the dimension the
#'              measurement is for. Here, \code{1} is for disc width,
#'              \code{2} is for disc length, and \code{3} is for
#'              cranial width.}
#' 
#'   \item{\code{measurement}}{The observed measurement value in
#'                             centimetres.}
#' 
#' }
"manta"
