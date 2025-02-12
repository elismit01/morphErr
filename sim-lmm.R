## In this file, we simulate one data set and use the full range of
## inference tools.
funs.file <- "lmm-funs.R"
source(funs.file)

set.seed(1234)
## Some parameter values, roughly aligning with Setyawan et al (2022).
mus <- c(290, 130, 75)
sigmas <- c(45, 25, 15)
rhos <- c(0.85, 0.90, 0.95)
psis <- c(2.00, 1.50, 1.00)
phis <- c(0.40, 0.50, 0.60)
pars <- c(mus, sigmas, rhos, psis, phis)
par.list <- list(mus = mus, sigmas = sigmas, rhos = rhos, psis = psis, phis = phis)

## Simulating data.
n.animals <- 10
n.photos <- rep(5, n.animals)
m <- 3
data <- sim.measurements(n.animals, n.photos, m, pars)

## Making RTMB object.
#obj <- MakeADFun(lmm.likelihood, pars)
## Parameter bounds.
#lower <- c(rep(-Inf, 3), rep(0, 3), rep(0, 3), rep(0, 3), rep(0, 3))
#upper <- c(rep(Inf, 3), rep(Inf, 3), rep(1, 3), rep(Inf, 3), rep(1, 3))
## Fitting model with TMB.
#fit.tmb <- nlminb(pars, obj$fn, obj$gr, control = list(iter.max = 10000),
#                  lower = lower, upper = upper)
#fit.tmb.rep <- sdreport(obj)
#summary(fit.tmb.rep)

## Same model, but with lme().
fit <- fit.morph(data)

## Estimates and standard errors.
summary(fit)
## Estimates and variance-covariance matrix.
vcov(fit)
## Estimated animal sizes.
coef(fit)

## Plotting data.
plot(fit, dims = c(1, 2))

## Coefficients to predict dimension 2 from dimension 1.
betas.2.1 <- summary(fit, type = "betas", y.dim = 2, x.dim = 1)
betas.2.1
## You can also get a PCA line.
betas.2.1.pca <- summary(fit, type = "betas-pca", y.dim = 2, x.dim = 1)
betas.2.1.pca
## Coefficients to predict dimension 1 from dimension 2.
betas.1.2 <- summary(fit, type = "betas", y.dim = 1, x.dim = 2)
betas.1.2
## Just checking: should get the same PCA line either way.
betas.1.2.pca <- summary(fit, type = "betas-pca", y.dim = 1, x.dim = 2)
c(-betas.1.2.pca[1, 1]/betas.1.2.pca[2, 1], 1/betas.1.2.pca[2, 1])
betas.2.1.pca[, 1]
## You can get coefficients to predict one dimension based on any
## number of other dimensions. Here we get coefficients to predict
## dimension 3 from dimensions 1 and 2.
summary(fit, type = "betas", y.dim = 3, x.dim = c(1, 2))

## You can plot a fitted line with confidence intervals like this.
plot(fit, dims = c(1, 2), line.type = "lm")
## You can do the same thing for a PCA line.
plot(fit, dims = c(1, 2), line.type = "pca")

## You can plot multiple lines on the same plot.
plot(fit, dims = c(1, 2), confints = FALSE, line.type = "lm")
plot(fit, dims = c(1, 2), line.type = "pca", add = TRUE, lty = "dotted")
## Here we're also adding the regression line of y vs x, and using
## reverse.axes to get the right line for a plot with the axes this
## way around.
plot(fit, dims = c(2, 1), line.type = "lm", add = TRUE, reverse.axes = TRUE)

## Estimates of mean ratios.
summary(fit, type = "ratios")
## Plotting PDFs of ratios.
plot(fit, dims = c(3, 1), type = "ratio-pdf")

## Testing for isometric growth using a normal approximation for the
## test statistic.
summary(fit, type = "isometric-pca")
## A bootstrap instead.
summary(fit, type = "isometric-pca-boot")

## Plotting how a ratio changes depending on one of the dimensions. We
## can use fitted lines of either the lm or pca type.
plot(fit, dims = c(1, 2), type = "ratio", line.type = "pca")
plot(fit, dims = c(1, 2), type = "ratio", line.type = "lm")

## Testing out predictions based on true and/or observed
## measurements. Here are some observed measurements subject to error.
obs <- matrix(data$measurement[data$animal.id == 1], ncol = 3, byrow = TRUE)
## Predicting all three dimensions using these measureents. Almost
## equal to the averages, which is what we'd expect.
predict.from.obs(fit, true = c(NA, NA, NA), obs = obs)
## What if we only had one photo?
predict.from.obs(fit, true = c(NA, NA, NA), obs = obs[1, ])
## Note that these predictions are similar to the average dimension
## sizes provided in obs, but are (very) slightly shrinking downwards
## towards the population means.
summary(fit)[1:3, ]

## Let's say we somehow know that the true value for the second
## dimension is 140. Here are updated predictions.
predict.from.obs(fit, true = c(NA, 140, NA), obs = obs)
## Let's say that somehow we didn't get measurements for the first
## dimension from any of the photographs. We can still make
## predictions.
obs[, 1] <- NA
predict.from.obs(fit, true = c(NA, NA, NA), obs = obs)
predict.from.obs(fit, true = c(NA, 140, NA), obs = obs)

## We can make predictions even if we don't have observed
## measurements, but then they just default to the population means.
predict.from.obs(fit, true = c(NA, NA, NA), obs = NULL)
summary(fit)[1:3, ]

## Here we're predicting the first and third dimensions based on a
## true measurement of the second, without any additional photograph
## measurements.
predict.from.obs(fit, true = c(NA, 140, NA), obs = NULL)

## Interestingly, we can do this in a different way using
## predict.lme.morph() and we get a prediction that is equal to the
## linear-model interpretation. I expect the two are mathematically
## equivalent somehow. The prediction using an RMA interpretation is a
## bit different.
predict(fit, y.dim = 1, newdata = data.frame(dim2 = 140), type = "lm")
predict(fit, y.dim = 1, newdata = data.frame(dim2 = 140), type = "pca")
