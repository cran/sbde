\name{coef.sbde}
\Rdversion{0.0-6}
\alias{coef.sbde}
\title{Coefficient Extraction from sbde Model Fit}
\description{Post process MCMC output from \code{\link{sbde}} to create summaries of parameter and quantile estimates.}
\usage{ 
\method{coef}{sbde}(object, burn.perc = 0.5, nmc = 200,
    prob = c(.001,.01,.1,1:99,99.9,99.99,99.999)/100, ...)
} 
\arguments{
\item{object}{a fitted model of the class \code{sbde}.}
\item{burn.perc}{a positive fraction indicating what fraction of the saved draws are to be discarded as burn-in}
\item{nmc}{integer giving the number of samples, post burn-in, to be used in Monte Carlo averaging}
\item{prob}{a numeric vector of probabiities at which quantiles are to be estimated.}
\item{...}{not currently implemented}
}
\value{
Extracts posterior summary of model parameters, as well as estimated quantiles. A list is returned invisibly with the following fields.

\item{psamp}{a matrix with 3 columns and \code{nmc} rows storing the posterior draws of the parameters of base distribution used in transformation}
\item{parametric}{a matrix with posterior median, 2.5th and 97.5th percentiles of the parameters of the base distribution.}
\item{prob}{numeric vector of probabilities at which quantiles have been estimated. Could differ slightly from the input vector \code{prob}, by removing repetitions, as well as values that are not strictly between zero and one.}
\item{qsamp}{a matrix with \code{nmc} columns giving the posterior draws of the quantile values at levels given by \code{prob}.}
\item{qest}{a summary of \code{qsamp} given by the posterior median and 95 precent credible interval end points.}
\item{ss}{a vector of integers giving the indices of the mcmc samples that were used in posterior summary calculations.}
}

\seealso{\code{\link{sbde}}, \code{\link{summary.sbde}} and \code{\link{predict.sbde}} for model fitting under sbde.}
\examples{
y <- abs(rt(n=1000, df=4))
fit <- sbde(y, blocking="all", fbase="gpd", verbose=FALSE)
coef(fit)
}
\keyword{programming}
