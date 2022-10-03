\name{predict.sbde}
\Rdversion{0.0-6}
\alias{predict.sbde}
\title{Posterior predictive Summary for Semiparametric Density Estimation}
\description{Extract posterior predictive density estimate for \code{\link{sbde}}}
\usage{ 
\method{predict}{sbde}(object, burn.perc = 0.5, nmc = 200, yRange = range(object$y), yLength = 401, ...)
} 
\arguments{
    \item{object}{a fitted model of the class 'sbde'.}
    \item{burn.perc}{a positive fraction indicating what fraction of the saved draws are to be discarded as burn-in}
    \item{nmc}{integer giving the number of samples, post burn-in, to be used in Monte Carlo averaging}
    \item{yRange}{Range of values over which posterior predictive density is to be evaluated.}
    \item{yLength}{Number of grid points spanning yRange for posterior predictive density evaluation.}
    \item{...}{no additional parameters are used.}
}
\value{
Returns a list with three items:
\item{y}{vector giving the grid over which the posterior predictive density is evaluated.}
\item{fsamp}{a matrix with \code{yLength} many rows and \code{nmc} many columns. Each column corresponds to a draw of the response density from the posterior predictive.}
\item{fest}{summary of the posterior predictive density given by point-wise median, 2.5th and 97.5th percentiles.}
}

\seealso{\code{\link{sbde}}, \code{\link{coef.sbde}} and \code{\link{summary.sbde}}.}
\examples{
y <- abs(rt(n=1000, df=4))
fit <- sbde(y, blocking="all", fbase="gpd", verbose=FALSE)
pp <- predict(fit)
hist(y, 50, freq=FALSE)
with(pp, for(j in 1:3) lines(y, fest[,j], lty=1+(j!=2)))
}
\keyword{programming}
