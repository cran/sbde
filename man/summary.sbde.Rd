\name{summary.sbde}
\Rdversion{0.0-6}
\alias{summary.sbde}
\title{Summary Method for Semiparametric Density Estimation}
\description{Summarize model fit for \code{\link{sbde}}}
\usage{ 
\method{summary}{sbde}(object, ntrace = 1000, burn.perc = 0.5, plot.dev = TRUE,
        more.details = FALSE, ...)
} 
\arguments{
    \item{object}{a fitted model of the class 'sbde'.}
    \item{ntrace}{number of draws to be included in trace plots}
    \item{burn.perc}{fraction of MCMC draws to be discarded as burn-in.}
    \item{plot.dev}{logical indicator of whether to show trace plot of deviance}
    \item{more.details}{logical indicating whether other details from MCMC are to be plotted}
    \item{...}{a limited number of plotting controls that are passed onto the deviance plot}
}
\value{
Displays the trace of the deviance statistic. More details include trace plots of of the proximity parameter of each GP, a plot of Geweke p-values for (from \code{\link{geweke.diag}}) convergence of each model parameter and an image plot of parameter correlation. 

The following quantities are returned invisibly.
\item{deviance}{vector deviance statistic of the samples parameter draws} 
\item{pg}{a matrix with \code{nsamp} number of columns. Each column gives the conditional posterior weights on the lambda grid values for the corresponding GP function.}       
\item{prox}{posterior draws of proximity parameter.}
\item{ll}{a matrix of \code{n*nsamp} containing observation level log-likelihood contributions. Used to calculate \var{waic}, and could be used for other AIC calculations.}       
\item{waic}{Two versions of Watanabe AIC from Gelman, Hwang and Vehtari (2014).}
}

\references{Gelman, A., Hwang, J., and Vehtari, A. (2014). Understanding predictive information criterion for Bayesian models. \emph{Stat Comput}, 24, 997-1016.}
\seealso{\code{\link{sbde}} and \code{\link{coef.sbde}}.}
\examples{
y <- abs(rt(n=1000, df=4))
fit <- sbde(y, blocking="all", fbase="gpd", verbose=FALSE)
sm <- summary(fit, more=TRUE)
print(sm$waic)
}
\keyword{programming}
