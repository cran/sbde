\name{sbde}
\Rdversion{0.0-6}
\alias{sbde}
\alias{update.sbde}
\title{Bayesian Semiparametric Density Estimation}
\description{Provides a semiparametric estimation of the density function of independent univariate data.
}
\usage{
sbde(y, nsamp = 1e3, thin = 10, cens = rep(0,length(y)), 
     wt = rep(1,length(y)), incr = list(knot=0.2, grid=0.01),
     par = c("Hill-kde", "pmean", "rand")[1], tail.warp = c(0,0), 
     hyper = list(sig = c(.1,.1), lam = c(6,4), kap = c(1.5,1.5,1)), 
     prox.range = c(.2,.95), acpt.target = 0.15, ref.size = 3, 
     blocking = c("all", "gp", "loc+scale+tail"), temp = 1, expo = 2, 
     blocks.mu, blocks.S, fix.nu = FALSE, 
     fbase = c("t", "t+", "gpd", "gpd-rescaled", "unif"), 
     spacing=list(knot="regular", grid="regular"),
     verbose = TRUE)

\method{update}{sbde}(object, nadd, append = TRUE, ...)
}
\arguments{
\item{y}{numeric vector of response data.}
\item{nsamp}{number of posterior samples to be saved; defaults to 1000.}
\item{thin}{thinning rate for the Markov chain sampler -- one posterior sample is saved per \code{thin} iterations. Defaults to 10. The Markov chain sampler runs for a total of \code{nsamp * thin} many iterations.}
\item{cens}{censoring status of response. Must be a vector of 0s and 1s of length same as length(y), with 1 indicating right censoring, and, 0 indicating no censoring. Defaults to all zeros.}
\item{wt}{weights attached to the observation units, expected to be non-negative numbers, and defaults to a vector of ones.}
\item{incr}{a list with two named elements, 'knot' and 'grid', giving the increment sizes for the knots in the predictive process approximation and the grid to be used for logistic Gaussian process likelihood evaluation. Defaults to 0.2 and 0.01 respectively}
\item{par}{either a numeric vector giving the parameter initialization or a character string indicating how the parameter should be initialized. If input numeric vector length is smaller than required parameter count, then supplied values are appended with zeros to create a full initialization. If input equals "pmean" then the mcmc is initialized at the prior center given by a vector of zeros, or if it equals "rand" then intialization is done by drawing randomly from the prior, or if it equals "Hill-kde" then the Hill estimate is used to estimate the shape parameter, the location and scale parameters are set based on data median and 95th percentile, and the initialization of the Gaussian process is done based on a kernel density estimate of the transformed data.}
\item{tail.warp}{a non-negative 2-vector giving the degrees of tail warping to be done at each tail. Larger values will allow more variation of the non-parametric density at the corresponding tail. }
\item{hyper}{hyperparameters of the prior distribution. Must be a list with one or both of the following two fields: \code{lam}: a two vector giving the parameters of the beta distribution on proximity = \eqn{\exp(-0.01* \lambda^2)}{exp(-0.01 * lambda^2)}, and \code{kap}:  a vector to be coerced into a \code{3 * nkap} matrix, with \code{nkap} being the number of components in the mixture of gamma prior on \code{kappa}, and each column of the matrix gives the shape, rate and mixing weight of a component.}
\item{prox.range}{for specifying the range of length-scale parameter of the Gaussian process prior.}
\item{acpt.target}{target acceptance rate of the adaptive Metropolis sampler; defaults to 0.15}
\item{ref.size}{adaptation rate of the adaptive Metropolis sampler. The proposal density is updated once every \code{ref.size} iterations. Could be a single number or a vector of length same as the number of blocks.}
\item{blocking}{type of blocking to be applied represented by a character vector with elements comprising of the strings: "gp", "loc", "scale", "tail" and their combinations separated by "+". Each of the basic string types will include the corresponding model parameters into the block. For example a valid input could be c("gp", "gp+loc+scale", "loc+scale+tail"), where the first block updates only the Gaussian process parameters, the second block jointly updates the GP parameters and the location and scale, and, the third block updates the location, scale and tail parameters. A combination of all four types can be represented as "all". }
\item{temp}{temperature of the log-likelihood function. The log-likelihood function is raised to the power of \code{temp}. Defaults to 1.}
\item{expo}{the exponent to be used in the covariance kernel of the Gaussian process priors. Defaults to 2, giving the standard squared-exponential covariance kernel.}
\item{blocks.mu}{initial block specific means in the form of  a list. If left unspecified then will be automatically generated as a list of vectors of zeros of appropriate lengths matching the corresponding block sizes.}
\item{blocks.S}{initial block specific covariance matrices in the form of a list. If left unspecified then will be automatically generated as a list of identity matrices of appropriate dimensions matching the corresponding block sizes.}
\item{fix.nu}{either the logical FALSE indicating that nu should be learned, or a positive real number giving the fixed value of nu, which is then excluded from MCMC updates}
\item{fbase}{either "t" (default) or "t+" (for half-t distributions on the positive real lines) or "gpd" (for generalized pareto distributions with location zero and parametrized by nu = 1 / shape) or "gpd-rescaled" (same as gpd, but scale parameter adjusted according to shape so that 90-th percentile matches that of gpd with shape=1/6 and scale=1) or  "unif" to indicate what base distribution is to be used.}
\item{spacing}{the type of spacing to be used for the predictive process knots and the likelihood evaluation grid. For either object, the default choice is "regular". Any other specification is taken to equal "irregular". A regular grid places points equally between 0 and 1 as given by the prespecified increment value. When the likelihood "grid" is chosen to be "irregular", the regular grid is appended with more points at both extremes by recursive bisection until 1/\code{n} or 1 - 1/\code{n} is reached. For predictive process knots, "irregular" applies only when \code{tail.warp} is different that \code{c(0,0)}, and more knots are appended at each extreme based on how much warping is done to it.}
\item{verbose}{logical indicating whether MCMC progress should be printed, defaults to TRUE}
\item{object}{a fitted model of the class 'qde'.}
\item{nadd}{number of additional MCMC samples.}
\item{append}{logical indicating whether new samples should be appended to old ones. If FALSE then old samples are discarded.}
\item{...}{no additional arguments are allowed}
}
\value{
\code{sbde(y, ...)} returns a `sbde' class object to be used by \code{\link{coef}}, \code{\link{summary}} and \code{\link{predict}}.

Returned object is a list containing the following variables.
\item{par}{latest draw of the parameter vector}
\item{y}{response vector}
\item{cens}{censoring status vector}
\item{wt}{vector of observation weights}
\item{hyper}{completed list of hyper-parameters}
\item{dim}{model dimension vector of the form c(n, length of tau grid, position of \eqn{\tau_0}{tau0} on the grid, nknots, length of lambda grid, nkap, total number of MCMC iterations, thin, nsamp)}
\item{gridmats}{details of covariance matrix factors etc, intended for internal use.}
\item{tau.g}{the tau grid}
\item{muV}{list of means for parameter blocks}
\item{SV}{list of covariance matrices for parameter blocks}
\item{blocks}{list of blocks}
\item{blocks.size}{vector of block lengths}
\item{dmcmcpar}{numeric vector containing details of adaptive MCMC runs, equals c(temp, decay rate of adaptation, vector of target acceptance rates for the blocks, vector of increment scales used in adaptation). Intended strictly for internal use.}
\item{imcmcpar}{numeric vector containing details of adaptive MCMC runs, equals c(number of parameter blocks, ref.size, indicator on whether details are to be printed during MCMC progress, rate of details printing, a vector of counters needed for printing). Intended strictly for internal use.}
\item{parsamp}{a long vector containing the parameter draws. Could be coerced into a matrix of dim \code{npar * nsamp}. Intended primarily for use by \code{\link{summary}} and \code{\link{coef}}.}
\item{acptsamp}{a long vector containing rates of acceptance statistics for parameter blocks. Could be coerced into a matrix of dim \code{nblocks * nsamp}. Not very informative, because thinning times and adaptation times may not be exactly synced.}
\item{lpsamp}{vector of log posterior values for the saved MCMC draws.}
\item{other.controls}{a vector of two integers, with the first storing the choice of the fbase, and the second storing the choice of the gridtype.}
\item{prox}{vector of proximity (exp(-0.01*lambda^2)) grid values}
\item{runtime}{run time of the MCMC}
\item{base.bundle}{a list of densty, distribution, quantile etc functions associated with the base distribution.}
}

\details{
For positive valued data, it is recommended to use fbase as "gpd", which yields much faster computation than the choice of "t+". The difference is entirely due to difference in machine time needed to compute the CDF of the generalized Pareto versus that of the Student-t. 
}
\references{Tokdar, S.T., Jiang, S. and Cunningham, E.L. (2022). Heavy-tailed density estimation. \emph{Journal of the American Statistical Association}, (just-accepted) <https://doi.org/10.1080/01621459.2022.2104727>.}
\seealso{\code{\link{summary.sbde}}, \code{\link{coef.sbde}} and \code{\link{predict.sbde}}.}
\examples{
y <- abs(rt(n=1000, df=4))
fit <- sbde(y, blocking="all", fbase="gpd", verbose=FALSE)
coef(fit)
}
\keyword{programming}

