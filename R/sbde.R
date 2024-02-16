sbde <- function(y, nsamp = 1e3, thin = 10, cens = rep(0,length(y)), 
                 wt = rep(1,length(y)), incr = list(knot=0.2, grid=0.01),
                 par = c("Hill-kde", "pmean", "rand")[1], tail.warp = c(0,0), 
                 hyper = list(sig = c(.1,.1), lam = c(6,4), kap = c(1.5,1.5,1)), 
                 prox.range = c(.2,.95), acpt.target = 0.15, ref.size = 3, 
                 blocking = c("all", "gp", "loc+scale+tail"), temp = 1, expo = 2, 
                 blocks.mu, blocks.S, fix.nu = FALSE, 
                 fbase = c("t", "t+", "gpd", "gpd-rescaled", "unif"), 
                 spacing=list(knot="regular", grid="regular"), 
                 verbose = TRUE){
  
  if(length(tail.warp)==1) tail.warp <- rep(tail.warp,2)
  fbase.choice <- match(fbase[1], names(base.functions))
  if(is.na(fbase.choice)) stop(paste0("Only", names(base.functions), "are allowed for the choice of fbase"))
  base.bundle <- base.functions[[fbase.choice]]
  fix.nu <- base.bundle$fix.nu
  pos.half <- (base.bundle$support == "half-pos")
  if(pos.half) if(any(y < 0)) stop("Negative data not allowed when using a base distribution with positive support")
  
  incr.g <- incr$grid; if(is.null(incr.g)) incr.g <- 0.01
  incr.k <- incr$knot; if(is.null(incr.k)) incr.k <- 0.2
  grid.type <- spacing$grid; if(is.null(grid.type)) grid.type <- "regular"; if(!(grid.type=="regular")) grid.type <- "irregular"
  knot.type <- spacing$knot; if(is.null(knot.type)) knot.type <- "regular"; if(!(knot.type=="regular")) knot.type <- "irregular"
  
  n <- length(y)
  incr.g <- 1/ceiling(1/incr.g)
  tau.g <- seq(0, 1, incr.g)
  grid.choice <- (grid.type == "irregular")
  if(grid.choice){
    M <- ceiling(log(n * incr.g,2))
    tail1 <- incr.g * 0.5^(M:1)
    tail2 <- 1 - incr.g * 0.5^(1:M)
    tau.g <- sort(c(tail1, tau.g, tail2))
  }
  L <- length(tau.g)
  mid <- which.min(abs(tau.g - 0.5))
  
  if(knot.type == "regular"){
    Qfn <- Qfn.reg
    tau.k <- seq(0,1,incr.k)
    nknots <- length(tau.k)
  } else {
    Qfn <- Qfn.irreg
    end.pt <- Qfn(c(0,1), tail.warp)
    incr.k <- diff(range(end.pt))/ceiling(diff(range(end.pt))/incr.k)
    tau.kb <- seq(end.pt[1], end.pt[2], incr.k)
    tau.k <- Qfn.inv(tau.kb, tail.warp)
    nknots <- length(tau.k)
  }
  a.sig <- hyper$sig; if(is.null(a.sig)) a.sig <- c(.1, .1)
  a.lam <- hyper$lam; if(is.null(a.lam)) a.lam <- c(6, 4)
  a.kap <- hyper$kap; if(is.null(a.kap)) a.kap <- c(1.5,1.5,1); a.kap <- matrix(a.kap, nrow = 3); nkap <- ncol(a.kap); a.kap[3,] <- log(a.kap[3,])
  hyper.reduced <- c(a.sig, c(a.kap))
  
  prox.grid <- proxFn(max(prox.range), min(prox.range), 0.5, tau.k, tail.warp, expo, Qfn)
  #ell <- 0.16 / c(0.1, 0.5, 1:5)
  #prox.grid <- exp(-0.5 * 0.01 / ell^2)
  ngrid <- length(prox.grid)
  prior.grid <- -diff(pbeta(c(1, (prox.grid[-1] + prox.grid[-ngrid])/2, 0), a.lam[1], a.lam[2]))
  #prior.grid <- rep(1/ngrid, ngrid)
  lamsq.grid <- lamFn(prox.grid)^2
  lp.grid <- log(prior.grid)
  
  d.kg <- outer(tau.k, tau.g, dfn, expo=expo, tail.warp=tail.warp, Qfn=Qfn)
  d.kk <- outer(tau.k, tau.k, dfn, expo=expo, tail.warp=tail.warp, Qfn=Qfn)
  gridmats <- matrix(NA, nknots*(L + nknots)+2, ngrid)
  K0 <- 0
  t1 <- Sys.time()
  for(i in 1:ngrid){
    K.grid <- exp(-lamsq.grid[i] * d.kg); K.knot <- exp(-lamsq.grid[i] * d.kk);	diag(K.knot) <- 1 + 1e-10	
    R.knot <- chol(K.knot); A.knot <- solve(K.knot, K.grid)
    gridmats[,i] <- c(c(A.knot), c(R.knot), sum(log(diag(R.knot))), lp.grid[i])		
    K0 <- K0 + prior.grid[i] * K.knot
  }
  t2 <- Sys.time()
  
  niter <- nsamp * thin
  dimpars <- c(n, L, mid - 1, nknots, ngrid, ncol(a.kap), niter, thin, nsamp)
  
  if(is.character(par[1])){
    if(par[1] == "rand") {
      par <- rnorm(nknots+3)
    } else if(par[1] == "Hill-kde"){
      Q0 <- base.bundle$Q0
      F0 <- base.bundle$F0
      loc <- 0
      if(!pos.half) loc <- median(y)
      y.descend <- sort(abs(y-loc), decreasing=TRUE)
      threshold <- extremefit::hill.adapt(y-loc)$Xadapt
      tail.p <- min(max(mean(y-loc > threshold), 0.99), 5/n)
      tail.n <- ceiling(n * tail.p)
      nu.est <- pmax(1, 1/(mean(log(y.descend[1:(tail.n-1)])) - log(y.descend[tail.n])))
      #cat("nu.est =", round(nu.est, 2), "\n")
      sig.est <- pmax(1e-10, quantile(y-loc, prob=1-tail.p) / Q0(1-tail.p, nu=nu.est))
      ytrans <- F0((y-loc)/sig.est, nu=nu.est)
      log.h <- pmax(-10, log(density(ytrans, bw=2/(nknots-1), n=nknots, from=0, to=1)$y))
      log.h.adj <- log.h - median(log.h)
      w <- pmax(pmin(log.h.adj, 5), -5)
      par <- as.numeric(c(w, loc, sigFn.inv(sig.est), nuFn.inv(nu.est)))
    } else par <- rep(0, nknots+3)
  } 
  if(length(par) < nknots+3) par <- c(par, rep(0, nknots+3 - length(par)))
  if(pos.half) par[nknots+1] <- 0
  if(fix.nu) par[nknots+3] <- nuFn.inv(fix.nu)
  
  npar <- nknots+3
  par.position <- list(gp = 1:nknots, loc = nknots + 1, scale = nknots + 2, tail = nknots + 3, all = 1:npar)
  form.block <- function(block.string){
    comps <- strsplit(block.string, "\\+")[[1]]
    block <- rep(FALSE, npar)
    for(pn in comps) block[par.position[[pn]]] <- TRUE
    return(block)
  }
  blocks <- lapply(blocking, form.block)
  nblocks <- length(blocks)
  if(fix.nu) for(j in 1:nblocks) blocks[[j]][nknots+3] <- FALSE
  if(pos.half) for(j in 1:nblocks) blocks[[j]][nknots+1] <- FALSE
  non.empty.blocks <- which(sapply(blocks, any))
  blocks <- blocks[non.empty.blocks]
  nblocks <- length(blocks)
  if(nblocks == 0) {
    blocks <- lapply(1:npar, function(i) return(1:npar == i))
    nblocks <- length(blocks)
    if(fix.nu) for(j in 1:nblocks) blocks[[j]][nknots+3] <- FALSE
    if(pos.half) for(j in 1:nblocks) blocks[[j]][nknots+1] <- FALSE
    non.empty.blocks <- which(sapply(blocks, any))
    blocks <- blocks[non.empty.blocks]
    nblocks <- length(blocks)
  }
  
  blocks.ix <- c(unlist(lapply(blocks, which))) - 1
  blocks.size <- sapply(blocks, sum)
  if(missing(blocks.mu)) blocks.mu <- rep(0, sum(blocks.size))
  if(missing(blocks.S)){
    blocks.S <- lapply(blocks.size, function(q) diag(1e-2, q))
    ix.GP <- c(rep(TRUE, nknots), rep(FALSE, 3))
    for(i in 1:nblocks){
      ovrlp.GP <- which(ix.GP & blocks[[i]])
      novrlp <- length(ovrlp.GP)
      if(novrlp > 0) blocks.S[[i]][1:novrlp,1:novrlp] <- K0[ovrlp.GP,ovrlp.GP]
    }
    blocks.S <- unlist(blocks.S)
  }
  
  imcmc.par <- c(nblocks, ref.size, verbose, max(10, niter/1e4), rep(0, nblocks))
  dmcmc.par <- c(temp, 0.999, rep(acpt.target, nblocks), 2.38 / sqrt(blocks.size))
  
  
  tm.c <- system.time(oo <- .C("SBDE", par = as.double(par), y = as.double(y), cens = as.integer(cens), wt = as.double(wt), hyper = as.double(hyper.reduced), dim = as.integer(dimpars), gridmats = as.double(gridmats), tau.g = as.double(tau.g), muV = as.double(blocks.mu), SV = as.double(blocks.S), blocks = as.integer(blocks.ix), blocks.size = as.integer(blocks.size), dmcmcpar = as.double(dmcmc.par), imcmcpar = as.integer(imcmc.par), parsamp = double(nsamp * length(par)), acptsamp = double(nsamp * nblocks), lpsamp = double(nsamp), other.controls = as.integer(c(fbase.choice,grid.choice))))
  if(verbose) cat("elapsed time:", round(tm.c[3]), "seconds\n")
  
  oo$y <- y; oo$gridmats <- gridmats; oo$prox <- prox.grid; oo$runtime <- tm.c[3]
  oo$base.bundle <- base.bundle
  class(oo) <- "sbde"
  return(oo)
}


update.sbde <- function(object, nadd, append = TRUE, ...){
  niter <- object$dim[7]; thin <- object$dim[8]; nsamp <- object$dim[9]
  if(missing(nadd)) nadd <- nsamp
  par <- object$par; npar <- length(par)
  dimpars <- object$dim
  dimpars[7] <- nadd * thin
  dimpars[9] <- nadd
  nblocks <- object$imcmcpar[1]
  object$imcmcpar[4] <- max(10, nadd * thin/1e4)
  
  tm.c <- system.time(oo <- .C("SBDE", par = as.double(par), y = as.double(object$y), cens = as.integer(object$cens), wt = as.double(object$wt),
                               hyper = as.double(object$hyper), dim = as.integer(dimpars), gridmats = as.double(object$gridmats),
                               tau.g = as.double(object$tau.g), muV = as.double(object$muV), SV = as.double(object$SV), blocks = as.integer(object$blocks),
                               blocks.size = as.integer(object$blocks.size), dmcmcpar = as.double(object$dmcmcpar),
                               imcmcpar = as.integer(object$imcmcpar), parsamp = double(nadd * npar),
                               acptsamp = double(nadd * nblocks), lpsamp = double(nadd), other.controls = as.integer(object$other.controls)))
  cat("elapsed time:", round(tm.c[3]), "seconds\n")
  
  oo$y <- object$y; oo$gridmats <- object$gridmats; oo$prox <- object$prox; oo$runtime <- object$runtime+tm.c[3]
  oo$base.bundle <- object$base.bundle
  if(append){
    oo$dim[7] <- niter + nadd * thin
    oo$dim[9] <- nsamp + nadd
    oo$parsamp <- c(object$parsamp, oo$parsamp)
    oo$acptsamp <- c(object$acptsamp, oo$acptsamp)
    oo$lpsamp <- c(object$lpsamp, oo$lpsamp)
  }
  class(oo) <- "sbde"
  return(oo)
}

coef.sbde <- function(object, burn.perc = 0.5, nmc = 200,
                      prob = c(.001,.01,.1,1:99,99.9,99.99,99.999)/100, ...){
  niter <- object$dim[7]
  nsamp <- object$dim[9]
  pars <- matrix(object$parsamp, ncol = nsamp)
  ss <- unique(round(nsamp * seq(burn.perc, 1, len = nmc + 1)[-1]))
  
  dimpars <- object$dim
  dimpars[7] <- length(ss)
  
  n <- dimpars[1]; nknots <- dimpars[4];
  
  ## parametric coefficients
  psamp <- cbind(gam0 = pars[nknots + 1,ss], sigma = sigFn(pars[nknots + 2,ss]), nu = nuFn(pars[nknots + 3,ss]))
  gamsignu <- t(apply(psamp, 2, quantile, pr = c(0.5, 0.025, 0.975)))
  dimnames(gamsignu)[[2]] <- c("Estimate", "Lo95%", "Up95%")
  
  ## semiparametric quantiles
  if(any(prob <= 0.0)){
    prob <- prob[prob > 0.0]
    warning("Removing prob values less than or equal to zero")
  }
  if(any(prob >= 1.0)){
    prob <- prob[prob < 1.0]
    warning("Removing prob values greater than or equal to one")
  }
  prob0 <- c(0.01, 0.1, 0.5, 0.9, 0.99)
  plen0 <- length(prob0)
  prob <- c(prob0, prob)
  #print(prob)
  plen <- length(prob)
  keep.prob <- (1:plen > plen0)
  
  ord <- order(prob)
  pSorted <- prob[ord]
  pUnique <- unique(pSorted)
  #print(pUnique)
  
  nReps <- diff(sapply(c(0, pUnique), function(a) sum(pSorted <= a)))
  pLength <- length(pUnique)
  
  
  ## q[ord] <- rep(qUnique, nReps)
  n <- dimpars[1]; p <- 0; ngrid <- dimpars[5]
  dimpars[1] <- pLength
  
  qnts <- .C("QUANT", pars = as.double(pars[,ss]), uGrid = as.double(pUnique), hyper = as.double(object$hyper), dim = as.integer(dimpars), gridmats = as.double(object$gridmats), tau.g = as.double(object$tau.g), qvalsamp = double(length(ss)*pLength), other.controls = as.integer(object$other.controls))
  qvals.unique <- matrix(qnts$qvalsamp, ncol = length(ss))
  #print(qvals.unique)
  qvals.sorted <- apply(qvals.unique, 2, function(z) return(rep(z, nReps)))
  #print(qvals.sorted)
  
  qvals <- matrix(NA, plen, length(ss))
  qvals[ord,] <- qvals.sorted
  
  out <- list(parametric = gamsignu, psamp = psamp, prob = prob[keep.prob], qsamp = qvals[keep.prob,],
              qest = t(apply(qvals[keep.prob,,drop=FALSE], 1, quantile, pr = c(.025, .5, .975))),
              ss=ss)
  
  qval0 <- t(apply(qvals[!keep.prob,,drop=FALSE], 1, quantile, pr = c(.025, .5, .975)))
  dimnames(qval0)[[1]] <- sprintf("Q%02d", 100*prob0)
  print(rbind(gamsignu, qval0))
  invisible(out)
}

summary.sbde <- function(object, ntrace = 1000, burn.perc = 0.5, plot.dev = TRUE, more.details = FALSE, ...){
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))

  thin <- object$dim[8]
  nsamp <- object$dim[9]
  pars <- matrix(object$parsamp, ncol = nsamp)
  ss <- unique(pmax(1, round(nsamp * (1:ntrace/ntrace))))
  post.burn <- (ss > nsamp * burn.perc)
  dimpars <- object$dim
  dimpars[7] <- length(ss)
  
  
  n <- object$dim[1]; p <- 0; nknots <- dimpars[4]; ngrid <- object$dim[5]; 

  psamp <- cbind(gam0 = pars[nknots + 1,ss], sigma = sigFn(pars[nknots + 2,ss]), nu = nuFn(pars[nknots + 3,ss]))
  sm <- .C("DEV", pars = as.double(pars[,ss]), y = as.double(object$y), cens = as.integer(object$cens), wt = as.double(object$wt), hyper = as.double(object$hyper), dim = as.integer(dimpars), gridmats = as.double(object$gridmats), tau.g = as.double(object$tau.g), devsamp = double(length(ss)), llsamp = double(length(ss)*n), pgsamp = double(length(ss)*ngrid), lnpdensendsamp = double(length(ss)*2), other.controls = as.integer(object$other.controls))
  deviance <- sm$devsamp
  ll <- matrix(sm$llsamp, ncol = length(ss))
  fit.waic <- waic(ll[,post.burn,drop=FALSE], print=FALSE)
  pg <- matrix(sm$pgsamp, ncol = length(ss))
  phi.ends <- t(matrix(sm$lnpdensendsamp, ncol=length(ss)))
  tail.scale <- exp(phi.ends[,2])*(psamp[,2]*psamp[,3])^(psamp[,3])
  prox.samp <- object$prox[apply(pg[1:ngrid,,drop=FALSE], 2, function(pr) sample(length(pr), 1, prob = pr))]
  
  if(more.details) par(mfrow = c(2,2), mar = c(5,4,3,2)+.1)
  if(plot.dev){
    plot(thin * ss, deviance, ty = "l", xlab = "Markov chain iteration", ylab = "Deviance", bty = "n", main = "Fit trace plot", ...)
    grid(col = "gray")
  }
  
  if(more.details){
    ngrid <- length(object$prox)
    prior.grid <- exp(object$gridmats[nrow(object$gridmats),])
    lam.priorsamp <- lamFn(sample(object$prox, ntrace, replace = TRUE, prob = prior.grid))
    lam.prior.q <- quantile(lam.priorsamp, pr = c(.025, .5, .975))
    lam.samp <- lamFn(prox.samp)
    a <- min(lamFn(object$prox))
    b <- diff(range(lamFn(object$prox))) * 1.2
    plot(thin * ss, lam.samp, ty = "n", ylim = a + c(0, b), bty = "n", ann = FALSE, axes = FALSE)
    axis(1)
    for(i in 1:1){
      abline(h = b * (i-1) + lamFn(object$prox), col = "gray")
      abline(h = b * (i - 1) + lam.prior.q, col = "red", lty = c(2,1,2))
      lines(thin * ss, b * (i-1) + lam.samp, lwd = 1, col = 4)
      if(i %% 2) axis(2, at = b * (i-1) + lamFn(object$prox[c(1,ngrid)]), labels = round(object$prox[c(1,ngrid)],2), las = 1, cex.axis = 0.6)
      mtext(substitute(beta[index], list(index = i - 1)), side = 4, line = 0.5, at = a + b * (i - 1) + 0.4*b, las = 1)
    }
    title(xlab = "Markov chain iteration", ylab = "Proxmity posterior", main = "Mixing over GP scaling")
    
    theta <- coda::as.mcmc(t(matrix(object$parsamp, ncol = nsamp)[,ss[post.burn]]))
    gg <- coda::geweke.diag(theta, .1, .5)
    zvals <- gg$z
    
    pp <- 2 * (1 - pnorm(abs(zvals)))
    plot(sort(pp), ylab = "Geweke p-values", xlab = "Parameter index (reordered)", main = "Convergence diagnosis", ty = "h", col = 4, ylim = c(0, 0.3), lwd = 2)
    abline(h = 0.05, col = 2, lty = 2)
    abline(a = 0, b = 0.1 / length(pp), col = 2, lty = 2)
    mtext(c("BH-10%", "5%"), side = 4, at = c(0.1, 0.05), line = 0.1, las = 0, cex = 0.6)
    
    npar <- length(object$par)
    image(1:npar, 1:npar, suppressWarnings(cor(theta)), xlab = "Parameter index", ylab = "Parameter index", main = "Parameter correlation")
    
  }
  invisible(list(deviance = deviance, pg = pg, prox = prox.samp, ll = ll, waic = fit.waic, tail.scale=tail.scale))
}

predict.sbde <- function(object, burn.perc = 0.5, nmc = 200, yRange = range(object$y), yLength = 401, ...){
  thin <- object$dim[8]
  nsamp <- object$dim[9]
  pars <- matrix(object$parsamp, ncol = nsamp)
  ss <- unique(round(nsamp * seq(burn.perc, 1, len = nmc + 1)[-1]))
  dimpars <- object$dim
  dimpars[7] <- length(ss)
  
  yGrid <- seq(yRange[1], yRange[2], len = yLength)
  n <- object$dim[1]; p <- 0; ngrid <- object$dim[5]
  dimpars[1] <- yLength
  
  pred <- .C("PRED", pars = as.double(pars[,ss]), yGrid = as.double(yGrid), hyper = as.double(object$hyper), dim = as.integer(dimpars), gridmats = as.double(object$gridmats), tau.g = as.double(object$tau.g), ldenssamp = double(length(ss)*yLength), other.controls = as.integer(object$other.controls))
  dens <- matrix(exp(pred$ldenssamp), ncol = length(ss))
  return(list(y = yGrid, fsamp = dens, fest = t(apply(dens, 1, quantile, pr = c(.025, .5, .975))), ss=ss))
}


waic <- function(logliks, print = TRUE){
  lppd <- sum(apply(logliks, 1, logmean))
  p.waic.1 <- 2 * lppd - 2 * sum(apply(logliks, 1, mean))
  p.waic.2 <- sum(apply(logliks, 1, var))
  waic.1 <- -2 * lppd + 2 * p.waic.1
  waic.2 <- -2 * lppd + 2 * p.waic.2
  if(print) cat("WAIC.1 =", round(waic.1, 2), ", WAIC.2 =", round(waic.2, 2), "\n")
  invisible(c(WAIC1 = waic.1, WAIC2 = waic.2))
}

lamFn <- function(prox) return(sqrt(-100*log(prox)))
#nuFn <- function(z) return(0.5 + 1.0 * log1p(exp(z)))
#nuFn.inv <- function(nu) return(log(expm1((nu - 0.5)/1.0)))
#nuFn <- function(z) return(0.5 + 5.5*exp(z/2))
#nuFn.inv <- function(nu) return(2*log((nu - 0.5)/5.5))
#nuFn <- function(z) return(1/(0.01 + 1.99/(1 + exp(-z))))
#nuFn.inv <- function(nu) return(-log(1.99/(1/nu - 0.01) - 1))
nuFn <- function(z) return(0.5 + 1.5*exp(z/1.5))
nuFn.inv <- function(nu) return(1.5*log((nu - 0.5)/1.5))
sigFn <- function(z) return(exp(z/2)) 
sigFn.inv <- function(s) return(2 * log(s))
unitFn <- function(u) return(pmin(1 - 1e-10, pmax(1e-10, u)))

sum_sq <- function(x) return(sum(x^2))
extract <- function(lo, vn) return(lo[[vn]])
logmean <- function(lx) return(max(lx) + log(mean(exp(lx - max(lx)))))
logsum <- function(lx) return(logmean(lx) + log(length(lx)))
shrinkFn <- function(x) return(1) ##(1/(1 + log(x)))
trape <- function(x, h, len = length(x)) return(c(0, cumsum(.5 * (x[-1] + x[-len]) * (h[-1] - h[-len]))))


klGP <- function(lam1, lam2, tau = seq(0,1,.1), tail.warp=c(0,0), expo=2, Qfn=Qfn.reg){
  nknots <- length(tau)
  dd <- outer(tau, tau, dfn, expo=expo, tail.warp=tail.warp, Qfn=Qfn)
  K1 <- exp(-lam1^2 * dd); diag(K1) <- 1 + 1e-10; R1 <- chol(K1); log.detR1 <- sum(log(diag(R1)))
  K2 <- exp(-lam2^2 * dd); diag(K2) <- 1 + 1e-10; R2 <- chol(K2); log.detR2 <- sum(log(diag(R2)))
  return(log.detR2-log.detR1 - 0.5 * (nknots - sum(diag(solve(K2, K1)))))
}



proxFn <- function(prox.Max, prox.Min, kl.step = 1, tau, tail.warp=c(0,0), expo=2, Qfn=Qfn.reg){
  prox.grid <- prox.Max
  j <- 1
  while(prox.grid[j] > prox.Min){
    prox1 <- prox.grid[j]
    prox2 <- prox.Min
    kk <- klGP(lamFn(prox1), lamFn(prox2), tau, tail.warp, expo, Qfn)
    while(kk > kl.step){
      prox2 <- (prox1 + prox2)/2
      kk <- klGP(lamFn(prox1), lamFn(prox2), tau, tail.warp, expo, Qfn)
    }
    j <- j + 1
    prox.grid <- c(prox.grid, prox2)
  }
  return(prox.grid)
}

transform.grid <- function(w, ticks, dists){
  return((1-dists) * w[ticks] + dists * w[ticks+1])
}

Qfn0 <- function(x, tail.warp) return(qbeta(x, 1+tail.warp[1], 1+tail.warp[2]))
Qfn0.inv <- function(y, tail.warp) return(pbeta(y, 1+tail.warp[1], 1+tail.warp[2]))
Qfn.reg <- function(x, tail.warp, v=c(.1,.9)) return(Qfn0(x, tail.warp) + (x - v[1])*(v[2]-Qfn0(v[2], tail.warp)) / diff(v) - (x - v[2])*(v[1]-Qfn0(v[1], tail.warp)) / diff(v))
Qfn.irreg <- function(x, tail.warp, v=c(.1,.9)) return(((v[1]-v[2])*Qfn0(x, tail.warp) + (v[2]*Qfn0(v[1],tail.warp)-v[1]*Qfn0(v[2],tail.warp)))/(Qfn0(v[1],tail.warp)-Qfn0(v[2],tail.warp)))
dfn <- function(x, y, expo=2, tail.warp=c(0,0), Qfn=Qfn.reg) return(abs(Qfn(x, tail.warp) - Qfn(y, tail.warp))^expo)
Qfn.inv <- function(y, tail.warp, v=c(.1,.9)) return(Qfn0.inv((Qfn0(v[1],tail.warp)*(y-v[2]) + Qfn0(v[2],tail.warp)*(v[1]-y))/(v[1]-v[2]), tail.warp))

base.functions <- list(## all base distributions for parametric transformations
  't' = list(support = "full",
             fix.nu = FALSE,
             s0 = function(nu = Inf) return(qt(0.9, df = nu)/qt(0.9, df = 6)),
             log_f0 = function(x, nu = Inf){
               s0 <- qt(0.9, df = nu)/qt(0.9, df = 6)
               return(dt(x*s0, df = nu, log = TRUE) + log(s0))},
             f0 = function(x, nu = Inf){
               s0 <- qt(0.9, df = nu)/qt(0.9, df = 6)
               return(dt(x * s0, df = nu) * s0)},
             F0 = function(x, nu = Inf){
               s0 <- qt(0.9, df = nu)/qt(0.9, df = 6)
               return(pt(x * s0, df = nu))},
             Q0 = function(u, nu = Inf) {
               s0 <- qt(0.9, df = nu)/qt(0.9, df = 6)
               return(qt(u, df = nu)/s0)}),
  't+' = list(support = "half-pos",
              fix.nu = FALSE,
              s0 = function(nu = Inf) return(qt(0.95, df = nu)/qt(0.95, df = 6)),
              log_f0 = function(x, nu = Inf){
                s0 <- qt(0.95, df = nu)/qt(0.95, df = 6)
                return(log(x > 0) * (log(2.0) + dt(x * s0, df = nu, log = TRUE) + log(s0)))},
              f0 = function(x, nu = Inf){
                s0 <- qt(0.95, df = nu)/qt(0.95, df = 6)
                return((x > 0) * 2.0 * dt(x*s0, df = nu)*s0)},
              F0 = function(x, nu = Inf){
                s0 <- qt(0.95, df = nu)/qt(0.95, df = 6)
                return((x > 0) * 2.0 * (pt(x*s0, df = nu) - 0.5))},
              Q0 = function(u, nu = Inf) {
                s0 <- qt(0.95, df = nu)/qt(0.95, df = 6)
                return(qt((1+u)/2, df = nu)/s0)}),
  'gpd' = list(support = "half-pos",
               fix.nu = FALSE,
               s0 = function(nu = 1) return(1),
               log_f0 = function(x, nu = 1) return(dgpd.local(x, shape = 1/nu, log = TRUE)),
               f0 = function(x, nu = 1) return(dgpd.local(x, shape = 1/nu)),
               F0 = function(x, nu = 1) return(pgpd.local(x, shape = 1/nu)),
               Q0 = function(u, nu = 1) return(qgpd.local(u, shape = 1/nu))),
  'gpd-rescaled' = list(support = "half-pos",
               fix.nu = FALSE,
               s0 = function(nu = 1) return(qgpd.local(0.9, shape = 1/nu)/qgpd.local(0.9, shape = 1/6)),
               log_f0 = function(x, nu = 1) {
                 s0 <- qgpd.local(0.9, shape = 1/nu)/qgpd.local(0.9, shape = 1/6)
                 return(dgpd.local(x * s0, shape = 1/nu, log = TRUE) + log(s0))},
               f0 = function(x, nu = 1) {
                 s0 <- qgpd.local(0.9, shape = 1/nu)/qgpd.local(0.9, shape = 1/6)
                 return(dgpd.local(x * s0, shape = 1/nu) * s0)},
               F0 = function(x, nu = 1) {
                 s0 <- qgpd.local(0.9, shape = 1/nu)/qgpd.local(0.9, shape = 1/6)
                 return(pgpd.local(x * s0, shape = 1/nu))},
               Q0 = function(u, nu = 1) {
                 s0 <- qgpd.local(0.9, shape = 1/nu)/qgpd.local(0.9, shape = 1/6)
                 return(qgpd.local((1+u)/2, shape = 1/nu) / s0)}),
  'unif' = list(support = "full",
                fix.nu = TRUE,
                s0 = function(nu = Inf) return(1),
                log_f0 = function(x, nu = Inf) return(dunif(x, -1,1, log = TRUE)),
                f0 = function(x, nu = Inf) return(dunif(x, -1,1)),
                F0 = function(x, nu = Inf) return(punif(x, -1,1)),
                Q0 = function(u, nu = Inf) return(qunif(u, -1,1))))


## Generalized Pareto Distribution
## defined here to avoid dependency on package 'evd'

dgpd.local <- function(x, shape=1, log=FALSE){
    val <- -(1/shape + 1) * log1p(shape*x)
    if(!log) val <- exp(val)
    return(val)
}
pgpd.local <- function(x, shape=1) return(1 - exp(-log1p(shape*x)/shape))
qgpd.local <- function(p, shape=1) return(expm1(-log1p(-p)*shape)/shape)



