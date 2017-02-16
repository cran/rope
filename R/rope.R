#' FDR controlled model selection
#'
#' Estimates a model from bootstap counts. The objective is to maximize accuracy
#' while controlling the false discovery rate of selected variables. Developed
#' for high-dimensional models with number of variables in the order of at least
#' 10000.
#'
#' @param data Matrix of variable presence counts. One column for each variable,
#'   one row for each parameter value (e.g. levels of regularization).
#' @param B Number of bootstraps used to construct \code{data}. At least 21 are
#'   needed for u-shape test heuristic to work, but in general it is recommended
#'   to use many more.
#' @param fdr Vector of target false discovery rates to return selections for
#' @param mc.cores Number of threads to run in parallel (1 turns of
#'   parallelization)
#' @param only.first Skip second part of algorithm. Saves time but gives worse
#'   results.
#'
#' @return A list with components
#'   \item{selection}{matrix (one row for each fdr target, one column for
#'     each variable)}
#'   \item{q}{vector of q-values, one for each variable}
#'   \item{level}{index of most separating parameter value}
#'   \item{alt.prop}{estimated proportion of alternative variables}
#'
#' @examples
#' \dontrun{
#' data # a matrix of selection counts, for 100 bootstraps, with ncol(data)
#'   # potential variables counted for nrow(data) different penalization levels
#' fdr <- c(0.05, 0.1)
#' result <- rope(data, 100, fdr)
#' }
#'
#' @importFrom graphics abline hist lines par plot
#' @importFrom stats density isoreg ksmooth optim quantile rnorm
#'
#' @author Jonatan Kallus, \email{kallus@@chalmers.se}
#' @keywords htest models multivariate
#'
#' @export
rope <- function(data, B, fdr=0.1, mc.cores=getOption("mc.cores", 2L),
        only.first=FALSE) {
    hs <- as.hists(data, B)
    fl <- localfitall(hs, B, mc.cores)

    plot.data <- list()
    plot.data$fit1 <- fl
    plot.data$B <- B

    ## Check if no histogram was u-shaped, if so return empty selection
    if (length(fl) == 1 && is.na(fl)) {
        return(list(selection=matrix(FALSE, length(fdr), dim(data)[2]),
                    q=rep(1, dim(data)[2]), level=NA, alt.prop=0,
                    plot.data=plot.data))
    }

    if (only.first) {
        i <- which.max(fl$goods)
        pi <- fl$pitp[i]
        q <- q.values(data[i, ], fl$fits[[i]], hs[i, ])

        plot.data$q <- q

        q <- q$by.var
        selection <- matrix(NA, length(fdr), length(q))
        for (j in 1:length(fdr)) {
            selection[j, ] <- q < fdr[j]
        }

        return(list(selection=selection, q=q, level=i, alt.prop=pi,
            plot.data=plot.data))
    }

    ## Estimate proportion of alternative variables
    maxloc <- floor(max.conf(fl$goods, 0.5))
    pitp.monotone <- rev(isoreg(rev(fl$pitp[maxloc:length(fl$pitp)]))$yf)
    if (maxloc > 1) {
        pitp.monotone <- c(rep(pitp.monotone[1], maxloc-1), pitp.monotone)
    }
    maxloc <- max.conf(fl$goods, 0.95)
    i <- min(length(pitp.monotone), ceiling(maxloc))
    pi <- pitp.monotone[i]

    plot.data$fit1$goods.loc <- maxloc
    plot.data$fit1$pi <- pi
    plot.data$fit1$pitp.monotone <- pitp.monotone

    fl <- localfitall(hs, B, mc.cores, pi=pi)

    plot.data$fit2 <- fl

    ## Select best histogram for variable selection
    minloc <- max.conf(fl$goods, 0.05)
    i <- max(1, floor(minloc))

    plot.data$fit2$goods.loc <- minloc

    q <- q.values(data[i, ], fl$fits[[i]], hs[i, ])

    plot.data$q <- q

    q <- q$by.var
    selection <- matrix(NA, length(fdr), length(q))
    for (j in 1:length(fdr)) {
        selection[j, ] <- q < fdr[j]
    }

    res <- list(selection=selection, q=q, level=i, alt.prop=pi,
      plot.data=plot.data)
    class(res) <- 'rope'
    return(res)
}

#' Run first step of model fitting to find good penalization interval
#'
#' @param data Matrix of variable presence counts. One column for each variable,
#'   one row for each parameter value (e.g. levels of regularization).
#' @param B Number of bootstraps used to construct \code{data}. At least 21 are
#'   needed for u-shape test heuristic to work, but in general it is recommended
#'   to use many more.
#' @param mc.cores Number of threads to run in parallel (1 turns of
#'   parallelization)
#' @return A list with components
#'   \item{pop.sep}{vector of values saying how separated true and false
#'     variables are for each level of penalization}
#'
#' @export
explore <- function(data, B, mc.cores=getOption("mc.cores", 2L)) {
    hs <- as.hists(data, B)
    fl <- localfitall(hs, B, mc.cores)

    plot.data <- list()
    plot.data$fit1 <- fl
    plot.data$B <- B

    res <- list(pop.sep=fl$goods, plot.data=plot.data)
    class(res) <- 'rope'
    return(res)
}

#' Convenience wrapper for \code{rope} for adjacency matrices
#'
#' When modeling graphs it may be more convenient to store data as matrices
#' instead of row vectors.
#'
#' @param data List of symmetric matrices, one matrix for each penalization
#'   level
#' @param B Number of bootstraps used to construct \code{data}. At least 21 are
#'   needed for u-shape test heuristic to work, but in general it is recommended
#'   to use many more.
#' @param ... Additional arguments are passed on to \code{rope}.
#'
#' @return A list with components
#'   \item{selection}{list of symmetric matrices, one matrix for each fdr
#'     target}
#'   \item{q}{symmetric matrix of q-values}
#'   \item{level}{index of most separating parameter value}
#'   \item{alt.prop}{estimated proportion of alternative variables}
#'
#' @examples
#' \dontrun{
#' data # a list of symmetric matrices, one matrix for each penalization level,
#'   # each matrix containing selection counts for each edge over 100 bootstraps
#' fdr <- c(0.05, 0.1)
#' result <- rope(data, 100, fdr)
#' }
#'
#' @export
ropegraph <- function(data, B, ...) {
  data <- symmetric.matrix.list2matrix(data)
  res <- rope(data, B, ...)
  res$selection <- lapply(1:nrow(res$selection),
    function(i) 1 == vector2symmetric.matrix(res$selection[i, ]))
  res$q <- vector2symmetric.matrix(res$q)
  return(res)
}

#' Convenience wrapper for \code{explore} for adjacency matrices
#'
#' When modeling graphs it may be more convenient to store data as matrices
#' instead of row vectors.
#'
#' @param data List of symmetric matrices, one matrix for each penalization
#'   level
#' @param B Number of bootstraps used to construct \code{data}. At least 21 are
#'   needed for u-shape test heuristic to work, but in general it is recommended
#'   to use many more.
#' @param ... Additional arguments are passed on to \code{explore}.
#'
#' @return A list with components
#'   \item{pop.sep}{vector of values saying how separated true and false
#'     variables are for each level of penalization}
#'
#' @export
exploregraph <- function(data, B, ...) {
  data <- symmetric.matrix.list2matrix(data)
  explore(data, B, ...)
}

symmetric.matrix.list2matrix <- function(data) {
  do.call(rbind, lapply(1:length(data), function(i) {
    if (requireNamespace("Matrix", quietly = TRUE)) {
      return(Matrix::Matrix(symmetric.matrix2vector(data[[i]]), 1))
    } else {
      return(symmetric.matrix2vector(data[[i]]))
    }
  }))
}

#' Take upper half of matrix and convert it to a vector
#'
#' If variable selection counts are in a matrix this function converts them into
#' vector to input into rope. Can be useful when variables correspond to edges
#' in a graph.
#'
#' @param m A symmetric matrix
#'
#' @export
symmetric.matrix2vector <- function(m) {
    return(m[upper.tri(m)])
}

#' Convert vector that represents half of a symmetric matrix into a matrix
#'
#' This can be convenient for using output when rope is used for selection of
#' graph models.
#'
#' @param v A vector with length p*(p-1)/2 for some integer p
#'
#' @export
vector2symmetric.matrix <- function(v) {
    d <- length(v)
    p <- sqrt(2*d+1/4)+1/2 ## since d=p*(p-1)/2
    M <- matrix(0, p, p)
    M[upper.tri(M)] <- v
    M <- M + t(M)
    return(M)
}

# Bootstrap estimate of confidence interval for maximum index
max.conf <- function(y, probs=c(0.05, 0.95)) {
  if (length(y) == 1) {
    return(rep(y, length(probs)))
  }
  x <- 1:length(y)
  y <- y/sum(y)
  maxs <- sapply(1:10000, function(i) {
    ix <- x + rnorm(length(x))
    bw <- density(ix, weights=y)$bw
    kr <- ksmooth(ix, y, 'normal', bw)
    kr$x[which.max(kr$y)]
  })
  as.vector(quantile(maxs, probs))
}

#' Plot rope results
#'
#' @param result An object returned by \code{rope} or \code{explore}
#' @param data Matrix of variable presence counts. One column for each variable,
#'   one row for each parameter value (e.g. levels of regularization).
#' @param types List of names of plots to draw (alternatives \code{'global'},
#'   \code{'q-values'} or \code{'fits'})
#' @param ... Pass level=v for a vector v of indices when drawing the fits plot
#'   to only plot for penalization levels corresponding to v
#' @export
plotrope <- function(result, data, types=c('global'), ...) {
    for (type in types) {
        if (type == 'global') plot.global(result, data)
        if (type == 'q-values') plot.qvalues(result, data)
        if (type == 'fits') plot.fits(result, data, ...)
    }
}

# Plot histograms and fits
plot.fits <- function(results, counts, level=NULL) {
    fits <- results$plot.data$fit1$fits
    hix <- 1:length(fits)
    fitname <- rep('step 1', length(fits))
    if ('fit2' %in% names(results$plot.data)) {
        fits2 <- results$plot.data$fit2$fits
        fits <- c(fits, fits2)
        hix <- c(hix, 1:length(fits2))
        fitname <- c(fitname, rep('step 2', length(fits2)))
    }
    B <- results$plot.data$B
    l <- length(fits)
    if (l <= 15) {
        ploti <- 1:l
    } else {
        ploti <- round(seq(1, l, length=15))
    }
    if (!is.null(level)) {
        ploti <- level
    }
    mfrow(length(ploti)*2, TRUE)
    par(mar=rep(2, 4))
    if (is.list(counts)) {
      counts <- symmetric.matrix.list2matrix(counts)
    }
    hs <- as.hists(counts, B)
    for (i in ploti) {
        plotlocal(hs[hix[i], ], B, fits[[i]],
            main=paste('index', hix[i], fitname[i]))
    }
}

# Plot pi_tp vs lambda and "good area" vs lambda
plot.global <- function(result, counts) {
    fit1 <- result$plot.data$fit1
    par(mfrow=c(1, 2))
    if ('fit2' %in% names(result$plot.data)) {
        fit2 <- result$plot.data$fit2
        par(mfrow=c(2, 2))
    }
    xlab <- 'Regularization (index)'

    plot(fit1$pitp, xlab=xlab, ylab='Alternative proportion',
        main='Unconstrained fits')
    lines(fit1$pitp.monotone)
    abline(v=fit1$goods.loc, lty=3)

    plot(fit1$goods, xlab=xlab, ylab='Component separation',
        main='Unconstrained fits')
    lines(fit1$goods.convex$convex)
    abline(h=max(fit1$goods) - fit1$goods.convex$margin, lty=2)
    abline(v=fit1$goods.loc, lty=3)

    if ('fit2' %in% names(result$plot.data)) {
        plot(fit2$pitp, xlab=xlab, ylab='Alternative proportion',
            main='Constrained fits')
        abline(v=fit2$goods.loc, lty=3)

        plot(fit2$goods, xlab=xlab, ylab='Component separation',
            main='Constrained fits')
        lines(fit2$goods.convex$convex)
        abline(h=max(fit2$goods) - fit2$goods.convex$margin, lty=2)
        abline(v=fit2$goods.loc, lty=3)
    }
}

# Plot q-values vs k
plot.qvalues <- function(result, counts) {
  plot(result$plot.data$q$by.k, type='l', xlab='k', ylab='q')
}

# Create histograms
#
# @param data Matrix of variable presence counts. One column for each variable,
#   one row for each parameter value (e.g. levels of regularization).
# @param B Number of bootstraps used to construct \code{data}
#
# @return Matrix with same number of rows as \code{data}. Each row is a
#   histogram. Number of columns is equal to \code{B+1} (since count of zeros
#   is included).
as.hists <- function(data, B) {
    t(sapply(1:(dim(data)[1]), function(i) {
        hist(data[i, ], seq(-0.5, B+0.5), plot=FALSE)$count
    }))
}

test.u <- function(count, B) {
    h <- hist(count, seq(-0.5, B+0.5), plot=FALSE)$count
    b1 <- round(0.025*B)
    b2 <- ceiling(0.005*B)
    b3 <- ceiling(600/B)
    return(!(h[B] <= b3 || mean(h[(B-b1):(B-b2)]) >= mean(h[(B-b2+1):(B+1)])))
}

# Fit model locally for each lambda
#
# Also filter out histograms with too little data in right part. Also does
# plotting for manual examination.
#
# @param hs Histograms
# @param B Number of bootstraps used to construct \code{data}
# @param mc.cores Number of threads to run in parallel (1 turns of
#   parallelization)
# @param plot Create figures?
# @param hts Histograms of truly alternative variables (only used for plotting)
# @param pi Upper bound on proportion of variables in alternative population
#
# @return List of "good area", alternative population proportion and estimated
#   model for each lambda
localfitall <- function(hs, B, mc.cores, plot=FALSE, hts=NULL, pi=NULL) {
    l <- dim(hs)[1]

    ## filter out hists with too little data in right part
    b1 <- round(0.025*B)
    b2 <- ceiling(0.005*B)
    b3 <- ceiling(600/B)
    for (i in 1:l) {
        if (hs[i, B] <= b3 ||
            mean(hs[i, (B-b1):(B-b2)]) >= mean(hs[i, (B-b2+1):(B+1)])) {
            l <- i-1
            break
        }
    }
    if (l == 0) {
        warning('No histogram with u-shape')
        return(NA)
    }

    if (plot) {
        if (l <= 15) {
            ploti <- 1:l
        } else {
            ploti <- round(seq(1, l, length=15))
        }
        mfrow(length(ploti)*2, TRUE)
        par(mar=rep(2, 4))
    }

    if (plot) { ## remove this branch when plotting is separated
        fits <- vector('list', l)
        for (i in 1:l) {
            fit <- localfit(hs[i, ], B, i, pi)

            if (plot && i %in% ploti) {
                plotlocal(hs[i, ], B, fit, hts[i, ])
            }
            fits[[i]] <- fit
        }
    } else if (requireNamespace("parallel", quietly = TRUE)) {
        fits <- parallel::mclapply(1:l, function(i) {
            localfit(hs[i, ], B, i, pi)
        }, mc.cores=mc.cores, mc.preschedule=FALSE)
    } else {
        fits <- lapply(1:l, function(i) {
            localfit(hs[i, ], B, i, pi)
        })
    }

    goods <- rep(NA, l)
    pitp <- rep(NA, l)
    for (i in 1:l) {
        pitp[i] <- sum(fits[[i]]$comps[3:4, ])
        goods[i] <- goodarea(fits[[i]]$comps[3, ]+fits[[i]]$comps[4, ],
                             fits[[i]]$comps[1, ]+fits[[i]]$comps[2, ])
    }

    return(list(goods=goods, pitp=pitp, fits=fits))
}

# Set number of rows and columns for plotting
#
# @param size Number of figures to plot
# @param force.even Enforce that the number of columns is even
#
# @return None
mfrow <- function(size, force.even=FALSE) {
    cols <- ceiling(sqrt(size))
    if (force.even) cols <- cols + (cols %% 2)
    rows <- ceiling(size/cols)
    par(mfrow=c(rows, cols))
    return(c(rows, cols))
}

# Fit model for one histogram (i.e. locally at one lambda value)
#
# @param h Histogram (vector of length \code{B+1})
# @param B Number of bootstraps used to construct \code{data}
# @param ix Index of histogram \code{h}, used only for debug output
# @param pi Upper bound on proportion of variables in alternative population
#
# @return List of individual components (matrix), total distribution and
#   optimal parameters
localfit <- function(h, B, ix, pi=NULL) {
    zero <- 1e-6
    one <- 1-zero

    hcum <- cumsum(h)
    ## model 25% of data with highest counts as powered BB mixture
    ## (variables with maximal count will also be excluded)
    c <- which.min(abs(hcum/hcum[length(hcum)]-0.75))

    myoptim <- function(par, f, lower, upper) {
        sol <- optim(par, f, method='L-BFGS-B', lower=lower, upper=upper,
                     control=list(fnscale=1, parscale=rep(0.0001, length(par)),
                                  maxit=1000, factr=.Machine$double.eps))
        if (sol$convergence != 0) {
            sol2 <- optim(par, f, method='L-BFGS-B', lower=lower, upper=upper,
                          control=list(fnscale=1,
                                       parscale=rep(0.00001, length(par)),
                                       maxit=1000, factr=.Machine$double.eps))
            if (sol2$value < sol$value) sol <- sol2
            if (sol$convergence != 0) {
                ## last resort: naive local search
                nlls <- function(p) {
                    if (any(p < lower)) return(Inf)
                    if (any(p > upper)) return(Inf)
                    return(nll(p))
                }
                for (i in 1:100) {
                    cands <- matrix(rep(sol$par, dim(permutations)[2]),
                        length(par), dim(permutations)[2]) +
                        permutations*0.001*i
                    vs <- apply(cands, 2, nlls)
                    j <- which.min(vs)
                    if (sol$value <= vs[j]) {
                        return(sol$par)
                    } else {
                        sol$par <- cands[, j]
                        sol$value <- vs[j]
                    }
                }
                warning(paste('Optimization did not converge for histogram', ix,
                              ':', sol$message))
            }
        }
        return(sol$par)
    }

    logpmf <- function(x, B) {
        term <- lgamma(B+1)-lgamma(x+1)-lgamma(B-x+1)
        log.dbb <- function(mu, sigma) {
            term + lgamma(1/sigma) + lgamma(x+mu/sigma) +
                lgamma(B+(1-mu)/sigma-x) - lgamma(mu/sigma) -
                lgamma((1-mu)/sigma) - lgamma(B+(1/sigma))
        }
        
        function(theta) log.dbb(theta[1], theta[2])
    }

    bb <- logpmf((c+1):(B-1), B)
    zpbb <- function (musigma) {
        x <- exp(bb(c(musigma, musigma)))
        if (sum(x) == 0) x[] <- zero
        x <- x/sum(x)
        x <- x - min(x)
        x[x < 0] <- 0
        if (sum(x) == 0) x[] <- zero
        x/sum(x)
    }
    pbb <- function(theta) {
        mu <- theta[1]
        sigma <- theta[2]
        gamma <- theta[3]
        x <- exp(bb(c(mu, sigma)))^gamma
        if (sum(x) == 0) x[] <- zero
        x/sum(x)
    }
    f <- function(theta) {
        fp <- pbb(theta[1:3])
        tp <- zpbb(theta[4])
        pi <- theta[5]
        x <- pi*tp + (1-pi)*fp
        x[x < zero] <- zero
        log(x)-log(sum(x))
    }
    growfast <- function(x) min(exp(exp(x))/exp(1), 10)
    if (is.null(pi)) {
        pim <- one
    } else {
        pim <- min(one, (pi*sum(h)-h[B+1])/sum(h[(c+2):B]))
    }
    spim <- min(0.005, pim)
    nll <- function(theta) {
        fac <- 1
        ## fp is decreasing in its right-most part
        if (theta[1]+theta[2] > 1) fac <-growfast(theta[1]+theta[2]-1)
        -sum(h[(c+2):B] * f(theta))*fac/B
    }

    ##            1      2      3            4      5
    ##       FP: mu, sigma, gamma, TP: musigma,    pi
    start <- c(0.25,  0.25,     1,        0.99,  spim)
    lower <- c(zero,  zero,     0,         0.5,  zero)
    upper <- c( one,   one,  1000,     1-1e-10,   pim)

    par <- myoptim(start, nll, lower, upper)

    if (par[1]+par[2] > 1) {
        warning('bad solution')
    }

    comps <- matrix(0, 4, B+1)
    comps[1, 1:(c+1)] <- sum(h[1:(c+1)])/sum(h)/(c+1)
    comps[2, (c+2):B] <- (1-par[5])*pbb(par[1:3])*sum(h[(c+2):B])/sum(h)
    comps[3, (c+2):B] <- par[5]*zpbb(par[4])*sum(h[(c+2):B])/sum(h)
    comps[4, B+1] <- h[B+1]/sum(h)

    bb <- logpmf(0:B, B)
    zpbb <- function (musigma) {
        x <- exp(bb(c(musigma, musigma)))
        if (sum(x) == 0) x[] <- zero
        x <- x/sum(x[(c+2):B])
        x <- x - min(x)
        x[x < 0] <- 0
        if (sum(x) == 0) x[] <- zero
        x/sum(x[(c+2):B])
    }
    pbb <- function(theta) {
        mu <- theta[1]
        sigma <- theta[2]
        gamma <- theta[3]
        x <- exp(bb(c(mu, sigma)))^gamma
        if (sum(x) == 0) x[] <- zero
        x/sum(x[(c+2):B])
    }
    fp <- (1-par[5])*pbb(par[1:3])[B+1]*sum(h[(c+2):B])/sum(h)
    comps[2, B+1] <- fp
    comps[4, B+1] <- comps[4, B+1]-fp
    tp <- par[5]*zpbb(par[4])[1:(c+1)]*sum(h[(c+2):B])/sum(h)
    comps[3, 1:(c+1)] <- tp
    comps[1, 1:(c+1)] <- comps[1, 1:(c+1)]-mean(tp)

    list(comps=comps, dist=apply(comps, 2, sum), par=par)
}

# Plot model fit for one histogram
#
# @param h Histogram (vector of length \code{B+1})
# @param B Number of bootstraps used to construct \code{data}
# @param fit Result from \code{localfit}
# @param ... Arguments passed on top plot.default
#
# @return None
plotlocal <- function(h, B, fit, ...) {
    plot(h/sum(h), ylim=c(0, max(fit$comp[2, ])), ...)
    lines(fit$dist, col='red', lwd=3)
    for (j in 1:4) lines(fit$comps[j, ], col='blue')
    plot(h/sum(h), xlim=c(0.5*B, B), ylim=c(0, max(fit$comps[2:3, B])), ...)
    lines(fit$dist, col='red')
    for (j in 1:4) lines(fit$comps[j, ], col='blue')
}

# Calculate difference between true positive rate and false positive rate at a
# maximum accuracy threshold
#
# @param t Estimated distribution of alternative variables
# @param f Estimated distribution of null variables
#
# @return Size of "good area"
goodarea <- function(t, f) {
    good <- t-f
    good[good < 0] <- 0
    sum(good)
}

# Calculate q-value for each variable
#
# @param d Vector of variable presence counts
# @param fit Estimated model for \code{d}
# @param h Histogram for \code{d}
#
# @return Vector of q-values (same length as \code{d})
q.values <- function(d, fit, h) {
    revcumsum <- function(v) rev(cumsum(rev(v)))
    fp <- apply(fit$comps[1:2, ], 2, sum)
    fpc <- revcumsum(fp)
    pf <- fpc/revcumsum(fit$dist)
    hc <- revcumsum(sum(h)*fit$dist)
    pfu <- pf + 2 * sqrt(hc * pf * (1-pf)) / hc
    ## conf int can be approximated with
    ## revcumsum(sapply(1:length(h), function(i) rbinom(1, h[i],
    ##   fp[i]/fit$dist[i])))/revcumsum(h)
    return(list(by.var=pfu[d+1], by.k=pfu))
}

## 242 directions in 5D space generated with
##   library(gtools)
##   dput(t(permutations(3, 5, c(-1, 0, 1), repeats.allowed=TRUE)[-122, ]))
permutations <- structure(c(-1, -1, -1, -1, -1, -1, -1, -1, -1, 0, -1, -1, -1,
-1, 1, -1, -1, -1, 0, -1, -1, -1, -1, 0, 0, -1, -1, -1, 0, 1,
-1, -1, -1, 1, -1, -1, -1, -1, 1, 0, -1, -1, -1, 1, 1, -1, -1,
0, -1, -1, -1, -1, 0, -1, 0, -1, -1, 0, -1, 1, -1, -1, 0, 0,
-1, -1, -1, 0, 0, 0, -1, -1, 0, 0, 1, -1, -1, 0, 1, -1, -1, -1,
0, 1, 0, -1, -1, 0, 1, 1, -1, -1, 1, -1, -1, -1, -1, 1, -1, 0,
-1, -1, 1, -1, 1, -1, -1, 1, 0, -1, -1, -1, 1, 0, 0, -1, -1,
1, 0, 1, -1, -1, 1, 1, -1, -1, -1, 1, 1, 0, -1, -1, 1, 1, 1,
-1, 0, -1, -1, -1, -1, 0, -1, -1, 0, -1, 0, -1, -1, 1, -1, 0,
-1, 0, -1, -1, 0, -1, 0, 0, -1, 0, -1, 0, 1, -1, 0, -1, 1, -1,
-1, 0, -1, 1, 0, -1, 0, -1, 1, 1, -1, 0, 0, -1, -1, -1, 0, 0,
-1, 0, -1, 0, 0, -1, 1, -1, 0, 0, 0, -1, -1, 0, 0, 0, 0, -1,
0, 0, 0, 1, -1, 0, 0, 1, -1, -1, 0, 0, 1, 0, -1, 0, 0, 1, 1,
-1, 0, 1, -1, -1, -1, 0, 1, -1, 0, -1, 0, 1, -1, 1, -1, 0, 1,
0, -1, -1, 0, 1, 0, 0, -1, 0, 1, 0, 1, -1, 0, 1, 1, -1, -1, 0,
1, 1, 0, -1, 0, 1, 1, 1, -1, 1, -1, -1, -1, -1, 1, -1, -1, 0,
-1, 1, -1, -1, 1, -1, 1, -1, 0, -1, -1, 1, -1, 0, 0, -1, 1, -1,
0, 1, -1, 1, -1, 1, -1, -1, 1, -1, 1, 0, -1, 1, -1, 1, 1, -1,
1, 0, -1, -1, -1, 1, 0, -1, 0, -1, 1, 0, -1, 1, -1, 1, 0, 0,
-1, -1, 1, 0, 0, 0, -1, 1, 0, 0, 1, -1, 1, 0, 1, -1, -1, 1, 0,
1, 0, -1, 1, 0, 1, 1, -1, 1, 1, -1, -1, -1, 1, 1, -1, 0, -1,
1, 1, -1, 1, -1, 1, 1, 0, -1, -1, 1, 1, 0, 0, -1, 1, 1, 0, 1,
-1, 1, 1, 1, -1, -1, 1, 1, 1, 0, -1, 1, 1, 1, 1, 0, -1, -1, -1,
-1, 0, -1, -1, -1, 0, 0, -1, -1, -1, 1, 0, -1, -1, 0, -1, 0,
-1, -1, 0, 0, 0, -1, -1, 0, 1, 0, -1, -1, 1, -1, 0, -1, -1, 1,
0, 0, -1, -1, 1, 1, 0, -1, 0, -1, -1, 0, -1, 0, -1, 0, 0, -1,
0, -1, 1, 0, -1, 0, 0, -1, 0, -1, 0, 0, 0, 0, -1, 0, 0, 1, 0,
-1, 0, 1, -1, 0, -1, 0, 1, 0, 0, -1, 0, 1, 1, 0, -1, 1, -1, -1,
0, -1, 1, -1, 0, 0, -1, 1, -1, 1, 0, -1, 1, 0, -1, 0, -1, 1,
0, 0, 0, -1, 1, 0, 1, 0, -1, 1, 1, -1, 0, -1, 1, 1, 0, 0, -1,
1, 1, 1, 0, 0, -1, -1, -1, 0, 0, -1, -1, 0, 0, 0, -1, -1, 1,
0, 0, -1, 0, -1, 0, 0, -1, 0, 0, 0, 0, -1, 0, 1, 0, 0, -1, 1,
-1, 0, 0, -1, 1, 0, 0, 0, -1, 1, 1, 0, 0, 0, -1, -1, 0, 0, 0,
-1, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1, 0, 0, 0,
1, -1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, -1, -1, 0, 0, 1,
-1, 0, 0, 0, 1, -1, 1, 0, 0, 1, 0, -1, 0, 0, 1, 0, 0, 0, 0, 1,
0, 1, 0, 0, 1, 1, -1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, -1,
-1, -1, 0, 1, -1, -1, 0, 0, 1, -1, -1, 1, 0, 1, -1, 0, -1, 0,
1, -1, 0, 0, 0, 1, -1, 0, 1, 0, 1, -1, 1, -1, 0, 1, -1, 1, 0,
0, 1, -1, 1, 1, 0, 1, 0, -1, -1, 0, 1, 0, -1, 0, 0, 1, 0, -1,
1, 0, 1, 0, 0, -1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1,
-1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, -1, -1, 0, 1, 1, -1,
0, 0, 1, 1, -1, 1, 0, 1, 1, 0, -1, 0, 1, 1, 0, 0, 0, 1, 1, 0,
1, 0, 1, 1, 1, -1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, -1, -1, -1,
-1, 1, -1, -1, -1, 0, 1, -1, -1, -1, 1, 1, -1, -1, 0, -1, 1,
-1, -1, 0, 0, 1, -1, -1, 0, 1, 1, -1, -1, 1, -1, 1, -1, -1, 1,
0, 1, -1, -1, 1, 1, 1, -1, 0, -1, -1, 1, -1, 0, -1, 0, 1, -1,
0, -1, 1, 1, -1, 0, 0, -1, 1, -1, 0, 0, 0, 1, -1, 0, 0, 1, 1,
-1, 0, 1, -1, 1, -1, 0, 1, 0, 1, -1, 0, 1, 1, 1, -1, 1, -1, -1,
1, -1, 1, -1, 0, 1, -1, 1, -1, 1, 1, -1, 1, 0, -1, 1, -1, 1,
0, 0, 1, -1, 1, 0, 1, 1, -1, 1, 1, -1, 1, -1, 1, 1, 0, 1, -1,
1, 1, 1, 1, 0, -1, -1, -1, 1, 0, -1, -1, 0, 1, 0, -1, -1, 1,
1, 0, -1, 0, -1, 1, 0, -1, 0, 0, 1, 0, -1, 0, 1, 1, 0, -1, 1,
-1, 1, 0, -1, 1, 0, 1, 0, -1, 1, 1, 1, 0, 0, -1, -1, 1, 0, 0,
-1, 0, 1, 0, 0, -1, 1, 1, 0, 0, 0, -1, 1, 0, 0, 0, 0, 1, 0, 0,
0, 1, 1, 0, 0, 1, -1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1,
-1, -1, 1, 0, 1, -1, 0, 1, 0, 1, -1, 1, 1, 0, 1, 0, -1, 1, 0,
1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, -1, 1, 0, 1, 1, 0, 1, 0,
1, 1, 1, 1, 1, -1, -1, -1, 1, 1, -1, -1, 0, 1, 1, -1, -1, 1,
1, 1, -1, 0, -1, 1, 1, -1, 0, 0, 1, 1, -1, 0, 1, 1, 1, -1, 1,
-1, 1, 1, -1, 1, 0, 1, 1, -1, 1, 1, 1, 1, 0, -1, -1, 1, 1, 0,
-1, 0, 1, 1, 0, -1, 1, 1, 1, 0, 0, -1, 1, 1, 0, 0, 0, 1, 1, 0,
0, 1, 1, 1, 0, 1, -1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1,
-1, -1, 1, 1, 1, -1, 0, 1, 1, 1, -1, 1, 1, 1, 1, 0, -1, 1, 1,
1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, -1, 1, 1, 1, 1, 0, 1, 1,
1, 1, 1), .Dim = c(5L, 242L))

