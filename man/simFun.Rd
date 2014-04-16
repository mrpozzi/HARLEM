\name{dataGenHARLEM}
\alias{dataGenHARLEM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Data Generating Distribution for Accelerated Residual Lifetime Models.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dataGenHARLEM(delta = 3, tau = 100, lambda = Vectorize(function(x) 1 * ((x >= 10) & (x <= 70))), par = rep(1, 2), rate = 0.017, copula = frankCopula(-5), margins = c("beta", "beta"), paramMargins = list(list(shape1 = 20, shape2 = 5), list(shape1 = 5, shape2 = 7)), nGrid = 200)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{delta}{
%%     ~~Describe \code{delta} here~~
}
  \item{tau}{
%%     ~~Describe \code{tau} here~~
}
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
}
  \item{par}{
%%     ~~Describe \code{par} here~~
}
  \item{rate}{
%%     ~~Describe \code{rate} here~~
}
  \item{copula}{
%%     ~~Describe \code{copula} here~~
}
  \item{margins}{
%%     ~~Describe \code{margins} here~~
}
  \item{paramMargins}{
%%     ~~Describe \code{paramMargins} here~~
}
  \item{nGrid}{
%%     ~~Describe \code{nGrid} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (delta = 3, tau = 100, lambda = Vectorize(function(x) 1 * 
    ((x >= 10) & (x <= 70))), par = rep(1, 2), rate = 0.017, 
    copula = frankCopula(-5), margins = c("beta", "beta"), paramMargins = list(list(shape1 = 20, 
        shape2 = 5), list(shape1 = 5, shape2 = 7)), nGrid = 200) 
{
    require("copula")
    simCop <- mvdc(copula, margins, paramMargins)
    gridX <- seq(0, tau, length.out = nGrid)
    lambdaGrid <- lambda(gridX)
    Sz <- function(tau) {
        para <- paste(unlist(lapply(1:length(simCop@paramMargins[[1]]), 
            function(j) paste(names(simCop@paramMargins[[1]])[j], 
                simCop@paramMargins[[1]][j], sep = "="))), collapse = ",")
        function(u) {
            u <- u/tau
            1 - eval(parse(text = paste("p", simCop@margins[1], 
                "(u,", para, ")", sep = "")))
        }
    }
    Sz <- Sz(tau)
    Sz <- Vectorize(Sz)
    intSz <- function(x) integrate(Sz, lower = x, upper = tau)$value
    intSz <- Vectorize(intSz)
    num <- intSz(gridX)
    den <- Sz(gridX)
    R1 <- num/den
    R1[num < sqrt(.Machine$double.eps) & den < sqrt(.Machine$double.eps)] <- 1
    R2 <- unlist(lapply(gridX, function(x) {
        slice1 <- function(z) {
            res <- (z - x) * dMvdc(cbind(z/tau, x/tau), simCop)
            if (is.nan(res)) 
                return(0)
            return(res)
        }
        slice1 <- Vectorize(slice1)
        slice2 <- function(z) {
            res <- dMvdc(cbind(z/tau, x/tau), simCop)
            if (is.nan(res)) 
                return(0)
            return(res)
        }
        slice2 <- Vectorize(slice2)
        num <- integrate(slice1, x + delta, tau)$value
        den <- integrate(slice2, x + delta, tau)$value
        res <- num/den
        res[num < sqrt(.Machine$double.eps) & den < sqrt(.Machine$double.eps)] <- 1
        res
    }))
    R <- log(R1) - log(R2)
    R[is.infinite(R)] <- 0
    Rfun <- function(gridX) {
        fitR <- interpSpline(R[lambdaGrid != 0] ~ gridX[lambdaGrid != 
            0])
        function(x) {
            if (lambda(x) == 0) 
                return(0)
            predict(fitR, x)$y
        }
    }
    Rfun <- Rfun(gridX)
    Rfun <- Vectorize(Rfun)
    lamR <- function(x) lambda(x) * Rfun(x)
    lamRx <- function(x) lambda(x) * x * Rfun(x)
    lamR <- Vectorize(lamR)
    lamRx <- Vectorize(lamRx)
    intR <- c(integrate(lamR, 0, tau)$value, integrate(lamRx, 
        0, tau)$value)
    lambdax <- function(x) x * lambda(x)
    lambdax2 <- function(x) x^2 * lambda(x)
    i0lambda <- integrate(lambda, 0, tau)$value
    i1lambda <- integrate(lambdax, 0, tau)$value
    i2lambda <- integrate(lambdax2, 0, tau)$value
    Lambda <- matrix(c(i0lambda, i1lambda, i1lambda, i2lambda), 
        2, 2)
    theta <- solve(Lambda, intR)
    simHARLEM <- function(n, seed = sample(.Random.seed, 1)) {
        set.seed(seed)
        cat("seed:")
        print(seed)
        fullData <- data.frame(w = rep(NA, n), x = rep(NA, n), 
            v = rep(NA, n), delta0 = rep(NA, n), delta1 = rep(NA, 
                n))
        count <- 0
        while (count < n) {
            xv <- tau * rMvdc(n - count, simCop)
            delta0 <- (xv[, 1] > xv[, 2] + delta) * 1
            w <- tau * rbeta(n - count, par[1], par[2])
            valid <- (xv[, 1] > w)
            nValid <- sum(valid)
            if (nValid != 0) {
                ind <- (count + 1):(count + nValid)
                c <- rexp(nValid, rate)
                fullData$w[ind] <- w[valid]
                delta0[valid][(xv[valid, 2] >= w[valid])] <- 0L
                fullData$delta0[ind] <- delta0[valid]
                fullData$v[ind] <- apply(cbind(w[valid] + c, 
                  xv[valid, 1]), 1, min)
                fullData$x[ind] <- xv[valid, 2] * delta0[valid]
                fullData$delta1[ind] <- 1 * (xv[valid, 1] <= 
                  w[valid] + c)
            }
            count <- count + sum(valid)
        }
        attr(fullData, "seed") <- seed
        fullData
    }
    return(list(theta = theta, simHARLEM = simHARLEM, lambda = lambda))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
