\name{bootHARLEM}
\alias{bootHARLEM}

\title{
Bootstrap for Accelerated Residual Lifetime Models.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
bootHARLEM(fullData, deltaFun = function(data) max(diff(sort(data$w))), lambda, B = 1000, nGrid = 300, tau = 100, Q1 = NULL, Q2tmle = NULL, T0 = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fullData}{
%%     ~~Describe \code{fullData} here~~
}
  \item{deltaFun}{
%%     ~~Describe \code{deltaFun} here~~
}
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
}
  \item{B}{
%%     ~~Describe \code{B} here~~
}
  \item{nGrid}{
%%     ~~Describe \code{nGrid} here~~
}
  \item{tau}{
%%     ~~Describe \code{tau} here~~
}
  \item{Q1}{
%%     ~~Describe \code{Q1} here~~
}
  \item{Q2tmle}{
%%     ~~Describe \code{Q2tmle} here~~
}
  \item{T0}{
%%     ~~Describe \code{T0} here~~
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
function (fullData, deltaFun = function(data) max(diff(sort(data$w))), 
    lambda, B = 1000, nGrid = 300, tau = 100, Q1 = NULL, Q2tmle = NULL, 
    T0 = NULL) 
{
    xGrid <- seq(0, tau, length.out = nGrid)
    vGrid <- seq(0, tau, length.out = nGrid)
    if (is.null(Q2tmle)) {
        bootFun <- function(b) {
            cat(paste(100 * b/B, "\% ", sep = ""))
            bootData <- fullData[sample(1:nrow(fullData), replace = TRUE), 
                ]
            obsData <- bootData[which(bootData$delta0 * bootData$delta1 == 
                1), ]
            obsData <- obsData[obsData$x != 0, ]
            lambdaObs <- lambda(obsData$x)
            deltaStar <- deltaFun(bootData)
            bootQ0 <- initQ(fullData, xGrid, vGrid, deltaStar, 
                tau = tau)
            bootQ1 <- bootQ0[[1]]
            bootQ2 <- bootQ0[[2]]
            Sn <- bootQ1$Sn
            Qc <- bootQ1$Qc
            dGn <- bootQ1$dGn
            bootT0 <- initT(bootData, Qc, dGn, deltaStar, xGrid, 
                vGrid)
            bootT1 <- bootT0[[1]]
            bootT2 <- bootT0[[2]]
            bootTMLE <- harlemTMLE(bootData, bootQ0, bootT0, 
                stepLinear(bootData, deltaStar, tau, lambda, 
                  1/sqrt(nrow(obsData))), linearPsi(bootData, 
                  deltaStar, tau), lambda, deltaStar, tau, nGrid = nGrid, 
                FALSE, FALSE)
            bootTMLE$theta
        }
    }
    else {
        lambdax <- function(x) x * lambda(x)
        lambdax2 <- function(x) x^2 * lambda(x)
        i0lambda <- integrate(lambda, 0, tau)$value
        i1lambda <- integrate(lambdax, 0, tau)$value
        i2lambda <- integrate(lambdax2, 0, tau)$value
        Lambda <- matrix(c(i0lambda, i1lambda, i1lambda, i2lambda), 
            2, 2)
        invLambda <- solve(Lambda)
        lambdaGrid <- lambda(xGrid)
        T1 <- T0[[1]]
        T2 <- T0[[2]]
        slice1 <- lambdaGrid * apply(T1 * Q2tmle, 2, intFun, 
            vGrid)
        slice2 <- lambdaGrid * apply(T2 * Q2tmle, 2, intFun, 
            vGrid)
        grad <- t(apply(T1, 1, function(t1) {
            r <- t1/slice1
            r[(abs(t1) <= sqrt(.Machine$double.eps) & abs(slice1) <= 
                sqrt(.Machine$double.eps))] <- 1
            r[lambdaGrid == 0] <- 0
            return(r)
        }) - apply(T2, 1, function(t2) {
            r <- t2/slice2
            r[(abs(t2) <= sqrt(.Machine$double.eps) & abs(slice2) <= 
                sqrt(.Machine$double.eps))] <- 1
            r[lambdaGrid == 0] <- 0
            return(r)
        }))
        IC <- linInt(grad, xGrid, vGrid)
        bootFun <- function(b) {
            cat(paste(100 * b/B, "\% ", sep = ""))
            bootData <- fullData[sample(1:nrow(fullData), replace = TRUE), 
                ]
            obsData <- bootData[which(bootData$delta0 * bootData$delta1 == 
                1), ]
            lambdaObs <- lambda(obsData$x)
            deltaStar <- deltaFun(obsData)
            QcFit <- survfit(Surv(bootData$v - bootData$w, 1 - 
                bootData$delta1) ~ 1)
            Qc <- stepfun(QcFit$time, c(1, QcFit$surv), f = 0)
            Snfit <- survfit(Surv(bootData$w, bootData$v, bootData$delta1) ~ 
                1)
            Sn <- stepfun(Snfit$time, c(1, Snfit$surv), f = 0)
            dGn <- 1/Sn(bootData$w)
            dGn[Sn(bootData$w) == 0] <- 0
            dGn <- dGn/sum(dGn)
            bootQ1 <- list(Sn = Sn, dGn = dGn, Qc = Qc)
            bootT0 <- initT(bootData, Qc, dGn, deltaStar, xGrid, 
                vGrid)
            bootT1 <- bootT0[[1]]
            bootT2 <- bootT0[[2]]
            bootPsi <- linearPsi(bootData, deltaStar, tau)(bootT1, 
                bootT2, bootQ1, Q2tmle, lambda, xGrid, vGrid)
            bootPsi - invLambda \%*\% colSums(apply(cbind(1, obsData$x), 
                2, "*", IC(obsData$x, obsData$v)))/n
        }
    }
    bootRep <- lapply(1:B, bootFun)
    cat("\n")
    return(apply(do.call(cbind, bootRep), 1, var))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
