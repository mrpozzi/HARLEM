\name{harlemTMLE}
\alias{harlemTMLE}

\title{
TMLE for Accelerated Residual Lifetime Models.
}

\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
harlemTMLE(fullData, Q0 = NULL, T0 = NULL, targetStep, Psi, lambda, deltaStar, tau, nGrid, verbose = FALSE, display = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fullData}{
%%     ~~Describe \code{fullData} here~~
}
  \item{Q0}{
%%     ~~Describe \code{Q0} here~~
}
  \item{T0}{
%%     ~~Describe \code{T0} here~~
}
  \item{targetStep}{
%%     ~~Describe \code{targetStep} here~~
}
  \item{Psi}{
%%     ~~Describe \code{Psi} here~~
}
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
}
  \item{deltaStar}{
%%     ~~Describe \code{deltaStar} here~~
}
  \item{tau}{
%%     ~~Describe \code{tau} here~~
}
  \item{nGrid}{
%%     ~~Describe \code{nGrid} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
}
  \item{display}{
%%     ~~Describe \code{display} here~~
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
function (fullData, Q0 = NULL, T0 = NULL, targetStep, Psi, lambda, 
    deltaStar, tau, nGrid, verbose = FALSE, display = FALSE) 
{
    if (nGrid\%\%2 == 1) 
        nGrid <- nGrid + 1
    xGrid <- seq(0, tau, length.out = nGrid)
    vGrid <- seq(0, tau, length.out = nGrid)
    w <- fullData$w
    obsData <- fullData[which(fullData$delta0 * fullData$delta1 == 
        1), ]
    obsData <- obsData[obsData$x != 0, ]
    if (is.null(Q0)) {
        if (verbose) 
            cat("Initializing the Q's...\n")
        Q0 <- initQ(fullData, xGrid, vGrid, deltaStar = deltaStar, 
            tau = tau, verbose = verbose)
    }
    Q1 <- Q0[[1]]
    Q2 <- Q0[[2]]
    Sn <- Q1$Sn
    Qc <- Q1$Qc
    dGn <- Q1$dGn
    if (is.null(T0)) {
        if (verbose) 
            cat("Initializing the T's...\n")
        T0 <- initT(fullData, Qc, dGn, deltaStar = deltaStar, 
            xGrid, vGrid)
    }
    T1 <- T0[[1]]
    T2 <- T0[[2]]
    lambdaGrid <- lambda(xGrid)
    lambdaObs <- lambda(obsData$x)
    if (display) {
        dev.set(which = 1)
        dev.set(which = 1)
        dev.set(which = 1)
    }
    count <- 0
    cond <- FALSE
    targetStep <- targetStep(T1, T2, Q1, xGrid, vGrid, verbose, 
        display)
    if (display | verbose) 
        thetaHat <- Psi(T1, T2, Q1, Q2, lambda, xGrid, vGrid)
    epsOpt <- normGrad <- logL <- c()
    while (!cond) {
        count <- count + 1
        cat(count, " ")
        new <- targetStep(Q2)
        if (verbose) 
            cat("\n")
        Q2 <- new$Q2
        cond <- new$cond
        logL <- c(logL, new$logL)
        normGrad <- rbind(normGrad, new$normGrad)
        epsOpt <- rbind(epsOpt, new$epsOpt)
        if (display) {
            dev.set(which = 2)
            par(mfrow = c(1, 2))
            plot(1:count, normGrad[, 1], xlab = "iter", ylab = expression(D[1]/sigma[1]), 
                type = "l", col = "blue")
            abline(h = 0, col = "green")
            points(1:count, normGrad[, 1], col = "red", pch = 16)
            plot(1:count, normGrad[, 2], xlab = "iter", ylab = expression(D[2]/sigma[2]), 
                type = "l", col = "blue")
            points(1:count, normGrad[, 2], col = "red", pch = 16)
            abline(h = 0, col = "green")
            dev.set(which = 3)
            plot(1:count, logL, xlab = "iter", ylab = "logLikelihood", 
                type = "l", col = "blue")
            points(1:count, logL, col = "red", pch = 16)
        }
    }
    cat("\n")
    theta <- Psi(T1, T2, Q1, Q2, lambda, xGrid, vGrid)
    return(list(theta = theta, Q = list(Q1 = Q1, Q2 = Q2), logL = logL, 
        normGrad = normGrad, epsOpt = epsOpt))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
