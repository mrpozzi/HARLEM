\name{Linear Model}
\alias{stepLinear}
\alias{linearPsi}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
TMLE Step and Calculation of the Parameter in the Linear case
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
stepLinear(fullData, deltaStar, tau, lambda, tol)
linearPsi(fullData, deltaStar,tau)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fullData}{
%%     ~~Describe \code{fullData} here~~
}
  \item{deltaStar}{
%%     ~~Describe \code{deltaStar} here~~
}
  \item{tau}{
%%     ~~Describe \code{tau} here~~
}
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
}
  \item{tol}{
%%     ~~Describe \code{tol} here~~
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
function (fullData, deltaStar, tau, lambda, tol) 
{
    obsData <- fullData[which(fullData$delta0 * fullData$delta1 == 
        1), ]
    obsData <- obsData[obsData$x != 0, ]
    w <- fullData$w
    v <- obsData$v
    x <- obsData$x
    n <- nrow(obsData)
    function(T1, T2, Q1, xGrid, vGrid, verbose = FALSE, display = FALSE) {
        if (display) {
            dev.set(which = 4)
            plot(0, 0, type = "n", xlim = c(-0.3, 0.3), ylim = c(-0.3, 
                0.3), xlab = expression(D[1]/sigma[1]), ylab = expression(D[2]/sigma[2]))
            abline(h = 0, col = "blue")
            abline(v = 0, col = "blue")
        }
        lambdaGrid <- lambda(xGrid)
        nx <- length(xGrid)
        nv <- length(vGrid)
        function(Q2) {
            delta <- deltaStar
            storage.mode(lambdaGrid) <- "double"
            logL <- 0
            normGrad <- rep(0, 2)
            epsOpt <- rep(0, 2)
            z <- .C("tmleStep", a_delta = delta, a_tau = tau, 
                a_n = n, a_xx = x, a_vv = v, a_lambda = lambdaGrid, 
                a_xGrid = xGrid, a_vGrid = vGrid, a_nx = nx, 
                a_nv = nv, a_T1 = T1, a_T2 = T2, a_Q2 = Q2, a_normGrad = normGrad, 
                a_logL = logL, a_epsOpt = epsOpt)
            Q2 <- z$a_Q2
            normGrad <- z$a_normGrad
            logL <- z$a_logL
            epsOpt <- z$a_epsOpt
            if (any(Q2 < 0)) 
                stop("Negative Elements in Q2...\a")
            if (display) {
                dev.set(which = 4)
                points(normGrad[1], normGrad[2], col = "red")
            }
            if (verbose) 
                cat("D: ", normGrad[1], " ", normGrad[2], " |  logLik =", 
                  logL, " ")
            cond <- (max(abs(normGrad)) < tol)
            return(list(Q2 = Q2, cond = cond, normGrad = normGrad, 
                logL = logL, epsOpt = epsOpt))
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
