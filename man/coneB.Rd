\name{coneB}
\alias{coneB}
\title{Cone Projection -- Constraint Cone}
\description{
This routine implements the hinge algorithm for cone projection to minimize \eqn{||y - \theta||^2} over the cone
\eqn{C} of the form \eqn{\{\theta: \theta = v + \sum b_i\delta_i,   i = 1,\ldots,m,  b_1,\ldots, b_m \ge 0\}}, \eqn{v} is in \eqn{V}.}

\usage{
coneB(y, delta, vmat, w = NULL, lst = TRUE)
}

\arguments{
  \item{y}{A vector of length \eqn{n}.}
  \item{delta}{A matrix whose rows are the constraint cone edges. The rows of delta must be irreducible. Its column number must equal the length of \eqn{y}. No row of delta is contained in the column space of vmat.}
  \item{vmat}{A matrix whose columns are the basis of the linear space contained in the constraint cone. Its row number must equal the length of \eqn{y}. The columns of vmat must be linearly independent.}
  \item{w}{An optional nonnegative vector of weights of length \eqn{n}. If w is not given, all weights are taken to equal 1. Otherwise, the minimization of \eqn{(y - \theta)'w(y - \theta)} over \eqn{C} is returned. The default is w = NULL.}
  \item{lst}{A logical scalar. If lst == TRUE, a list of df, yhat, steps, coefs, message and convergence is returned; Otherwise, message and convergence are omitted. The default is lst = TRUE.}  }
\details{The routine coneB dynamically loads a C++ subroutine "coneBCpp". }

\value{
\item{df}{The dimension of the face of the constraint cone on which the projection lands.}
\item{yhat}{The projection of \eqn{y} on the constraint cone.}
\item{steps}{The number of iterations before the algorithm converges.}
\item{coefs}{The coefficients of the constraint cone edges and the basis of the linear space contained in the constraint cone.}
\item{message}{A warning message showing the convergence of the algorithm. If steps > \eqn{n^2}, \eqn{n} is the dimension of \eqn{y}, a warning message is given.}
\item{convergence}{A number showing the convergence of the algorithm. If convergence is 0, the algorithm converges; If convergence is 1, the algorithm fails to converge.}  }

\references{
Meyer, M. C. (1999) An extension of the mixed primal-dual bases algorithm to the case of more constraints than 
dimensions. \emph{Journal of Statistical Planning and Inference \bold{81}}, 13--31.

Meyer, M. C. (2013b) A simple new algorithm for quadratic programming with applications in statistics. \emph{Communications
in Statistics \bold{42(5)}}, 1126--1139.
}

\author{
Mary C. Meyer and Xiyue Liao
}

\seealso{
\code{\link{coneA}}, \code{\link{shapereg}}
}

\examples{
#generate y
    n <- 50
    y <- sort(rnorm(n)) + rnorm(n)

#create the edges of the constraint cone to make y monotonically increasing
    n <- 50
    delta <- matrix(0, n - 1 , n)
    for(i in 1:(n - 1)){delta[i, (i + 1):n] <- 1}

#make the basis of V
    vmat <- matrix(rep(1, n), ncol = 1)

#call coneB
    ans1 <- coneB(y, delta, vmat)
    ans2 <- coneB(y, delta, vmat, w = (1:n)/n)

#make a plot to compare the unweighted fit and weighted fit
    par(mar = c(4, 4, 1, 1))
    plot(y, cex = .7, ylab = "y")
    lines(ans1$yhat, col = 2, lty = 2)
    lines(ans2$yhat, col = 4, lty = 2)
    legend("topleft", c("unweighted fit", "weighted fit"), col = c(2, 4), lty = c(2, 2))
    title("ConeB Example Plot")
}

\keyword{cone projection routines}
