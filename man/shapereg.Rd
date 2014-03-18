\name{shapereg}
\alias{shapereg}
\title{Shape-Restricted Regression}
\description{The regression model \eqn{y_i = f(t_i) + x_i'\beta + \epsilon_i, i = 1,\ldots,n} is considered, where the only assumptions about \eqn{f} concern its shape. The vector expression for the model is \eqn{y = \theta + X\beta + \epsilon}. \eqn{X} represents a categorical covariate. The shapereg function allows eight shapes: increasing, decreasing, convex, concave, increasing-convex, increasing-concave, decreasing-convex, and decreasing-concave. This routine employs a single cone projection to find \eqn{\theta} and \eqn{\beta} simultaneously.  }

\usage{shapereg(y, t, shape, xmat = NULL, w = NULL, test = FALSE)}

\arguments{
  \item{y}{A vector of length \eqn{n}.}

  \item{t}{A continuous and constrained predictor of length \eqn{n}.}

  \item{shape}{A value ranging from 1 to 8 indicating the shape of \eqn{f}: 1 = increasing; 2 = decreasing; 3 = convex; 4 = concave; 5 = increasing convex; 6 = decreasing convex; 7 = increasing concave; 8 = decreasing concave. }

  \item{xmat}{A \eqn{n} by \eqn{k} full colunm rank matrix. If xmat is given as a vector, it will be transformed into a \eqn{n} by \eqn{1} matrix. It represents a categorical covariate. The user can choose to put a constant vector in xmat or not. If xmat is not given, xmat is created as a \eqn{n} by \eqn{1} matrix with all elements equal to 1. The default is xmat = NULL.}

  \item{w}{An optional nonnegative vector of weights of length \eqn{n}. If w is not given, all weights are taken to equal 1. Otherwise, the minimization of \eqn{(y - \phi)'w(y - \phi)} over \eqn{C} is returned. The default is w = NULL.} 
  \item{test}{A logical scalar. If test == TRUE, then the p-value for the hypothesis test \eqn{H_0: \phi} is in \eqn{V} versus \eqn{H_1: \phi} is in \eqn{C} is returned. \eqn{C} is the constraint cone, and \eqn{V} is the linear space contained in the constraint cone.}   
}

\details{
This routine constrains \eqn{\theta} in the equation \eqn{y = \theta + X\beta + \epsilon} by a shape parameter. 

The constraint cone \eqn{C} has the form \eqn{\{\phi: \phi = v + \sum b_i\delta_i, i = 1,\ldots,m, b_1,\ldots, b_m \ge 0 \}}, \eqn{v} is in \eqn{V}. The column vectors of \eqn{X} are in \eqn{V}, i.e., the linear space contained in the constraint cone. 

The hypothesis test \eqn{H_0: \phi} is in \eqn{V} versus \eqn{H_1: \phi} is in \eqn{C} is an exact one-sided test, and the test statistic is \eqn{E_{01} = (SSE_0 - SSE_1)/(SSE_0)}, which has a mixture-of-betas distribution when \eqn{H_0} is true and \eqn{\epsilon} is a vector following a standard multivariate normal distribution with mean 0. The mixing parameters are found through simulations. The number of simulations used to obtain the mixing distribution parameters for the test is 10,000. Such simulations usually take some time. For the "feet" data set used as an example in this section, whose sample size is 39, the time to get a p-value is roughly between 4 seconds. 

This routine calls coneB for the cone projection part. 
}

\value{
\item{pval}{The p-value for the hypothesis test \eqn{H_0: \phi} is in \eqn{V} versus \eqn{H_1: \phi} is in \eqn{C}. \eqn{C} is the constraint cone of the form \eqn{\{\phi: \phi = v + \sum b_i\delta_i, i = 1,\ldots,m, b_1,\ldots, b_m \ge 0 \}}, \eqn{v} is in \eqn{V}, and \eqn{V} is the linear space contained in the constraint cone. If test == TRUE, a p-value is returned. Otherwise, the test is skipped and no p-value is returned.}

\item{coefs}{The estimated coefficients for \eqn{X}, i.e., the estimation for the vector \eqn{\beta}. Note that even if the user does not provide a constant vector in \eqn{X}, the coefficient for the intercept will be returned.}

\item{constr.fit}{The shape-restricted fit over the constraint cone \eqn{C} of the form \eqn{\{\phi: \phi = v + \sum b_i\delta_i, 
i = 1,\ldots,m, b_1,\ldots, b_m \ge 0 \}}, \eqn{v} is in \eqn{V}.}

\item{linear.fit}{The least-squares regression of \eqn{y} on \eqn{V}, i.e., the linear space contained in the constraint cone. If shape is 3 or shape is 4, \eqn{V} is spanned by \eqn{X} and \eqn{t}. Otherwise, it is spanned by \eqn{X}. \eqn{X} must be full column rank, and the matrix formed by combining \eqn{X} and \eqn{t} must also be full column rank.}

\item{se.beta}{The standard errors for the estimation of the vector \eqn{\beta}. The degree of freedom is returned by coneB and is multiplied by 1.5. Note that even if the user does not provide a constant vector in \eqn{X}, the standard error for the intercept will be returned.}

\item{pvals.beta}{The approximate p-values for the estimation of the vector \eqn{\beta}. A t-distribution is used as the approximate distribution. Note that even if the user does not provide a constant vector in \eqn{X}, the approximate p-value for the intercept will be returned.}

\item{shape}{The shape given by the user.}

\item{test}{The test parameter given by the user.}

\item{SSE0}{The sum of squared residuals for the linear part.}

\item{SSE1}{The sum of squared residuals for the full model.}

\item{call}{The matched call.}
}

\references{Raubertas, R. F., C.-I. C. Lee, and E. V. Nordheim (1986) Hypothesis tests for normals
means constrained by linear inequalities. \emph{Communications in Statistics - Theory and 
Methods \bold{15 (9)}}, 2809--2833.

Robertson, T., F. Wright, and R. Dykstra (1988) \emph{Order Restricted Statistical Inference}
New York: John Wiley and Sons.

Fraser, D. A. S. and H. Massam (1989) A mixed primal-dual bases algorithm for regression
under inequality constraints application to concave regression. \emph{Scandinavian Journal of Statistics \bold{16}}, 65--74.

Meyer, M. C. (2003) A test for linear vs convex regression function using shape-restricted
regression. \emph{Biometrika \bold{90 (1)}}, 223--232.

Cheng, G.(2009) Semiparametric additive isotonic regression. \emph{Journal of Statistical Planning and Inference \bold{139}}, 1980--1991.

Meyer, M.C.(2013a) Semiparametric additive constrained regression. \emph{Journal of Nonparametric Statistics to appear}.  }

\author{Mary C. Meyer and Xiyue Liao}

\seealso{
\code{\link{coneB}}
}

\examples{
#load the feet data set
    data(feet)

#extract the continuous and constrained predictor
    l <- feet$length

#extract the continuous response
    w <- feet$width

#extract the categorical covariate
    s <- feet$sex

#create the dummy variable indicating sex: if sex = "G", x = 0; if sex = "B", x = 1 
    n <- length(s)
    x <- as.numeric(s != "G")

#create the xmat in the shapereg function
    xmat <- x

#make an increasing fit with test set as FALSE
    shape <- 1
    ans <- shapereg(w, l, shape, xmat)

#check the summary table 
    summary(ans)

#make an increasing fit with test set as TRUE
    ans <- shapereg(w, l, shape, xmat, test = TRUE)

#check the summary table 
    summary(ans)

#make a plot comparing the unconstrained fit and the constrained fit
    par(mar = c(4, 4, 1, 1))
    ord <- order(l)
    plot(sort(l), w[ord], type = "n", xlab = "foot length (cm)", ylab = "foot width (cm)")
    title("Shapereg Example Plot")

#sort l according to sex
    ord1 <- order(l[s == "G"])
    ord2 <- order(l[s == "B"])

#make the scatterplot of l vs w for boys and girls
    points(sort(l[s == "G"]), w[s == "G"][ord1], pch = 21, col = 1)
    points(sort(l[s == "B"]), w[s == "B"][ord2], pch = 24, col = 2)

#make an unconstrained fit to boys and girls
    fit <- lm(w ~ l + factor(s))

#plot the unconstrained fit 
    lines(sort(l), (coef(fit)[1] + coef(fit)[2] * l + coef(fit)[3])[ord], lty = 2)
    lines(sort(l), (coef(fit)[1] + coef(fit)[2] * l)[ord], lty = 2, col = 2)
    legend(21.5, 9.8, c("boy","girl"), pch = c(24, 21), col = c(2, 1)) 

#plot the constrained fit
    lines(sort(l), (ans$constr.fit - ans$linear.fit + ans$coefs[1])[ord], col = 1)
    lines(sort(l), (ans$constr.fit - ans$linear.fit + ans$coefs[1] + ans$coefs[2])[ord], col = 2)
}
\keyword{cone projection routines}

