\name{lzrq}
\alias{lzrq}
\alias{print.lzrq}
\alias{summary.lzrq}
\title{
Quantile regression for logarithmic relationships with non-positive outcome values
}
\description{
\code{lzrq} fits a quantile regression model where the raw dependent variable
is transformed with a calibrated extensive margin (CEM) transformation of the
form

\deqn{y^* = \ln(y) \quad \text{if } y > 0}
\deqn{y^* = -\psi \quad \text{if } y \leq 0}

\code{lzrq} is a wrapper for \code{\link[quantreg]{rq}}. If all fitted values
from a quantile regression on the CEM-transformed outcome are greater than
\eqn{-\psi}, results are returned. The resulting coefficients can be
interpreted as logarithmic intensive-margin relationships between the outcome
and the independent variables, even when the outcome contains non-positive
values. If the condition does not hold, \code{lzrq} iteratively increases
\eqn{\psi} and re-checks. After ten failed iterations, \code{lzrq} returns an
error and suppresses results. This is an automated adaptation of the algorithm
described by Liu & Kaplan (2025).

All methods supported by \code{\link[quantreg]{rq}} are available, including
\code{coef()}, \code{residuals()}, \code{fitted()}, \code{predict()},
\code{summary()}, \code{AIC()}, and \code{logLik()}.
}
\usage{
lzrq(formula, data, tau = 0.5, floor_val = -1e35, ...)
}
\arguments{
  \item{formula}{
A formula with the raw (untransformed) outcome on the left-hand side.
  }
  \item{data}{
Optional data frame containing the variables in the model.
  }
  \item{tau}{
Quantile level. Numeric scalar strictly between 0 and 1. Defaults to 0.5.
  }
  \item{psi_init}{
Initial lower bound constant used in the bisection algorithm.
Defaults to \code{-1e35} to mirror Stata's
\code{lzqreg} convention.
  }
  \item{...}{
Additional arguments passed to \code{\link[quantreg]{rq}}, such as
\code{weights} or \code{method}.
  }
}
\value{
An object of class \code{c("lzrq", "rq")}. This is the fitted \code{rq}
object with two additional fields:

\item{lb_constant}{The lower bound constant \eqn{-\psi} used to replace non-positive
outcomes in the final iteration.}
\item{n_nonpos}{Number of non-positive outcome values in the estimation
sample.}

All \code{rq} fields are preserved, so the full suite of
\code{\link[quantreg]{rq}} postestimation methods works automatically.
}
\examples{
library(MASS)

# The 'epil' dataset records epileptic seizure counts.
# The outcome 'y' is zero for about 10% of observations.
table(epil$y == 0)

# Because the median of 'y', conditional on 'trt', is above zero for
# both treatment groups, lzrq returns results at quantile 0.5.
result <- lzrq(y ~ trt + age, data = epil, tau = 0.5)
print(result)
summary(result)

# However, the tenth percentile of 'y' is zero for some treatment groups,
# so lzrq suppresses results and returns an error at quantile 0.1.
tryCatch(
  lzrq(y ~ trt + age, data = epil, tau = 0.1),
  error = function(e) message(e$message)
)
}
\references{
Fitzgerald, J., Adema, J., Fiala, L., Kujansuu, E., & Valenta, D. (2026).
Non-Robustness in Log-Like Specifications. MetaArXiv.
\doi{10.31222/osf.io/juda7_v1}

Liu, X., & Kaplan, D. M. (2025). Quantile Regression with Log(0) Outcomes.
\url{https://drive.google.com/file/d/1F3dnhm8MrlO5aRrGt48rBWAEaBqdCBH-/view}
}
\author{
David Valenta (dvalenta@uottawa.ca) and Jack Fitzgerald (jackfitzgeraldresearch@gmail.com)
}
\seealso{
\code{\link[quantreg]{rq}}, \code{\link[quantreg]{summary.rq}}
}