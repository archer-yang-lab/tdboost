\name{relative.influence}
\alias{relative.influence}
\alias{permutation.test.NPtweedie}
\alias{NPtweedie.loss}
\title{ Methods for estimating relative influence }
\description{
Helper functions for computing the relative influence of each variable in the NPtweedie object.
}
\usage{
relative.influence(object, n.trees)
permutation.test.NPtweedie(object, n.trees)
NPtweedie.loss(y,f,w,offset,dist,baseline)
}
\arguments{
\item{object}{a \code{NPtweedie} object created from an initial call to \code{\link{NPtweedie}}.}
\item{n.trees}{ the number of trees to use for computations.}
\item{y,f,w,offset,dist,baseline}{For \code{NPtweedie.loss}: These components are the
outcome, predicted value, observation weight, offset, distribution, and comparison
loss function, respectively.}
}
\details{
This is not intended for end-user use. These functions offer the different
methods for computing the relative influence in \code{\link{summary.NPtweedie}}.
\code{NPtweedie.loss} is a helper function for \code{permutation.test.NPtweedie}.
}
\value{
Returns an unprocessed vector of estimated relative influences.
}
\references{
Yang, Y. and Zou, H. (2013), \dQuote{Nonparametric Multiple Expectile Regression via ER-Boost,} \emph{Journal of Statistical Computation and Simulation}. Accepted.

BugReport: \url{http://code.google.com/p/NPtweedie/}\cr

G. Ridgeway (1999). \dQuote{The state of boosting,} \emph{Computing Science and
Statistics} 31:172-181.

\url{http://cran.r-project.org/web/packages/gbm/}\cr


J.H. Friedman (2001). "Greedy Function Approximation: A Gradient Boosting
Machine," Annals of Statistics 29(5):1189-1232.

L. Breiman (2001). "Random Forests," Available at \url{ftp://ftp.stat.berkeley.edu/pub/users/breiman/randomforest2001.pdf}.

}
\author{Yi Yang \email{yiyang@umn.edu} and Hui Zou \email{hzou@stat.umn.edu}}

\seealso{ \code{\link{summary.NPtweedie}} }

\keyword{ hplot }
