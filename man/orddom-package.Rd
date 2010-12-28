\name{orddom-package}
\alias{orddom-package}
\alias{cliff's delta}
\docType{package}
\title{Ordinal Dominance Statistics}
\description{
This package provides ordinal statistics as an alternative to independent or paired group mean comparisons, especially for Cliff's delta statistics.
}
\details{
\tabular{ll}{
Package: \tab orddpm-package\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2010-12-29\cr
License: \tab GPL-2\cr
}
}
\author{
Jens J. Rogmann, University of Hamburg, Dept of Psychology, Germany
Maintainer: Jens J. Rogmann <Jens.Rogmann@uni-hamburg.de>
}
\references{
Cliff, N. (1993). Dominance statistics: Ordinal analyses to answer ordinal questions. \emph{Psychological Bulletin, 114}, 494-509.\cr
Cliff, N. (1996a). \emph{Ordinal Methods for Behavioral Data Analysis}. Mahwah, NJ: Lawrence Erlbaum.\cr
Cliff, N. (1996b). Answering ordinal questions with ordinal data using ordinal statistics. \emph{Multivariate Behavioral Research, 31}, 331-350.\cr
Long, J. D., Feng, D., & Cliff, N. (2003). Ordinal analysis of behavioral data. In J. Schinka & W. F. Velicer (eds.), \emph{Research Methods in Psychology. Volume 2 of Handbook of Psychology} (I. B. Weiner, Editor-in-Chief). New York: John Wiley & Sons.\cr
Feng, D., & Cliff, N. (2004). Monte Carlo Evaluation of Ordinal d with Improved Confidence Interval. \emph{Journal of Modern Applied Statistical Methods, 3}(2), 322-332.\cr
Feng, D. (2007). Robustness and Power of Ordinal d for Paired Data. In Shlomo S. Sawilowsky (Ed.), \emph{Real Data Analysis} (pp. 163-183). Greenwich, CT : Information Age Publishing.\cr
}
\examples{\dontrun{orddom (x,y,paired=FALSE) provides delta statistics for independent groups x and y\cr
orddom (x,y,paired=TRUE) provides delta statistics for paired data (e.g. x:Pre and y:Post)}}
\keyword{nonparametric}
\keyword{robust}
\keyword{htest}
