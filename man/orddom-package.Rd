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
Package: \tab anRpackage\cr
Type: \tab Package\cr
Version: \tab 1.5\cr
Date: \tab 2012-01-02\cr
License: \tab GPL-2\cr
}
}
\author{
Jens J. Rogmann, University of Hamburg, Dept of Psychology, Germany
Maintainer: Jens J. Rogmann <Jens.Rogmann@uni-hamburg.de>
}
\references{
Cliff, N. (1993). Dominance statistics: Ordinal analyses to answer ordinal questions. \emph{Psychological Bulletin, 114}, 494-509.\cr \cr
Cliff, N. (1996a). \emph{Ordinal Methods for Behavioral Data Analysis}. Mahwah, NJ: Lawrence Erlbaum.\cr \cr 
Cliff, N. (1996b). Answering ordinal questions with ordinal data using ordinal statistics. \emph{Multivariate Behavioral Research, 31}, 331-350.\cr \cr
Cohen, J. (1988). \emph{Statistical power analysis for the behavioral sciences (2nd edition)}. New York: Academic Press.\cr \cr
Delaney, H.D. & Vargha, A. (2002). Comparing Several Robust Tests of Stochastic Equality With Ordinally Scaled Variables and Small to Moderate Sized Samples. \emph{Psychological Methods, 7}, 485-503.\cr \cr
Feng, D., & Cliff, N. (2004). Monte Carlo Evaluation of Ordinal d with Improved Confidence Interval. \emph{Journal of Modern Applied Statistical Methods, 3}(2), 322-332.\cr \cr
Feng, D. (2007). Robustness and Power of Ordinal d for Paired Data. In Shlomo S. Sawilowsky (Ed.), \emph{Real Data Analysis} (pp. 163-183). Greenwich, CT : Information Age Publishing.\cr \cr
Grissom, R.J. (1994). Probability of the superior outcome of one treatment over another. \emph{Journal of Applied Psychology, 79}, 314-316.\cr \cr
Grissom, R.J. & Kim, J.J. (2005). \emph{Effect sizes for research. A broad practical approach}. Mahwah, NJ, USA: Erlbaum.\cr \cr
McGraw, K.O. & Wong, S.P. (1992). A common language effect size statistic. \emph{Psychological Bulletin, 111}, 361-365.\cr \cr
Long, J. D., Feng, D., & Cliff, N. (2003). Ordinal analysis of behavioral data. In J. Schinka & W. F. Velicer (eds.), \emph{Research Methods in Psychology. Volume 2 of Handbook of Psychology} (I. B. Weiner, Editor-in-Chief). New York: John Wiley & Sons.\cr \cr
Romano, J., Kromrey, J. D., Coraggio, J., & Skowronek, J. (2006). \emph{Appropriate statistics for ordinal level data: Should we really be using t-test and Cohen's d for evaluating group differences on the NSSE and other surveys?}. Paper presented at the annual meeting of the Florida Association of Institutional Research, Feb. 1-3, 2006, Cocoa Beach, Florida. Last retrieved January 2, 2012 from www.florida-air.org/romano06.pdf \cr \cr
Vargha, A., & Delaney, H. D. (1998). The Kruskal-Wallis test and stochastic homogeneity. \emph{Journal of Educational and Behavioral Statistics, 23}, 170-192.\cr \cr
Vargha, A., & Delaney, H. D. (2000). A critique and improvement of the CL common language effect size statistic of McGraw and Wong. \emph{Journal of Educational and Behavioral Statistics, 25}, 101-132.\cr \cr
}
\examples{\dontrun{
#orddom (x,y,paired=FALSE) provides delta statistics for independent groups x and y
#orddom (x,y,paired=TRUE) provides delta statistics for paired data (e.g. x:Pre and y:Post)
}}
\keyword{nonparametric}
\keyword{robust}
\keyword{htest}
