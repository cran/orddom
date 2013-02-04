\name{orddom-package}
\alias{orddom-package}
\alias{cliff's delta}
\docType{package}
\title{Ordinal Dominance Statistics}
\description{
This package provides ordinal, nonparametric statistics and effect sizes as an alternative to independent or paired group mean comparisons, with special reference to Cliff's delta statistics (or success rate difference, SRD), but also providing McGraw and Wong's common language effect size for the discrete case (i.e. Grissom and Kim's Probability of Superiority), Vargha and Delaney's A (or the Area Under a Receiver Operating Characteristic Curve AUC), and Cook & Sackett's number needed to treat (NNT) effect size (cf. Kraemer & Kupfer, 2006). For the nonparametric effect sizes, various bootstrap CI estimates may also be obtained. Nonparametric effect sizes are also expressed as Cohen's d based on percentages of group non-overlap (cf. Huberty & Lowman, 2000).
}
\details{
\tabular{ll}{
Package: \tab anRpackage\cr
Type: \tab Package\cr
Version: \tab 3.1\cr
Date: \tab 2013-02-07\cr
License: \tab GPL-2\cr
}}
\seealso{\code{\link{orddom}}, \code{\link{dmes}}, \code{\link{dmes.boot}} and \code{\link{orddom_f}}.}
\author{
Dr. Jens J. Rogmann, University of Hamburg, Dept of Psychology, Germany
Maintainer: Jens J. Rogmann <Jens.Rogmann@uni-hamburg.de>
}
\note{
Please cite as: \cr \cr
Rogmann, J. J. (2013). \emph{Ordinal Dominance Statistics (orddom): An R Project for Statistical Computing package to compute ordinal, nonparametric alternatives to mean comparison} (Version 3.1). Available online from the CRAN website \url{http://cran.r-project.org/}. \cr \cr
\emph{Major changes from orddom version 3.0 to 3.1}
\itemize{
 \item{Correction for dmes list names}
 }
\emph{Major changes from orddom version 2.0 to 3.0}
\itemize{
 \item{New function dmes to easily calculate nonparametric effect size measures independently from orddom}
 \item{Easier and more reliable input possible (vectors, lists, arrays, data frames) (by means of new function return1colmatrix)}
 \item{Individual variable label and test descriptions can now be assigned}
 \item{Outputs now also contain Number Needed to Treat (NNT) effect size}
 \item{(Para)metric Common Language effect size McGraw & Wong (1992) added}
 \item{Metric d CI in orddom now based on Hedges & Olkin (1985)}
 \item{Elimination of negative delta-between variance estimates for paired comparisons}
 \item{Correction of symmetric CI for independent Cliff's delta statistics}
 \item{New function dmes.boot to calculate bootstrap-based CI for nonparametric effect size measures and Cohen's d,}
 \item{dmes.boot was largely based on R code provided by J. Ruscio and T. Mullen (2011) reused with kind permission}
 \item{New function delta_gr now yields a graphical and interpretational output for Cliff's delta statistic}
 \item{New options for one- and two-tailed CI in orddom and orddom_f, resulting in changes of rows 21 and 22 of independent and rows 18 and 19 of paired orddom result matrix}
 \item{New Metric_t function (t, p and df can now be calculated and embedded in orddom as standard or Welch approximated) \cr \cr}}
\emph{Major changes from orddom Version 1.5 to 2.0}
\itemize{
\item{orddom now also accepts simple vectors as x or y.}
\item{New orddom_f() function file allows for file output of statistics for multiple sample comparisons (e.g. csv or analyses in MS Excel or Open Office Calc).}
\item{New orddom_p() function file allows for detailed tab-formatted output for single sample comparisons.}
\item{Package dependency on compute.es package was suspended (tes-Function for metric Cohen's d in orddom).}
\item{New metric_t() function for additional information on metric t-test results.}
\item{The dm() function can now also return difference matrices.}
\item{Improved stability of orddom function as well as minor corrections in orddom output and manuals. \cr \cr }}
\emph{Major changes from orddom Version 1.0 to 1.5}
\itemize{
\item{Calculation of CI and delta z-score-estimates can now be based on Students t-distribution rather than using fixed normal distribution z-scores.}
\item{Symmetric CIs can now be obtained to increase power of the delta statistics in certain cases.}
\item{Formulas used for calculation added in orddom manual.}
\item{Probability of Superiority statistic as well as variance estimates for delta in the independent groups analyses were corrected.}
\item{Minor changes were implemented to allow for calculation of \eqn{d = \pm 1}{d=+-1} extreme cases without error abort.}
\item{Output of raw y-dataset in independent group analysis was corrected.}
\item{Dependencies on packages \emph{psych} and \emph{compute.es} declared in DESCRIPTION and NAMESPACE files.}
}}
\references{
Cliff, N. (1993). Dominance statistics: Ordinal analyses to answer ordinal questions. \emph{Psychological Bulletin, 114}, 494-509.\cr \cr
Cliff, N. (1996a). \emph{Ordinal Methods for Behavioral Data Analysis}. Mahwah, NJ: Lawrence Erlbaum.\cr \cr 
Cliff, N. (1996b). Answering ordinal questions with ordinal data using ordinal statistics. \emph{Multivariate Behavioral Research, 31}, 331-350.\cr \cr
Cohen, J. (1988). \emph{Statistical power analysis for the behavioral sciences (2nd edition)}. New York: Academic Press.\cr \cr
Cook, R.J. & Sackett, D. L. (1995). The number needed to treat: A clinically useful measure of treatment effect. \emph{British Medical Journal, 310}, 452-454. \cr \cr
Delaney, H.D. & Vargha, A. (2002). Comparing Several Robust Tests of Stochastic Equality With Ordinally Scaled Variables and Small to Moderate Sized Samples. \emph{Psychological Methods, 7}, 485-503.\cr \cr
Feng, D., & Cliff, N. (2004). Monte Carlo Evaluation of Ordinal d with Improved Confidence Interval. \emph{Journal of Modern Applied Statistical Methods, 3}(2), 322-332.\cr \cr
Feng, D. (2007). Robustness and Power of Ordinal d for Paired Data. In Shlomo S. Sawilowsky (Ed.), \emph{Real Data Analysis} (pp. 163-183). Greenwich, CT : Information Age Publishing.\cr \cr
Grissom, R.J. (1994). Probability of the superior outcome of one treatment over another. \emph{Journal of Applied Psychology, 79}, 314-316.\cr \cr
Grissom, R.J. & Kim, J.J. (2005). \emph{Effect sizes for research. A broad practical approach}. Mahwah, NJ, USA: Erlbaum.\cr \cr
Huberty, C. J. & Lowman, L. L. (2000). Group overlap as a basis for effect size. \emph{Educational and Psychological Measurement, 60}, 543-563.\cr\cr
Kraemer, H.C. &  Kupfer, D.J. (2006). Size of Treatment Effects and Their Importance to Clinical Research and Practice. \emph{Biological Psychiatry, 59}, 990-996. \cr \cr
McGraw, K.O. & Wong, S.P. (1992). A common language effect size statistic. \emph{Psychological Bulletin, 111}, 361-365.\cr \cr
Long, J. D., Feng, D., & Cliff, N. (2003). Ordinal analysis of behavioral data. In J. Schinka & W. F. Velicer (eds.), \emph{Research Methods in Psychology. Volume 2 of Handbook of Psychology} (I. B. Weiner, Editor-in-Chief). New York: John Wiley & Sons.\cr \cr
Romano, J., Kromrey, J. D., Coraggio, J., & Skowronek, J. (2006). \emph{Appropriate statistics for ordinal level data: Should we really be using t-test and Cohen's d for evaluating group differences on the NSSE and other surveys?}. Paper presented at the annual meeting of the Florida Association of Institutional Research, Feb. 1-3, 2006, Cocoa Beach, Florida. Last retrieved January 2, 2012, from www.florida-air.org/romano06.pdf \cr \cr
Ruscio, J. & Mullen, T. (2012). Confidence Intervals for the Probability of Superiority Effect Size Measure and the Area Under a Receiver Operating Characteristic Curve. \emph{Multivariate Behavioral Research, 47}, 221-223.
Vargha, A., & Delaney, H. D. (1998). The Kruskal-Wallis test and stochastic homogeneity. \emph{Journal of Educational and Behavioral Statistics, 23}, 170-192.\cr \cr
Vargha, A., & Delaney, H. D. (2000). A critique and improvement of the CL common language effect size statistic of McGraw and Wong. \emph{Journal of Educational and Behavioral Statistics, 25}, 101-132.\cr \cr
}
\examples{\dontrun{
#ordinal comparison and delta statistics for independent groups x and y 
#(e.g. x:comparison/control group and y:treatment/experimental group)
orddom(x,y,paired=FALSE) 
#
#ordinal comparison and delta statistics for paired data 
#(e.g. x:Pretest/Baseline and y:Posttest)
orddom(x,y,paired=TRUE)  
#
#Dominance Matrix production
dms(x,y,paired=T)
#
#Print dominance matrix 
orddom_p(x,y,sections="4a")
#
#Graphic output and interpretational text for Cliff's delta statistics
delta_gr(x,y)
#
#nonparametric effect sizes (SRD/delta, A/AUC, CL/PS, NNT)
#(e.g. C:control group scores, T:treatment group scores)
dmes(C,T)
#
#Convert Cliff's delta value to Cohen's d (as distributional non-overlap)
delta2cohd(dmes(C,T)$dc)
#
#Confidence Interval estimate of AUC (by bootstrap)
#cf. Ruscio, J. & Mullen, T. (2012)
#(e.g. C:control group scores, T:treatment group scores)
dmes.boot(C,T,theta.es="Ac")
}}
\keyword{nonparametric}
\keyword{robust}
\keyword{htest}
