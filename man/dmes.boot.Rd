\name{dmes.boot}
\alias{dmes.boot}
\title{Dominance Matrix Effect Sizes}
\description{Bootstrap-based calculation of standard error and CI constructs for Cohen's d and the statistics used in the Dominance Matrix Effect Size (dmes) function}
\usage{dmes.boot(x,y,theta.es="dc",ci.meth="BCA",B=1999,alpha=.05,seed=1)}
\arguments{
  \item{x}{A vector or 1 column matrix with \eqn{n_x} values from (control or pre-test or comparison) group X}
  \item{y}{A vector or 1 column matrix with \eqn{n_y} values from (treatment or post-test) group Y}
  \item{theta.es}{Specification of the nonparametric effect size for which the SE and CI is to be constructed. All output values of the \link{dmes} function can be used, e.g. "PSc", "Ac", "dc", "NNTc", "PSw", "Aw", "dw", "NNTw", "PSb", "Ab", "db"  or "NNTb".}
  \item{ci.meth}{Specify type of method used for bootstrap confidence interval construction: "BSE", "BP" or "BCA". \cr \cr
\emph{"BSE"} uses the bootstrap standard error estimate of the respective nonparametric effect size to construct a confidence interval with \eqn{\hat{\theta} \pm z_{\alpha/2} \cdot \widehat{SE_{\theta}}}{theta +- z(alpha/2) * SE(theta)}, where \eqn{\hat{\theta}} ist the observed effect size, \eqn{z_{\alpha/2}}{z(alpha/2)} the z value of the standard normal table at the given (two-tailed) significance level (e.g. z=1.96 when alpha=5%) and \eqn{\widehat{SE_{\theta}}}{SE(theta)} is the standard error estimate of the respective nonparametric effect size, based upon B bootstrap sample estimates of theta. \cr \cr
\emph{"BP"} calculates confidence intervals based on bootstrap percentiles. B bootstrap sample estimates of the respective nonparametric effect size \eqn{\theta}{theta} are generated and ordered, and the \eqn{(B \cdot 100 \cdot \alpha)}{(B*100*alpha)}th as well as the \eqn{(B \cdot 100 \cdot (1-\alpha)}{(B*100*(1-alpha))}th of these ordered estimates are used to determine the confidence intervals. For example, if B=2000 bootstrap samples are calculated and \eqn{\alpha=.05}{alpha=.05}, then the 100th and 1900th of the ordered values are selected as lower and upper CI limits.\cr \cr
\emph{"BCA"} calculates bias-corrected and accelerated confidence intervals (also based on bootstrap percentiles). Here, however, the \eqn{\alpha} levels (or percentiles) are corrected depending on the bias and the rate of change of the standard error with formulas suggested by Efron & Tibshirani (1993, Chapter 14).}
   \item{B}{Number of bootstrap samples to be used for the estimates.}
   \item{alpha}{Significance level.}
   \item{seed}{Integer argument to set random number generation seeds, see \link{Random}.}
}
\details{Returns an associative list with the following values:}
\value{
\item{\emph{$theta}}{Type and observed value of the respective nonparametric effect size estimate for samples Y and X.} 
\item{\emph{$theta.SE}}{The bootstrap-based estimated standard error of the respective nonparametric effect size estimate.}
\item{\emph{$bci.meth}}{String indicating which type of bootstrap (BSE, BP or BCA) was used to construct the confidence interval for the respective nonparametric effect size estimate and Cohen's d.}
\item{\emph{$theta.bci.lo}}{Lower end of the confidence interval for the respective nonparametric effect size estimate as determined by type of bootstrap used (BSE, BP or BCA).}
\item{\emph{$theta.bci.up}}{Upper end of the confidence interval for the respective nonparametric effect size estimate as determined by type of bootstrap used (BSE, BP or BCA).}
\item{\emph{$Coh.d}}{Effect size estimate of Cohen's d based on student's t and assuming pooled variance. For details, see \link{metric_t}.}
\item{\emph{$Coh.d.bSE}}{The bootstrap-based estimated standard error of Cohen's d.}
\item{\emph{$Coh.d.bci.lo}}{Lower end of the confidence interval for the Cohen's d estimated through bootstrapping (type BSE, BP or BCA).}
\item{\emph{$Coh.d.bci.up}}{Upper end of the confidence interval for the Cohen's d estimated through bootstrapping (type BSE, BP or BCA).}
}
\references{
Efron, B. & Tibshirani (1993). \emph{An Introduction to the Bootstrap}. New York/London: Chapman & Hall. \cr \cr
Ruscio, J. & Mullen, T. (2011). \emph{Bootstrap CI for A} (R program code, last updated April 11,2011). Retrieved from \url{http://www.tcnj.edu/~ruscio/Bootstrap\%20CI\%20for\%20A.R} . \cr \cr
Ruscio, J. & Mullen, T. (2012). Confidence Intervals for the Probability of Superiority Effect Size Measure and the Area Under a Receiver Operating Characteristic Curve. \emph{Multivariate Behavioral Research, 47}, 221-223. \cr
}
\author{Jens J. Rogmann}
\seealso{dmes}
\note{\emph{dmes.boot} was largely based on R code provided by John Ruscio and Tara Mullen (2011) which was reused with kind permission from the authors.
}
\examples{\dontrun{
> # cf. Efron & Tibshirani (1993, Ch. 14) 
> # Spatial Test Data (Table 14.1, p.180)
> A<-c(48,36,20,29,42,42,20,42,22,41,45,14,6,0,33,28,34,4,32,24,47,41,24,26,30,41)
> B<-c(42,33,16,39,38,36,15,33,20,43,34,22,7,15,34,29,41,13,38,25,27,41,28,14,28,40)
> dmes.boot(A,B)
$theta
         dc 
-0.08136095 

$theta.SE
[1] 0.1656658

$bci.meth
[1] "BCA"

$theta.bci.lo
[1] -0.4008876

$theta.bci.up
[1] 0.2440828

$Coh.d
[1] -0.06364221

$Coh.d.bSE
[1] 0.2895718

$Coh.d.bci.lo
[1] -0.6106167

$Coh.d.bci.up
[1] 0.5031792
}
\dontrun{
> ############################################################################
> #Example from Ruscio & Mullen (2012, p. 202)
> x <- c(6,7,8,7,9,6,5,4,7,8,7,6,9,5,4) # Treatment Group
> y <- c(4,3,5,3,6,2,2,1,6,7,4,3,2,4,3) # Control Group
> dmes.boot(y,x,theta.es="Ac")          #AUC
$theta
       Ac 
0.8844444 

$theta.SE
[1] 0.05910963

$bci.meth
[1] "BCA"

$theta.bci.lo
[1] 0.7022222

$theta.bci.up
[1] 0.9644444

$Coh.d
[1] 1.727917

$Coh.d.bSE
[1] 0.4932543

$Coh.d.bci.lo
[1] 0.7753663

$Coh.d.bci.up
[1] 2.573305
}}
\keyword{array}
\keyword{nonparametric}
\keyword{robust}
