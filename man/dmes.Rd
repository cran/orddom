\name{dmes}
\alias{dmes}
\title{Dominance Matrix Effect Sizes}
\description{Generates simple list of nonparametric ordinal effect size measures such as \cr
-the Probability of Superiority (or discrete case Common Language) effect size, \cr
-the Vargha and Delaney's A (or area under the receiver operating characteristic curve, AUC) \cr
-Cliff's delta (or success rate difference, SRD), and
-the number needed to treat (NNT) effect size (based on Cliff's delta value).}
\usage{dmes(x,y)}
\arguments{
  \item{x}{A vector or 1 column matrix with \eqn{n_x} values from (control or pre-test or comparison) group X}
  \item{y}{A vector or 1 column matrix with \eqn{n_y} values from (treatment or post-test) group Y}}
\details{Based on the dominance matrix created by direct ordinal comparison of values of Y with values of X, an associative list is returned.}
\value{
\item{\emph{$nx}}{Vector or sample size of x, \eqn{n_x}.} 
\item{\emph{$ny}}{Vector or sample size of y, \eqn{n_y}} 
\item{\emph{$PSc}}{Discrete case Common Language CL effect size or Probability of Superiority (PS) of all values of Y over all values of X: \deqn{PS_c(Y>X)=\frac{\#(y_i>x_j)}{n_y n_x}}{PS(Y>X)=\#(y_i>x_j)/(n_y n_x)}, \cr
where \eqn{i=\{1, 2, ... , n_y\}} and \eqn{j=\{1, 2, ... , n_x\}}. See \link{orddom} \emph{PS Y>X} for details.)}
\item{\emph{$Ac}}{Vargha & Delaney's A or Area under the receiver operating characteristics curve (AUC) for all possible comparisons: \deqn{A(Y>X)=[\#(y_i>x_j) + .5(\#(y_i=x_j)] (n_y n_x)^{-1}}{A(Y>X)=(\#(y_i>x_j)+.5(\#(y_i=x_j)))/(n_y n_x)}, \cr 
where \eqn{i=\{1, 2, ... , n_y\}} and \eqn{j=\{1, 2, ... , n_x\}}. See \link{orddom} \emph{A Y>X} for details.)} 
\item{\emph{$dc}}{Success rate difference when comparing all values of Y with all values of X: \deqn{d_c(Y>X)=\frac{\#(y_i>x_j)-\#(y_i<x_j)}{n_y n_x}}{delta(Y>X)=(\#(y_i>x_j)-\#(y_i<x_j))/(n_y n_x)}, \cr
where \eqn{i=\{1, 2, ... , n_y\}} and \eqn{j=\{1, 2, ... , n_x\}}. See \link{orddom} \emph{Cliff's delta} for independent groups for details. \cr
\emph{Note that in the paired samples case with \eqn{n_y=n_x}, $dc does not return the combined estimate, i.e. \eqn{\$dc \neq \$dw+\$db}{$dc<>dw+db}!}}
\item{\emph{$NNTc}}{Number needed to treat, based on the success rate difference or \eqn{\$dc^{-1}}. See \link{orddom} "NNT" for details.} 
\item{\emph{$PSw}}{When sample sizes are equal, this value returns the Probability of Superiority (PS) for within-changes, i.e. alle paired values: \eqn{PS_c(Y>X)=\frac{\#(y_i>x_i)}{n_y n_x}}{PS(Y>X)=\#(y_i>x_i)/(n_y n_x)}, limited to the \eqn{n_x=n_y} paired cases where \eqn{i=\{1,2,...,n_x=n_y\}}. (For unequal sample sizes, this equals $PSc.)} 
\item{\emph{$Aw}}{When sample sizes are equal, this value returns A for the paired subsample values, i.e. limited to the \eqn{n_x=n_y} paired cases where \eqn{i=j=\{1,2,...,n_x=n_y\}}. (For unequal sample sizes, this equals $Ac.)} 
\item{\emph{$dw}}{When \eqn{n_x=n_y}, this value returns \emph{Cliff's delta-within}, i.e. paired comparisons limited to the diagonal of the dominance matrix or those cases where \eqn{i=j}. (For unequal sample sizes, this equals $dc.)} 
\item{\emph{$NNTw}}{Number needed to treat, based on the within-case-success rate difference or \eqn{\$dw^{-1}}. See \link{orddom} \emph{NNT within} for dependent groups for details.} 
\item{\emph{$PSb}}{When sample sizes are equal, this gives the Probability of Superiority (PS) for all cases but within-pair changes, i.e.: \deqn{PS_b(Y>X)=\frac{\#(y_i>x_j)}{n_y n_x}}{PS(Y>X)=\#(y_i>x_j)/(n_y n_x)}, \cr
limited to those cases where \eqn{i \neq j}{i<>j}. (For unequal sample sizes, this equals $PSc and $PSw.)} 
\item{\emph{$Ab}}{When sample sizes are equal, this value returns A for all cases where \eqn{i \neq j}{i<>j}. (For unequal sample sizes, this equals $Ac.)} 
\item{\emph{$db}}{When \eqn{n_x=n_y}, this value returns Cliff's delta-between, i.e. all but the paired comparisons or excepting the diagonal of the dominance matrix. The parameter is calculated by taking only those ordinal comparisons into account where \eqn{i \neq j}{i<>j}. (For unequal sample sizes, this equals $dc.)} 
\item{\emph{$NNTb}}{Number needed to treat, based on Cliff's delta-between or \eqn{\$db^{-1}}. See \link{orddom} \emph{NNT between} for dependent groups for details.}
}
\references{
Delaney, H.D. & Vargha, A. (2002). Comparing Several Robust Tests of Stochastic Equality With Ordinally Scaled Variables and Small to Moderate Sized Samples. \emph{Psychological Methods, 7}, 485-503.\cr \cr
Kraemer, H.C. &  Kupfer, D.J. (2006). Size of Treatment Effects and Their Importance to Clinical Research and Practice. \emph{Biological Psychiatry, 59}, 990-996. \cr \cr
Ruscio, J. & Mullen, T. (2012). Confidence Intervals for the Probability of Superiority Effect Size Measure and the Area Under a Receiver Operating Characteristic Curve. \emph{Multivariate Behavioral Research, 47}, 221-223.
Vargha, A., & Delaney, H. D. (1998). The Kruskal-Wallis test and stochastic homogeneity. \emph{Journal of Educational and Behavioral Statistics, 23}, 170-192.\cr \cr
Vargha, A., & Delaney, H. D. (2000). A critique and improvement of the CL common language effect size statistic of McGraw and Wong. \emph{Journal of Educational and Behavioral Statistics, 25}, 101-132.\cr \cr
}
\author{Jens J. Rogmann}
\seealso{\link{dm}, \link{orddom}}
\examples{\dontrun{
> #Example from Efron & Tibshirani (1993, Table 2.1, p. 11)
> #cf. Efron, B. & Tibshirani (1993). An Introduction to the Bootstrap. New York/London: Chapman&Hall.
> y<-c(94,197,16,38,99,141,23) # Treatment Group
> x<-c(52,104,146,10,50,31,40,27,46) # Control Group
> dmes(x,y)
$nx
[1] 9

$ny
[1] 7

$PSc
[1] 0.5714286

$Ac
[1] 0.5714286

$dc
[1] 0.1428571

$NNTc
[1] 7

$PSw
[1] 0.5714286

$Aw
[1] 0.5714286

$dw
[1] 0.1428571

$NNTw
[1] 7

$PSb
[1] 0.5714286

$Ab
[1] 0.5714286

$db
[1] 0.1428571

$NNTb
[1] 7

> ############################################################################
> #Example from Ruscio & Mullen (2012, p. 202)
> #Ruscio, J. & Mullen, T. (2012). Confidence Intervals for the Probability of Superiority Effect Size Measure and the Area Under a Receiver Operating Characteristic Curve, Multivariate Behavioral Research, 47, 201-223.
> x <- c(6,7,8,7,9,6,5,4,7,8,7,6,9,5,4) # Treatment Group
> y <- c(4,3,5,3,6,2,2,1,6,7,4,3,2,4,3) # Control Group
> dmes(y,x)
$nx
[1] 15

$ny
[1] 15

$PSc
[1] 0.8444444

$Ac
[1] 0.8844444

$dc
[1] 0.7688889

$NNTc
[1] 1.300578

$PSw
[1] 1

$Aw
[1] 1

$dw
[1] 1

$NNTw
[1] 1

$PSb
[1] 0.8333333

$Ab
[1] 0.8761905

$db
[1] 0.752381

$NNTb
[1] 1.329114
}}
\keyword{array}
\keyword{nonparametric}
\keyword{robust}

