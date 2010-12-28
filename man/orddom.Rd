\name{orddom}
\alias{orddom}
\title{Ordinal Dominance Statistics}
\description{Returns an array of ordinal dominance statistics based on the input of two 1-column matrices as an alternative to independent or paired group mean comparisons (especially for Cliff's delta statistics).}
\usage{orddom(x, y, alpha = 0.05, paired = FALSE, outputfile = "")}
\arguments{
  \item{x}{A 1-column matrix with optional column name containing all values of group 1 (e.g. Pretest group.)}
  \item{y}{A 1-column matrix with optional column name containing all values of group 2. (e.g. Posttest group)}
  \item{alpha}{Alpha level used for the calculation of the confidence interval. Thus, for multiple dominance comparisons, a Bonferroni procedure may be implemented. Default value is \eqn{\alpha} = 5 Percent, giving a 95 Percent CI.}
  \item{paired}{As a default, independence of the two groups is assumed. If the number of cases in x and y are equal and paired (e.g. pre-post) comparisons, this should be set to TRUE to return the full array of within, between, combined and metric delta statistics.}
  \item{outputfile}{If a a detailed report of the ordinal dominance analysis is wanted, a filename should be given here. The report as standard text file is written to the current working directory.}
}
\value{\strong{INDEPENDENT GROUPS}\cr
In the case of independent groups, a 2-column-matrix containing 28 rows with values is returned.\cr
The ordinal statistics can be retrieved from the first column (named "ordinal") while the second column (named "metric") contains the metric comparison data.\cr
\item{[1  \emph{or} ["var1_X", \emph{col#}]}{Original column name of the x input matrix.}
\item{[2  \emph{or} ["var2_Y", \emph{col#}]}{Original column name of the y input matrix.}
\item{[3  \emph{or} ["type", \emph{col#}]}{Returns type of the comparison, in this case "indep".}
\item{[4  \emph{or} ["n in X", \emph{col#}]}{Number of cases in x.}
\item{[5  \emph{or} ["n in Y", \emph{col#}]}{Number of cases in y.}
\item{[6  \emph{or} ["N #Y>X", \emph{col#}]}{Number of occurences of an observation from group y having a higher value than an observation from group x.}
\item{[7  \emph{or} ["N #Y=X", \emph{col#}]}{Number of occurences of an observation from group y having the same value as an observation from group x.}
\item{[8  \emph{or} ["N #Y<X", \emph{col#}]}{Number of occurences of an observation from group y having a smaller value than an observation from group x.}
\item{[9  \emph{or} ["PS X>Y", \emph{col#}]}{Common Language CL effect size or Probability of Superiority (PS) of X over Y (Grissom, 1994; Grissom & Kim, 2005).}
\item{[10  \emph{or} ["PS Y>X", \emph{col#}]}{Common Language CL effect size or Probability of Superiority (PS) of Y over X (Grissom, 1994; Grissom & Kim, 2005).}
\item{[11  \emph{or} ["A X>Y", \emph{col#}]}{Vargha and Delaney's A as stochastic superiority of X over Y, calculated as \eqn{A_{X>Y} = PS(X>Y)+.5 PS(X=Y)} (Vargha & Delaney, 1998, 2000; Delaney & Vargha, 2002).}
\item{[12  \emph{or} ["A Y>X", \emph{col#}]}{Vargha and Delaney's A as stochastic superiority of Y over X.}
\item{[13  \emph{or} ["delta", \emph{col#}]}{For column 1 ("ordinal"): Cliff's delta for independent groups (Cliff, 1996; Long et al., 2003). For column 2 ("metric"): Difference between means of x and y.}
\item{[14  \emph{or} ["1-alpha", \emph{col#}]}{Alpha level for CI estimation.}
\item{[15  \emph{or} ["CI low", \emph{col#}]}{Confidence interval lower boundary estimate of delta or mean difference. For Cliff's delta, an improved formula was used (Feng & Cliff, 2004).}
\item{[16  \emph{or} ["CI high", \emph{col#}]}{Confidence interval upper boundary estimate of delta or mean difference.}
\item{[17  \emph{or} ["s delta", \emph{col#}]}{Standard deviation of delta statistics.}
\item{[18  \emph{or} ["var delta", \emph{col#}]}{Variance of delta.}
\item{[19  \emph{or} ["se delta", \emph{col#}]}{Only in column 2 ("Metric"): Standard error of mean differences.}
\item{[20  \emph{or} ["z/t score", \emph{col#}]}{z score of delta. In column 2 ("metric") equal to the t-test score.}
\item{[21  \emph{or} ["p (1-tail)", \emph{col#}]}{Probability of z score (1 tailed comparison).}
\item{[22  \emph{or} ["p (2-tail)", \emph{col#}]}{Probability of z score (2-tailed comparison).}
\item{[23  \emph{or} ["Cohen's d", \emph{col#}]}{Cohen's d estimate of delta. For Cliff's delta inferred from distributional non-overlap.}
\item{[24  \emph{or} ["d CI low", \emph{col#}]}{Cohen's d estimate of lower boundary of confidence interval.}
\item{[25  \emph{or} ["d CI high", \emph{col#}]}{Cohen's d estimate of upper boundary of confidence interval.}
\item{[26  \emph{or} ["var d.i", \emph{col#}]}{Row variance of dominance/difference matrix.}
\item{[27  \emph{or} ["var dj.", \emph{col#}]}{Column variance of dominance/difference matrix.}
\item{[28  \emph{or} ["var dij", \emph{col#}]}{Variance of dominance/difference matrix.}
\cr
\strong{DEPENDENT/PAIRED GROUPS}\cr
In the case of paired data, a 4-column-matrix containing 27 rows with values is returned.\cr
The ordinal statistics can be retrieved from the first three columns (named\cr
-"within" [...,1] for the within pair changes \eqn{d_{ii}},\cr
-"between" [...,2] for the overall distribution changes \eqn{d_{ij}, i \ne j}, and\cr
-"combined" [...,3] for combined inferences \eqn{d_w + d_b}.\cr
Here, the fourth column (named "metric") contains the metric comparison data.\cr
\item{[1  \emph{or} ["var1_X_pre", \emph{col#}]}{Original column name of the x or pretest input matrix.}
\item{[2  \emph{or} ["var2_Y_post", \emph{col#}]}{Original column name of the y or posttest input matrix.}
\item{[3  \emph{or} ["type", \emph{col#}]}{Returns type of the comparison, in this case "paired".}
\item{[4  \emph{or} ["N #Y>X", \emph{col#}]}{Number of occurences of a posttest observation having a higher value than a pretest observation, limited to the respective pairs under observation in "within", "between" or "combined".}
\item{[5  \emph{or} ["N #Y=X", \emph{col#}]}{Number of occurences of a posttest observation having the same value as a pretest observation, limited to the respective pairs under observation in "within", "between" or "combined".}
\item{[6  \emph{or} ["N #Y<X", \emph{col#}]}{Number of occurences of a posttest observation having a smaller value than a pretest observation, limited to the respective pairs under observation in "within", "between" or "combined".}
\item{[7  \emph{or} ["PS X>Y", \emph{col#}]}{Common Language CL effect size or Probability of Superiority (PS) of X over Y (Grissom, 1994; Grissom & Kim, 2005) (limited to the respective pairs under observation in "within", "between" or "combined").}
\item{[8  \emph{or} ["PS Y>X", \emph{col#}]}{Common Language CL effect size or Probability of Superiority (PS) of Y over X (Grissom, 1994; Grissom & Kim, 2005) (limited to the respective pairs under observation in "within", "between" or "combined").}
\item{[9  \emph{or} ["A X>Y", \emph{col#}]}{Vargha and Delaney's A as stochastic superiority of X over Y, calculated as \eqn{A_{X>Y} = PS(X>Y)+.5 PS(X=Y)} (Vargha & Delaney, 1998, 2000; Delaney & Vargha, 2002), limited to the respective pairs under observation in "within", "between" or "combined".}
\item{[10  \emph{or} ["A Y>X", \emph{col#}]}{Vargha and Delaney's A as stochastic superiority of X over Y, limited to the respective pairs under observation in "within", "between" or "combined".}
\item{[11  \emph{or} ["delta", \emph{col#}]}{For column 1 ("ordinal"): Cliff's delta for independent groups (Cliff, 1996; Long et al., 2003). For column 2 ("metric"): Difference between means of x and y.}
\item{[12  \emph{or} ["1-alpha", \emph{col#}]}{Alpha level for CI estimation.}
\item{[13  \emph{or} ["CI low", \emph{col#}]}{Confidence interval lower boundary estimate of delta or mean difference. For Cliff's delta, an improved formula was used (Feng & Cliff, 2004; Feng, 2007).}
\item{[14  \emph{or} ["CI high", \emph{col#}]}{Confidence interval upper boundary estimate.}
\item{[15  \emph{or} ["s delta", \emph{col#}]}{Standard deviation of the respective delta statistic.}
\item{[16  \emph{or} ["var delta", \emph{col#}]}{Variance of the respective delta statistic.}
\item{[17  \emph{or} ["z/t score", \emph{col#}]}{z score of delta. In column 4 ("metric") equal to the t-test score.}
\item{[18  \emph{or} ["p (1-tail)", \emph{col#}]}{Probability f z-score (1-tailed comparison).}
\item{[19  \emph{or} ["p (2-tail)", \emph{col#}]}{Probability f z-score (2-tailed comparison).}
\item{[20  \emph{or} ["Cohen's d", \emph{col#}]}{Cohen's d estimate of the respective delta. For Cliff's delta inferred from distributional non-overlap.}
\item{[21  \emph{or} ["d CI low", \emph{col#}]}{Cohen's d estimate of lower boundary of the respective confidence interval.}
\item{[22  \emph{or} ["d CI high", \emph{col#}]}{Cohen's d estimate of lower boundary of the respective confidence interval.}
\item{[23,3] \emph{or} ["var d.i","combined"]}{Component of \eqn{s_{d_b}^2}: \eqn{s_{di.}^2} (Third column only.)}
\item{[24,3] \emph{or} ["var dj.","combined"]}{Component of \eqn{s_{d_b}^2}: \eqn{s_{d.i}^2} (Third column only.)}
\item{[25,3] \emph{or} ["cov(di,dj)","combined"]}{Component of \eqn{s_{d_b}^2}: \eqn{cov(d_{i.},d_{.j})} (Third column only.)}
\item{[26,3] \emph{or} ["var dij","combined"]}{Component of \eqn{s_{d_b}^2}: \eqn{s_{d_{ij}}^2} (Third column only.)}
\item{[27,3] \emph{or} ["cov(dih,dhi)","combined"]}{Component of \eqn{s_{d_b}^2}: \eqn{cov(d_{ih},d_{hi})} (Third column only.)}
\item{[28,3] \emph{or} ["cov(db,dw)","combined"]}{Covariance between \eqn{d_b} and \eqn{d_w}: \eqn{cov(d_b,d_w)} (for purposes of combined inferences). (Third column only.)}
}
\references{Cliff, N. (1996). \emph{Ordinal Methods for Behavioral Data Analysis}. Mahwah, NJ: Lawrence Erlbaum.\cr
Feng, D. (2007). Robustness and Power of Ordinal d for Paired Data. In Shlomo S. Sawilowsky (Ed.), \emph{Real Data Analysis} (pp. 163-183). Greenwich, CT : Information Age Publishing.\cr
Feng, D., & Cliff, N. (2004). Monte Carlo Evaluation of Ordinal d with Improved Confidence Interval. \emph{Journal of Modern Applied Statistical Methods, 3}(2), 322-332.\cr
Long, J. D., Feng, D., & Cliff, N. (2003). Ordinal analysis of behavioral data. In J. Schinka & W. F. Velicer (eds.), \emph{Research Methods in Psychology. Volume 2 of Handbook of Psychology} (I. B. Weiner, Editor-in-Chief). New York: John Wiley & Sons.\cr
Delaney, H.D. & Vargha, A. (2002). Comparing Several Robust Tests of Stochastic Equality With Ordinally Scaled Variables and Small to Moderate Sized Samples. \emph{Psychological Methods, 7}, 485-503.\cr
Grissom, R.J. (1994). Probability of the superior outcome of one treatment over another. \emph{Journal of Applied Psychology, 79}, 314-316.\cr
Grissom, R.J. & Kim, J.J. (2005). \emph{Effect sizes for research. A broad practical approach}. Mahwah, NJ, USA: Erlbaum.\cr
McGraw, K.O. & Wong, S.P. (1992). A common language effect size statistic. \emph{Psychological Bulletin, 111}, 361-365.\cr
Vargha, A., & Delaney, H. D. (1998). The Kruskal-Wallis test and stochastic homogeneity. \emph{Journal of Educational and Behavioral Statistics, 23}, 170-192.\cr
Vargha, A., & Delaney, H. D. (2000). A critique and improvement of the CL common language effect size statistic of McGraw and Wong. \emph{Journal of Educational and Behavioral Statistics, 25}, 101-132.\cr
}
\author{Jens Rogmann}
\examples{\dontrun{
#Independent Groups}
x<-t(matrix(c(1,4,6,7,7,14,14,18,19,20,21,24,25,26,26,26,27,28,28,30,33,33,44,45,50),1))
colnames(x)<-c("Group 1")
y<-t(matrix(c(3,3,3,4,5,6,12,12,13,14,15,15,15,15,15,16,18,18,18,23,23,27,28,28,43),1))
colnames(y)<-c("Group 2")
orddom(x,y,,paired=FALSE,"tmp_r.txt")
\dontrun{
#Paired Comparison with CSV-Output}
x<-t(matrix(c(2,6,6,7,7,8,8,9,9,9,10,10,10,11,11,12,13,14,15,16),1))
colnames(x)<-c("Incidental")
y<-t(matrix(c(4,11,8,9,10,11,11,5,14,12,13,10,14,16,14,13,15,15,16,10),1))
colnames(y)<-c("Intentional")
write.table(orddom(x,y,.05,TRUE),"tmp_r.txt",quote=TRUE,sep=",",eol=";\\n")
}
\keyword{nonparametric}
\keyword{robust}
\keyword{htest}
