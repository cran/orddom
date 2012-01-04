\name{orddom}
\alias{orddom}
\title{Ordinal Dominance Statistics}
\description{Returns an array of ordinal dominance statistics based on the input of two 1-column matrices as an alternative to independent or paired group mean comparisons (especially for Cliff's delta statistics).}
\usage{orddom(x,y,alpha=.05,paired=FALSE,outputfile="",studdist=TRUE,symmetric=FALSE)}
\arguments{
  \item{x}{A 1-column matrix with optional column name containing all \eqn{n_x} values or scores of group X or 1 (e.g. experimental or pretest group.), e.g. declared in R as\cr
  \emph{x<-t(matrix(c(\eqn{x_1,x_2,x_3,x_4, \ldots , x_{n_x}}),1))\cr
colnames(x)<-c("[label of group X]").}}
  \item{y}{A 1-column matrix with optional column name containing all \eqn{n_y} values of group Y or 2 (e.g. comparison or post-test group). For paired comparisons (e.g. pre-post), \eqn{n_x = n_y} is required.}
  \item{alpha}{Significance or \eqn{\alpha}-level used for the calculation of the confidence intervals. Default value is \eqn{\alpha = .05} or 5 Percent, giving a 95 Percent CI. For multiple dominance comparisons, a Bonferroni procedure may be implemented: Cliff (1996, p.150) suggested dividing \eqn{\alpha} by the number of possible comparisons, i.e. \eqn{\alpha (\frac{1}{2}k(k-1))^{-1}}{alpha / (.5*k(k-1))} for comparisons beteen \eqn{k} data sets.}
  \item{paired}{By default, independence of the two groups or data sets is assumed. If the number of cases in x and y are equal and paired (e.g. pre-post) comparisons, this should be set to TRUE to return the full array of within, between, combined and metric delta statistics.}
  \item{outputfile}{If a a detailed report of the ordinal dominance analysis is wanted, a filename should be given here. The report as standard text file is written to the current working directory.}
  \item{studdist}{By default, it is assumed that small samples are being examined. In this case, z-values based on Student's t-distribution are used for estimating upper and lower limits of the confidence intervals (CI) as well as z-probabilities. If larger sample sizes  are used, these values approximate estimates based on normally distributed z-values. In this case or if comparing with estimates calculated with orddom versions <1.5 (where z-values based on the Standard Normal Distributions were used), this parameter may be set to FALSE.}
  \item{symmetric}{By default, asymmetric confidence intervals (CI) are being calculated to compensate for positive correlations between the samples as generally recommended by the literature on the delta statistics. To increase power in certain cases, however, (e.g. in small paired samples, cf. Cliff 1996, p. 165), symmetric CIs may also be obtained by setting this argument to TRUE.\cr \cr \cr}
} 
\value{\strong{INDEPENDENT GROUPS} (\emph{paired} argument set to FALSE)\cr \cr
In the case of independent groups or data sets X and Y (e.g. treatment group vs. comparison group), a 2-column-matrix containing 29 rows with values is returned.\cr \cr
The ordinal statistics can be retrieved from the first column (named "ordinal") while the second column (named "metric") contains metric comparison data where appropriate.\cr \cr
\item{[1  \emph{or} ["var1_X", \emph{col#}]}{Original column name of the x input matrix.}
\item{[2  \emph{or} ["var2_Y", \emph{col#}]}{Original column name of the y input matrix.}
\item{[3  \emph{or} ["type", \emph{col#}]}{Returns type of the comparison, in this case "indep".}
\item{[4  \emph{or} ["n in X", \emph{col#}]}{Number of cases in x (i.e. group X sample size).}
\item{[5  \emph{or} ["n in Y", \emph{col#}]}{Number of cases in y (i.e. group Y sample size).}
\item{[6  \emph{or} ["N #Y>X", \emph{col#}]}{Number of occurences of an observation from group y having a higher value than an observation from group x when comparing all x scores with all y scores: \eqn{N_{\#Y>X}=\#(y_i>x_j) ,} where \# denotes "the number of times" whilst comparing each \eqn{i=1, 2, 3, \ldots n_y} score in sample Y with each \eqn{j=1, 2, 3, \ldots n_x} score in sample X (resulting in \eqn{n_x \cdotp n_y}{(n_x)(n_y)} comparisons).}
\item{[7  \emph{or} ["N #Y=X", \emph{col#}]}{Number of occurences of an observation from group y having the same value as an observation from group x: \eqn{N_{\#Y=X}=\#(y_i=x_j)}.}
\item{[8  \emph{or} ["N #Y<X", \emph{col#}]}{Number of occurences of an observation from group y having a smaller value than an observation from group x: \eqn{N_{\#Y<X}=\#(y_i<x_j)}.}
\item{[9  \emph{or} ["PS X>Y", \emph{col#}]}{Common Language CL effect size or Probability of Superiority (PS) of X over Y; \deqn{PS(Y>X)=\frac{\#(y_i>x_j)}{n_y \cdotp n_x}}{PS(Y>X)=\#(y_i>x_j)/(n_y n_x)} (cf. Grissom, 1994; Grissom & Kim, 2005).}
\item{[10  \emph{or} ["PS Y>X", \emph{col#}]}{Common Language CL effect size or Probability of Superiority (PS) of Y over X.}
\item{[11  \emph{or} ["A X>Y", \emph{col#}]}{Vargha and Delaney's A as stochastic superiority of X over Y, calculated as \deqn{A_{X>Y} = PS(X>Y)+.5 PS(X=Y)}{A(X>Y) = PS(X>Y)+.5 PS(X=Y)} (cf. Vargha & Delaney, 1998, 2000; Delaney & Vargha, 2002).}
\item{[12  \emph{or} ["A Y>X", \emph{col#}]}{Vargha and Delaney's A as stochastic superiority of Y over X.}
\item{[13  \emph{or} ["delta", \emph{col#}]}{For column 1 ("ordinal"): Cliff's delta for independent groups (Cliff, 1996; Long et al., 2003):\deqn{d=\frac{\#(y_i>x_j)-\#(y_i<x_j)}{n_y \cdotp n_x}=\frac{\sum_{i} \sum_{j} d_{ij}}{n_y \cdotp n_x}}{d=SUM(SUM(d_ij))/(n_x*n_y)} where \eqn{d_{ij}=sign(y_i-x_j)}{d_ij=sign(y_i-x_j)} across all score comparisons.\cr \cr
In column 2, the metric differences between the means are given: \eqn{\bar{y}-\bar{x}=\bar{d_{ij}}}{My-Mx = M(d_ij)} between all comparable x and y scores with \eqn{d_{ij}=y_i-x_j}{d_ij=y_i-x_j}.}
\item{[14  \emph{or} ["1-alpha", \emph{col#}]}{Significance or \eqn{\alpha}-level for CI estimation, given as percentage between 0 and 100.}
\item{[15  \emph{or} ["CI low", \emph{col#}]}{Unless the default \emph{symmetric} parameter is explicitly set to \emph{TRUE}, improved formulas are used (Feng & Cliff, 2004) to caculate asymmetric confidence interval (CI) boundary estimates of delta or mean difference: \deqn{CI_{lower/upper}=\frac{d - d^{3} \pm t_{\alpha/2} s_{d} \sqrt{1 - 2 d^{2} + d^{4} + t_{\alpha/2}^{2} s_{d}^{2} )}}{1 - d^{2} + t_{\alpha/2}^{2} s_{d}^{2}},}{CI(lower/upper)=(d-d^3+-t s_d ((1-2d^2+d^4+t^2 s_d^2)^(-1/2)))/(1-d^2+t^2 s_d^2)} with t-values at the given \eqn{\alpha}{alpha}-level taken from Student's t distribution by default (unless the \emph{studdist} is set FALSE, in which case t-values are based on z-values from the Standard Normal Distribution ).\cr \cr
In case the \emph{symmetric} argument is explicitly set to \emph{TRUE}, however, ordinary CIs are being calculated with \eqn{CI_{lower/upper}=d \pm t_{\alpha/2} s_{d}}{CI_{lower/upper}=d +- t s_d}.\cr \cr
In any case, if Cliffs' \eqn{d = \pm 1}{d=1 or d=-1}, one CI is assumed being equal to d, the respective other is calculated as \deqn{CI_{lower/upper}=((n_{b}-t_{\alpha/2}^{2}))(n_{b}+t_{\alpha/2}^{2})^{-1},}{CI_lower/upper=(n-t^2)/(n+t^2),} where \eqn{t_{\alpha/2}}{t} is the t-value or z-score at the selected \eqn{\alpha} level (2-tailed) of the respective \emph{studdist}-controlled distribution, and \eqn{n_{b}}{n} the number of observations or cases in the smaller of the two samples.}
\item{[16  \emph{or} ["CI high", \emph{col#}]}{Confidence interval upper boundary estimate of delta or mean difference.}
\item{[17  \emph{or} ["s delta", \emph{col#}]}{Unbiased sample estimate of the delta standard deviation in column 1.\cr
In column 2 ("metric"): Pooled standard deviation of metric mean difference with \eqn{s_{xy}=[((n_x-1)s_x + (n_y-1)s_y)/(n_x+n_y-2)]^{1/2}}{sd_xy=sqrt((((n_x-1)*var(x))+((n_y-1)*var(y)))/(n_x-1+n_y-1))} .}
\item{[18  \emph{or} ["var delta", \emph{col#}]}{Column 1: Variance of delta (unbiased sample estimate), calculated as
\deqn{s_{d}^{2} = \frac{n_y^{2} \sum(d_{i\cdot} - d)^{2} + n_{x}^{2} \sum(d_{\cdot j} - d)^{2} - \sum \sum (d_{ij} - d)^{2}}{n_{x} n_{y} (n_{x}-1) (n_{y}-1)},}{s_d^2 = (n_y^2*SUM((d_i.-d)^2) + n_x^2*SUM((d_.j-d)^2)) - SUM(SUM((dij-d)^2))) / (n_x*n_y*(n_x-1)*(n_y-1)),} or, using the partial variances
\deqn{s_{d}^{2}= \frac{n_y^{2}(n_x-1)s_{d_{i \cdot}}^2 + n_x^{2} (n_y-1) s_{d_{\cdot j}}^{2} - (n_x n_y -1) s_{d_{ij}}^{2}}{n_x n_y (n_x-1) (n_y-1)} ,}{s_d^2=(n_y^2(n_x-1)s_di.^2 + n_x^2(n_y-1)s_d.j^2 - (n_x n_y -1)s_dij^2)/(n_x n_y (n_x-1) (n_y-1)) ,} which can also alternatively be put as
\deqn{s_{d}^{2} = \frac{n_{y} s_{d_{i\cdot}}^{2}}{n_x (n_y-1)} + \frac {n_{x} s_{d_{\cdot j}}^{2}}{n_y (n_x-1)} - \frac{(n_x n_y -1) s_{d_{ij}}^{2}} {n_{x} n_{y}(n_{x}-1) (n_{y}-1)}.}{s_d^2 = (n_y s_di.^2)/(n_x(n_y-1)) + (n_x s_d.j^2)/(n_y(n_x-1)) - (n_x n_y -1)s_dij^2) / (n_x n_y (n_x-1) (n_y-1)).}
(For differences to Cliff's (1996, p. 138) formula see notes to Row 28 ("var dij") below.)\cr \cr
In case this calculation of \eqn{s_{d}^{2}}{s_d^2} yields values of less than \eqn{(1 - d^{2})/(n_x n_y -1)}{(1-d^2)/((n_x n_y)-1)}, this latter formula is used for calculating the variance of delta.\cr \cr
Column 2 contains the pooled \eqn{s_{xy}^{2}}{s_xy^2}.\cr}
\item{[19  \emph{or} ["se delta", \emph{col#}]}{Column 2 only: metric Standard error of mean difference:\cr
\eqn{SE_{xy}=s_{xy} \sqrt{1/n_x+1/n_y}}{SE_xy=s_xy * sqrt(1/n_x+1/n_y)} .}
\item{[20  \emph{or} ["z/t score", \emph{col#}]}{Column 1: z score of delta on the of the respective \emph{studdist}-controlled distribution (Student's t or standard normal).\cr \cr
Column 2: Metric z/t-score (= \eqn{\bar{d_{ij}}/SE_{xy}}{M(d_ij) / SE_xy}).}
\item{[21  \emph{or} ["p (1-tail)", \emph{col#}]}{Probability of z/t score (1 tailed comparison).}
\item{[22  \emph{or} ["p (2-tail)", \emph{col#}]}{Probability of z/t score (2-tailed comparison).}
\item{[23  \emph{or} ["Cohen's \emph{d}", \emph{col#}]}{Cohen's \emph{d} effect size estimate of delta. For Cliff's delta inferred from distributional non-overlap as suggested by Grissom & Kim (2005, p. 106 f.) as well as Romano, Kromrey, Coraggio, & Skowronek (2006, p. 14-15), relating to the relative positions of the distributions of X and Y. When Cliff's delta equals 0, there is no effect, and the Y and X distributions overlap completely. If there are effects, a certain percentage of non-overlap between X and Y is created, and the relative positions of the X and Y distribtions shift. The degree of non-overlap thus is a measure of effect size and is expressed as Cohen's \emph{d} in terms of non-overlap between two normal distributions (based on U1 in Table 2.2.1, Cohen, 1988, p.22). See delta2cohd manual of orddom package.}
\item{[24  \emph{or} ["d CI low", \emph{col#}]}{Cohen's \emph{d} effect size estimate of lower boundary of confidence interval.}
\item{[25  \emph{or} ["d CI high", \emph{col#}]}{Cohen's \emph{d} estimate of upper boundary of confidence interval.}
\item{[26  \emph{or} ["var d.i", \emph{col#}]}{Row variance of dominance/difference matrix, calculated as \cr \eqn{(n_x-1)^{-1} \sum (d_{i\cdot}-d)^2}{SUM((di.-d)^2)/(n_x-1)}.}
\item{[27  \emph{or} ["var dj.", \emph{col#}]}{Column variance of dominance/difference matrix, calculated as \cr \eqn{(n_y-1)^{-1} \sum (d_{\cdot j}-d)^2}{SUM((d.j-d)^2)/(n_y-1)}.}
\item{[28  \emph{or} ["var dij", \emph{col#}]}{Variance of dominance/difference matrix as sample estimate according to Long et al. (2003, section 3.3 before eqn. 67):\deqn{s_{d_{ij}}^{2}=\frac{\sum \sum(d_{ij}-d)^{2}}{n_x n_y -1}=\frac{\sum d_{ij}^{2} - \frac{(\sum d_{ij})^{2}}{n_x n_y}}{n_x n_y -1},}{s_d_ij=SUM(SUM((d_ij-d)^2))/(n_x*n_y-1),} thus avoiding Cliff's original (1996, p. 138) suggestion to use \eqn{(n_x-1)(n_y-1)} as the denominator).}
\item{[29  \emph{or} ["df", \emph{col#}]}{If the \emph{studdist} parameter is not set to FALSE, column 1 returns the degrees of freedom (\eqn{df}) used for CI as well as z/t-score and z-probability estimates. \cr
In column 2 ("metric") df as used for metric t-test.
\cr \cr \cr \cr \cr \cr \cr \cr}
\strong{DEPENDENT/PAIRED GROUPS} (\emph{paired} argument set to TRUE)\cr \cr
In the case of paired data (e.g. pretest-posttest comparisons of the \eqn{n_x=n_y} same subjects), a 4-column-matrix containing 29 rows with values is returned.\cr \cr
The ordinal statistics for \eqn{d_{ij}} can be retrieved from the first three columns (named \cr
\item{within   [.,1]}{for the \eqn{n_x=n_y} within-pair changes (where \eqn{i=j} in all cases);}
\item{between  [.,2]}{for the overall distribution changes, based on all \eqn{n^2-n = n(n-1)} comparisons where \eqn{i \ne j}; \cr
and}
\item{combined [.,3]}{for combined inferences \eqn{d_w + d_b}.}
Here, the fourth column (named "metric") contains metric comparison data.\cr \cr
\item{[1  \emph{or} ["var1_X_pre", \emph{col#}]}{Original column name of the x (or pretest) input matrix.}
\item{[2  \emph{or} ["var2_Y_post", \emph{col#}]}{Original column name of the y (or posttest) input matrix.}
\item{[3  \emph{or} ["type", \emph{col#}]}{Returns type of the comparison, in this case \emph{paired}.}
\item{[4  \emph{or} ["N #Y>X", \emph{col#}]}{Number of occurences (\#) of a posttest observation \eqn{y_i} having a higher value than a pretest observation \eqn{x_j}: \eqn{N_{\#Y>X}=\#(y_i>x_j) ,} limited to the respective pairs under observation in \emph{within}, \emph{between} or \emph{combined}.\cr Column 4 equals column 3.}
\item{[5  \emph{or} ["N #Y=X", \emph{col#}]}{Number of occurences of a posttest observation having the same value as a pretest observation, limited to the respective pairs under observation in \emph{within}, \emph{between} or \emph{combined}.\cr
Column 4 equals column 3.}
\item{[6  \emph{or} ["N #Y<X", \emph{col#}]}{Number of occurences of a posttest observation having a smaller value than a pretest observation, limited to the respective pairs under observation in \emph{within}, \emph{between} or \emph{combined}.\cr
Column 4 equals column 3.}
\item{[7  \emph{or} ["PS X>Y", \emph{col#}]}{Common Language CL effect size or Probability of Superiority (PS) of Y over X (Grissom, 1994; Grissom & Kim, 2005) (limited to the respective pairs under observation in \emph{within}, \emph{between} or \emph{combined}): \deqn{PS(Y>X)=\frac{\#(y_i>x_j)}{n_y \cdotp n_x} .}{PS(Y>X)= \#(y_i>x_j)/(n_y n_x) .}
Column 4 equals column 3.}
\item{[8  \emph{or} ["PS Y>X", \emph{col#}]}{Common Language CL effect size or Probability of Superiority (PS) of X over Y (Grissom, 1994; Grissom & Kim, 2005) (limited to the respective pairs under observation in \emph{within}, \emph{between} or \emph{combined}).\cr
Column 4 equals column 3.}
\item{[9  \emph{or} ["A X>Y", \emph{col#}]}{Vargha and Delaney's A as stochastic superiority of X over Y, calculated as \eqn{A_{X>Y} = PS(X>Y)+.5 PS(X=Y)} (Vargha & Delaney, 1998, 2000; Delaney & Vargha, 2002), limited to the respective pairs under observation in \emph{within}, \emph{between} or \emph{combined}.\cr
Column 4 equals column 3.}
\item{[10  \emph{or} ["A Y>X", \emph{col#}]}{Vargha and Delaney's A as stochastic superiority of Y over X, limited to the respective pairs under observation in \emph{within}, \emph{between} or \emph{combined}.\cr
Column 4 equals column 3.}
\item{[11  \emph{or} ["delta", \emph{col#}]}{For columns 1 to 3 ("ordinal"), the respective delta for dependent groups (Cliff, 1996; Long et al., 2003; Feng, 2007) is reported. With \eqn{d_{ij}=sign(y_i-x_j)}{d_ij=sign(y_i-x_j)},\cr \cr
Column 1  reports the (\emph{within}) value, which is the "difference between the proportion of individual subjects who change in one direction and the proportion of individuals who change in the other" (Cliff, 1996, p. 159), calculated as
\deqn{d_w=({\sum_{i}}_{=}\sum_{j} d_{ii})/ n ,}{d_w=SUM(SUM(d_ii))/n,}
where \eqn{i=j} in the \eqn{n = n_x = n_y} possible paired comparisons.\cr \cr
"The extent to which the overall distribution has moved, except for the self-comparisons" (Cliff, 1996, p. 160) is given in column 2, the delta-(\emph{between}) statistic. It is estimated by the average between-subject dominance, calculated as
\deqn{d_b=({\sum_{i}}_{\not=}\sum_{j} d_{ij})/ (n (n-1)) ,}{d_b=SUM(SUM(d_ij))/(n(n-1)),}
where \eqn{i \not= j}{i <> j}.\cr \cr
Column 3 reports the combination effect \eqn{d_w + d_b}.
\cr \cr
In column 4 ("metric"), the differences between subsample means are reported: \eqn{\bar{y}-\bar{x}}{My-Mx}. \cr}
\item{[12  \emph{or} ["1-alpha", \emph{col#}]}{Significance or \eqn{\alpha}-level for CI estimation, given as percentage between 0 and 100.}
\item{[13  \emph{or} ["CI low", \emph{col#}]}{Confidence interval (CI) lower boundary estimate. Unless the default \emph{symmetric} parameter is explicitly set to \emph{TRUE}, asymmetric Confidence interval (CI) boundary estimates for ordinal differences are calculated (Feng & Cliff, 2004; Feng, 2007) as 
\deqn{CI_{lower/upper}=\frac{d - d^{3} \pm t_{\alpha/2} s_{d} \sqrt{1 - 2 d^{2} + d^{4} + t_{\alpha/2}^{2} s_{d}^{2} )}}{1 - d^{2} + t_{\alpha/2}^{2} s_{d}^{2}} ,}{CI(lower/upper)=(d-d^3+-t s_d ((1-2d^2+d^4+t^2 s_d^2)^(-1/2)))/(1-d^2+t^2 s_d^2) ,}
with t-values at the respective significance level based on either Student's t or on z-values from the Standard Normal Distribution, depending on the \emph{studdist} argument. \cr \cr 
However, using an asymmetric CI may reduce power of \eqn{d_b} value given in \samp{[10,2]}, especially in small paired samples (Cliff, 1996, p. 165). To obtain symmetric CI estimates with \eqn{CI_{lower/upper}=d \pm t_{\alpha/2} s_{d}}{CI_{lower/upper}=d +- t s_d}, the default \emph{symmetric} argument must be set to \emph{TRUE}.\cr \cr
In any case, if \eqn{d = \pm 1}{d=1 or d=-1}, one CI is set as equal to d, the other is calculated as \deqn{CI_{lower/upper}=((n_{b}-t_{\alpha/2}^{2}))(n_{b}+t_{\alpha/2}^{2})^{-1},}{CI_lower/upper=(n-t^2)/(n+t^2),} where \eqn{t_{\alpha/2}}{t} is the t-value or z-score at the selected \eqn{\alpha} level (2-tailed) of the respective \emph{studdist}-controlled distribution, and \eqn{n_{b}}{n} the number of observations or cases in the smaller of the two samples.}
\item{[14  \emph{or} ["CI high", \emph{col#}]}{Confidence interval upper boundary estimate.}
\item{[15  \emph{or} ["s delta", \emph{col#}]}{Estimated standard deviation of the respective delta statistic. Column 4 reports the metric standard deviation of the paired (within) differences.}
\item{[16  \emph{or} ["var delta", \emph{col#}]}{Unbiased estimates of the variances of the respective delta statistic.\cr \cr
Column 1 reports the \emph{within} value, calculated as 
\deqn{s_{d_w}^{2}=(n(n-1))^{-1}(\sum (d_{ii}-d_w)^{2}).}{s_dw^2=SUM((d_ii-d_w)^2)/(n(n-1)).}
Please note that in various pieces of the available research literature (e.g. Cliff, 1996, eq. 6.8, p. 161), \eqn{s_{d_w}^{2}} is erroneously reported to be calculated as \eqn{s_{d_w}^{2}=(n-1)^{-1}(\sum (d_{ii}-d_w)^{2}).}{s_dw^2=SUM((d_ii-d_w)^2)/(n-1).} The denominator, however must read \eqn{n(n-1)} as "using just \eqn{(n-1)} would give the variance of the individual \eqn{d_{ii}} whereas we want the variance of \eqn{d_w}, which is a kind of mean" (Feng, 07.02.2011, personal communication).\cr \cr
The (\emph{between}) unbiased estimate in column 2 is calculated as \cr
\eqn{s_{d_b}^{2}=[(n-1)^{2}(\sum (d_{i \cdot}-d_b)^{2} + \sum (d_{\cdot j}-d_b)^{2} + 2 \sum (d_{i \cdot}-d_b)(d_{\cdot j}-d_b)) - \sum \sum (d_{ij}-d_b)^{2} - \sum \sum (d_{ij}-d_b)(d_{ji}-d_b)][n(n-1)(n-2)(n-3)]^{-1}
.}{s_db^2=[(n-1)^2 * (SUM((d_i.-d_b)^2) + SUM((d_.j-d_b)^2) + 2 SUM((d_i.-d_b)(d_.j-d_b))) - SUM(SUM((d_ij-d_b)^2)) - SUM(SUM((d_ij-d_b)(d_ji-d_b)))]/(n(n-1)(n-2)(n-3)).} \cr \cr
Since \eqn{d_w} and \eqn{d_b} are interdependent, the \emph{combined} effect involves taking into account their estimated covariance when calculating the unbiased estimate for the variance for the sum of \eqn{d_w} and \eqn{d_b}, which is reported in column 3 as\cr
\eqn{s_{d_w+d_b}^{2}=s_{d_w}^{2} + s_{d_b}^{2} + 2 \widehat{cov}(d_b,d_w),}{var(d_w+d_b)= s_dw^2 + s_db^2 + 2 Est[cov(d_b,d_w)],} with \cr
\eqn{\widehat{cov}(d_b,d_w)=}{Est[cov(d_b,d_w)=} \eqn{\left( \sum_{i} [d_{ii}(\sum_{j}(d_{ij}) + \sum_{j}(d_{ji})] - 2n(n-1) d_b d_w \right) (n(n-1)(n-2))^{-1}}{(SUM_i(d_ii(SUM_j(d_ij)+SUM_j(d_ji))) - 2n(n-1)(d_b)(d_w)) / (n(n-1)(n-2))}.\cr \cr
Column 4 reports the metric variance of the paired (within) differences. \cr \cr}
\item{[17  \emph{or} ["z/t score", \emph{col#}]}{z score of delta. In column 4 ("metric") equal to the t-test score.}
\item{[18  \emph{or} ["p (1-tail)", \emph{col#}]}{Probability of z-score (1-tailed comparison).}
\item{[19  \emph{or} ["p (2-tail)", \emph{col#}]}{Probability of z-score (2-tailed comparison).}
\item{[20  \emph{or} ["Cohen's \emph{d}", \emph{col#}]}{Cohen's \emph{d} estimate of the respective delta value (see above). \cr Not available for the combined delta in column 3.}
\item{[21  \emph{or} ["d CI low", \emph{col#}]}{Cohen's \emph{d} estimate of lower boundary of the respective confidence interval. \cr Not available for the combined delta in column 3.}
\item{[22  \emph{or} ["d CI high", \emph{col#}]}{Cohen's \emph{d} estimate of lower boundary of the respective confidence interval. \cr Not available for the combined delta in column 3.}
\item{[23,3] \emph{or} ["var d.i",\emph{combined}]}{Component of \eqn{s_{d_w+d_b}^2}: \eqn{s_{di.}^2} (Available for the combined analyses in column 3 only.)}
\item{[24,3] \emph{or} ["var dj.",\emph{combined}]}{Component of \eqn{s_{d_w+d_b}^2}: \eqn{s_{d.i}^2} (Third column only.)}
\item{[25,3] \emph{or} ["cov(di,dj)",\emph{combined}]}{Component of \eqn{s_{d_w+d_b}^2}: \eqn{cov(d_{i.},d_{.j})} (Third column only.)}
\item{[26,3] \emph{or} ["var dij",\emph{combined}]}{Component of \eqn{s_{d_w+d_b}^2}: \eqn{s_{d_{ij}}^2} (Third column only.)}
\item{[27,3] \emph{or} ["cov(dih,dhi)",\emph{combined}]}{Component of \eqn{s_{d_w+d_b}^2}: \eqn{cov(d_{ih},d_{hi})} (Third column only.)}
\item{[28,3] \emph{or} ["cov(db,dw)",\emph{combined}]}{Estimated covariance between \eqn{d_b} and \eqn{d_w}: \eqn{\widehat{cov}(d_b,d_w)}{Est[cov(d_b,d_w)]} (for purposes of combined inferences). (Third column only.)}
\item{[29  \emph{or} ["df", \emph{col#}]}{Unless the \emph{studdist} argument is not set to FALSE, the degrees of Freedom \emph{df} used for the CI and z-score calculations are reported in column 1.\cr \cr
Column 2 returns Welch's \emph{df} for the metric t-test for dependent samples.}
}
\references{Cliff, N. (1996). \emph{Ordinal Methods for Behavioral Data Analysis}. Mahwah, NJ: Lawrence Erlbaum.\cr \cr
Cohen, J. (1988). \emph{Statistical power analysis for the behavioral sciences (2nd edition)}. New York: Academic Press.\cr \cr
Feng, D. (2007). Robustness and Power of Ordinal d for Paired Data. In Shlomo S. Sawilowsky (Ed.), \emph{Real Data Analysis} (pp. 163-183). Greenwich, CT : Information Age Publishing.\cr \cr
Feng, D., & Cliff, N. (2004). Monte Carlo Evaluation of Ordinal d with Improved Confidence Interval. \emph{Journal of Modern Applied Statistical Methods, 3}(2), 322-332.\cr \cr
Long, J. D., Feng, D., & Cliff, N. (2003). Ordinal analysis of behavioral data. In J. Schinka & W. F. Velicer (eds.), \emph{Research Methods in Psychology. Volume 2 of Handbook of Psychology} (I. B. Weiner, Editor-in-Chief). New York: John Wiley & Sons.\cr \cr
Delaney, H.D. & Vargha, A. (2002). Comparing Several Robust Tests of Stochastic Equality With Ordinally Scaled Variables and Small to Moderate Sized Samples. \emph{Psychological Methods, 7}, 485-503.\cr \cr
Grissom, R.J. (1994). Probability of the superior outcome of one treatment over another. \emph{Journal of Applied Psychology, 79}, 314-316.\cr \cr
Grissom, R.J. & Kim, J.J. (2005). \emph{Effect sizes for research. A broad practical approach}. Mahwah, NJ, USA: Erlbaum.\cr \cr
McGraw, K.O. & Wong, S.P. (1992). A common language effect size statistic. \emph{Psychological Bulletin, 111}, 361-365.\cr \cr
Romano, J., Kromrey, J. D., Coraggio, J., & Skowronek, J. (2006). \emph{Appropriate statistics for ordinal level data: Should we really be using t-test and Cohen's d for evaluating group differences on the NSSE and other surveys?}. Paper presented at the annual meeting of the Florida Association of Institutional Research, Feb. 1-3, 2006, Cocoa Beach, Florida. Last retrieved January 2, 2012 from www.florida-air.org/romano06.pdf \cr \cr
Vargha, A., & Delaney, H. D. (1998). The Kruskal-Wallis test and stochastic homogeneity. \emph{Journal of Educational and Behavioral Statistics, 23}, 170-192.\cr \cr
Vargha, A., & Delaney, H. D. (2000). A critique and improvement of the CL common language effect size statistic of McGraw and Wong. \emph{Journal of Educational and Behavioral Statistics, 25}, 101-132.\cr \cr
}
\note{\emph{Major changes from orddom Version 1.0 to 1.5}
\itemize{
\item{Calculation of CI and delta z-score-estimates can now be based on Students t-distribution rather than using fixed normal distribution z-scores.}
\item{Symmetric CIs can now be obtained to increase power of the delta statistics in certain cases.}
\item{Formulas used for calculation added in orddom manual.}
\item{Probability of Superiority statistic as well as variance estimates for delta in the independent groups analyses were corrected.}
\item{Minor changes were implemented to allow for calculation of \eqn{d = \pm 1}{d=+-1} extreme cases without error abort.}
\item{Output of raw y-dataset in independent group analysis was corrected.}
\item{Dependencies on packages \emph{psych} and \emph{compute.es} declared in DESCRIPTION and NAMESPACE files.}
}}
\author{Jens J. Rogmann, University of Hamburg, Department of Psychology, \cr Hamburg, Germany (Jens.Rogmann@uni-hamburg.de)}
\examples{\dontrun{
#Independent Samples (Data taken from Long et al. (2003), Table 4}
x<-t(matrix(c(3,3,3,4,5,6,12,12,13,14,15,15,15,15,15,16,18,18,18,23,23,27,28,28,43),1))
colnames(x)<-c("Nonalcohol.")
y<-t(matrix(c(1,4,6,7,7,14,14,18,19,20,21,24,25,26,26,26,27,28,28,30,33,33,44,45,50),1))
colnames(y)<-c("Alcoholic")
orddom(x,y,,paired=FALSE,"tmp_r.txt",,)
\dontrun{
#Paired Comparison with CSV-Output (Data taken from Long et al. (2003), Table 4}
x<-t(matrix(c(2,6,6,7,7,8,8,9,9,9,10,10,10,11,11,12,13,14,15,16),1))
colnames(x)<-c("Incidental")
y<-t(matrix(c(4,11,8,9,10,11,11,5,14,12,13,10,14,16,14,13,15,15,16,10),1))
colnames(y)<-c("Intentional")
write.table(orddom(y,x,,TRUE,,FALSE,TRUE),"tmp_r.txt",quote=TRUE,sep=",",dec=",",eol=";\n")
\dontrun{
#Directly returns d_b of the paired comparison}
orddom(x,y,,TRUE,,,)[11,2]}
\keyword{nonparametric}
\keyword{robust}
\keyword{htest}

