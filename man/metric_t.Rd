\name{metric_t}
\alias{metric_t}
\title{Metric t-test parameters matrix}
\description{Returns a matrix of independent or paired t-test data for comparison to ordinal alternatives}
\usage{metric_t(a,b,alpha=0.05,paired=FALSE,t.welch=TRUE)}
\arguments{
  \item{a}{First dataset (vector or matrix).}
  \item{b}{Second dataset (vector or matrix).}
  \item{alpha}{Significance or \eqn{\alpha}-level used for the calculation of the confidence intervals. Default value is \eqn{\alpha = .05} or 5 Percent.}
  \item{paired}{By default, independence of the two groups or data sets is assumed. If the number of cases in x and y are equal and paired (e.g. pre-post) comparisons, this should be set to TRUE.}
  \item{t.welch}{By default, the variances of the two datasets are not assumed equal. If the pooled variance is needed for t, p, and df this should be set to FALSE. This setting has no effect on the calculation of Cohens'd.}
}
\value{
\item{[1,1] \emph{or} ["Diff M" ,1]}{Mean Difference \eqn{\bar{y}-\bar{x}}{My-Mx} or estimate (in the paired case) See \link{t.test} for details.}
\item{[2,1] \emph{or} ["t value" ,1]  \emph{or} ["t(dep.)" ,1]}{Value of the t-statistic for the independent or the paired case. See \link{t.test} for details.}
\item{[3,1] \emph{or} ["df" ,1]  \emph{or} ["df" ,1]}{Degrees of freedom for the t-statistic. For independent samples, the Welch approximation of degrees of freedom is returned unless t.welch is set to FALSE. See \link{t.test} for details.}
\item{[4,1] \emph{or} ["p value" ,1]}{The p-value of the test. See \link{t.test} for details. For independent samples, the Welch approximation of degrees of freedom is returned unless t.welch is set to FALSE.}
\item{[5,1] \emph{or} ["Cohen's d" ,1]}{Cohen's d effect size for both the independent and the paired case calculated using student's t (i.e. assuming pooled variance) as
\deqn{{d_{Cohen}} = {t_{(pooled var)}} \sqrt{\frac{n_{y} + n_{x}}{n_{y} n_{x}}},}{Cohen's d = t ((n_y+n_x)/(n_y n_x))^(1/2),}
following the advice of Dunlap, Cortina, Vaslow and Burke (1996) who suggested using the independent group t-value and the original standard deviations also for the paired case to avoid overestimation of the effect size.}
}
\author{Jens J. Rogmann, University of Hamburg, Department of Psychology, \cr Hamburg, Germany (Jens.Rogmann@uni-hamburg.de)}
\references{Dunlap, W. P., Cortina, J. M., Vaslow, J. B., & Burke, M. J. (1996). Meta-analysis of experiments with matched groups or repeated measures designs. \emph{Psychological Methods, 1}, 170-177.}
\seealso{\link{t.test}}
\examples{\dontrun{
> #Example from Dunlap et al. (1996), Table 1
> y<-c(27,25,30,29,30,33,31,35)
> x<-c(21,25,23,26,27,26,29,31)
> metric_t(x,y)
                 [,1]
Diff M     4.00000000
t value    2.52982213
df Welch  14.00000000
p value    0.02403926
Cohen's d  1.26491106
> metric_t(x,y,paired=TRUE)
                 [,1]
Diff M    4.000000000
t(dep.)   4.512608599
df        7.000000000
p value   0.002756406
Cohen's d 1.264911064
}}
