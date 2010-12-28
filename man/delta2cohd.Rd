\name{delta2cohd}
\alias{delta2cohd}
\title{Cliff's delta to Cohen's d}
\description{Converts Cliff's delta estimate to Cohen's d effect size as non-overlap between two standard normal distributions}
\usage{delta2cohd(d)}
\arguments{\item{d}{Cliff's delta estimate \eqn{\delta}.}}
\details{Returns Cohen's d (or non-overlap, based on U1 in Table 2.2.1, Cohen, 1988, p.22).}
\value{\eqn{d(\delta) = 2 z_{\frac{-1}{\delta-2}}}, where \eqn{z_p \equiv \Phi^{-1}(p) = AUC^{-1}(p)}}
\references{Cohen, J. (1988). \emph{Statistical Power Analysis for the Behavioral Sciences} (2nd ed.). Hillsdale, NJ, USA: Lawrence Erlbaum Associates.}
\author{Jens Rogmann}
\seealso{cohd2delta}
\examples{\dontrun{> delta2cohd(-.10)
[1] -0.1194342
> delta2cohd(-.86)
[1] -0.7725292
> delta2cohd(.10)
[1] 0.1320236
> delta2cohd(.774)
[1] 1.797902
}}
\keyword{distribution}
