\name{cohd2delta}
\alias{cohd2delta}
\title{Cohen's d to Cliff's delta}
\description{Converts Cohen's d effect size to Cliff's delta as non-overlap between two standard normal distributions}
\usage{cohd2delta(d)}
\arguments{
  \item{d}{Cohen's d value}
}
\details{Returns delta (or non-overlap, see Table 2.2.1 in Cohen, 1988, p.22).}
\value{\deqn{\delta(d) = \frac{2AUC(\frac{d}{2})-1}{AUC(\frac{d}{2})}},
where \eqn{AUC(x) = \frac{1}{\sqrt{2\pi}} \int_{-\infty}^x e^{-t^2/2} \, dt}}
\references{
Cohen, J. (1988). \emph{Statistical Power Analysis for the Behavioral Sciences} (2nd ed.). Hillsdale, NJ, USA: Lawrence Erlbaum Associates.}

\author{Jens Rogmann}
\seealso{delta2cohd}
\examples{\dontrun{> cohd2delta(1.1)
[1] 0.589245
> cohd2delta(2.1)
[1] 0.8278607
> cohd2delta(2.2)
[1] 0.8430398
> cohd2delta(4.0)
[1] 0.9767203}}
\keyword{distribution}
