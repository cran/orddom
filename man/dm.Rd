\name{dm}
\alias{dm}
\title{Dominance Matrix Creation}
\description{Returns a dominance matrix based on the comparison of all values of two 1-column matrices x and y}
\usage{dm(x, y)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{1 column matrix with \eqn{n_1} values (e.g. from group X)}
  \item{y}{1 column matrix with \eqn{n_2} values (e.g. from group Y)}
}
\details{
Each dominance matrix' cell value \eqn{d_{ij}} is calculated as \eqn{sign(y_j-x_i)} across all \eqn{i=1,2,3,...,n_1} values (=rows) of x and \eqn{i=1,2,3,...,n_2} values (=rows) of y.}
\value{Returns dominance matrix with X values as rownames and with Y values as columnnames}
\references{Cliff, N. (1996). \emph{Ordinal Methods for Behavioral Data Analysis}. Mahwah, NJ: Lawrence Erlbaum.}
\author{Jens Rogmann}
\seealso{dms}
\examples{\dontrun{
> x<-t(matrix(c(1,1,2,2,2,3,3,3,4,5),1))
> y<-t(matrix(c(1,2,3,4,4,5),1))
> dm(x,y)
  1  2  3  4  4  5
1 0 -1 -1 -1 -1 -1
1 0 -1 -1 -1 -1 -1
2 1  0 -1 -1 -1 -1
2 1  0 -1 -1 -1 -1
2 1  0 -1 -1 -1 -1
3 1  1  0 -1 -1 -1
3 1  1  0 -1 -1 -1
3 1  1  0 -1 -1 -1
4 1  1  1  0  0 -1
5 1  1  1  1  1  0
}}
\keyword{array}
