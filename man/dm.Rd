\name{dm}
\alias{dm}
\title{Dominance or Difference Matrix Creation}
\description{Returns a dominance or difference matrix based on the comparison of all values of two 1-column matrices x and y}
\usage{dm(x, y, diff=FALSE)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{1 column matrix with \eqn{n_1} values (e.g. from group X)}
  \item{y}{1 column matrix with \eqn{n_2} values (e.g. from group Y)}
  \item{diff}{If argument is set to true, the function will return a difference matrix. Otherwise, a dominance matrix is produced.}
}
\details{
Each difference matrix cell value \eqn{d_{ij}} is calculated as \eqn{y_j-x_i} across all \eqn{i=1,2,3,...,n_1} values (=rows) of x and \eqn{i=1,2,3,...,n_2} values (=rows) of y.
Dominance matrix cell values are calculated as \eqn{sign(y_j-x_i)}.}
\value{Returns difference or dominance matrix with X values as rownames and with Y values as columnnames}
\references{Cliff, N. (1996). \emph{Ordinal Methods for Behavioral Data Analysis}. Mahwah, NJ: Lawrence Erlbaum.}
\author{Jens Rogmann}
\seealso{dms}
\examples{\dontrun{
> x<-t(matrix(c(1,1,2,2,2,3,3,3,4,5),1))
> y<-t(matrix(c(1,2,3,4,4,5),1))
> dm(x,y,diff=TRUE)
  1  2  3  4  4  5
1 0 -1 -2 -3 -3 -4
1 0 -1 -2 -3 -3 -4
2 1  0 -1 -2 -2 -3
2 1  0 -1 -2 -2 -3
2 1  0 -1 -2 -2 -3
3 2  1  0 -1 -1 -2
3 2  1  0 -1 -1 -2
3 2  1  0 -1 -1 -2
4 3  2  1  0  0 -1
5 4  3  2  1  1  0
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
