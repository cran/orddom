\name{dms}
\alias{dms}
\title{Dominance Matrix in Symbols}
\description{Returns a character-based dominance matrix based on the signs of all cell values of a given matrix}
\usage{dms(dom, paired = FALSE)}
\arguments{
  \item{dom}{Input matrix, typically raw difference or dominance matrix}
  \item{paired}{Should only be set to TRUE if the number of rows equal the number of columns and if the difference data in the matrix diagonal are to be given different symbols.}}
\details{According to the sign of each input matrix' cell value (\eqn{sign(d_{ij})}), a respective symbol is written to the output matrix ("-" for -1, "O" for 0" and "+" for 1).\cr
If paired==TRUE, the diagonal vector of the output matrix receives different symbols (i.e. "<" for -1, "=" for 0, "<" for 1).}
\references{Cliff, N. (1996). \emph{Ordinal Methods for Behavioral Data Analysis}. Mahwah, NJ: Lawrence Erlbaum.}
\author{Jens Rogmann}
\seealso{dm}
\examples{\dontrun{> x<-t(matrix(c(1,1,2,2,2,3,3,3,4,5),1))
> y<-t(matrix(c(1,2,3,4,4,5),1))
> write.table(dms(dm(x,y)),quote=FALSE,row.names=FALSE,col.names=FALSE,sep="")
O-----
O-----
+O----
+O----
+O----
++O---
++O---
++O---
+++OO-
+++++O}}
\keyword{array}
