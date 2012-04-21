\name{return1colmatrix}
\alias{return1colmatrix}
\title{Convert vectors, data frames, lists, or arrays to 1-column matrix for use in \link{orddom}}
\description{Converts vectors, data frames, lists, and arrays to 1-column matrix with optional column name and sorting option for use in various orddom functions}
\usage{return1colmatrix(x,grp.name="",sortx=FALSE)}
\arguments{
  \item{x}{Vector, data frame, list or array with \eqn{n_x} values and an optional header or name}
  \item{grp.name}{A name or column title for x may be assigned. By default, the variable name is returned as \emph{var(x)}.}
  \item{sortx}{If argument is set to \emph{TRUE}, the function will return a matrix with sorted scores.}
}
\value{Returns a 1-column matrix with n scores in n rows with X columnname.}
\author{Jens Rogmann}
\seealso{\link{orddom}}
\keyword{array}
