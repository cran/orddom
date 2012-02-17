\name{orddom_p}
\alias{orddom_p}
\title{Ordinal Dominance Matrices and Statistics: Printer-friendly Tab-Delimited Report Output File}
\description{Generates a sectioned report file with ordinal dominance matrices and statistics.}
\usage{orddom_p(x,y,alpha=.05,paired=FALSE,sections="1234a4b5a5b",header="Y",sorted="XY",
       outfile="orddom_csv_tab.txt",appendfile=FALSE,show=1)}
\arguments{
  \item{x}{A 1-column matrix with optional column name containing all \eqn{n_x} values or scores of group X or 1 (e.g. control or pretest group.); see \code{\link{orddom}} for details.}
  \item{y}{A 1-column matrix with optional column name containing all \eqn{n_y} values of group Y or 2 (e.g. experimental or post-test group); see \code{\link{orddom}} for details.}
  \item{alpha}{Significance or \eqn{\alpha}-level used for the calculation of the confidence intervals; see \code{\link{orddom}} for details.}
  \item{paired}{By default, independence of the two groups or data sets is assumed. For paired comparisons, set to TRUE; see \code{\link{orddom}} for details.}
  \item{sections}{By default all of the following report sections are written to the file. If only a selection of all sections is needed, a string should be given containing all section numbers needed in the output, e.g. \emph{...,sections="135a",...} for sections 1, 3 and 5a.\cr
The following sections are available for output:\cr
\emph{"1"} - Raw data of the x and y data sets \cr
\emph{"2"} - Metric descriptives for x and y  \cr
\emph{"3"} - Metric difference tests \cr
\emph{"4a"} - Metric difference matrix with x in rows and y in columns \cr
\emph{"4b"} - Metric difference matrix with y in rows and x in columns \cr
\emph{"5a"} - Ordinal dominance matrix with x in rows and y in columns \cr
\emph{"5b"} - Ordinal dominance matrix with y in rows and x in columns \cr}
  \item{header}{By default, section headers are part of the output. If headers are to be omitted, this argument should be set to FALSE.}
  \item{sorted}{All outputs in sections 1,4a,4b,5a and 5b may be automatically sorted ascendingly for the x data set (string is to contain "X") and/or for the y data set (string is to contain "Y"). This is the default option.}
  \item{outfile}{A filename for the report should be given here. The report as standard text file is written to the current working directory.}
  \item{appendfile}{By default, new report files are created. If a given report file ist to be appended, set to TRUE.}
  \item{show}{By default, the generated file is displayed. Set to FALSE to avoid the resulting file to be shown.\cr \cr \cr}
} 
\seealso{
       \code{\link{orddom}}.
     }
\author{Jens J. Rogmann, University of Hamburg, Department of Psychology, \cr Hamburg, Germany (Jens.Rogmann@uni-hamburg.de)}
\examples{\dontrun{
#Independent Samples (Data taken from Long et al. (2003), Table 4}
x<-t(matrix(c(3,3,3,4,5,6,12,12,13,14,15,15,15,15,15,16,18,18,18,23,23,27,28,28,43),1))
colnames(x)<-c("Nonalcohol.")
y<-t(matrix(c(1,4,6,7,7,14,14,18,19,20,21,24,25,26,26,26,27,28,28,30,33,33,44,45,50),1))
colnames(y)<-c("Alcoholic")
orddom_p(x,y,,paired=FALSE,outfile="orddom_csv_tab.txt")}
\keyword{nonparametric}
\keyword{robust}
\keyword{htest}