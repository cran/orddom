\name{orddom_f}
\alias{orddom_f}
\title{Ordinal Dominance Statistics: File output of statistics for multiple comparisons}
\description{Writes ordinal dominance statistics to tailored target output file, e.g. for purposes of multiple comparisons.}
\usage{orddom_f(x,y,alpha=0.05,paired=FALSE,outputfile="orddom_csv.txt",studdist=TRUE,
       symmetric=FALSE,quotechar=TRUE,decimalpt=".",separator="\t",
       notavailable="NA",endofline="\n",description="NA")}
\arguments{
  \item{x}{A 1-column matrix with optional column name containing all \eqn{n_x} values or scores of group X or 1 (e.g. control or pretest group.); see \code{\link{orddom}} for details.}
  \item{y}{A 1-column matrix with optional column name containing all \eqn{n_y} values of group Y or 2 (e.g. treatment or post-test group); see \code{\link{orddom}} for details.}
  \item{alpha}{Significance or \eqn{\alpha}-level used for the calculation of the confidence intervals; see \code{\link{orddom}} for details.}
  \item{paired}{By default, independence of the two groups or data sets is assumed. For paired comparisons, set to TRUE; see \code{\link{orddom}} for details.}
  \item{outputfile}{A filename for the report should be given here. The report as standard text file is written to the current working directory. All data are appended to this file. 
  If the file does not exist initially, row headers are produced.}
  \item{studdist}{By default, z-values are based on Student's t-distribution. Set to FALSE if z-values should be based on the Standard Normal Distribution; see \code{\link{orddom}} for details.}
  \item{symmetric}{By default, asymmetric confidence intervals (CI) are being calculated for delta statistics. Set to TRUE to have symmetric CIs returned; see \code{\link{orddom}} for details.}
  \item{quotechar}{By default, string outputs are quoted.}
  \item{decimalpt}{By default, numeric outputs use the colon as decimal point. Where commas are used instead, this argument should be set to \emph{...,decimalpt=",",...} .}
  \item{separator}{By default, field entries are separated by tabulators (\emph{...,separator="\\tab",...}). 
  If, for example, .csv files are to be produced using the semicolon as the field separator, this argument should be set to \emph{...,separator=";",...}.}
  \item{notavailable}{By default, if field entries are ot available, \emph{"NA"} is printed to the file. 
  Other values to be printed can be given, e.g. \emph{...,notavailable="",...} or \emph{...,notavailable="NULL",...}.} 
  \item{endofline}{By default, a carriage return denotes the end of the single output line. 
  Other values may be given, such as the IETF standard for csv files (\emph{...,endofline="\\r\\n",...}).}
  \item{description}{A string description can be given to give even more details about the very comparison and statistics written to the file. 
  This information is then also written to the file.\cr \cr \cr}
} 
\seealso{\code{\link{orddom}}}
\author{Jens J. Rogmann, University of Hamburg, Department of Psychology, \cr Hamburg, Germany (Jens.Rogmann@uni-hamburg.de)}
\examples{\dontrun{
# Example: Experiment with experimental group "ex" and control group "con"
# Data sets:
ex_pre<-c(52,53,55,59,57)
con_pre<-c(51,56,54,60,56)
ex_post<-c(58,62,63,64,69)
con_post<-c(48,58,57,62,55)
# Two independent and two paired comparisons are possible
# These are to be written to a csv-file
# Alpha-level = 10% / 4 comparisons (Bonferroni) = .025
orddom_f(con_pre,ex_pre,alpha=0.025,decimalpt=",",description="EXP 01: Between groups at time 01")   
# result delta=-.04
orddom_f(con_post,ex_post,alpha=0.025,decimalpt=",",description="EXP 01: Between groups at time 02") 
# result delta=.84
orddom_f(ex_pre,ex_post,alpha=0.025,paired=TRUE,decimalpt=",",description="EXP 01: Within exp 01 to 02") 
# result delta_b=.9
orddom_f(con_pre,con_post,alpha=0.025,paired=TRUE,decimalpt=",",description="EXP 01: Within con 01 to 02") 
# result delta_b=.2
file.show(file.path(getwd()),"orddom_csv.txt")
}}
\keyword{nonparametric}
\keyword{robust}
\keyword{htest}

