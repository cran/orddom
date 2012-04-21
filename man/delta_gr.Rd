\name{delta_gr}
\alias{delta_gr}
\title{Cliff's delta Graphics and Interpretation}
\description{Returns a graphical representation and interpretation of Cliff's delta} 
\usage{delta_gr(x,y, ... ,dv=2)}
\arguments{
 \item{x}{A 1-column matrix with optional column name containing all \eqn{n_x} values or scores of group X or 1 (e.g. control or pretest group.). \cr }
 \item{y}{A 1-column matrix with optional column name containing all \eqn{n_y} values of group Y or 2 (e.g. experimental or post-test group). For paired comparisons (e.g. pre-post), \eqn{n_x = n_y} is required. \cr See \code{\link{orddom}} for details. \cr }
 \item{...}{Other arguments to be passed on to the \code{\link{orddom}} function, such as (for example):\cr
  - \emph{paired}: to compare dependent data (e.g. pre-post) set to \emph{paired=TRUE} ,\cr
  - \emph{alpha} for the respective significance level to be used, e.g.\emph{alpha=.01} for 1%;\cr
  - \emph{onetailed} to generate one-sided testing p and confidence interval (CI) values set to \emph{onetailed=TRUE} ,\cr
  - \emph{studdist} to obtain  CI based on normal distribution z values (instead of Student distribution t) set to \emph{studdist=FALSE} for 1%;\cr
  - \emph{symmetric} to obtain symmetric rather than asymmetric CIs (see \code{\link{orddom}} for details) set to \emph{symmetric=TRUE} .\cr 
  - \emph{onetailed} for one-sided rather than the default two-tailed testing.\cr
  - \emph{x.name} to assign an individual label to group x (i.e. 1st or control or pretest group). \cr
  - \emph{y.name} to label the y input matrix or group y (i.e. 2nd or experimental or posttest group). \cr
  - \emph{description} This argument allows for assigning a string (as title or description) for the ordinal comparison outputs. \cr \cr}
 \item{dv}{(For paired comparisons (\emph{dv=3}) only.) Determines which ordinal \eqn{\delta} statistics are to be returned. Set to: \cr
  - \emph{dv=1} [within] to return an analysis for the \eqn{n_x=n_y} within-pair changes,\cr 
  - \emph{dv=2} [between] to return an analysis for the overall distribution changes, based on all \eqn{n^2-n = n(n-1)} score comparisons between y and x where \eqn{i \ne j},\cr
  - \emph{dv=3} [combined] to return an analysis for the combined inference \eqn{d_w + d_b}. It is advisable to use \emph{dv=3} in combination with \emph{symmetric=TRUE}.\cr
}}
\value{Returns a graphical representation and text interpretation of Cliff's delta.}
\author{Jens Rogmann}
\seealso{\code{\link{orddom}}}
\examples{\dontrun{
#Paired comparison combined inference (Data taken from Long et al. (2003), Table 4)
x2<-t(matrix(c(2,6,6,7,7,8,8,9,9,9,10,10,10,11,11,12,13,14,15,16),1))
colnames(x2)<-c("Incidental")
y2<-t(matrix(c(4,11,8,9,10,11,11,5,14,12,13,10,14,16,14,13,15,15,16,10),1))
colnames(y2)<-c("Intentional")
delta_gr(x2,y2,paired=TRUE,studdist=FALSE,dv=3)
#
#Journal of Statistics Education Dataset: Oral Contraceptive Drug Interaction Study
#Journal of Statistics Education, Volume 12, Number 1 (March 2004).
columns<-c("SubjectNo","Seq","Period","Treatment","EEAUC","EECmax","NETAUC","NETCmax")
data<-read.table("http://www.amstat.org/publications/jse/datasets/ocdrug.dat",col.names=columns)
#
#returns delta (between) and 95% confidence interval 1-sided
x<-subset(data,data$Treatment==0)[6] #Placebo EECmax
colnames(x)<-"Placebo Phase"
y<-subset(data,data$Treatment==1)[6] #Treatment EECmax
colnames(y)<-"Treatment Phase"
delta_gr(x,y,paired=TRUE,onetailed=TRUE,dv=2)
#
#checks treatment groups delta equivalence in placebo phase
#returns delta and 95% confidence interval
plac<-subset(data,data$Treatment==0)
x<-subset(plac,plac$Period==1)[6] #control (placebo before drug)
colnames(x)<-"Control (before Drug)"
y<-subset(plac,plac$Period==2)[6] #experimental (placebo after drug)
colnames(y)<-"Exp (Placebo after Drug)"
delta_gr(x,y)
#
}}
\keyword{distribution}
