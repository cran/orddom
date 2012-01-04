orddom <- function (x,y,alpha=0.05,paired=FALSE,outputfile="",studdist=TRUE,symmetric=FALSE) {
library(psych)
library(compute.es)
if (paired==FALSE) {
#Intro
n_x=length(x)
n_y=length(y)
prepost <- matrix (nrow=n_x, ncol=n_y)
dom <- matrix (nrow=n_x, ncol=n_y)
di_<- matrix (nrow =n_x)
dj_<- matrix (nrow =n_y)
ki_<- matrix (nrow =n_x)
kj_<- matrix (nrow =n_y)
di<- matrix (nrow =n_x)
dj<- matrix (nrow =n_y)
da <- matrix (nrow =29,ncol=2)
rownames(da)<-c("var1_X","var2_Y","type","n in X","n in Y","N #Y>X","N #Y=X","N #Y<X","PS Y>X","PS X>Y","A X>Y","A Y>X","delta","1-alpha","CI low","CI high","s delta","var delta","se delta","z/t score","p (1-tail)","p (2-tail)","Cohen's d","d CI low","d CI high", "var d.i","var dj.","var dij","df")
colnames(da)<-c("ordinal","metric")
#Dominance Matrix creation
dom<-dm(x,y)
rownames(prepost)<-x
colnames(prepost)<-y
da[4,1:2]<-c(n_x,n_x)
da[5,1:2]<-c(n_y,n_y)
for (i in 1:n_x) {
 for (j in 1:n_y) {
  prepost[i,j] <- -c(y[j]-x[i])
  k <- (5+sign(y[j]-x[i]))
  } }
#Number of Differences Analysis
     da[6,1] <- sum(as.numeric(-dom==1))
 da[7,1] <- sum(as.numeric(-dom==0))
     da[8,1] <- sum(as.numeric(-dom==-1))
 for (i in 6:8) { 
  da[i,2]<- da[i,1]
  }
#Berechnungen Probability of Superiority/Common Language ES
#cf. Grissom, RJ & Kim, JJ (2005). Effect sizes for reseach: A broad practical approach. Mahwah, NJ, USA: LEA.
 for (i in 1:2) { 
da[9,i]<- data.matrix(da[8,i])/sum(da[6:8,i]) #(eq. 5.1 p. 98)
da[10,i]<- data.matrix(da[6,i])/sum(da[6:8,i])
#cf. Vargha, A. & Delaney, H. (2000). A critique and improvement of the CL common language effect size statistics of McGraw and Wong. Journal of Educational and Behavioral Statistics 25(2), 101-132.
da[11,i]<-(da[8,i]+(.5*da[7,i]))/sum(da[6:8,i]) #(eq. 51, p. 127)
da[12,i]<-(da[6,i]+(.5*da[7,i]))/sum(da[6:8,i])
  }
#t-test
t_t <- matrix (nrow =5,ncol=1)
rownames(t_t)<-c("Diff M", "t value", "df Welch", "p value", "Cohen's d")
t_t[1]<-mean(y)-mean(x)
if ((t_t[1]!=0)&&(da[9,2]!=1)&&(da[10,2]!=1)) {
t_t[2]<-t.test(y,x,paired=FALSE)$statistic
t_t[3]<-t.test(y,x,paired=FALSE)$parameter
t_t[4]<-t.test(y,x,paired=FALSE)$p.value
t_t[5]<-tes(t.test(y,x,paired=FALSE)$statistic,n_x,n_y)$MeanDifference[1]
} else {for (i in 2:5) {t_t[i]<-0}} #special case to avoid errors
#Berechnungen Dominance Analysis
dw<-mean(dom)
kw<-mean(prepost)
for (i in 1:n_x) {
  di_[i]<-mean(dom[i,]) 
  ki_[i]<-mean(prepost[i,]) }
for (j in 1:n_y) {
  dj_[j]<-mean(dom[,j]) 
  kj_[j]<-mean(prepost[,j]) }
s2di_ <- sum((di_-dw)^2)/(n_x-1)
s2ki_ <- sum((ki_-kw)^2)/(n_x-1)
s2dj_ <- sum((dj_-dw)^2)/(n_y-1)
s2kj_ <- sum((kj_-kw)^2)/(n_y-1)
s2dij <- sum((dom-dw)^2)/((n_x*n_y)-1) #as sample estimate according to Long et al. (2003) eqn. 67, not Cliff's (1996, p. 138) suggestion of /((n_x-1)*(n_y-1))
s2kij <- sum((prepost-kw)^2)/((n_x*n_y)-1)
# up to version 1.0 orddom: s2dw <- ((n_x*s2di_)+(n_y*s2dj_)+s2dij)/(n_x*n_y)
s2dw <- ((((n_y*n_y)*sum((di_-dw)^2))+(((n_x*n_x)*sum((dj_-dw)^2)))-(sum((dom-dw)^2)))/(n_x*n_y*(n_x-1)*(n_y-1)))
# up to version 1.0 orddom: ((n_x*s2ki_)+(n_y*s2kj_)+s2kij)/(n_x*n_y)
s2kw <- ((n_y*n_y*(n_x-1)*s2ki_)+(n_x*n_x*(n_y-1)*s2kj_)+ (((n_x*n_y)-1)*s2kij))/(n_x*n_y*(n_x-1)*(n_y-1)) 
if(s2dw<((1-dw^2)/((n_x*n_y)-1))) {s2dw<-((1-dw^2)/((n_x*n_y)-1))} #Minimum allowable value for s2d, cf. Cliff (1996, p. 138), Feng & Cliff (2004, p. 323) or Long et al. (2003, par after eqn. 66).
dw<--mean(dom)
kw<--mean(prepost)
#df
da[29,1]<-((n_x)+(n_y)-2)
da[29,2]<-t_t[3]
#Confidence Intervals
#cf. Feng & Cliff (2004), Journal of Modern Applied Statistical Methods, 3(2), 322-332 and
#cf. Feng (2007), in:  Shlomo S. Sawilowsky (Ed.), Real Data Analysis (pp. 163-183).
#t_level<-qnorm(1-(alpha/2)) when approx NormDistrib t_level=1.96 as in Du Feng's Fortran Code
#here: Derived fm t-test: qt((1-(alpha/2)),df)) 
if (studdist==TRUE) {t_level<-qt((1-(alpha/2)),da[29,1])} #Student instead of Normal Distribution
else {t_level<-qnorm(1-(alpha/2)) #Normal Distribution instead of Student's t
da[29,1]<-NaN
da[29,2]<-t_t[3]}
if (symmetric==TRUE) {#symmetric CI
ci_dw_lo<-(dw-(dw^3)-(t_level*sqrt(s2dw)*sqrt(1-(2*dw*dw)+(dw^4)+(t_level*t_level*s2dw)))) / (1-(dw*dw)+(t_level*t_level*s2dw))
ci_dw_hi<-(dw-(dw^3)+(t_level*sqrt(s2dw)*sqrt(1-(2*dw*dw)+(dw^4)+(t_level*t_level*s2dw)))) / (1-(dw*dw)+(t_level*t_level*s2dw))
} else {#asymmetric CI
ci_dw_lo<-(dw-(dw^3)-(t_level*sqrt(s2dw)*sqrt(1-(2*dw*dw)+(dw^4)+(t_level*t_level*s2dw)))) / (1-(dw*dw)+(t_level*t_level*s2dw))
ci_dw_hi<-(dw-(dw^3)+(t_level*sqrt(s2dw)*sqrt(1-(2*dw*dw)+(dw^4)+(t_level*t_level*s2dw)))) / (1-(dw*dw)+(t_level*t_level*s2dw))
}
#special case Feng & Cliff (2004), eq. 8, p.324: d=+-1
if ((dw*dw)==1) {n_b<-n_x
 if (n_y<n_x){n_b<-n_y}}
if (dw==1)  {ci_dw_lo <-((n_x-(t_level*t_level))/(n_x+(t_level*t_level)))
             if (s2dw==0){ci_dw_hi<-Inf}}
if (dw==-1) {ci_dw_hi <--((n_x-(t_level*t_level))/(n_x+(t_level*t_level)))
             if (s2dw==0){ci_dw_lo<--Inf}}
if (ci_dw_lo < -1) {ci_dw_lo<--1}
if (ci_dw_hi > 1) {ci_dw_hi<-1}
#z and p levels
zdw<-dw/sqrt(s2dw)
if (studdist==TRUE) {pdw<-1-pt((sign(zdw)*zdw),da[29,1])} #Student 
 else {pdw<-1-pnorm((sign(zdw)*zdw))} #Normal Distribution
#Ausgabe Intro
da[13,1]<-dw #delta
da[13,2]<--mean(prepost)
da[14,1]<-(1-alpha)*100
da[14,2]<-da[14,1]
da[15,1]<-ci_dw_lo # CI low
da[16,1]<-ci_dw_hi # CI high
da[17,1]<-sqrt(s2dw) # s delta
da[17,2]<-sqrt((((n_x-1)*var(x))+((n_y-1)*var(y)))/(n_x-1+n_y-1)) #Pooled SD M Diff
da[19,2]<-da[17,2]*sqrt((1/n_x)+(1/n_y)) #STD Error M Diff
da[15,2]<--mean(prepost)-((qt((1-(alpha/2)),da[29,2])*da[19,2]))
da[16,2]<--mean(prepost)+((qt((1-(alpha/2)),da[29,2])*da[19,2]))
da[18,1]<-s2dw # var delta
da[18,2]<-da[17,2]^2
da[20,1]<-zdw # z value
da[20,2]<-(da[13,2]/da[19,2])
da[21,1]<-pdw # p value
da[21,2]<-t_t[4]/2
da[22,1]<-pdw*2
da[22,2]<-t_t[4]
da[23,1]<-delta2cohd(dw) #Cohen's d as overlap
da[23,2]<-t_t[5]
da[24,1]<-delta2cohd(ci_dw_lo) #Cohen's d as overlap
da[24,2]<-tes(da[15,2],n_x,n_y)$MeanDifference[1]
da[25,1]<-delta2cohd(ci_dw_hi) #Cohen's d as overlap
da[25,2]<-tes(da[16,2],n_x,n_y)$MeanDifference[1]
da[26,1]<-s2di_
da[26,2]<-s2ki_
da[27,1]<-s2dj_
da[27,2]<-s2kj_
da[28,1]<-s2dij
da[28,2]<-s2kij
for (i in 1:2){
da[1,i]<-colnames(x)[1]
da[2,i]<-colnames(y)[1]
da[3,i]<-c("indep")
}
#Ausgabe
    if(outputfile!="") {
    err<-getOption("show.error.messages")
options("show.error.messages"=FALSE)
sink(outputfile, append=FALSE, type="output", split=FALSE)
 cat("\n\n=== ORDINAL DOMINANCE ANALYSIS FOR INDEPENDENT GROUPS ===\n")
 cat("1. Raw Data(X-Grp1: ",append=TRUE)
 cat(colnames(x),append=TRUE)
 cat(")\n",append=TRUE)
     write.table(data.frame(x),,append=TRUE,quote=FALSE,sep="\t",,,,,col.names=FALSE) #Ausgabe Rohdaten Gruppe 1
 cat("           (Y-Grp2:",append=TRUE)
 cat(colnames(y),append=TRUE)
 cat(")\n",append=TRUE)
     write.table(data.frame(y),,append=TRUE,quote=FALSE,sep="\t",,,,,col.names=FALSE) #Ausgabe Rohdaten Gruppe 2
 cat("\n2. Descriptives\n")
 cat("\t X (Grp1)\n",append=TRUE)
 write.table(format(t(describe(x)),digits=4,width=12),append=TRUE,quote=FALSE,sep="\t",col.names=FALSE) #Ausgabe Descriptives
 cat("\t Y (Grp2)\n",append=TRUE)
 write.table(format(t(describe(y)),digits=4,width=12),append=TRUE,quote=FALSE,sep="\t",col.names=FALSE) #Ausgabe Descriptives
 cat("\n\n3. t-Test Outcome\n",append=TRUE)
 write.table(t_t,append=TRUE,quote=FALSE,sep="\t",col.names=FALSE) #Ausgabe Paired t-Test
 cat("\n4. Difference Matrix (Rows:X-Grp1/ Cols:Y-Grp2)\nX / Y  ",append=TRUE)
 cat(colnames(prepost),append=TRUE)
 cat("\n")
 write.table(format(-prepost,width=2),quote=FALSE,append=TRUE,col.names=FALSE) #Ausgabe Differenzmatrix (Zeilen=Pr<e4>/X, Spalten=Post/Y)
 cat("\n4. Dominance Matrix (Rows:X-Grp1/ Cols:Y-Grp2)\n   (+ indicate higher values in the Y group)\n",append=TRUE)
 write.table(dms(-dom,FALSE),quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE,sep="") #Ausgabe Dominanzmatrix (Zeilen=Pr<e4>/X, Spalten=Post/Y)
 cat("\n5. Dominance Analysis\n",append=TRUE)#Ausgabe DA
 write.table(da,,append=TRUE,quote=FALSE,sep="\t") 
     sink()
file.show(outputfile)
closeAllConnections() 
options("show.error.messages"=err)}
#if (xlsxfile!="") {write.xlsx(data.frame(da), file_path_as_absolute(file.path("tmp_r.txt")), "depCliffsDelta", NULL, TRUE, TRUE, file.access(xlsxfile,0))}
return (da)
} else
{#Paired Data **************************************************************************************************************************************
n_x=length(x)
n_y=length(y)
if ((n_x==n_y)&&(n_x>3)) {
    prepost <- matrix (nrow=n_x, ncol=n_y)
dom <- matrix (nrow=n_x, ncol=n_y)
di_<- matrix (nrow =n_x)
dj_<- matrix (nrow =n_y)
di<- matrix (nrow =n_x)
dj<- matrix (nrow =n_y)
da <- matrix (nrow =29,ncol=4)
rownames(da)<-c("var1_X_pre","var2_Y_post","type","N #Y>X","N #Y=X","N #Y<X","PS X>Y","PS Y>X","A X>Y","A Y>X","delta","1-alpha","CI low","CI high","s delta","var delta","z/t score","p (1-tail)","p (2-tail)","Cohen's d","d CI low","d CI high","var d.i","var dj.","cov(di,dj)","var dij","cov(dih,dhi)","cov(db,dw)","df")
colnames(da)<-c("within","between","combined","metric")
#Dominance Matrix creation
dom<-dm(x,y)
rownames(prepost)<-x
colnames(prepost)<-y
for (i in 1:n_x) {
 for (j in 1:n_y) {
  prepost[i,j] <- -c(y[j]-x[i])
  k <- (5+sign(y[j]-x[i]))
  } }
#Number of Differences Analysis
     da[4,1] <- sum(as.numeric(diag(-dom)==1))
 da[5,1] <- sum(as.numeric(diag(-dom)==0))
     da[6,1] <- sum(as.numeric(diag(-dom)==-1))
     da[4,3] <- sum(as.numeric(-dom==1))
 da[5,3] <- sum(as.numeric(-dom==0))
     da[6,3] <- sum(as.numeric(-dom==-1))
 for (i in 4:6) { 
  da[i,2]<- sum(c(sum(as.numeric(-dom==(5-i))),-sum(as.numeric(diag(-dom)==(5-i)))))
      da[i,4]<- da[i,3]   
  }
#Berechnungen Probability of Superiority/Common Language ES
#cf. Grissom, RJ & Kim, JJ (2005). Effect sizes for reseach: A broad practical approach. Mahwah, NJ, USA: LEA.
 for (i in 1:4) { 
da[7,i]<- data.matrix(da[6,i])/sum(da[4:6,i]) #(eq. 5.9 p. 115)
da[8,i]<- data.matrix(da[4,i])/sum(da[4:6,i])
#cf. Vargha, A. & Delaney, H. (2000). A critique and improvement of the CL common language effect size statistics of McGraw and Wong. Journal of Educational and Behavioral Statistics 25(2), 101<96>132.
da[9,i]<-(da[6,i]+(.5*da[5,i]))/sum(da[4:6,i]) #(eq. 51, p. 127)
da[10,i]<-(da[4,i]+(.5*da[5,i]))/sum(da[4:6,i])
  }
#t-test
t_t <- matrix (nrow =5,ncol=1)
rownames(t_t)<-c("Diff M", "t (dep.)", "deg fdm df", "p value", "Cohen's d")
if ((da[7,4]!=1)&&(da[8,4]!=1)) {
t_t[1]<-t.test(y,x,paired=TRUE)$estimate
t_t[2]<-t.test(y,x,paired=TRUE)$statistic
t_t[3]<-t.test(y,x,paired=TRUE)$parameter
t_t[4]<-t.test(y,x,paired=TRUE)$p.value
t_t[5]<-tes(t.test(y,x,paired=TRUE)$statistic,n_x,n_y)$MeanDifference[1]
} else {for (i in 2:5) {t_t[i]<-0}} #special case to avoid errors
#Berechnungen Dominance Analysis
dw<--mean(diag(dom))
s2dw<-sum((diag(dom)--dw)^2)/((n_x-1)*n_x) #correct! (Personal communication with Du Feng 07-Feb-2011: /(N-1) (as in Eq 6.8 in Cliff (1993), p. 168 or Eq 4. in Feng(2007), p.165) are errors!)
db<--((sum(dom)-sum(diag(dom)))/(n_x*(n_y-1)))
for (i in 1:n_x) {
  di_[i]<--((sum(dom[i,])-diag(dom)[i]))/(n_x-1)
  di[i]<--((sum(dom[i,])-diag(dom)[i]))
  dj_[i]<--((sum(dom[,i])-diag(dom)[i]))/(n_y-1)
  dj[i]<--((sum(dom[,i])-diag(dom)[i])) 
  }
s2di_ <- sum((di_-db)^2)/(n_x-1)
s2dj_ <- sum((dj_-db)^2)/(n_y-1)
covdidj <- sum((di_-db)*(dj_-db))/(n_x-1)
s2dij <- (sum((-dom-db)^2)-sum((diag(-dom)-db)^2))/(n_x*(n_y-1))
covdihhi <- (sum((-dom-db)*(-t(dom)-db))-sum((diag(-dom)-db)^2))/(n_x*(n_y-1))
s2db <- ((((n_x-1)*(n_y-1))*(sum((di_-db)^2)+sum((dj_-db)^2)+(2*(sum((di_-db)*(dj_-db))))))-(sum((-dom-db)^2)-sum((diag(-dom)-db)^2))-(sum((-dom-db)*(-t(dom)-db))-sum((diag(-dom)-db)^2)))/(n_x*(n_y-1)*(n_x-2)*(n_y-3))
dwdb <- dw+db
covdwdb <- (sum(((di+dj)*diag(-dom))) -(2*n_x*(n_y-1)*db*dw))/((n_x*(n_y-1)*(n_x-2)))#Fengs fortran: /(N-1)^2*(N-2) is incorrect!
sdwdb <- sqrt(s2dw + s2db + (2*covdwdb))
#df
da[29,1]<-((n_x-1))#df dw
da[29,2]<-da[29,1] #df db
da[29,3]<-da[29,1] #df dw+db
da[29,4]<-t_t[3] #df metric t-test 
#Confidence Intervals
#cf. Feng & Cliff (2004), Journal of Modern Applied Statistical Methods, 3(2), 322-332 and
#cf. Feng (2007), in:  Shlomo S. Sawilowsky (Ed.), Real Data Analysis (pp. 163-183).
#t_level<-qnorm(1-(alpha/2)) approx NormDistrib t_level=1.96 as in Du Feng's Fortran Code
#Here: Derived fm t-test: qt((1-(alpha/2)),df) :
if (studdist==TRUE) {#Student instead of Normal Distribution
t_level1<-qt((1-(alpha/2)),da[29,1])
t_level2<-qt((1-(alpha/2)),da[29,2])
t_level3<-qt((1-(alpha/2)),da[29,3])} 
else {#Normal Distribution instead of Student's t
t_level1<-qnorm(1-(alpha/2))
t_level2<-qnorm(1-(alpha/2))
t_level3<-qnorm(1-(alpha/2))
da[29,1]<-NaN
da[29,2]<-NaN
da[29,3]<-NaN
}
if (symmetric==TRUE) {#symmetric CI 
ci_dw_lo<-dw-(t_level1*sqrt(s2dw))
ci_dw_hi<-dw+(t_level1*sqrt(s2dw))
ci_db_lo<-(db-(t_level2*sqrt(s2db)))
ci_db_hi<-(db+(t_level2*sqrt(s2db)))
ci_dwdb_lo<-(dwdb-(t_level3*sdwdb))
ci_dwdb_hi<-(dwdb+(t_level3*sdwdb))
} else {#asymmetric CI
ci_dw_lo<-(dw-(dw^3)-(t_level1*sqrt(s2dw)*sqrt(1-(2*dw*dw)+(dw^4)+(t_level1*t_level1*s2dw)))) / (1-(dw*dw)+(t_level1*t_level1*s2dw))
ci_dw_hi<-(dw-(dw^3)+(t_level1*sqrt(s2dw)*sqrt(1-(2*dw*dw)+(dw^4)+(t_level1*t_level1*s2dw)))) / (1-(dw*dw)+(t_level1*t_level1*s2dw))
ci_db_lo<-(db-(db^3)-(t_level2*sqrt(s2db)*sqrt(1-(2*db*db)+(db^4)+(t_level2*t_level2*s2db)))) / (1-(db*db)+(t_level2*t_level2*s2db))
ci_db_hi<-(db-(db^3)+(t_level2*sqrt(s2db)*sqrt(1-(2*db*db)+(db^4)+(t_level2*t_level2*s2db)))) / (1-(db*db)+(t_level2*t_level2*s2db))
ci_dwdb_lo<-(dwdb-(dwdb^3)-(t_level3*sdwdb*sqrt(1-(2*dwdb*dwdb)+(dwdb^4)+(t_level3*t_level3*sdwdb*sdwdb)))) / (1-(dwdb*dwdb)+(t_level3*t_level3*sdwdb*sdwdb))
ci_dwdb_hi<-(dwdb-(dwdb^3)+(t_level3*sdwdb*sqrt(1-(2*dwdb*dwdb)+(dwdb^4)+(t_level3*t_level3*sdwdb*sdwdb)))) / (1-(dwdb*dwdb)+(t_level3*t_level3*sdwdb*sdwdb))
}
if (db==1) {ci_db_lo <-((n_x-(t_level2*t_level2))/(n_x+(t_level2*t_level2)))
            if (s2db==0){ci_db_hi<-Inf}}
if (db==-1) {ci_db_hi <--((n_x-(t_level2*t_level2))/(n_x+(t_level2*t_level2)))
            if (s2db==0){ci_db_lo<--Inf}}
if (dw==1) {ci_dw_lo <-((n_x-(t_level1*t_level1))/(n_x+(t_level1*t_level1)))
            if (s2dw==0){ci_dw_hi<-Inf}}
if (dw==-1) {ci_dw_hi <--((n_x-(t_level1*t_level1))/(n_x+(t_level1*t_level1)))
            if (s2dw==0){ci_dw_lo<--Inf}}
if (dwdb==1){if (sdwdb==0){ci_dwdb_hi<-Inf}}
if (dwdb==-1){if (sdwdb==0){ci_dwdb_lo<--Inf}}
if (ci_dw_lo < -1) {ci_dw_lo<--1}
if (ci_dw_hi > 1) {ci_dw_hi<-1}
if (ci_db_lo < -1) {ci_db_lo<--1}
if (ci_db_hi > 1) {ci_db_hi<-1}
#z and p levels
zdw<-dw/sqrt(s2dw)
zdb<-db/sqrt(s2db)
zdwdb<-dwdb/sdwdb
if (studdist==TRUE) {#Student instead of Normal Distribution
pdw<-1-pt((sign(zdw)*zdw),da[29,1]) #instead of Normal Distribution 1-pnorm((sign(zdw)*zdw))
pdb<-1-pt((sign(zdb)*zdb),da[29,2])#instead of Normal Distribution 1-pnorm((sign(zdb)*zdb))
pdwdb<-1-pt((sign(zdwdb)*zdwdb),da[29,3])}
else {#Normal Distribution
pdw<-1-pnorm((sign(zdw)*zdw)) 
pdb<-1-pnorm((sign(zdb)*zdb))
pdwdb<-1-pnorm((sign(zdwdb)*zdwdb))}
#Ausgabe Intro
da[11,1]<-dw #delta
da[11,2]<-db
da[11,3]<-dwdb
da[11,4]<--mean(diag(prepost))
da[12,1]<-(1-alpha)*100
da[12,2]<-da[12,1]
da[12,3]<-da[12,1]
da[12,4]<-da[12,1]
da[13,1]<-ci_dw_lo # CI 95% low
da[13,2]<-ci_db_lo
da[13,3]<-ci_dwdb_lo
da[13,4]<--mean(diag(prepost))-((t_level1*sqrt(var(diag(prepost))))/sqrt(n_x))
da[14,1]<-ci_dw_hi # CI 95% high
da[14,2]<-ci_db_hi
da[14,3]<-ci_dwdb_hi
da[14,4]<--mean(diag(prepost))+((t_level1*sqrt(var(diag(prepost))))/sqrt(n_x))
da[15,1]<-sqrt(s2dw) # s delta
da[15,2]<-sqrt(s2db)
da[15,3]<-sdwdb
da[15,4]<-sqrt(var(diag(prepost)))
da[16,1]<-s2dw # var delta
da[16,2]<-s2db
da[16,3]<-(sdwdb*sdwdb)
da[16,4]<-var(diag(prepost))
da[17,1]<-zdw # z value
da[17,2]<-zdb
da[17,3]<-zdwdb
da[17,4]<--mean(diag(prepost))/sqrt(var(diag(prepost)))*sqrt(n_x)
da[18,1]<-pdw # p value
da[18,2]<-pdb
da[18,3]<-pdwdb
da[18,4]<-t_t[4]/2
da[19,1]<-pdw*2
da[19,2]<-pdb*2
da[19,3]<-pdwdb*2
da[19,4]<-t_t[4]
da[20,1]<-delta2cohd(dw) #Cohen's d as overlap
da[20,2]<-delta2cohd(db)
da[20,4]<-t_t[5]
da[21,1]<-delta2cohd(ci_dw_lo) #Cohen's d as overlap
da[21,2]<-delta2cohd(ci_db_lo)
da[21,4]<-tes(da[13,4],n_x,n_y)$MeanDifference[1]
da[22,1]<-delta2cohd(ci_dw_hi) #Cohen's d as overlap
da[22,2]<-delta2cohd(ci_db_hi)
da[22,4]<-tes(da[14,4],n_x,n_y)$MeanDifference[1]
da[23,3]<-s2di_
da[24,3]<-s2dj_
da[25,3]<-covdidj
da[26,3]<-s2dij
da[27,3]<-covdihhi
da[28,3]<-covdwdb
for (i in 1:4){
da[1,i]<-colnames(x)[1]
da[2,i]<-colnames(y)[1]
da[3,i]<-c("paired")
}
#Ausgabe
   if(outputfile!="") {
    err<-getOption("show.error.messages")
options("show.error.messages"=FALSE)
sink(outputfile, append=FALSE, type="output", split=FALSE)
cat("\n\n=== ORDINAL DOMINANCE ANALYSIS FOR PAIRED DATA ===\n")
cat("1. Raw Data(X-Pre/Y-Post)\n",append=TRUE)
write.table(data.frame(x,y),,append=TRUE,quote=FALSE,sep="\t") #Ausgabe Rohdaten
cat("\n2. Descriptives\n      \t X (Pre)    \t Y (Post)\n",append=TRUE)
write.table(format(t(describe(cbind(x,y))),digits=4,width=12),append=TRUE,sep="\t",quote=FALSE,col.names=FALSE) #Ausgabe Descriptives
cat("\n3. Paired T-Test Outcome\n",append=TRUE)
write.table(t_t,append=TRUE,sep="\t",quote=FALSE,col.names=FALSE) #Ausgabe Paired t-Test
cat("\n4. Difference Matrix (Rows:X-Pre/ Cols:Y-Post)\nX / Y  ",append=TRUE)
write.table(format(-prepost,width=2),,append=TRUE,quote=FALSE,sep="\t") #Ausgabe Differenzmatrix (Zeilen=Pr<e4>/X, Spalten=Post/Y)
cat("\n4. Dominance Matrix (Rows:X-Pre/ Cols:Y-Post)\n   (+/> indicate higher values in the Y group)\n",append=TRUE)
write.table(dms(-dom,paired=TRUE),quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE,sep="") #Ausgabe Dominanzmatrix (Zeilen=Pr<e4>/X, Spalten=Post/Y)
cat("\n5. Dominance Analysis\n",append=TRUE)#Ausgabe DA
write.table(da,,append=TRUE,quote=FALSE,sep="\t",)
    sink()
file.show(outputfile)
closeAllConnections() 
options("show.error.messages"=err)}
#if (xlsxfile!="") {write.xlsx(data.frame(da), file_path_as_absolute(file.path("tmp_r.txt")), "depCliffsDelta", NULL, TRUE, TRUE, file.access(xlsxfile,0))}
return (da)
} else {return (warnings("Error. Unequal number of cases or number of pairs <3."))}}}
