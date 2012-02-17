metric_t<-function (a,b,alpha=0.05,paired=FALSE) {
#returns matrix with metric t-test parameters for output
 mat<-matrix (nrow =5,ncol=1)
 if(paired==FALSE){
 rownames(mat)<-c("Diff M", "t value", "df Welch", "p value", "Cohen's d")}
 else {
 rownames(mat)<-c("Diff M", "t(dep.)", "df", "p value", "Cohen's d")}
 if (paired==TRUE) {
 mat[1]<-t.test(b,a,paired=TRUE,conf.level=(1-alpha))$estimate
 mat[2]<-t.test(b,a,paired=TRUE,conf.level=(1-alpha))$statistic
 mat[3]<-t.test(b,a,paired=TRUE,conf.level=(1-alpha))$parameter
 mat[4]<-t.test(b,a,paired=TRUE,conf.level=(1-alpha))$p.value 
 } else {
 mat[1]<-mean(b)-mean(a)
 mat[2]<-t.test(b,a,conf.level=(1-alpha))$statistic
 mat[3]<-t.test(b,a,conf.level=(1-alpha))$parameter
 mat[4]<-t.test(b,a,conf.level=(1-alpha))$p.value }
 #Cohen's d in both cases according to Dunlop, Cortina, Vaslow, & Burke (1996)
 mat[5]<-(t.test(b,a,conf.level=(1-alpha))$statistic)*sqrt((length(b)+length(a))/(length(b)*length(a)))
 return (mat)
}
 
