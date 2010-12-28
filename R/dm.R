dm <-
function (x,y) { #produces dominance matrix
n_x=length(x)
n_y=length(y)
dx <- matrix(nrow=n_x, ncol=n_y)
for (i in 1:n_x)
  {for (j in 1:n_y) 
  {dx[i,j]<--sign(y[j]-x[i])}
  } 
 rownames(dx)<-x
colnames(dx)<-y
return(dx)}

