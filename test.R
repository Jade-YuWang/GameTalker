
a1 <- read.csv("FT-LN18.csv")[1:11,-1]
a2 <- read.csv("ST-LN18.csv")[1:11,-1]



b1 <- read.csv("FT-3T3.csv")[1:11,-1]
b2 <- read.csv("ST-3T3.csv")[1:11,-1]


a <- (a1 + a2)/2
b <- (b1 + b2)/2

aa1 <-  matrix(as.numeric(unlist(c(a))),nrow=11)
aa2 <-  matrix(rep(as.numeric(a[1,]),11),nrow=11,byrow=T)
aa <- aa1 -aa2+1

bb1 <-  matrix(as.numeric(unlist(c(b))),nrow=11)
bb2 <-  matrix(rep(as.numeric(b[1,]),11),nrow=11,byrow=T)
bb <- bb1 -bb2+1

par.curve <- c(6,0.32)
rr <- optim_BFGS.ind(par.curve,max.loop=5,z=as.numeric(aa[,6]),times=0:10)

plot(1:11,aa[,6],pch=16)
lines(1:11,LC_get_mu1(rr$par,0:10),lwd=2,col="darkred")

dat1 <- list(x=aa[,1],y=aa[,2],z=aa[,3],xc=aa[,4],yc=aa[,5],zc=aa[,6],times=1:11)

dat2 <- list(x=bb[,1],y=bb[,2],z=bb[,3],xc=bb[,4],yc=bb[,5],zc=bb[,6],times=1:11)


