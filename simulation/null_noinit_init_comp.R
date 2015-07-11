# standard normal without initialization

tmp<-proc.time()

sink(file("null_standard_normal_log_noninit.txt",open="wt"),type='message')
options(warn=1)

t<-c()

N<-1000

for (i in 1:N){
	beta1<-rnorm(1000,0,1)
	beta2<-rnorm(1000,0,1)
	s<-rep(1,1000)
	betahat1<-rnorm(1000,beta1,s)
	betahat2<-rnorm(1000,beta2,s)
	
	r12<-ash(c(betahat1,betahat2),c(s,s),mixcompdist="normal",method="shrink")
	r1<-ash(betahat1,s,mixcompdist="normal",method="shrink")
	r2<-ash(betahat2,s,mixcompdist="normal",method="shrink")
	
	t<-c(t,2*(r1$loglik+r2$loglik-r12$loglik))
}

write.table(t,file="null_standard_normal_stat_noninit.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)

tmp<-proc.time()-tmp
write(tmp,file="null_standard_normal_time_noninit.txt")

# mixture normal without initialization

tmp<-proc.time()

sink(file("null_mixture_normal_log_noninit.txt",open="wt"),type='message')
options(warn=1)

betagen<-function(n,sigma,pi,theta){
	k<-length(sigma)
	I<-as.numeric(rmultinom(1,n,c((1-pi)*theta,pi)))
	beta<-c()
	for (i in 1:k){
		beta<-c(beta,rnorm(I[i],0,sigma[i]))
	}
	beta<-c(beta,rep(0,I[k+1]))
	return(beta)
}


n<-10
t<-c()
for (i in 1:n){
M<-sample(1000:2000,1)
N<-sample(1000:2000,1)
K<-sample(1:10,1)
pi<-as.numeric(rdirichlet(1,rep(1,K)))
sigma0<-runif(K,exp(-2),exp(3))
beta1<-betagen(M,sigma0,0,pi)
beta2<-betagen(N,sigma0,0,pi)

betahat1<-rnorm(M,beta1,1)
betahat2<-rnorm(N,beta2,1)

r12<-ash(c(betahat1,betahat2),1,mixcompdist="normal",method="shrink")
r1<-ash(betahat1,1,mixcompdist="normal",method="shrink")
r2<-ash(betahat2,1,mixcompdist="normal",method="shrink")
t<-c(t,2*(r1$loglik+r2$loglik-r12$loglik))
}

write.table(t,file="null_mixture_normal_stat_noninit.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)

tmp<-proc.time()-tmp
write(tmp,file="null_mixture_normal_time_noninit.txt")

# standard normal with initialization

tmp<-proc.time()

sink(file("null_standard_normal_log.txt",open="wt"),type='message')
options(warn=1)

t<-c()

N<-1000

for (i in 1:N){
	beta1<-rnorm(1000,0,1)
	beta2<-rnorm(1000,0,1)
	s<-rep(1,1000)
	betahat1<-rnorm(1000,beta1,s)
	betahat2<-rnorm(1000,beta2,s)
	
	r12<-ash(c(betahat1,betahat2),c(s,s),mixcompdist="normal",method="shrink")
	r1<-ash(betahat1,s,mixcompdist="normal",method="shrink",g=r12$fitted.g)
	r2<-ash(betahat2,s,mixcompdist="normal",method="shrink",g=r12$fitted.g)
	
	t<-c(t,2*(r1$loglik+r2$loglik-r12$loglik))
}

write.table(t,file="null_standard_normal_stat.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)

tmp<-proc.time()-tmp
write(tmp,file="null_standard_normal_time.txt")

# mixture normal with initialization

tmp<-proc.time()

sink(file("null_mixture_normal_log.txt",open="wt"),type='message')
options(warn=1)

betagen<-function(n,sigma,pi,theta){
	k<-length(sigma)
	I<-as.numeric(rmultinom(1,n,c((1-pi)*theta,pi)))
	beta<-c()
	for (i in 1:k){
		beta<-c(beta,rnorm(I[i],0,sigma[i]))
	}
	beta<-c(beta,rep(0,I[k+1]))
	return(beta)
}


n<-10
t<-c()
for (i in 1:n){
M<-sample(1000:2000,1)
N<-sample(1000:2000,1)
K<-sample(1:10,1)
pi<-as.numeric(rdirichlet(1,rep(1,K)))
sigma0<-runif(K,exp(-2),exp(3))
beta1<-betagen(M,sigma0,0,pi)
beta2<-betagen(N,sigma0,0,pi)

betahat1<-rnorm(M,beta1,1)
betahat2<-rnorm(N,beta2,1)

r12<-ash(c(betahat1,betahat2),1,mixcompdist="normal",method="shrink")
r1<-ash(betahat1,1,mixcompdist="normal",method="shrink",g=r12$fitted.g)
r2<-ash(betahat2,1,mixcompdist="normal",method="shrink",g=r12$fitted.g)
t<-c(t,2*(r1$loglik+r2$loglik-r12$loglik))
}

write.table(t,file="null_mixture_normal_stat.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)

tmp<-proc.time()-tmp
write(tmp,file="null_mixture_normal_time.txt")