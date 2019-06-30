library(Rcpp)
library(RcppArmadillo)
library(bayesm)
library(ggplot2)
library(reshape2)
library(clues)

sourceCpp("dmrpfunctions.cpp")

# ------------------------------------------------------------------ #
# functions
# ------------------------------------------------------------------ #

# compute trace plot from set of draws
pitrace = function(draws){
  levels = unique(draws)
  nlevels = nrow(levels)
  counts = match(data.frame(t(draws)),data.frame(t(levels)))
  freq = table(counts)
  trace = match(data.frame(t(draws)),data.frame(t(levels)))
  return(trace)
}

# ------------------------------------------------------------------ #
# simulate data
# ------------------------------------------------------------------ #

# generate data
set.seed(1)
pi = rep(1:5,each=4)
n = length(pi)
K = max(pi)
Pi = matrix(0,n,K)
Pi[cbind(1:n,pi)] = 1
Bmult = Pi %*% diag(K) %*% t(Pi)
nbeta = sum(Bmult)
B = matrix(runif(n*n,0,2),n,n)*Bmult
diag(B) = -abs(diag(B))

# values of exogenous variables
T = 1000
P = matrix(runif(n*T,0,1),nrow=T)
Sigma = matrix(0.1,n,n) + 2*diag(n)
error = matrix(rnorm(T*n),ncol=n)%*%chol(Sigma)
Zlist = NULL
psimat = NULL
Zpsi = NULL
p = 2
nvar = 0
for(i in 1:n){
  Z = cbind(rep(1,T),matrix(runif(T*(p-1),-1,1),nrow=T))
  nvar = nvar + ncol(Z)
  psi = runif(p,-2,2)
  Zlist[[i]] = Z
  Zpsi = c(Zpsi,Z%*%psi)
  psimat = cbind(psimat,psi)
}
Zpsi = matrix(Zpsi,T,n)
Y = P%*%B + Zpsi + error

# true values
true.pi = pi
true.Pi = Pi
true.B = B
true.beta = B[Bmult==1]
true.psi = as.vector(psimat)
true.Sigma = Sigma
true.loglike = mvregll(Y,P%*%B+Zpsi,solve(Sigma))

# ------------------------------------------------------------------ #
# MCMC settings
# ------------------------------------------------------------------ #

R = 100000
keep = 10
nprint = 10000
end = R/keep
burn = end*.5

# ------------------------------------------------------------------ #
# unrestricted demand model
# ------------------------------------------------------------------ #

# prior
betab = 0
abeta = .1
psibar = double(nvar)
Apsi = .01*diag(nvar)
nu = 3+n
V = nu*diag(n)

# run sampler
out = mvregmcmc(Y, P, Zlist, 
                abeta, betab, Apsi, psibar, nu, V,
                R, keep, nprint)

# log likelihood trace plot
plot(out$loglikedraw,type="l")
abline(h=true.loglike,col=2)

# recovery of beta parameters 
plot(true.B,apply(out$betadraw[burn:end,],2,mean),xlab="true beta",ylab="estimated beta")
abline(0,1)
segments(
  true.B,  
  apply(out$betadraw[burn:end,],2,quantile,.025),
  true.B,  
  apply(out$betadraw[burn:end,],2,quantile,.975)
)

# recovery of psi parameters 
plot(true.psi,apply(out$psidraw[burn:end,],2,mean),xlab="true psi",ylab="estimated psi")
abline(0,1)
segments(
  true.psi,  
  apply(out$psidraw[burn:end,],2,quantile,.025),
  true.psi,  
  apply(out$psidraw[burn:end,],2,quantile,.975)
)

# ------------------------------------------------------------------ #
# isolated demand model with random partition and DP prior
# ------------------------------------------------------------------ #

# prior
betab = 0
abeta = .1
psibar = double(nvar)
Apsi = .01*diag(nvar)
nu = 3+n
V = nu*diag(n)
alpha = 1

# run sampler
outdp = isomvregdpmcmc(Y, P, Zlist, 
                       abeta, betab, Apsi, psibar, nu, V, alpha,
                       pi_step=1/(n*log(n)),
                       R, keep, nprint)

# log likelihood trace plot
plot(outdp$loglikedraw,type="l")
abline(h=true.loglike,col=2)

# partition trace plot
plot(pitrace(outdp$pidraw),ylab="partition",xlab="iteration",cex=.5)
abline(h=match(data.frame(true.pi),data.frame(t(unique(outdp$pidraw)))),col=2,lwd=1)

# recovery of beta parameters 
plot(true.B,apply(outdp$betadraw[burn:end,],2,mean),xlab="true beta",ylab="estimated beta")
abline(0,1)
segments(
  true.B,  
  apply(outdp$betadraw[burn:end,],2,quantile,.025),
  true.B,  
  apply(outdp$betadraw[burn:end,],2,quantile,.975)
)

# recovery of psi parameters 
plot(true.psi,apply(outdp$psidraw[burn:end,],2,mean),xlab="true psi",ylab="estimated psi")
abline(0,1)
segments(
  true.psi,  
  apply(outdp$psidraw[burn:end,],2,quantile,.025),
  true.psi,  
  apply(outdp$psidraw[burn:end,],2,quantile,.975)
)

# ------------------------------------------------------------------ #
# isolated demand model with random partition and LSP prior
# ------------------------------------------------------------------ #

# prior
betab = 0
abeta = .1
psibar = double(nvar)
Apsi = .01*diag(nvar)
nu = 3+n
V = nu*diag(n)
rho = rep(1,n)
tau = 0.1/(n*log(n))

# run sampler
outlsp = isomvreglspmcmc(Y, P, Zlist, 
                         abeta, betab, Apsi, psibar, nu, V, rho, tau,
                         pi_step=1/(n*log(n)),
                         R, keep, nprint)

# log likelihood trace plot
plot(outlsp$loglikedraw,type="l")
abline(h=true.loglike,col=2)

# partition trace plot
plot(pitrace(outlsp$pidraw),ylab="partition",xlab="iteration",cex=.5)
abline(h=match(data.frame(true.pi),data.frame(t(unique(outlsp$pidraw)))),col=2,lwd=1)

# recovery of beta parameters 
plot(true.B,apply(outlsp$betadraw[burn:end,],2,mean),xlab="true beta",ylab="estimated beta")
abline(0,1)
segments(
  true.B,  
  apply(outlsp$betadraw[burn:end,],2,quantile,.025),
  true.B,  
  apply(outlsp$betadraw[burn:end,],2,quantile,.975)
)

# recovery of psi parameters 
plot(true.psi,apply(outlsp$psidraw[burn:end,],2,mean),xlab="true psi",ylab="estimated psi")
abline(0,1)
segments(
  true.psi,  
  apply(outlsp$psidraw[burn:end,],2,quantile,.025),
  true.psi,  
  apply(outlsp$psidraw[burn:end,],2,quantile,.975)
)

# ------------------------------------------------------------------ #
# plot results
# ------------------------------------------------------------------ #

# posterior similarity matrix w/ dp prior
psm = compsm(outdp$pidraw[burn:end,])
psm[lower.tri(psm)] = .5*compsm(matrix(true.pi,1))[lower.tri(psm)]
psm = melt(psm)
ggplot(data=psm,aes(x=Var1,y=Var2,fill=value)) +
  geom_tile(color="white") +
  guides(fill=guide_legend(label=TRUE,keywidth=1,keyheight=0,reverse=TRUE,
                           title="true partition",
                           title.position="left",
                           title.theme = element_text(angle=270,size=12),
                           title.hjust=1)) +
  labs(x="",y="",title="posterior similarity",fill="") +
  scale_x_discrete(expand=c(0.001,0.001)) +
  scale_y_discrete(expand=c(0.001,0.001)) +
  scale_fill_gradient2(low="white",high="indianred2") +
  theme(axis.text.x=element_text(angle=90,hjust=1,size=5),
        axis.text.y=element_text(size=5),
        plot.title = element_text(hjust = 0.5)) +
  coord_fixed()

# posterior similarity matrix w/ lsp prior
psm = compsm(outlsp$pidraw[burn:end,])
psm[lower.tri(psm)] = .5*compsm(matrix(true.pi,1))[lower.tri(psm)]
psm = melt(psm)
ggplot(data=psm,aes(x=Var1,y=Var2,fill=value)) +
  geom_tile(color="white") +
  guides(fill=guide_legend(label=TRUE,keywidth=1,keyheight=0,reverse=TRUE,
                           title="true partition",
                           title.position="left",
                           title.theme = element_text(angle=270,size=12),
                           title.hjust=1)) +
  labs(x="",y="",title="posterior similarity",fill="") +
  scale_x_discrete(expand=c(0.001,0.001)) +
  scale_y_discrete(expand=c(0.001,0.001)) +
  scale_fill_gradient2(low="white",high="indianred2") +
  theme(axis.text.x=element_text(angle=90,hjust=1,size=5),
        axis.text.y=element_text(size=5),
        plot.title = element_text(hjust = 0.5)) +
  coord_fixed()

# compare mean and SD of own elasticities
wchdiag = diag(matrix(1:n^2,n,n))
plot(apply(out$betadraw[burn:end,wchdiag],2,mean),apply(outlsp$betadraw[burn:end,wchdiag],2,mean),
     xlab="unrestricted",ylab="restricted w/ LSP prior",main="posterior means",
     xlim=c(min(diag(B))-1,max(diag(B))+1),ylim=c(min(diag(B))-1,max(diag(B))+1))
abline(0,1)
plot(apply(out$betadraw[burn:end,wchdiag],2,sd),apply(outlsp$betadraw[burn:end,wchdiag],2,sd),
     xlab="unrestricted",ylab="restricted w/ LSP prior",main="posterior standard deviations",
     xlim=c(0,1),ylim=c(0,1))
abline(0,1)

# compare distribution over groups
hist(apply(outdp$pidraw[burn:end,],1,max),col=rgb(0,0,1,1/4),breaks=seq(-0.5,n+0.5, by=1),prob=TRUE,ylim=c(0,0.7),xlab="number of groups",ylab="density",main="")
hist(apply(outlsp$pidraw[burn:end,],1,max),breaks=seq(-0.5,n+0.5, by=1),col=rgb(1,0,0,1/4),border=1,prob=TRUE,add=TRUE)
legend("topright",c("DP","LSP"),y.intersp=1,bty="n",pch=16,col=c(rgb(0,0,1,1/2),rgb(1,0,0,1/2)))
abline(v=max(true.pi),lty=2)

# ------------------------------------------------------------------ #
# compute in-sample and predictive fit
# ------------------------------------------------------------------ #

# generate test data
Tout = 25
Pout = matrix(runif(n*Tout,0,1),nrow=Tout)
errorout = matrix(rnorm(Tout*n),ncol=n)%*%chol(Sigma)
Zlistout = NULL
Zpsiout = NULL
for(i in 1:n){
  Zout = cbind(rep(1,Tout),matrix(runif(Tout*(p-1),-1,1),nrow=Tout))
  Zlistout[[i]] = Zout
  Zpsiout = c(Zpsiout,Zout%*%psimat[,i])
}
Zpsiout = matrix(Zpsi,Tout,n)
Yout = Pout%*%B + Zpsiout + errorout

# RMSE and predicted RMSE
rmse = double(end-burn+1)
rmsedp = double(end-burn+1)
rmselsp = double(end-burn+1)
prmse = double(end-burn+1)
prmsedp = double(end-burn+1)
prmselsp = double(end-burn+1)
i = 0
for(r in burn:end){
  
  i = i + 1
  
  # in-sample RSME ----------------------------------------- #
  
  logYhat = matrix(bdiag(Zlist)%*%out$psidraw[r,],T,n) + P%*%matrix(out$betadraw[r,],n,n)
  rmse[i] = sqrt(mean((logYhat-Y)^2))
  
  logYhat = matrix(bdiag(Zlist)%*%outdp$psidraw[r,],T,n) + P%*%matrix(outdp$betadraw[r,],n,n)
  rmsedp[i] = sqrt(mean((logYhat-Y)^2))
  
  logYhat = matrix(bdiag(Zlist)%*%outlsp$psidraw[r,],T,n) + P%*%matrix(outlsp$betadraw[r,],n,n)
  rmselsp[i] = sqrt(mean((logYhat-Y)^2))
  
  # predicted RMSE ----------------------------------------- #
  
  logYhat = matrix(bdiag(Zlistout)%*%out$psidraw[r,],Tout,n) + Pout%*%matrix(out$betadraw[r,],n,n)
  prmse[i] = sqrt(mean((logYhat-Yout)^2))
  
  logYhat = matrix(bdiag(Zlistout)%*%outdp$psidraw[r,],Tout,n) + Pout%*%matrix(outdp$betadraw[r,],n,n)
  prmsedp[i] = sqrt(mean((logYhat-Yout)^2))
  
  logYhat = matrix(bdiag(Zlistout)%*%outlsp$psidraw[r,],Tout,n) + Pout%*%matrix(outlsp$betadraw[r,],n,n)
  prmselsp[i] = sqrt(mean((logYhat-Yout)^2))
  
}

# in-sample RMSE
round(apply(cbind(rmse,rmsedp,rmselsp),2,mean),3)
round(apply(cbind(rmse,rmsedp,rmselsp),2,sd),3)
boxplot(cbind(rmse,rmsedp,rmselsp),names=c("unrestricted","DP","LSP"),
        las=2,cex.axis=.75,main="RMSE")

# predictive RMSE
round(apply(cbind(prmse,prmsedp,prmselsp),2,mean),3)
round(apply(cbind(prmse,prmsedp,prmselsp),2,sd),3)
boxplot(cbind(prmse,prmsedp,prmselsp),names=c("unrestricted","DP","LSP"),
        las=2,cex.axis=.75,main="Predicted RMSE")

