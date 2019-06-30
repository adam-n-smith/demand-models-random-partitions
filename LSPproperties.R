library(Rcpp)
library(RcppArmadillo)
library(ggplot2)
library(reshape2)
library(shallot)
library(clues)

sourceCpp("dmrpfunctions.cpp")

# ------------------------------------------------------------------ #
# functions
# ------------------------------------------------------------------ #

# sample partitions from distant dependent CRP
rddcrp = function(alpha,D,a){
  
  a = ifelse(a>0,a,1e-5)
  n = nrow(D)
  pi = double(n)
  for(i in 1:n){
    pvec = double(n)
    for(j in 1:n){
      if(i==j){
        pvec[j] = alpha
      }
      else{
        pvec[j] = exp(-D[i,j]/a)
      }
    }
    pi[i] = which.max(rmultinom(1,1,pvec/sum(pvec)))
  }
  return(pi)
}

# ------------------------------------------------------------- #
# LSP pairwise similarity (FIGURE 1)
# ------------------------------------------------------------- #

set.seed(1)
R = 10000
tau = .05
rho = rep(1:5,each=20)
n = length(rho)
draws = rlsp(R,rho,tau)
psm = compsm(draws)
upsm = melt(psm)
labs = pretty(1:n)
labs[1] = 1
ggplot(data=upsm,aes(x=Var1,y=Var2,fill=value)) + 
  geom_tile(color="white") + 
  guides(fill=FALSE) +
  labs(x="",y="",fill="") +
  scale_x_discrete(expand=c(0.001,0.001),limits=labs) +
  scale_y_discrete(expand=c(0.001,0.001),limits=labs) +
  scale_fill_gradient2(low="white",high="grey24") +
  coord_fixed()

# ------------------------------------------------------------- #
# LSP pairwise similarity w/ covariates (FIGURE 2)
# ------------------------------------------------------------- #

set.seed(1)
rho = rep(1:5,each=20)
n = length(rho)
tau = 0.05
X = as.matrix(c(rnorm(n/2,0,.1),rnorm(n/2,1,.1)))
X = apply(X,2,function(x)(x-mean(x))/sd(x))
p = ncol(X)
lambda = 0.05

# plot similarity matrix
R = 10000
draws = rlspx(R,rho,tau,X,lambda)
psm = compsm(draws)
upsm = melt(psm)
labs = pretty(1:n)
labs[1] = 1
ggplot(data=upsm,aes(x=Var1,y=Var2,fill=value)) + 
  geom_tile(color="white") + 
  guides(fill=FALSE) +
  labs(x="",y="",fill="") +
  scale_x_discrete(expand=c(0.001,0.001),limits=labs) +
  scale_y_discrete(expand=c(0.001,0.001),limits=labs) +
  scale_fill_gradient2(low="white",high="grey24") +
  coord_fixed()

# ------------------------------------------------------------- #
# compare LSP to other partition distributions (FIGURE 3)
# ------------------------------------------------------------- #

set.seed(1)
rho = rep(1:2,each=5)
n = length(rho)
D = 1-compsm(matrix(rho,1))
R = 10000
tau = seq(0,5,by=.1)
out = matrix(0,length(tau),4)
for(i in 1:length(tau)){
  
  # LSP
  lspdraws = rlsp(R,rho,tau[i])
  
  # EPA
  permutation = permutation(n.items=n,fixed=FALSE)
  decay = decay.exponential(temperature(1/tau[i]),as.dist(D+1e-5))
  x = sample.partitions(ewens.pitman.attraction(mass(max(rho)),discount(0),attraction(permutation,decay)),R, parallel = TRUE)
  x = process.samples(x)
  epadraws = x$labels
  
  # ddCRP
  ddcrpdraws = matrix(0,nrow=R,ncol=n)
  for(r in 1:R){
    ddcrpdraws[r,] = rddcrp(max(rho),D,tau[i])
  }
  
  # DP
  x = sample.partitions(ewens(mass(max(rho)),n),R, parallel = TRUE)
  x = process.samples(x)
  dpdraws = x$labels
  
  # similarity
  sim = matrix(0,R,4)
  for(r in 1:R){
    sim[r,1] = adjustedRand(lspdraws[r,],rho,randMethod = "HA")
    sim[r,2] = adjustedRand(epadraws[r,],rho,randMethod = "HA")
    sim[r,3] = adjustedRand(ddcrpdraws[r,],rho,randMethod = "HA")
    sim[r,4] = adjustedRand(dpdraws[r,],rho,randMethod = "HA")
  }
  out[i,] = apply(sim,2,mean)
}
matplot(tau,out,type="l",xlab="scaling parameter",ylab="similarity to location partition",lwd=2)
legend("topright",c("LSP","EPA","ddCRP","DP"),lt=1:4,col=1:4,lwd=2,bt="n",y.intersp=2)

# ------------------------------------------------------------- #
# 3D plots of ddCRP and EPA (FIGURE B.1)
# ------------------------------------------------------------- #

rho = rep(1:2,each=5)
n = length(rho)
D = 1-compsm(matrix(rho,1))
R = 10000

# ddCRP: scaling vs. mass
set.seed(1)
tau = seq(0,5,by=.1)
alpha = seq(0.1,5,by=.1)
out = matrix(0,length(tau),length(alpha))
for(i in 1:length(tau)){
  
  for(j in 1:length(alpha)){
    
    # ddCRP
    ddcrpdraws = matrix(0,nrow=R,ncol=n)
    for(r in 1:R){
      ddcrpdraws[r,] = rddcrp(alpha[j],D,tau[i])
    }
    
    # similarity
    sim = double(R)
    for(r in 1:R){
      sim[r] = adjustedRand(ddcrpdraws[r,],rho,randMethod = "HA")
    }
    
    out[i,j] = mean(sim)
    
  }
  
  print(tau[i])
  
}
persp(tau,alpha,out,zlim=c(0,1),ticktype="detailed",theta=40,phi=10,
      xlab="scaling parameter",ylab="mass parameter",
      zlab="similarity to location partition")

# EPA: scaling vs. mass
set.seed(1)
tau = seq(0,5,by=.1)
alpha = seq(0.1,5,by=.1)
out = matrix(0,length(tau),length(alpha))
for(i in 1:length(tau)){
  
  for(j in 1:length(alpha)){
    
    # EPA
    permutation = permutation(n.items=n,fixed=FALSE)
    decay = decay.exponential(temperature(1/tau[i]),as.dist(D+1e-5))
    x = sample.partitions(ewens.pitman.attraction(mass(alpha[j]),discount(0),attraction(permutation,decay)),R, parallel = TRUE)
    x = process.samples(x)
    epadraws = x$labels
    
    # similarity
    sim = double(R)
    for(r in 1:R){
      sim[r] = adjustedRand(epadraws[r,],rho,randMethod = "HA")
    }
    
    out[i,j] = mean(sim)
    
  }
  
  print(tau[i])
  
}
persp(tau,alpha,out,zlim=c(0,1),ticktype="detailed",theta=40,phi=10,
      xlab="scaling parameter",ylab="mass parameter",
      zlab="similarity to location partition")

# EPA: scaling vs. discount
set.seed(1)
tau = seq(0,5,by=.1)
delta = seq(0,0.9,by=.1)
out = matrix(0,length(tau),length(delta))
for(i in 1:length(tau)){
  
  for(j in 1:length(delta)){
    
    # EPA
    permutation = permutation(n.items=n,fixed=FALSE)
    decay = decay.exponential(temperature(1/tau[i]),as.dist(D+1e-5))
    x = sample.partitions(ewens.pitman.attraction(mass(max(rho)),discount(delta[j]),attraction(permutation,decay)),R, parallel = TRUE)
    x = process.samples(x)
    epadraws = x$labels
    
    # similarity
    sim = double(R)
    for(r in 1:R){
      sim[r] = adjustedRand(epadraws[r,],rho,randMethod = "HA")
    }
    
    out[i,j] = mean(sim)
    
  }
  
  cat(tau[i])
  
}
persp(tau,delta,out,zlim=c(0,1),ticktype="detailed",theta=40,phi=10,
      xlab="scaling parameter",ylab="discount parameter",
      zlab="similarity to location partition")

# EPA: mass vs. discount
set.seed(1)
delta = seq(0,0.9,by=.1)
alpha = seq(0.1,5,by=.1)
out = matrix(0,length(alpha),length(delta))
for(i in 1:length(alpha)){
  
  for(j in 1:length(delta)){
    
    # EPA
    permutation = permutation(n.items=n,fixed=FALSE)
    decay = decay.exponential(temperature(1/0),as.dist(D+1e-5))
    x = sample.partitions(ewens.pitman.attraction(mass(alpha[i]),discount(delta[j]),attraction(permutation,decay)),R, parallel = TRUE)
    x = process.samples(x)
    epadraws = x$labels
    
    # similarity
    sim = double(R)
    for(r in 1:R){
      sim[r] = adjustedRand(epadraws[r,],rho,randMethod = "HA")
    }
    
    out[i,j] = mean(sim)
    
  }
  
  cat(alpha[i],sep=" ")
  
}
persp(alpha,delta,out,zlim=c(0,1),ticktype="detailed",theta=40,phi=10,
      xlab="mass parameter",ylab="discount parameter",
      zlab="similarity to location partition")
