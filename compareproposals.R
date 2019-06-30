library(Rcpp)
library(RcppArmadillo)
library(doParallel)
library(clues)
library(ggplot2)

sourceCpp("dmrpfunctions.cpp")

#----------------------------------------------------------------------------------#
# functions
#----------------------------------------------------------------------------------#

# simulate data set
simdatafn = function(n,T,sigma,Pi,minx,maxx){
  beta = rep(1,n)
  X = matrix(runif(n*T,minx,maxx),nrow=T)
  y = apply(((X*matrix(beta,T,n,byrow=TRUE))%*%Pi)^2,1,sum) + rnorm(T,0,sigma)
  return(list(y=y,X=X,beta=beta,sigma=sigma))
}

# split-merge proposal algorithm
splitmerge = function(pi){
  n = length(pi)
  index = 1:n
  gg = pi
  ij = sample(1:n,2)
  i = ij[1]
  j = ij[2]
  S = which(pi%in%pi[ij] & !(index%in%ij))
  shift = n + 1000
  
  # split
  if(pi[i]==pi[j]){
    gi.split = max(pi)+1
    gj.split = pi[j]
    gg[i] = gi.split
    gg[j] = gj.split
    gg[S] = sample(c(gi.split,gj.split),length(S),replace=TRUE)

    # counts
    ni.split = length(which(gg==gi.split))
    nj.split = length(which(gg==gj.split))
    
    # reorder
    gnew = gg + shift
    for(k in 1:max(gnew-shift)){
      gnew = replace(gnew,which(gnew==(unique(gnew-shift)[k]+shift)),k)
    }
    
    # transition probabilities
    logqnum = 0
    logqdenom = (ni.split+nj.split-2)*log(0.5)
  }
  
  # merge
  else{
    gi.merge = pi[j]
    gj.merge = pi[j]
    gg[i] = gi.merge
    gg[j] = gj.merge
    gg[S] = pi[j]
    
    # counts
    ni = length(which(gg==pi[i]))
    nj = length(which(gg==pi[j]))
    
    # reorder
    gnew = gg + shift
    for(k in 1:max(gnew-shift)){
      gnew = replace(gnew,which(gnew==(unique(gnew-shift)[k]+shift)),k)
    }
    
    # transition probabilities
    logqnum = (ni+nj-2)*log(0.5)
    logqdenom = 0
  }
  
  return(list(gnew=gnew,logqnum=logqnum,logqdenom=logqdenom))
}

# mcmc sampler for regression with split-merge proposals
regsm = function(Data, Mcmc) {
  
  # set data 
  y = Data$y
  X = Data$X
  n = Data$n
  sigma = Data$sigma
  beta = Data$beta
  true.pi = Data$true.pi
  
  # set mcmc
  R = Mcmc$R
  keep = Mcmc$keep
  step = Mcmc$step
  v = Mcmc$v
  pi0 = Mcmc$pi0
  
  # dimensions of storage matrices
  pidraw = matrix(nrow=floor(R/keep),ncol=n)
  
  # print time
  itime = proc.time()[3]
  cat("MCMC Iteration (estimated time to end in minutes)",fill=TRUE)
  
  # initial conditions
  if(is.null(pi0)){
    oldpi = as.vector(rlsp(1,1:n,n))
  }
  else{
    oldpi = pi0
  }
  naccept = 0
  
  # begin MH algorithm
  for (rep in 1:R){
    
    # set old values
    pi = oldpi
    K = max(pi)
    
    # generate new values
    sm = splitmerge(pi)
    pi.c = sm$gnew
    K.c = max(pi.c)
    Pi.c = matrix(0,nrow=n,ncol=K.c)
    Pi.c[cbind(1:n,pi.c)] = 1
    
    # proposal
    lognew.prop = sm$logqdenom
    logold.prop = sm$logqnum
    
    # prior
    logold.prior = 0
    lognew.prior = 0
    
    # compute the log likelihood of the old and new pi
    logold = regll(y,X,pi,sigma)
    lognew = regll(y,X,pi.c,sigma)
    
    # compute acceptance probability alpha
    ldiff = lognew + lognew.prior + logold.prop - logold - logold.prior - lognew.prop
    alpha = min(1,exp(ldiff))
    if(is.na(alpha)){alpha=-1}
    u = runif(1)
    if (u < alpha){
      oldpi = pi.c
      naccept = naccept+1
    }
    
    # save posterior draws
    mkeep = rep/keep
    if (mkeep*keep == (floor(mkeep)*keep)){      
      pidraw[mkeep,] = oldpi
    }
    
    # print run time and acceptance rate
    if (rep%%5000 == 0) {
      ctime = proc.time()[3]
      timetoend = ((ctime - itime)/rep)*(R - rep)
      cat(" ",rep," (",round(timetoend/60,1),")", " [",round(naccept/rep,4),"] ",fill = TRUE)
    }
  }
  
  # print total run time
  ctime = proc.time()[3]
  cat(" Total Time Elapsed: ",round((ctime - itime)/60,2),fill = TRUE)
  
  return(list(pidraw=pidraw))
}

# mcmc sampler for regression with Gibbs proposals
reggibbs = function(Data, Mcmc) {
  
  # set data 
  y = Data$y
  X = Data$X
  n = Data$n
  sigma = Data$sigma
  beta = Data$beta
  true.pi = Data$true.pi
  K = max(true.pi)
  
  # set mcmc
  R = Mcmc$R
  keep = Mcmc$keep
  step = Mcmc$step
  v = Mcmc$v
  pi0 = Mcmc$pi0
  
  # dimensions of storage matrices
  pidraw = matrix(nrow=floor(R/keep),ncol=n)
  
  # print time
  itime = proc.time()[3]
  cat("MCMC Iteration (estimated time to end in minutes)",fill=TRUE)
  
  # initial conditions
  if(is.null(pi0)){
    oldpi = as.vector(rlsp(1,1:n,n))
  }
  else{
    oldpi = pi0
  }
  
  # begin MH algorithm
  for (rep in 1:R){
    
    for(i in 1:n){
      lprobs = double(K)
      pi.c = oldpi
      for(k in 1:K){
        pi.c[i] = k
        lprobs[k] = regll(y,X,pi.c,sigma)
      }
      probs = exp(lprobs-max(lprobs))/sum(exp(lprobs-max(lprobs)))
      oldpi[i] = which.max(rmultinom(1,1,prob=probs))
    }
    
    # save posterior draws
    mkeep = rep/keep
    if (mkeep*keep == (floor(mkeep)*keep)){      
      pidraw[mkeep,] = oldpi
    }
    
    # print run time and acceptance rate
    if (rep%%5000 == 0) {
      ctime = proc.time()[3]
      timetoend = ((ctime - itime)/rep)*(R - rep)
      cat(" ",rep," (",round(timetoend/60,1),")",fill = TRUE)
    }
  }
  
  # print total run time
  ctime = proc.time()[3]
  cat(" Total Time Elapsed: ",round((ctime - itime)/60,2),fill = TRUE)
  
  return(list(pidraw=pidraw))
}

# compute the adjusted Rand similarity between draws and single partition pi
randfn = function(draws,pi){
  out = double(nrow(draws))
  for(i in 1:nrow(draws)){
    out[i] = adjustedRand(draws[i,],pi,"HA")
  }
  return(out)
}

#----------------------------------------------------------------------------------#
# simulation study in parallel (small n)
#----------------------------------------------------------------------------------#

n = 6
B = 25
R = 2000
keep = 5
end = R/keep
burn = .5*end
ncores = detectCores() - 1

# simulation w/ covariates in (-1,1)
registerDoParallel(cores=ncores)
cl = makeCluster(ncores)  
itime = proc.time()[3]
out1 = foreach(b=1:B,.combine=rbind) %dopar% {
  
  # simulate data
  pi = as.vector(rlsp(1,rep(1,n),1))
  Pi = matrix(0,nrow=n,ncol=max(pi))
  Pi[cbind(1:n,pi)] = 1
  simdata = simdatafn(n,T=20,sigma=1,Pi,minx=-1,maxx=1)
  Data = list(y=simdata$y,X=simdata$X,n=n,sigma=simdata$sigma,true.pi=pi,beta=simdata$beta)
  
  # initial partition
  pi0 = as.vector(rlsp(1,rep(1,n),n))
  
  # run samplers
  Mcmc = list(R=R,keep=keep,pi0=pi0)
  outgibbs = reggibbs(Data,Mcmc)
  outsm = regsm(Data,Mcmc)
  Mcmc = list(R=R,keep=keep,v=1/(n*log(n)),pi0=pi0,nprint=0)
  outlsp = reglsp(Data,Mcmc)
  Mcmc = list(R=R,keep=keep,v=1/(n*log(n)),blength=2,p=0,pi0=pi0,nprint=0)
  outblocklsp = regblocklsp(Data,Mcmc)
  
  cbind(randfn(outgibbs$pidraw[burn:end,],pi),
        randfn(outsm$pidraw[burn:end,],pi),
        randfn(outlsp$pidraw[burn:end,],pi),
        randfn(outblocklsp$pidraw[burn:end,],pi))
}
stopCluster(cl)  
ctime = proc.time()[3]
round((ctime - itime)/60,2)

# simulation w/ covariates in (0,2)
registerDoParallel(cores=ncores)
cl = makeCluster(ncores)  
itime = proc.time()[3]
out2 = foreach(b=1:B,.combine=rbind) %dopar% {
  
  # simulate data
  pi = as.vector(rlsp(1,rep(1,n),1))
  Pi = matrix(0,nrow=n,ncol=max(pi))
  Pi[cbind(1:n,pi)] = 1
  simdata = simdatafn(n,T=20,sigma=1,Pi,minx=0,maxx=2)
  Data = list(y=simdata$y,X=simdata$X,n=n,sigma=simdata$sigma,true.pi=pi,beta=simdata$beta)
  
  # initial partition
  pi0 = as.vector(rlsp(1,rep(1,n),n))
  
  # run samplers
  Mcmc = list(R=R,keep=keep,pi0=pi0)
  outgibbs = reggibbs(Data,Mcmc)
  outsm = regsm(Data,Mcmc)
  Mcmc = list(R=R,keep=keep,v=1/(n*log(n)),pi0=pi0,nprint=0)
  outlsp = reglsp(Data,Mcmc)
  Mcmc = list(R=R,keep=keep,v=1/(n*log(n)),blength=2,p=0,pi0=pi0,nprint=0)
  outblocklsp = regblocklsp(Data,Mcmc)
  
  cbind(randfn(outgibbs$pidraw[burn:end,],pi),
        randfn(outsm$pidraw[burn:end,],pi),
        randfn(outlsp$pidraw[burn:end,],pi),
        randfn(outblocklsp$pidraw[burn:end,],pi))
}
stopCluster(cl)  
ctime = proc.time()[3]
round((ctime - itime)/60,2)

# compare boxplots (FIGURE 6)
names = factor(c("GIBBS","SM","LSP","BLSP"),levels=c("GIBBS","SM","LSP","BLSP"))
boxdf = data.frame(y=as.vector(rbind(out1,out2)),
                   x=rep(c("Unif(-1,1)","Unif(0,2)"),each=B*(end-burn+1)),
                   method=rep(names,each=2*B*(end-burn+1)))
ggplot(boxdf,aes(method,y)) + 
  geom_boxplot(aes(method,y),outlier.size=.1) + 
  facet_grid(.~x) + 
  theme_classic() +
  labs(x = "",y="adjusted Rand index") + 
  theme(strip.background=element_rect(fill=alpha("grey",0),colour=alpha("grey",0)),
        strip.text.x=element_text(size=10),
        panel.border = element_rect(fill=NA,colour="black",size=.5))

#----------------------------------------------------------------------------------#
# simulation study in parallel (large n)
#----------------------------------------------------------------------------------#

# set dimensions
B = 5
nindex = c(5,10,15,20)
nindex = c(25,50,75,100)
Rindex = nindex*500

# simulation w/ many observations (T=10n)
ncores = detectCores() - 1
registerDoParallel(cores=ncores)
cl = makeCluster(ncores)
itime = proc.time()[3]
outhigh = foreach(n=nindex,i=icount()) %:%
  foreach(b=1:B,.combine=rbind) %dopar% {
    
    # mcmc settings
    R = Rindex[i]
    keep = 10
    end = R/keep
    burn = round(.5*end)
    
    # generate data
    pi = as.vector(rlsp(1,rep(1,n),1))
    Pi = matrix(0,nrow=n,ncol=max(pi))
    Pi[cbind(1:n,pi)] = 1
    simdata = simdatafn(n,T=10*n,sigma=1,Pi,minx=-1,maxx=1)
    Data = list(y=simdata$y,X=simdata$X,n=n,sigma=simdata$sigma,beta=simdata$beta)
    
    # initial partition
    pi0 = as.vector(rlsp(1,rep(1,n),n))
    
    # run
    Mcmc = list(R=R,keep=keep,v=1/(n*log(n)),pi0=pi0,nprint=0)
    outlsp = reglsp(Data,Mcmc)
    Mcmc = list(R=R,keep=keep,v=1/(n*log(n)),blength=2,p=0.5,pi0=pi0,nprint=0)
    outblocklsp = regblocklsp(Data,Mcmc)

    cbind(randfn(outlsp$pidraw[burn:end,],pi),
          randfn(outblocklsp$pidraw[burn:end,],pi),
          apply(outlsp$pidraw[burn:end,],1,max),
          apply(outblocklsp$pidraw[burn:end,],1,max),
          max(pi))
  }
stopCluster(cl)  
ctime = proc.time()[3]
round((ctime - itime)/60,2)

# simulation w/ fewer observations (T=5n)
ncores = detectCores() - 1
registerDoParallel(cores=ncores)
cl = makeCluster(ncores)
itime = proc.time()[3]
outlow = foreach(n=nindex,i=icount()) %:%
  foreach(b=1:B,.combine=rbind) %dopar% {
    
    # mcmc settings
    R = Rindex[i]
    keep = 10
    end = R/keep
    burn = round(.5*end)
    
    # generate data
    pi = as.vector(rlsp(1,rep(1,n),1))
    Pi = matrix(0,nrow=n,ncol=max(pi))
    Pi[cbind(1:n,pi)] = 1
    simdata = simdatafn(n,T=5*n,sigma=1,Pi,minx=-1,maxx=1)
    Data = list(y=simdata$y,X=simdata$X,n=n,sigma=simdata$sigma,beta=simdata$beta)
    
    # initial partition
    pi0 = as.vector(rlsp(1,rep(1,n),n))
    
    # run
    Mcmc = list(R=R,keep=keep,v=1/(n*log(n)),pi0=pi0,nprint=0)
    outlsp = reglsp(Data,Mcmc)
    Mcmc = list(R=R,keep=keep,v=1/(n*log(n)),blength=2,p=0.5,pi0=pi0,nprint=0)
    outblocklsp = regblocklsp(Data,Mcmc)
    
    cbind(randfn(outlsp$pidraw[burn:end,],pi),
          randfn(outblocklsp$pidraw[burn:end,],pi),
          apply(outlsp$pidraw[burn:end,],1,max),
          apply(outblocklsp$pidraw[burn:end,],1,max),
          max(pi))
  }
stopCluster(cl)  
ctime = proc.time()[3]
round((ctime - itime)/60,2)

# names of samplers and data dimensions
snames = factor(c("LSP","BLSP"),levels=c("LSP","BLSP"))
nnames = ordered(paste0("$n=",nindex,"$"),levels=paste0("$n=",nindex,"$"))
onames = factor(c("High","Low"),levels=c("High","Low"))

# Rand similarity
randdf = data.frame(y=c(unlist(lapply(outhigh,function(x)x[,1])),
                        unlist(lapply(outhigh,function(x)x[,2])),
                        unlist(lapply(outlow,function(x)x[,1])),
                        unlist(lapply(outlow,function(x)x[,2]))),
                    n=rep(rep(nnames,Rindex/10*B/2+B),4),
                    obs=rep(onames,each=2*sum(Rindex/10*B/2+B)),
                    sampler=rep(rep(snames,each=sum(Rindex/10*B/2+B)),2))
means = aggregate(randdf$y,list(n=randdf$n,obs=randdf$obs,sampler=randdf$sampler),FUN=mean)
sds = aggregate(randdf$y,list(n=randdf$n,obs=randdf$obs,sampler=randdf$sampler),FUN=sd)
means = melt(means,c("obs","n","sampler","x"))
sds = melt(sds,c("obs","n","sampler","x"))

# groups 
Kdf = data.frame(y=c(unlist(lapply(outhigh,function(x)x[,3]-x[,5])),
                     unlist(lapply(outhigh,function(x)x[,4]-x[,5])),
                     unlist(lapply(outlow,function(x)x[,3]-x[,5])),
                     unlist(lapply(outlow,function(x)x[,4]-x[,5]))),
                 n=rep(rep(nnames,Rindex/10*B/2+B),4),
                 obs=rep(onames,each=2*sum(Rindex/10*B/2+B)),
                 sampler=rep(rep(snames,each=sum(Rindex/10*B/2+B)),2))
Kmeans = aggregate(Kdf$y,list(n=Kdf$n,obs=Kdf$obs,sampler=Kdf$sampler),FUN=mean)
Ksds = aggregate(Kdf$y,list(n=Kdf$n,obs=Kdf$obs,sampler=Kdf$sampler),FUN=sd)
Kmeans = melt(Kmeans,c("obs","n","sampler","x"))
Ksds = melt(Ksds,c("obs","n","sampler","x"))

# create table summarizing results (TABLE 1)
tmp = means[means$sampler=="LSP",-3]
colnames(tmp)[3] = "LSPmean"
tmp$LSPsd = sds[sds$sampler=="LSP",]$x
tmp$BLSPmean = means[means$sampler=="BLSP",]$x
tmp$BLSPsd = sds[sds$sampler=="BLSP",]$x
tmp$blank = rep("",nrow(tmp))
tmp$LSPKmean = Kmeans[Kmeans$sampler=="LSP",]$x
tmp$LSPKsd = Ksds[Ksds$sampler=="LSP",]$x
tmp$BLSPKmean = Kmeans[Kmeans$sampler=="BLSP",]$x
tmp$BLSPKsd = Ksds[Ksds$sampler=="BLSP",]$x
