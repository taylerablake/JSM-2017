
cl <- makeCluster(detectCores()-1)
registerDoParallel(cl)
clusterCall(cl,function() {
      .libPaths("~/Rlibs/lib")
      library(doBy)
      library(doBy)
      library(lattice)
      library(MASS)
      library(magrittr)
      library(rlist)
      library(plyr)
      library(stringr)
      library(dplyr)
      library(doParallel)
})
setwd("/Users/taylerblake/GitRepos/JSM-2017")

source("/Users/taylerblake/GitRepos/pspline-mixed-models/lib/entropy-loss.R")
source("/Users/taylerblake/GitRepos/pspline-mixed-models/lib/quadratic-loss.R")
source("/Users/taylerblake/GitRepos/pspline-mixed-models/lib/help-functions.R")
source("/Users/taylerblake/GitRepos/pspline-mixed-models/lib/build-grid.r")



N <- 30
M <- m <- 20
set.seed(1985)

# Simulation data parameters

myGrid <- build_grid(20)
myGrid <- myGrid %>% transform(.,l_s=l/max(myGrid$l),
                               m_s=m/max(myGrid$m))

sig2 <- 0.3^2
Sigma <- diag(rep(sig2,M))
Omega <- solve(Sigma)


bPars <- rbind(c(0,1,length(unique(grid$l))-5,3,100,3),
               c(0,1,length(unique(grid$l))-5,3,100,3))

Bl <- bsplbase(as.vector(grid$l)/max(grid$l), bPars[1,])$base
Bm <- bsplbase(as.vector(grid$m)/max(grid$m), bPars[2,])$base

n1 <- ncol(Bl)
n2 <- ncol(Bm)

knots.l <- bsplbase(as.vector(grid$l)/max(grid$l), bPars[1,])$knots
interior.knots.l <- knots.l[1:(length(knots.l)-(bPars[1,4]+1))]

knots.m <- bsplbase(as.vector(grid$m)/max(grid$m), bPars[2,])$knots
interior.knots.m <- knots.m[1:(length(knots.m)-(bPars[1,4]+1))]

B. <- kronecker(Bm,
                t(as.vector(rep(1,ncol(Bl))))) * kronecker(t(as.vector(rep(1,ncol(Bm)))),
                                                           Bl)


knot_grid <- expand.grid(m=interior.knots.m,l=interior.knots.l)[,2:1]
keep <- (knot_grid$m >= 0.5) & (knot_grid$m <= max(knot_grid$l)-(0.5*knot_grid$l)) |
      (knot_grid$m <= 0.5) & (knot_grid$m >= min(knot_grid$l)+(0.5*knot_grid$l))
knot_grid <- transform(knot_grid,keep=factor(keep))
knot_grid <- data.frame(knot_grid,
                        expand.grid(m_index=1:length(interior.knots.m),
                                    l_index=(1:length(interior.knots.l)))[,2:1])

basisKeepIndex <-( knot_grid$keep[!duplicated(knot_grid[,c("l_index","m_index")])]=="TRUE") %>% which
B. <- B.[,basisKeepIndex]
knot_grid <- subset(knot_grid,keep==TRUE)


dl <- 0
if(dl>0){
      if(sum(knot_grid$m==unique(knot_grid$m)[1])>dl){
            Pl <- matrix(data=0,nrow=sum(knot_grid$m==unique(knot_grid$m)[1])-dl,
                         ncol=nrow(knot_grid))
            Pl[,which(knot_grid$m==unique(knot_grid$m)[1])] <- diff(diag(sum(knot_grid$m==unique(knot_grid$m)[1])),
                                                                    differences = dl)            
            for(i in 2:length(unique(knot_grid$m))){
                  if(sum(knot_grid$m==unique(knot_grid$m)[i])>dl){
                        pl <- matrix(data=0,nrow=sum(knot_grid$m==unique(knot_grid$m)[i])-dl,
                                     ncol=nrow(knot_grid))
                        pl[,which(knot_grid$m==unique(knot_grid$m)[i])] <- diff(diag(sum(knot_grid$m==unique(knot_grid$m)[i])),
                                                                                differences = dl)      
                        Pl <- rbind(Pl,pl)
                  }
            }
      }
      
      if((sum(knot_grid$m==unique(knot_grid$m)[1])<=dl) &(sum(knot_grid$m==unique(knot_grid$m)[2])>dl)){
            Pl <- matrix(data=0,nrow=sum(knot_grid$m==unique(knot_grid$m)[2])-dl,
                         ncol=nrow(knot_grid))
            Pl[,which(knot_grid$m==unique(knot_grid$m)[2])] <- diff(diag(sum(knot_grid$m==unique(knot_grid$m)[2])),
                                                                    differences = dl)            
            for(i in 3:length(unique(knot_grid$m))){
                  if(sum(knot_grid$m==unique(knot_grid$m)[i])>dl){
                        pl <- matrix(data=0,nrow=sum(knot_grid$m==unique(knot_grid$m)[i])-dl,
                                     ncol=nrow(knot_grid))
                        pl[,which(knot_grid$m==unique(knot_grid$m)[i])] <- diff(diag(sum(knot_grid$m==unique(knot_grid$m)[i])),
                                                                                differences = dl)      
                        Pl <- rbind(Pl,pl)      
                  }
            }
      }
      
}
if(dl==0){
      Pl <- diag(nrow(knot_grid))
}


dm <- 0
if(dm>0){
      if((sum(knot_grid$m==unique(knot_grid$m)[1])>dl)){
            Pm <- matrix(data=0,nrow=sum(knot_grid$l==unique(knot_grid$l)[1])-dm,
                         ncol=nrow(knot_grid))
            Pm[,which(knot_grid$l==unique(knot_grid$l)[1])] <- diff(diag(sum(knot_grid$l==unique(knot_grid$l)[1])),
                                                                    differences = dm)
            for(i in 2:length(unique(knot_grid$l))){
                  if(sum(knot_grid$m==unique(knot_grid$m)[i])>dl){
                        pm <- matrix(data=0,nrow=sum(knot_grid$l==unique(knot_grid$l)[i])-dm,
                                     ncol=nrow(knot_grid))
                        pm[,which(knot_grid$l==unique(knot_grid$l)[i])] <- diff(diag(sum(knot_grid$l==unique(knot_grid$l)[i])),
                                                                                differences = dm)      
                        Pm <- rbind(Pm,pm)
                  }
            }  
      }
      
      if((sum(knot_grid$m==unique(knot_grid$m)[1])<=dl) &(sum(knot_grid$m==unique(knot_grid$m)[2])>dl)){
            Pm <- matrix(data=0,nrow=sum(knot_grid$l==unique(knot_grid$l)[2])-dm,
                         ncol=nrow(knot_grid))
            Pm[,which(knot_grid$l==unique(knot_grid$l)[2])] <- diff(diag(sum(knot_grid$l==unique(knot_grid$l)[2])),
                                                                    differences = dm)
            for(i in 3:length(unique(knot_grid$l))){
                  if(sum(knot_grid$m==unique(knot_grid$m)[i])>dl){
                        pm <- matrix(data=0,nrow=sum(knot_grid$l==unique(knot_grid$l)[i])-dm,
                                     ncol=nrow(knot_grid))
                        pm[,which(knot_grid$l==unique(knot_grid$l)[i])] <- diff(diag(sum(knot_grid$l==unique(knot_grid$l)[i])),
                                                                                differences = dm)      
                        Pm <- rbind(Pm,pm)
                  }
            }  
      }
}
if(dm==0){
      Pm <- diag(nrow(knot_grid))
}


lambdas <- expand.grid(lam_l=exp(seq(-5,10,length.out=30)),
                       lam_m=exp(seq(-5,10,length.out=30)))

y <- mvrnorm(n=N,mu=rep(0,M),Sigma=Sigma)
y_vec <- as.vector(t(y[,-1]))

X <- matrix(data=0,nrow=N*(M-1),ncol=choose(M,2))
no.skip <- 0
for (t in 2:M){
      X[((0:(N-1))*(M-1)) + t-1,(no.skip+1):(no.skip+t-1)] <- y[,1:(t-1)]
      no.skip <- no.skip + t - 1
}

U. <- X%*%B.

clusterExport(cl,c("fit_cholesky_PS",
                   "Sigma",
                   "N",
                   "m",
                   "M",
                   "Bl",
                   "Bm",
                   "B.",
                   "Pl",
                   "Pm",
                   "lambdas",
                   "dl",
                   "dm",
                   "C",
                   "U."))

startTS <- Sys.time()
PS_fit_sim <- foreach(l=iter(lambdas,by="row")) %dopar% {
      
      Fit <- fit_cholesky_PS(y_vec,
                             N=N,
                             U.,
                             D=diag(diag(C)),
                             Pl,l$lam_l,
                             Pm,l$lam_m,
                             0.0001)    
      list(fit = Fit,
           lambda_l = l$lam_l,
           lambda_m = l$lam_m)
      
}
endTS <- Sys.time()
endTS-startTS


lam <- data.frame(lambda_l=unlist(lapply(PS_fit_sim,function(l){
      l$lambda_l
})),lambda_m=unlist(lapply(PS_fit_sim,function(l){
      l$lambda_m
})))
fit_list <- lapply(PS_fit_sim,function(l){
      Phi <- B. %*% l$fit$coef      
      T_mat <- diag(M)
      T_mat[lower.tri(T_mat)] <- -Phi
      list(fit=Phi,
           T_mod=T_mat)
})


omega_list <- lapply(fit_list,function(l){
      
      Omega_hat <- t(l$T_mod)%*%diag(1/diag(D)^2)%*%l$T_mod
      
})


quad_loss <- lapply(omega_list, function(omega_hat){
      quadratic_loss(omega_hat,Sigma)
}) %>%
      unlist

entrpy_loss <- lapply(omega_list, function(omega_hat){
      entropy_loss(omega_hat,Sigma)
}) %>%
      unlist


sse <- lapply(fit_list,function(l){
      t(y_vec - X %*% l$fit) %*% diag(1/rep(diag(C^2)[-1],N)) %*% (y_vec - X %*% l$fit)    
}) %>% unlist

ed <- lapply(PS_fit_sim,function(l){
      l$fit$eff.dim      
}) %>% unlist


loglik <- lapply(fit_list,function(l){
      sum(apply(y, MARGIN=1, FUN=function(x){
            log_lik( as.numeric(x), l$T_mod, diag(diag(C^2)))
      }))
}) %>% unlist
aic <- -2*loglik + 2*ed





## TODO : LEGEND NEEDS REFORMATTING, CHOOSE NON-DIVERGENT COLOR PALETTE
image.plot(y=log(sort(unique(lambdas$lam_l))),
           x=log(sort(unique(lambdas$lam_m))[-1]),
           z=matrix(aic[lambdas$lam_m>min(lambdas$lam_m)],
                    ncol=length(unique(lambdas$lam_l)),
                    nrow=length(unique(lambdas$lam_m))-1,
                    byrow=TRUE),
           ylab=expression("log "~lambda["l"]),
           xlab=expression("log "~lambda["m"]))



persp(y=log(sort(unique(lambdas$lam_l))),
      x=log(sort(unique(lambdas$lam_m))),
      z=matrix(aic[lambdas$lam_m>min(lambdas$lam_m)],
               ncol=length(unique(lambdas$lam_l)),
               nrow=length(unique(lambdas$lam_m)),
               byrow=TRUE),
      phi=15,
      theta=90,
      ylab="lambda l",
      xlab="lambda m",
      col = "light blue",
      main="",
      ticktype = "detailed",
      expand = 0.7,
      cex.axis=0.6,
      cex.lab=0.6,
      zlab="aic")


persp(y=log(sort(unique(lambdas$lam_l))),
      x=log(sort(unique(lambdas$lam_m))),
      z=matrix(ed[lambdas$lam_m>min(lambdas$lam_m)],
               ncol=length(unique(lambdas$lam_l)),
               nrow=length(unique(lambdas$lam_m)),
               byrow=TRUE),
      phi=15,
      theta=90,
      ylab="lambda l",
      xlab="lambda m",
      col = "light blue",
      main="effective model dimension",
      ticktype = "detailed",
      expand = 0.7,
      cex.axis=0.6,
      cex.lab=0.6,
      zlab="ed")

persp(z=diag(M)-T_mod,
      phi=15,
      theta=10,
      xlab="s",
      ylab="t",
      col = "light blue",
      main="true Cholesky surface",
      ticktype = "detailed",
      expand = 0.7,
      cex.axis=0.6,
      cex.lab=0.6,
      zlab="")
png(file.path(getwd(),"TeX","img","identity-cov-estimated-cholesky.png"))
persp(z=diag(M)-fit_list[[which.min(aic)]]$T_mod,
      phi=15,
      theta=10,
      xlab="s",
      ylab="t",
      col = "light blue",
      main="Estimated Cholesky surface",
      ticktype = "detailed",
      expand = 0.7,
      cex.axis=0.6,
      cex.lab=0.6,
      zlab="",
      zlim=c(-0.001,0.001))
dev.off()
Omega <- solve(Sigma)

persp(x=seq(0,1,length.out=M),
      y=seq(0,1,length.out=M),
      z=Omega,
      phi=15,
      theta=15,
      xlab="s",
      ylab="t",
      col = "light blue",
      zlab="",
      main=expression(paste(Sigma^{-1})),
      ticktype = "detailed",
      expand = 0.7,
      cex.axis=0.6,
      cex.lab=0.6)



persp(x=seq(0,1,length.out=M),
      y=seq(0,1,length.out=M),
      z=omega_list[[which.min(aic)]],
      phi=15,
      theta=15,
      xlab="s",
      ylab="t",
      col = "light blue",
      zlab="",
      main=expression(paste(hat(Sigma)^{-1})),
      ticktype = "detailed",
      expand = 0.7,
      cex.axis=0.6,
      cex.lab=0.6)


persp(x=seq(0,1,length.out=M),
      y=seq(0,1,length.out=M),
      z=omega_list[[which.min(entrpy_loss)]],
      phi=15,
      theta=15,
      xlab="s",
      ylab="t",
      col = "light blue",
      zlab="",
      main=expression(paste(hat(Sigma)^{-1})),
      ticktype = "detailed",
      expand = 0.7,
      cex.axis=0.6,
      cex.lab=0.6)



persp(x=seq(0,1,length.out=M),
      y=seq(0,1,length.out=M),
      z=omega_list[[which.min(aic)]],
      phi=15,
      theta=15,
      xlab="s",
      ylab="t",
      col = "light blue",
      zlab="",
      main=expression(paste(hat(Sigma)^{-1})),
      ticktype = "detailed",
      expand = 0.7,
      cex.axis=0.6,
      cex.lab=0.6)





png(file.path(getwd(),"TeX","img","identity-cov-true-covariance.png"))
persp(x=seq(0,1,length.out=M),
      y=seq(0,1,length.out=M),
      z=Sigma,
      phi=15,
      theta=15,
      xlab="s",
      ylab="t",
      col = "light blue",
      zlab="",
      main=expression(paste(Sigma)),
      ticktype = "detailed",
      expand = 0.7,
      cex.axis=0.6,
      cex.lab=0.6)
dev.off()


png(file.path(getwd(),"TeX","img","identity-cov-estimated-covariance.png"))
persp(x=seq(0,1,length.out=M),
      y=seq(0,1,length.out=M),
      z=solve(omega_list[[which.min(aic)]]),
      phi=15,
      theta=15,
      xlab="s",
      ylab="t",
      col = "light blue",
      zlab="",
      main=expression(paste(hat(Sigma))),
      ticktype = "detailed",
      expand = 0.7,
      cex.axis=0.6,
      cex.lab=0.6)
dev.off()
S <- outer(y[1,],y[1,])
for(subject in 2:nrow(y)){
      S <- S + outer(y[subject,],y[subject,])
}
S <- S/N

persp(x=seq(0,1,length.out=M),
      y=seq(0,1,length.out=M),
      z=S,
      phi=15,
      theta=15,
      xlab="s",
      ylab="t",
      col = "light blue",
      zlab="",
      main=expression(paste(tilde(Sigma))),
      ticktype = "detailed",
      expand = 0.7,
      cex.axis=0.6,
      cex.lab=0.6)

entropy_loss(solve(S),Sigma)
min(entrpy_loss)









fileName <- paste0("identity-covariance","-N-",N,"-M-",M, "-pordl-", dl, "-pordm-", dm)
save(list=ls(),file=file.path(getwd(),"data",paste0(fileName,".Rdata")))









Pen <- rbind(cbind(rep(0,nrow(Dl)),lambda_pair$laml*Dl,
                   matrix(data=0,nrow=nrow(Dl),ncol=ncol(Dm))),
             cbind(rep(0,nrow(Dl)),
                   matrix(data=0,nrow=nrow(Dm),ncol=ncol(Dl)),
                   lambda_pair$lamm*Dm))


Pen <- rbind(Pen,
             ridge.adj*cbind(rep(0,nl+nm),diag(nl+nm)))





