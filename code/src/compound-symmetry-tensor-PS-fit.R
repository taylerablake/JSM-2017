
# Load functions

source("/Users/taylerblake/GitRepos/pspline-mixed-models/lib/entropy-loss.R")
source("/Users/taylerblake/GitRepos/pspline-mixed-models/lib/quadratic-loss.R")
source("/Users/taylerblake/GitRepos/pspline-mixed-models/lib/help-functions.R")
source("/Users/taylerblake/GitRepos/pspline-mixed-models/lib/build-grid.r")
source(file.path(getwd(),"code/lib/helper-functions-2.R"))
source(file.path(getwd(),"code/lib/fit_cholesky_PS.R"))

cl <- makeCluster(detectCores()-1)
registerDoParallel(cl)
clusterCall(cl,function() {
      .libPaths("~/Rlibs/lib")
      library(doBy)
      library(splines)
      library(lattice)
      library(MASS)
      library(grDevices)
      library(magrittr)
      library(rlist)
      library(plyr)
      library(stringr)
      library(dplyr)
      library(doParallel)
})

set.seed(1985)

# Simulation data parameters
M = 20  # within-subject sample size

myGrid <- build_grid(20)
myGrid <- myGrid %>% transform(.,l_s=l/(max(myGrid$l)+1),
                               m_s=m/(max(myGrid$m)+1))

sigma. = 0.05  # error noise

flag. = TRUE  # to monitor Schall algorithm iterations


# B-spline parameters
bdeg = 3  # B-spline degree

nsegl = 14  # number of inner knots to build B-spline bases for x1 and x2
nsegm = 14

div = 1  # divisor to build nested bases for the interaction terms
# = 1 (no nesting) div>1 (nested bases for interaction with
# nseg/div inner knots) To ensure full nesting, nseg/div
# has to be an integer

N <- 30

Sigma <- matrix(0.7,nrow=M,ncol=M) + diag(rep(0.3),M)
Omega <- solve(Sigma)


C <- t(chol(Sigma))
D <- diag(diag(C))
L <- C%*%solve(D)
T_mod <- solve(L)
phi <- as.vector(T_mod[lower.tri(T_mod)])



Bl <- bspline(myGrid$l_s, 0, 1, nsegl, bdeg)$B
Bm <- bspline(myGrid$m_s, 0, 1, nsegl, bdeg)$B

n1 <- ncol(Bl)
n2 <- ncol(Bm)

knots.l <- bspline(as.vector(grid$l)/max(grid$l), 0, 1, nsegl, bdeg)$knots
interior.knots.l <- knots.l[1:(length(knots.l)-(bdeg+1))]

knots.m <- bspline(as.vector(grid$m)/max(grid$m), 0, 1, nsegl, bdeg)$knots
interior.knots.m <- knots.m[1:(length(knots.m)-(bdeg+1))]

B. <- kronecker(Bm,
                t(as.vector(rep(1,ncol(Bl))))) * kronecker(t(as.vector(rep(1,ncol(Bm)))),
                                                           Bl)


ggplot(myGrid,aes(l_s,m_s)) +
      geom_point(size=0.8) +
      theme_minimal() +
      xlab("l") +
      ylab("m") +
      geom_point(data=knot_grid[which(colSums(B. != 0)==0),],
                 aes(x=l,y=m),size=0.5,colour="red")



# Create a function interpolating colors in the range of specified colors
jet.colors <- colorRampPalette( c("blue", "green") )
# Generate the desired number of colors from this palette
nbcol <- 100
color <- jet.colors(nbcol)
# Compute the z-value at the facet centres
z <- matrix(Phi,
            nrow=100,
            ncol=100,
            byrow = FALSE)
if(Type.!="f05"){
      zfacet <- z[-1, -1] + z[-1, -100] + z[-100, -1] + z[-100, -100]
      # Recode facet z-values into color indices
      facetcol <- cut(zfacet, nbcol)
      plotColor <- color[facetcol]
} else {plotColor <- color}


persp(x=seq(0,1,length.out=100),
      y=seq(0,1,length.out=100),
      z=z,
      phi=15,
      theta=25,
      xlab="l",
      ylab="m",
      col = plotColor,
      #      zlab=expression(phi^"\u2217"),
      #zlab=expression(paste(phi,"^",*)),
      ticktype = "detailed",
      main="True Cholesky surface",
      expand = 0.5,
      zlab="")


N <- 50
y <- matrix(data=NA,nrow=N,M)
y[,1] <- rnorm(N,sd=sigma.)
for(this_t in 2:M){
      y[,this_t] <- phi[myGrid$t==this_t & myGrid$s == 1]*y[,1] + rnorm(N,sd=sigma.)
      for(this_s in 1:(this_t-1)){
            y[,this_t] <- y[,this_t] + phi[myGrid$t==this_t & myGrid$s == this_s]*y[,this_s]
      }
}
matplot(t(y),type="l")


T_mat <- diag(M)
T_mat[lower.tri(T_mat)] <- -phi
if(Type. != "f05"){
      Tfacet <- (T_mat-diag(M))[-1, -1] + (T_mat-diag(M))[-1, -M] + (T_mat-diag(M))[-M, -1] + (T_mat-diag(M))[-M, -M]
      # Recode facet z-values into color indices
      facetcol <- cut(Tfacet, nbcol)
      persp(z=T_mat-diag(M),
            phi=15,
            theta=20,
            xlab="l",
            ylab="m",
            col = color[facetcol],
            zlab="phi",
            ticktype = "detailed",
            expand = 0.5)
}

W <- matrix(data=0,nrow=N*(M-1),ncol=choose(M,2))
no.skip <- 0
for (t in 2:M){
      W[((0:(N-1))*(M-1)) + t-1,(no.skip+1):(no.skip+t-1)] <- y[,1:(t-1)]
      no.skip <- no.skip + t - 1
}
yVec <- as.vector(y[,-1])





start1 <- proc.time()[3]

# Build Mixed Model Bases For nested bases for the interaction terms, create
# new mixed model basis with number of segments as a integer divisor of nseg1,
# and nseg2 for div=1, G1=G1n and G2=G2n

pordl <- 3
pordm <- 2


no_support <- which(colSums(B. != 0)==0)

dm <- 2
if ( dm > 0 ){
      Dm <- diff(diag(n2),differences = dm)    
} else {
      Dm <- diag(n2)
}
Pm <- kronecker(t(Dm)%*%Dm, diag(n1))


dl <- 3
if ( dl > 0 ){
      Dl <- diff(diag(n1),differences = dl)     
} else {
      Dl <- diag(n1)
}
Pl <- kronecker(diag(n2),t(Dl)%*%Dl)
difference_rows_remove <- which(apply((abs(Pl) + abs(Pm))[,no_support],
                                MARGIN=1,
                                FUN=max)>0)
Pl <- Pl[-difference_rows_remove,-no_support]
Pm <- Pm[-difference_rows_remove,-no_support]



y <- mvrnorm(n=N,mu=rep(0,M),Sigma=Sigma)
y_vec <- as.vector(t(y[,-1]))

X <- matrix(data=0,nrow=N*(M-1),ncol=choose(M,2))
no.skip <- 0
for (t in 2:M){
      X[((0:(N-1))*(M-1)) + t-1,(no.skip+1):(no.skip+t-1)] <- y[,1:(t-1)]
      no.skip <- no.skip + t - 1
}

U. <- X%*%B.[,-no_support]

lambdas <- expand.grid(lam_l=exp(seq(-5,11,length.out=30)),
                       lam_m=exp(seq(-5,11,length.out=30)))

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


PS_fit_sim <- foreach(l=iter(lambdas,by="row")) %dopar% {
      
      Fit <- fit_cholesky_PS(y_vec,
                             N=N,
                             U.,
                             D=diag(diag(C)),
                             Pl,l$lam_l,
                             Pm,l$lam_m,
                             0.00000000001)    
      list(fit = Fit,
           lambda_l = l$lam_l,
           lambda_m = l$lam_m)
      
}

Fit <- MM$M %*% b
ll <- ldply(ll,data.frame)
sigmasq <- sigmasq[-1]
AIC. <- AIC.[-1]
SSR. <- SSR.[-1]
MSE. <- MSE.[-1]
tt <- ldply(tt,data.frame)
effdim <- ldply(effdim,data.frame)

rm(V)
rm(Vu)
rm(WXtWX.)
rm(WXtWZ.)
rm(WXty.)
rm(WZty.)
rm(WZtWZ.)



ggplot(data.frame(it=(1:length(SSR.))*2,ssr=SSR.),
       aes(x=it,y=ssr)) +
      geom_line() +
      theme_minimal() +
      #      ylim(0,0.02) +
      xlab("iteration") +
      ylab("ssr")

ggplot(data.frame(it=(1:length(MSE.))*2,mse=MSE.),
       aes(x=it,y=mse)) +
      geom_line() +
      theme_minimal() +
      #      ylim(0,0.02) +
      xlab("iteration") +
      ylab("mse")


ggplot(data.frame(it=(1:length(AIC.))*2,aic=AIC.),
       aes(x=it,y=aic)) +
      geom_line() +
      theme_minimal() +
      #      ylim(0,0.02) +
      xlab("iteration") +
      ylab("aic")


ggplot(data.frame(it=(1:length(sigmasq))*2,sigma2=sigmasq),
       aes(x=it,y=sigma2)) +
      geom_line() +
      theme_minimal() +
      #      ylim(0,0.02) +
      xlab("iteration") +
      ylab(expression(hat(sigma)^2)) +
      annotate(geom="text",
               x=125,
               y=(sigma.^2)+0.00001,
               label= paste("sigma^2 == ", sigma.^2),parse=TRUE)
ggsave(file=file.path(getwd(),"reports-1-simulations",paste0(Type.,"-sigma2-convergence.png")))
ggplot(ll,aes(x=it,y=log(lambda))) +
      geom_line(aes(colour=component)) +
      theme_minimal() +
      xlab("iteration") +
      ylab(expression("log("~lambda~")"))
ggsave(file=file.path(getwd(),"reports-1-simulations","f6-lambdas-convergence.png"))
# ggplot(tt,aes(x=it,y=tau)) +
#       geom_line(aes(colour=component)) +
#       theme_minimal() +
#       xlab("iteration") +
#       ylab(expression(tau)) +
#       ylim(0,.2)
ggplot(effdim,aes(x=it,y=ed)) +
      geom_line(aes(colour=component)) +
      theme_minimal() +
      xlab("iteration") +
      ylab("ED")
ggsave(file=file.path(getwd(),"reports-1-simulations","f6-ED-convergence.png"))

#####################################################################################
T_mat <- diag(M)
T_mat[lower.tri(T_mat)] <- -phi

if (Type. != "f05"){
      Tfacet <- (-T_mat+diag(M))[-1, -1] - (T_mat-diag(M))[-1, -M] - (T_mat-diag(M))[-M, -1] - (T_mat-diag(M))[-M, -M]
      # Recode facet z-values into color indices
      facetcol <- cut(Tfacet, nbcol)
      # png(file=file.path(getwd(),"reports-1-simulations",
      #                    "img",
      #                    paste0(Type.,"-persp-true-vs-est-cholesky-surface.png")))
      persp(z=diag(M)-T_mat,
            phi=15,
            theta=10,
            xlab="s",
            ylab="t",
            col = color[facetcol],
            #main=expression(paste(phi," = ", e^{-6}*plain(cos),"(6",pi,"l)")),
            ticktype = "detailed",
            expand = 0.7,
            cex.axis=0.6,
            cex.lab=0.6,
            zlab="")
      #main="True Cholesky surface")
      # dev.off()
}

T_mat_hat <- diag(M)
T_mat_hat[lower.tri(T_mat_hat)] <- -Fit
# Recode facet z-values into color indices

persp(z=T_mat_hat-diag(M),
      phi=15,
      theta=10,
      xlab="s",
      ylab="t",
      col = plotColor,
      zlab="",
      main=expression(paste(hat(phi))),
      ticktype = "detailed",
      expand = 0.7,
      cex.axis=0.6,
      cex.lab=0.6)



Omega <- t(T_mat)%*%diag(rep(1/sigma.^2, M))%*%T_mat
Omegafacet <- Omega[-1, -1] + Omega[-1, -M] + Omega[-M, -1] + Omega[-M, -M]
# Recode facet z-values into color indices
facetcol <- cut(Omegafacet, nbcol)
#par(mfrow=c(1,2))
if (Type. !="f05") {
      plotColor <- color[facetcol]
} else{
      plotColor <- "light blue"
}
persp(x=seq(0,1,length.out=M),
      y=seq(0,1,length.out=M),
      z=Omega,
      phi=15,
      theta=15,
      xlab="s",
      ylab="t",
      col = color[facetcol],
      zlab="",
      main=expression(paste(Sigma^{-1})),
      ticktype = "detailed",
      expand = 0.7,
      cex.axis=0.6,
      cex.lab=0.6)

#Omega_hat <- t(T_mat_hat)%*%diag(rep(1/sig2, M))%*%T_mat_hat  
Omega_hat <- t(T_mat_hat)%*%diag(rep(1/sigma.^2, M))%*%T_mat_hat  
facetcol <- cut(Omegafacet, nbcol)


# png(file=file.path(getwd(),"reports-1-simulations","img",
#                    paste0(Type.,"-persp-true-vs-estimated-precision-matrix.png")))
persp(x=seq(0,1,length.out=M),
      y=seq(0,1,length.out=M),
      z=Omega_hat,
      phi=15,
      theta=15,
      xlab="s",
      ylab="t",
      col = plotColor,
      zlab="",
      main=expression(paste(hat(Sigma)^{-1})),
      ticktype = "detailed",
      expand = 0.7,
      cex.axis=0.6,
      cex.lab=0.6)
#dev.off()


Sigma <- solve(Omega)
Sigmafacet <- Sigma[-1, -1] + Sigma[-1, -M] + Sigma[-M, -1] + Sigma[-M, -M]
Sigma_hat <- solve(Omega_hat)
#par(mfrow=c(3,1))
# Recode facet z-values into color indices
if (Type. !="f05") {
      facetcol <- cut(Sigmafacet, nbcol)
      plotColor <- color[facetcol]
} else{
      plotColor <- "light blue"
}


# png(file=file.path(getwd(),"reports-1-simulations","img",
#                    paste0(Type.,"-persp-true-vs-estimated-covariance-matrix.png")))
persp(x=seq(0,1,length.out=M),
      y=seq(0,1,length.out=M),
      z=Sigma,
      phi=10,
      theta=10,
      xlab="s",
      ylab="t",
      col = plotColor,
      zlab="",
      main=expression(paste(Sigma)),
      ticktype = "detailed",
      expand = 0.7,
      cex.axis=0.6,
      cex.lab=0.6)#,
#      zlim=c(min(min(Sigma),min(Sigma_hat)),
#             max(max(Sigma),max(Sigma_hat))))

S <- outer(y[1,],y[1,])
for(subject in 2:nrow(y)){
      S <- S + outer(y[subject,],y[subject,])
}
S <- S/N

persp(x=seq(0,1,length.out=M),
      y=seq(0,1,length.out=M),
      z=S,
      xlab="s",
      ylab="t",
      col = plotColor,
      #      zlab=expression(phi^"\u2217"),
      #zlab=expression(paste(phi,"^",*)),
      ticktype = "detailed",
      main="S",
      phi=10,
      theta=10,
      zlab="",
      expand = 0.7,
      cex.axis=0.6,
      cex.lab=0.6)


persp(x=seq(0,1,length.out=M),
      y=seq(0,1,length.out=M),
      z=Sigma_hat,
      phi=10,
      theta=10,
      xlab="s",
      ylab="t",
      col = plotColor,
      zlab="",
      main=expression(paste(hat(Sigma))),
      ticktype = "detailed",
      expand = 0.7,
      cex.axis=0.6,
      cex.lab=0.6)
#dev.off()
#####################################################################################


omega_list <- list.zip(b=bb,sig2=as.list(sigmasq)) %>%
      lapply(.,function(l){
            phiHat <- MM$M%*%l$b
            T_hat <- diag(M)
            T_hat[lower.tri(T_hat)] <- Fit
            Omega_hat <- t(T_hat)%*%diag(rep(1/l$sig2, M))%*%T_hat
      })
quad_loss <- lapply(omega_list,quadratic_loss,Sig=Sigma)
entrpy_loss <- lapply(omega_list,entropy_loss,Sig=Sigma)

data.frame(it=2*(1:length(quad_loss)),
           quadLoss=unlist(quad_loss),
           entropyLoss=unlist(entrpy_loss)) %>%
      ggplot(.,aes(x=it)) +
      geom_line(aes(y=entropyLoss)) +
      theme_minimal() +
      xlab("iteration") +
      ylab("") +
      ggtitle(expression(paste(L[1],"(",Sigma,",",hat(Omega),")")))
ggsave(file=file.path(getwd(),"reports-1-simulations",
                      "img",
                      paste0(Type.,"-entropy-loss-vs-iteration.png")))

data.frame(it=2*(1:length(quad_loss)),
           quadLoss=unlist(quad_loss),
           entropyLoss=unlist(entrpy_loss)) %>%
      ggplot(.,aes(x=it)) +
      geom_line(aes(y=quadLoss)) +
      theme_minimal() +
      xlab("iteration") +
      ylab("") +
      ggtitle(expression(paste(L[2],"(",Sigma,",",hat(Omega),")")))

which.min(unlist(quad_loss))

b <- bb[[which.min(unlist(quad_loss))]]
sig2 <- sigmasq[which.min(unlist(quad_loss))]
bestFit <- MM$M %*% b

T_mat_hat <- diag(M)
T_mat_hat[lower.tri(T_mat_hat)] <- -bestFit
# Recode facet z-values into color indices

if(Type.!="f05"){
      facetcol <- cut(Tfacet, nbcol)     
      plotColor <- color[facetcol]
} else if (Type.=="f05") {
      plotColor <- "light blue"
}
persp(z=T_mat_hat-diag(M),
      phi=15,
      theta=10,
      xlab="s",
      ylab="t",
      col = plotColor,
      zlab="",
      main=expression(paste(hat(phi))),
      ticktype = "detailed",
      expand = 0.7,
      cex.axis=0.6,
      cex.lab=0.6)
dev.off()


Omega <- t(T_mat)%*%diag(rep(1/sigma.^2, M))%*%T_mat
Omegafacet <- Omega[-1, -1] + Omega[-1, -M] + Omega[-M, -1] + Omega[-M, -M]
# Recode facet z-values into color indices
facetcol <- cut(Omegafacet, nbcol)
if (Type. !="f05") {
      plotColor <- color[facetcol]
} else{
      plotColor <- "light blue"
}

Omega_hat <- t(T_mat_hat)%*%diag(rep(1/sig2, M))%*%T_mat_hat  
facetcol <- cut(Omegafacet, nbcol)
persp(z=Omega_hat,
      phi=15,
      theta=10,
      xlab="s",
      ylab="t",
      col = plotColor,
      zlab="",
      main=expression(paste(hat(Omega))),
      ticktype = "detailed",
      expand = 0.7,
      cex.axis=0.6,
      cex.lab=0.6)

entropy_loss(Omega_hat,Sigma)      
Sigma <- solve(Omega)
Sigmafacet <- Sigma[-1, -1] + Sigma[-1, -M] + Sigma[-M, -1] + Sigma[-M, -M]
Sigma_hat <- solve(Omega_hat)
persp(z=Sigma_hat,
      phi=15,
      theta=10,
      xlab="s",
      ylab="t",
      col = plotColor,
      zlab="",
      main=expression(paste(hat(Sigma))),
      ticktype = "detailed",
      expand = 0.7,
      cex.axis=0.6,
      cex.lab=0.6)

persp(z=Sigma,
      phi=15,
      theta=10,
      xlab="s",
      ylab="t",
      col = plotColor,
      zlab="",
      main=expression(paste(Sigma)),
      ticktype = "detailed",
      expand = 0.7,
      cex.axis=0.6,
      cex.lab=0.6)

entropy_loss(solve(S),Sigma)

data.frame(it=2*(1:length(quad_loss)),
           quadLoss=unlist(quad_loss),
           entropyLoss=unlist(entrpy_loss)) %>%
      ggplot(.,aes(x=it)) +
      geom_line(aes(y=quadLoss)) +
      theme_minimal() +
      xlab("iteration") +
      ylab("") +
      ggtitle(expression(paste(L[2],"(",Sigma,",",hat(Omega),")")))
ggsave(file=file.path(getwd(),"reports-1-simulations",
                      "img",
                      paste0(Type.,"-quadratic-loss-vs-iteration.png")))
#####################################################################################
par(mfrow=c(1,1))
## Select spatials points for prediction
MMgrid <- mmodel.pred.frame2d(myGrid$l_s, 
                              myGrid$m_s,
                              nsegl, nsegm,
                              nsegl/2, nsegm/2,
                              bdeg,
                              pordl, pordm,
                              ngrid=100)

Fit.grid <- matrix(MM$M %*% b[[which.min(unlist(quad_loss))]],
                   100, 100,byrow=TRUE)
persp(x=seq(0,1,length.out=100),
      y=seq(0,1,length.out=100),
      z=Fit.grid,
      phi=15,
      theta=35,
      xlab="l",
      ylab="m",
      col = "light blue",
      zlab="",
      main=expression(paste(hat(phi))),
      ticktype = "detailed",
      expand = 0.7,
      cex.axis=0.6)
persp(x=seq(0,1,length.out=100),
      y=seq(0,1,length.out=100),
      z=matrix(Phi,nrow=100,ncol=100,byrow=TRUE),
      phi=15,
      theta=35,
      xlab="l",
      ylab="m",
      col = "light blue",
      zlab="",
      main=expression(paste(phi)),
      ticktype = "detailed",
      expand = 0.7,
      cex.axis=0.6)

indexl <- 1 * (idx == 2)
if (pordl > 0) indexl[2:pordl] <- 1
Fitl <- matrix(MMgrid$M %*% (indexl * b), 100, 100)
Phi1.grid <- matrix(Phi1,nrow=100,ncol=100)
l_facet <- Phi1.grid[-1, -1] + Phi1.grid[-1, -M] + Phi1.grid[-M, -1] + Phi1.grid[-M, -M]

# Recode facet z-values into color indices
par(mfrow=c(1,2))
facetcol <- cut(l_facet, nbcol)
persp(x=seq(0,1,length.out=100),
      y=seq(0,1,length.out=100),
      z=Fitl,
      phi=15,
      theta=15,
      xlab="l",
      ylab="m",
      col = color[facetcol],
      zlab="",
      ticktype = "detailed",
      expand = 0.7,
      main=expression(paste(hat(phi)[1])))
persp(x=seq(0,1,length.out=100),
      y=seq(0,1,length.out=100),
      z=Phi1.grid,
      phi=15,
      theta=15,
      xlab="l",
      ylab="m",
      col = color[facetcol],
      zlab="",
      ticktype = "detailed",
      expand = 0.7,
      main=expression(paste(phi[1])))

data.frame(l=rep(seq(0,1,length.out = 100),2),
           phi=c(Phi1.grid[,1],Fitl[,1]),
           func=c(rep("true",100),
                  rep("estimated",100))) %>%
      ggplot(data=.,aes(x=l,y=phi)) +
      geom_line(aes(colour=func)) +
      guides(colour=guide_legend(title="")) +
      theme_minimal() +
      ylab("") +
      ggtitle(expression("estimated vs true "~phi[1]))


indexm <- 1 * (idx == 3)
if ( pordm > 0 ) {
      indexm[pordl*(1:(pordm-1))+1] <- 1      
}
Fitm <- matrix(MMgrid$M %*% (indexm * b), 100, 100)
data.frame(m=rep(seq(0,1,length.out = 100),2),
           phi=c(rep(0,nrow(Fitm)),Fitm[1,]),
           func=c(rep("true",100),
                  rep("estimated",100))) %>%
      ggplot(data=.,aes(x=m,y=phi)) +
      geom_line(aes(colour=func)) +
      guides(colour=guide_legend(title="")) +
      theme_minimal() +
      ylab("") +
      ylim(-0.0001,0.000001) +
      ggtitle(expression("estimated vs true "~phi[2]))


indexlm <- 1 - (indexl + indexm)
indexlm[1] <- 0
Fitlm <- matrix(MMgrid$M %*% (indexlm * b), 100, 100)

# indexlma <- 1 * (idx == 4)
# #indexlma[4] <- 1
# Fitlma <- matrix(MMgrid$M %*% (indexlma * b), 100,100)
# 
# indexlmb <- 1 * (idx == 5)
# Fitlmb <- matrix(MMgrid$M %*% (indexlmb * b), 100,100)
# 
# indexlmc <- 1 * (idx == 6)
# indexlmc[4] <- 1
# Fitlmc <- matrix(MMgrid$M %*% (indexlmc * b), 100, 100)



cat("----------- PS-ANOVA-Schall -------------\n")
cat("Sample size           ", length(yVec), "\n")
cat("sigma                 ", sigma., "\n")
cat("sigma (estimated)     ", round(sig2^0.5, 3), "\n")
if (div == 1) {
      cat("No nested bases \n")
} else {
      cat("Nested bases for interactions, with nseg1=", nseg1, ", nseg2=", nseg2, "and div=", 
          div, "\n")
}
cat("------------------------------------------\n")

cat("-----------------------------------------\n")
cat("   AIC", AIC., "\n")
cat("          fixed    fx1       fx2     fx1:x2     x1:fx2    fx1:fx2\n")
cat("    ed", sprintf("%8.3f", ed), "\n")
cat("sum ed", sum(ed), "\n")
cat("MSE", MSE., "\n")
cat("-----------------------------------------\n")



ZLIM. <- range(c(-0.5, .5))

subset(fullGrid, (s<t) & (s >= 0) & (s < 1) & (t > 0) & (t < 1) )$l
image.plot(matrix(data=truePhi,nrow=100,ncol=100,byrow=TRUE),xlab="l",ylab="m",col=color)
points(subset(fullGrid, !((s<t) & (s >= 0) & (s < 1) & (t > 0) & (t < 1)) )$l,
       subset(fullGrid, !((s<t) & (s >= 0) & (s < 1) & (t > 0) & (t < 1)) )$m,
       pch=20,cex=2.5)

windows()
jet.colors <- colorRampPalette( c("blue", "green") )
# Generate the desired number of colors from this palette
nbcol <- 100
color <- jet.colors(nbcol)
par(mfrow = c(1, 1))
image.plot(Fit.grid, xlab = "l", ylab = "m", zlim = ZLIM.,col=color)
points(subset(fullGrid, !((s<t) & (s >= 0) & (s < 1) & (t > 0) & (t < 1)) )$l,
       subset(fullGrid, !((s<t) & (s >= 0) & (s < 1) & (t > 0) & (t < 1)) )$m,
       pch=20,cex=2.5)
title("Fitted surface")



drape.plot(seq(0,1,length.out=100),
           seq(0,1,length.out=100), Fitl,
           phi=25,
           theta=25,
           expand = 0.5, 
           border = NA, shade = 0.1, add.legend = T, xlab = "l", ylab = "m", 
           zlab = "Phi1", main = "(a) main effect of l",
           col = color,
           ticktype = "detailed")

drape.plot(seq(0,1,length.out=100),
           seq(0,1,length.out=100),
           Fitl,
           phi=25,
           theta=25,
           expand = 0.5, 
           xlab = "l", ylab = "m", 
           zlab = expression("phi[l]"), main = "(a) main effect of l: nonparametric part",
           ticktype = "detailed",
           col=color)

plot(seq(0,100,length.out=100),
     matrix(Phi1,nrow=100,ncol=100)[,1],
     type="l",
     xlab="l",
     ylab=expression("phi[l]"))
plot(seq(0,100,length.out=100),
     matrix(c(seq(0,1,length.out=100),
              seq(0,1,length.out=100)),nrow=100,ncol=2,byrow=FALSE) %*% b[2:3],
     type="l",
     xlab="l",
     ylab=expression("phi[l]"))

drape.plot(seq(0,1,length.out=100),
           seq(0,1,length.out=100),
           Fitlnp,
           phi=25,
           theta=25,
           expand = 0.5, 
           border = NA, shade = 0.1, add.legend = T, xlab = "l", ylab = "m", 
           zlab = "Phi1", main = "(a) main effect of l: nonparametric part",
           col = color,
           ticktype = "detailed")



drape.plot(seq(0,1,length.out=100), seq(0,1,length.out=100), Fitm, theta = -30, phi = 30, expand = 0.85, 
           border = NA, shade = 0.1, add.legend = T, xlab = "l", ylab = "m", 
           zlab = "Y", main = "(b) main effect of m", zlim = ZLIM., col=color)
drape.plot(ml, mm, Fitlma, theta = -30, phi = 30, expand = 0.85, 
           border = NA, shade = 0.1, add.legend = T, xlab = "l", ylab = "m", 
           zlab = "Y", main = "(c) f(x1):x2 interaction", zlim = ZLIM.)
drape.plot(ml, mm, Fitlmb, theta = -30, phi = 30, expand = 0.85, 
           border = NA, shade = 0.1, add.legend = T, xlab = "l", ylab = "m", 
           zlab = "Y", main = "(d) x1:f(x2) interaction", zlim = ZLIM.)
drape.plot(ml, mm, Fitlmc, theta = -30, phi = 30, expand = 0.85, 
           border = NA, shade = 0.1, add.legend = T, xlab = "l", ylab = "m", 
           zlab = "Y", main = "(e) f(x1,x2) interaction", zlim = ZLIM.)
drape.plot(ml, mm, Fitlma + Fitlmb + Fitlmc, theta = -30, phi = 30, 
           expand = 0.85, border = NA, shade = 0.1, add.legend = T, 
           xlab = "x1", ylab = "x2", zlab = "True", main = "(f) sum of all interactions", 
           zlim = ZLIM.)


drape.plot(ml, mm, Fit.grid, theta = -30, phi = 30, 
           expand = 0.85, border = NA, shade = 0.1, add.legend = T, 
           xlab = "l", ylab = "m", # main = "(f) sum of all interactions", 
           zlim = ZLIM.)

drape.plot(seq(0,1,length.out=100), seq(0,1,length.out=100),
           matrix(truePhi,nrow=100,ncol=100,byrow=TRUE),
           theta = -30, phi = 30, 
           expand = 0.85, border = NA, shade = 0.1, add.legend = T, 
           xlab = "l", ylab = "m", # main = "(f) sum of all interactions", 
           zlim = ZLIM.)
