
# Load functions

source(file.path(getwd(),"code/lib/helper-functions-2.R"))
source(file.path(getwd(),"code/lib/fit_cholesky_PS.R"))
source(file.path(getwd(),"code/lib/help-functions.R"))
source(file.path(getwd(),"code/lib/quadratic-loss.R"))
source(file.path(getwd(),"code/lib/entropy-loss.R"))
source(file.path(getwd(),"code/lib/build-grid.r"))


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




fit_list <- lapply(PS_fit_sim,function(l){
      Phi <- B.[,-no_support] %*% l$fit$coef      
      T_mat <- diag(M)
      T_mat[lower.tri(T_mat)] <- -Phi
      list(fit=Phi,
           T_mod=T_mat)
})


omega_list <- lapply(fit_list,function(l){
      
      Omega_hat <- t(l$T_mod)%*%diag(1/diag(C)^2)%*%l$T_mod
      
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


image.plot(y=log(sort(unique(lambdas$lam_l))),
           x=log(sort(unique(lambdas$lam_m))[-1]),
           z=matrix(aic[lambdas$lam_m>min(lambdas$lam_m)],
                    ncol=length(unique(lambdas$lam_l)),
                    nrow=length(unique(lambdas$lam_m))-1,
                    byrow=TRUE),
           ylab=expression("log "~lambda["l"]),
           xlab=expression("log "~lambda["m"]))



persp(y=log(sort(unique(lambdas$lam_l))),
      x=log(sort(unique(lambdas$lam_m))[-1]),
      z=matrix(aic[lambdas$lam_m>min(lambdas$lam_m)],
               ncol=length(unique(lambdas$lam_l)),
               nrow=length(unique(lambdas$lam_m))-1,
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
      x=log(sort(unique(lambdas$lam_m))[-1]),
      z=matrix(ed[lambdas$lam_m>min(lambdas$lam_m)],
               ncol=length(unique(lambdas$lam_l)),
               nrow=length(unique(lambdas$lam_m))-1,
               byrow=TRUE),
      phi=15,
      theta=90,
      ylab="log lambda l",
      xlab="log lambda m",
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

#png(file.path(getwd(),"TeX","img","compound-symmetry-estimated-cholesky.png"))
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
      zlab="")



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
      z=omega_list[[which.min(entrpy_loss)]]%*%Sigma,
      phi=15,
      theta=15,
      xlab="s",
      ylab="t",
      col = "light blue",
      zlab="",
      main=expression(paste(hat(Sigma)^{-1},Sigma)),
      ticktype = "detailed",
      expand = 0.7,
      cex.axis=0.6,
      cex.lab=0.6)

persp(x=seq(0,1,length.out=M),
      y=seq(0,1,length.out=M),
      z=omega_list[[900]],
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

png(file.path(getwd(),"TeX","img","compound-symmetry-true-covariance.png"))
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
      cex.lab=0.6,
      zlim=c(0.4,1.1))
dev.off()

png(file.path(getwd(),"TeX","img","compound-symmetry-estimated-covariance.png"))
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



persp(x=seq(0,1,length.out=M),
      y=seq(0,1,length.out=M),
      z=solve(omega_list[[900]]),
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

