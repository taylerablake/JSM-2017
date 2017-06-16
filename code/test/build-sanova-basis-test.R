

for (helper.file in list.files(file.path(getwd(),"code","lib"))[list.files(file.path(getwd(),"code","lib"))!="build-bd-penalty.R"]){
  source(file.path(getwd(),"code","lib",helper.file))
}


cl <- makeCluster(detectCores()-3)
registerDoParallel(cl)
clusterCall(cl,function() {
  .libPaths("~/Rlibs/lib")
  library(doBy)
  library(lattice)
  library(MASS)
  library(magrittr)
  library(rlist)
  library(Matrix)
  library(plyr)
  library(stringr)
  library(dplyr)
  library(doParallel)
  library(splines)
})


m <- 20
N <- 50
grid <- build_grid(20)




bPars <- rbind(c(0,1,15,3,100,3),
               c(0,1,10,3,100,3))

Bl <- bsplbase(grid$l/max(grid$l),
               bPars[1,],outer.okay = FALSE)$base
Bm <- bsplbase(grid$m/max(grid$m),
               bPars[2,],outer.okay = FALSE)$base

n1 <- ncol(Bl)
n2 <- ncol(Bm)

Blm <- kronecker(Bm,
                t(as.vector(rep(1,n1)))) * kronecker(t(as.vector(rep(1,n2))),
                                                           Bl)
B <- cbind(as.vector(rep(1,nrow(Bl))),Bl,Bm,Blm)


## validity testing ##################################################################

if(ncol(Blm) != n1*n2) print("Error: interaction basis Blm is not of dimension n1xn2")
if(ncol(B) != (1 + n1 + n2 + (n1*n2))) print("Error: interaction basis Blm is not of dimension 1 + n1 + n2 + n1xn2")


rankMatrix(Bl)
rankMatrix(Bm)
rankmatrix(kronecker(Bm,t(as.vector(rep(1,n1)))))
rankMatrix(kronecker(t(as.vector(rep(1,n2))),Bl))
rankMatrix(kronecker(t(as.vector(rep(1,n2))),Bl)*kronecker(Bm,t(as.vector(rep(1,n1)))))
rankMatrix(Blm)
rankMatrix(B)


######################################################################################



dl=3
Dl <- diff(diag(ncol(Bl)),
           differences=dl)
Ul <- eigen(t(Dl)%*%Dl)$vectors
Ul1 <- Ul[,1:rankMatrix(t(Dl)%*%Dl)]
Ul0 <- Ul[,-(1:rankMatrix(t(Dl)%*%Dl))]

Xl <- Bl%*%Ul0
Zl <- Bl%*%Ul1


dm=3
Dm <- diff(diag(ncol(Bm)),
           differences=dm)
Um <- eigen(t(Dm)%*%Dm)$vectors
Um1 <- Um[,1:rankMatrix(t(Dm)%*%Dm)]
Um0 <- Um[,-(1:rankMatrix(t(Dm)%*%Dm))]

Xm <- Bm%*%Um0
Xm <- cbind(rep(1,),
            poly(grid$m/max(grid$m), degree=dm-1)
Zm <- Bm%*%Um1
rankMatrix(bdiag(0,t(Dl)%*%Dl,2.2323*t(Dm)%*%Dm, kronecker(3.231*t(Dm)%*%Dm,4.6*t(Dl)%*%Dl)))


rankMatrix(bdiag(0,t(Dl)%*%Dl, kronecker(3.231*t(Dm)%*%Dm,4.6*t(Dl)%*%Dl)))
rankMatrix(bdiag(0,kronecker(3.231*t(Dm)%*%Dm,4.6*t(Dl)%*%Dl)))
rankMatrix(t(Dl)%*%Dl)
