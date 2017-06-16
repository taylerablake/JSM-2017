fit_cholesky_PS <- function(Y,
                            y_aug,
                            U,
                            D,
                            P_l,
                            lambda_l,
                            P_m,
                            lambda_m,
                            lambda_ridge,
                            P_l_shape, lambda_l_shape){
  N.subjects <- nrow(Y)
  if (!is.null(y_aug)) {
    yVec <- as.vector(c(y_aug,
                        t(Y[,-1])))
  }
  
  if (is.null(y_aug)) {
    yVec <- as.vector(t(Y[,-1]))
  }
  Pen <- rbind(lambda_l*P_l,
               lambda_m*P_m)
  if (lambda_ridge > 0) {
    Pen <- rbind(Pen,lambda_ridge*diag(ncol(P_l)))
  }

  if (lambda_l_shape > 0) {
    Pen <- rbind(Pen,lambda_l_shape*P_l_shape)
  }
  n.col <- ncol(U)
  nix <- rep(0,nrow(Pen))
  
  nix.ridge <- rep(0, ncol(U))
  coef.est <- rep(1, ncol(U))
  
  mu <- rep(mean(yVec), length(yVec))
  it <- 0
  repeat {
    if(it == 0) {
      eta <- mu
    }
    
    it <- it + 1
    if(it > 25)
      break
    
    mu <- eta
    h.prime <- 1
    if (is.null(y_aug)) {
      w <- rep(1, length(yVec))
    }
        
    if (!is.null(y_aug)) {
      w <- c(rep(0,length(y_aug)),
             rep(1/diag(D)[-1]^2,N.subjects))      
    }

    u <- (yVec - mu)/h.prime + eta
    
    startTS <- Sys.time()
    f <- lsfit(rbind(U,Pen), c(u, nix), 
               wt = c(w, nix + 1) *c(w, nix + 1),
               intercept = F)
    endTS <- Sys.time()
    endTS-startTS
    
    coef.old <- coef.est
    coef.est <- as.vector(f$coef)
    
    d.coef <- max(abs((coef.est[coef.old>0] - coef.old[coef.old>0])/coef.old[coef.old>0]))
    if(d.coef < 1e-008)
      break
    print(c(it, d.coef))
    eta <- U %*% coef.est
  }
  
  if(it > 24) {
    warning(paste("parameter estimates did NOT converge in 25 iterations"
    ))
  }      
  
  H <- hat(f$qr, intercept = F)[1:(M-1)]
  trace <- eff.dim <- sum(H)
  
  list(coef=f$coefficients,
       H=H,
       eff.dim=eff.dim,
       trace=trace)
}