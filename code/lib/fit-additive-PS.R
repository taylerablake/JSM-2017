

fit_additive_PS <- function(yVec,
                            N,
                            U,
                            D,
                            P_l,
                            lambda_l,
                            P_m,
                            lambda_m,
                            lambda_ridge){
      
      construct.block2 <- function(A1, A2, A4) {
            block <- rbind(cbind(A1, A2), cbind(t(A2), A4))
            return(block)
      }
      Pen <- construct.block2(lambda_l*Pl,
                       matrix(data=0,nrow=nrow(Pl),ncol=ncol(Pm)),
                       lambda_m*Pm) %>% cbind(rep(0,nrow(Pl)+nrow(Pm)),.) %>%
            rbind(rep(0,1+ncol(Pl)+ncol(Pm)),.)
      ridge_pen <- diag(ncol(Pl)+ncol(Pm))
      ridge_pen[1,1] <- 0
      Pen <- rbind(Pen, lambda_ridge*ridge_pen)
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
            #w <- rep(1, length(y_vec))
            w <- rep(1/diag(D)[-1],N)
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
