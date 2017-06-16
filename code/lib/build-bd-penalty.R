

build_bd_pen <- function (lambda_list,
                          basis_list,
                          diff_ord_list) {
  D1 <- diff(diag(ncol(B1)),
             differences=dl)
  D2 <- diff(diag(ncol(Bm)),
             differences=dm)
  D3 <- diff(diag(ncol(B3)),
             differences=dl)
  D4 <- diff(diag(ncol(B4)),
             differences=dl)
  D5 <- diff(diag(ncol(B5)),
             differences=dm)
  D6 <- diff(diag(ncol(B6)),
             differences=dm)
  lam_list <- as.list(c(lambda_l,lambda_m,lambda_l, lambda_l, lambda_m, lambda_m))
  if (length(diff_ord_list) != length(basis_list)) {
    stop( "you must provide a differencing order for each additive basis" )
  }
  Diff_mat_list <- list()
  for (this_basis in 1:length(basis_list) ) {
    Diff_mat_list[[this_basis]] <- diff(diag(ncol(basis_list[[this_basis]])),
                                        differences=diff_ord_list[[this_basis]])
  }
  if (length(lambda_list) != length(basis_list)) {
    stop( "you must provide a penalty parameter for each additive basis" )
  }
  for(  )
  P1 <- cbind(lambda1*D1, matrix(data=0, 
                                 nrow=nrow(D1),
                                 ncol=))
  Pen <- rbind(## P1 ####################################################
               cbind(rep(0,nrow(Dl)),
                     lambda_pair$lam_l1*Dl,
                     matrix(data=0,nrow=nrow(Dl),
                            ncol=((2*ncol(Dm))+ncol(Dl)))),
               ## P2 ####################################################
               cbind(rep(0,nrow(Dl)),
                     matrix(data=0,nrow=nrow(Dm),
                            ncol=ncol(Dl)),
                     lambda_pair$lam_m1*Dm,
                     matrix(data=0,nrow=nrow(Dm),ncol=ncol(Dl)+ncol(Dm))),
               ## P3 ####################################################
               cbind(rep(0,nrow(Dl)),
                     matrix(data=0,nrow=nrow(Dm),
                            ncol=ncol(Dl)+ncol(Dm)),
                     lambda_pair$lam_l2*Dl,
                     matrix(data=0,nrow=nrow(Dm),ncol=ncol(Dm))),
               ## P4 ####################################################
               cbind(rep(0,nrow(Dl)),
                     matrix(data=0,nrow=nrow(Dm),
                            ncol=((2*ncol(Dl))+ncol(Dm))),
                     lambda_pair$lam_m2*Dm))
  
}
