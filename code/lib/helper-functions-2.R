

## D = diag(sigsq1,...,sigsqM)

log_lik <- function( y_arg, T_arg, D_sigma_arg ){
      D_arg_inv <- diag(1/diag(D_sigma_arg))
      -0.5*( - log( det(t(T_arg) %*% (D_arg_inv) %*% T_arg) ) +
                   t(y_arg) %*% t(T_arg) %*% (D_arg_inv^2) %*% T_arg %*% y_arg)

}



