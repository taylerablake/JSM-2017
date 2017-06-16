bsplbase<- function(x, bpars, outer.okay=TRUE){
      # Compute a B-spline basis
      # Input:
      #   x: abcissae
      #   bpars: B-spline parameters: xmin, xmax, nseg, degree (= one less than "order")
      # Output:
      #   base: matrix with nrow = length(x) and nseg + degree columns
      #
      # Paul Eilers, 2000
      require(splines)
      dx <- (bpars[2] - bpars[1])/bpars[3]
      knots <- seq(bpars[1] - bpars[4] * dx, bpars[2] + bpars[4] * dx, by = 
                         dx)
      base <- as.matrix(spline.des(knots, x, bpars[4] + 1, 0 * x,outer.ok=outer.okay)$design)
      list(base=base, knots=knots)
}