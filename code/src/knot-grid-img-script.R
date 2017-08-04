


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


# B-spline parameters
bdeg = 3  # B-spline degree

nsegl = 12  # number of inner knots to build B-spline bases for x1 and x2
nsegm = 12

Bl <- bspline(myGrid$l_s, min(myGrid$l_s),
              max(myGrid$l_s),
              nsegl, bdeg)$B
lKnots <- bspline(myGrid$l_s, min(myGrid$l_s),
                 max(myGrid$l_s),
                 nsegl, bdeg)$knots



Bm <- bspline(myGrid$m_s, min(myGrid$m_s),
              max(myGrid$m_s),
              nsegm, bdeg)$B
mKnots <- bspline(myGrid$m_s, min(myGrid$m_s)-0.01,
                  max(myGrid$m_s)+0.01,
                  nsegm, bdeg)$knots

B. <- kronecker(Bm,
                t(as.vector(rep(1,ncol(Bl))))) * kronecker(t(as.vector(rep(1,ncol(Bm)))),
                                                           Bl)




      
expand.grid(t=seq(min(lKnots),max(lKnots),length.out=length(lKnots)),
            s=seq(min(lKnots),max(lKnots),length.out=length(lKnots))) %>%
      ggplot(.,aes(t,s)) +
      geom_point(size=0.5, alpha=0.1) +
      theme_minimal() +
      xlab("t") +
      ylab("s") +
      scale_x_continuous(breaks = seq(0,1.2,by=0.2)) +
      scale_y_continuous(breaks = seq(0,1.2,by=0.2)) +
      geom_polygon(data=data.frame(s=c(0,0,1,1),t=c(0,0,1,0)),
                   aes(s,t),
                   colour="light grey",
                   fill="light grey",
                   alpha=0.4) +
      theme(panel.grid.minor=element_blank()) 
      

which(colSums(B. != 0) == 0) 
knot_grid <- expand.grid(m=mKnots,l=lKnots)[,2:1]
ggplot(knot_grid,aes(l,m)) +
      geom_point(size=0.5) +
      theme_minimal() +
      xlab("l") +
      ylab("m") +
      scale_x_continuous(breaks = seq(0,1.2,by=0.2)) +
      scale_y_continuous(breaks = seq(0,1.2,by=0.2)) +
      geom_polygon(data=data.frame(l=c(0,1,0),m=c(0,0.5,1)),
                   aes(l,m),
                   colour="light grey",
                   fill="light grey",
                   alpha=0.4) +
      theme(panel.grid.minor=element_blank())



ggplot(knot_grid,aes(l,m)) +
      geom_point(size=0.5) +
      theme_minimal() +
      xlab("l") +
      ylab("m") +
      geom_segment(aes(x = 0, y = 0, xend = 1, yend = 0.5), colour = "forestgreen") +
      geom_segment(aes(x = 0, y = 1, xend = 1, yend = 0.5), colour = "forestgreen") +
      geom_segment(aes(x = 0, y = 0, xend = 0, yend = 1), colour = "forestgreen") 

keep <- (knot_grid$m >= 0.5) & (knot_grid$m <= max(knot_grid$l)-(0.5*knot_grid$l)) |
      (knot_grid$m <= 0.5) & (knot_grid$m >= min(knot_grid$l)+(0.5*knot_grid$l))
knot_grid <- transform(knot_grid,keep=factor(keep))
knot_grid <- data.frame(knot_grid,
                        expand.grid(m_index=1:length(interior.knots.m),
                                    l_index=(1:length(interior.knots.l)))[,2:1])