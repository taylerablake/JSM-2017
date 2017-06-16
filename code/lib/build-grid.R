
## arguments: M - number of repeated measurements on each subject
## value:     grid - dataframe containing the coordinates of observations
##                   in terms of s,t and l,m

build_grid <- function(m) {
   require(magrittr)
   obs_grid <- expand.grid(t=1:m,s=1:m) %>%
     subset(.,t>s) %>%
     transform(.,l=t-s,
               m=0.5*(t+s))
   obs_grid
}


