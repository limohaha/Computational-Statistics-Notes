Par1Sim <- function(mu=NULL, phi=NULL, sig2=NULL, alpha=NULL, n=NULL, T=NULL){
  
  y   <- rep(0, n)
  et  <- rep(0, n)
  Z   <- rnorm(n)
  
  y[1] <- mu[1]+alpha*et[1]
  
  for(i in 2:n){
    nu    <- i - T*floor((i-1)/T)
    et[i] <- phi[nu]*et[i-1] + sig2[nu]*Z[i]
    y[i]  <- mu[nu] + alpha*i + et[i]
  }
  
  return(y)
}




# simulation parameters
mu0   <- c(-0.61,  0.99,  2.35,  4.91, 
            8.74, 12.15, 15.51, 15.47, 
           12.79,  7.82,  2.32, -0.25)
phi0  <- c(0.272, 0.284, 0.478, 0.286, 
           0.335, 0.279, 0.245, 0.137, 
          -0.127, 0.082, 0.196, 0.214)
sig20 <- c(2.712, 2.748, 1.871, 1.717, 
           2.474, 2.403, 2.569, 1.910, 
           2.826, 2.488, 2.394, 2.256)
alpha0 <- 0
T      <- 12
d      <- 100
K      <- 3
d      <- d+K*2 # for edge effect
n      <- T*d

mydat <- Par1Sim(mu=mu0, phi=phi0, sig2=sig20, alpha=alpha0, n=n, T=T)
plot(as.ts(mydat))