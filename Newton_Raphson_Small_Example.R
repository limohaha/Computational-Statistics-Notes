# Newton-Raphson Algorithm Small Example
# Censored exponentially distributed observations

C <- 10  # censored constant
theta_true <- 1/5

set.seed(12345)

# simulate censored data with exp(1/5)
Y <- rexp(10, theta_true)
# detecting which observation is cencored
R <- ifelse(Y>C, 0, 1)
Y[R==0]=C

# loglikelihood
loglikelihood <- function(theta, Y, R){
  # m: number of uncencored osbervations
  m <- sum(R==1)
  loglihd <- m*log(theta)-theta*sum(Y)
  attr(loglihd, "gradient") <- m/theta-sum(Y)
  attr(loglihd, "hession") <- -m/(theta^2)
  return(loglihd)
}

# Newton-Raphson
Newton_Rp <- function(theta, Y, R){
  l <- loglikelihood(theta, Y, R)
  temp <- theta - attr(l,"gradient")/attr(l, "hession")
  return(temp)
}


# theta0 <- runif(1,0,1)
# if choose the theta value far away from the true value, may produce Nan in log from log(theta)
theta0 <- 0.4
theta0
theta <- Newton_Rp(theta0, Y, R)
theta 
theta <- Newton_Rp(theta, Y, R)
theta 
theta <- Newton_Rp(theta, Y, R)
theta
theta <- Newton_Rp(theta, Y, R)
theta