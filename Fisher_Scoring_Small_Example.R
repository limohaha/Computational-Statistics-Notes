fisher_scoring <- function(y=NULL,x=NULL,initial_value=NULL){
  n <- nrow(x)
  p <- ncol(x)
  theta <- initial_value
  score <- rep(0, p+1)
  ii <- 0
  diff <- 0.01
  while(ii<1000 & diff>1e-05){
    # derivatives w.r.t. betas
    score[1:p] <- (1/theta[p+1])*t((y-x%*%theta[1:p]))%*%x
    # derivatives w.r.t. sigma^2
    score[p+1] <- -(n/(2*theta[p+1]))+(1/(2*theta[p+1]^2))*crossprod(y-x%*%theta[1:p])
    # new
    hessMat <- matrix(0,ncol=p+1,nrow=p+1)
    for(j in 1:n)
    {
      # Estimate derivative of likelihood for each observation
      estVec <- c((1/theta[p+1]) * 
                    t((y[j] - x[j,]%*%theta[1:p])) %*% x[j,],
                  -(n/(2*theta[p+1])) + (1/(2*theta[p+1]^2)) *
                    crossprod(y[j] - x[j,] %*% theta[1:p]))
      # Add them up as suggested to get an estimate of the Hessian.
      hessMat <- hessMat + estVec%*%t(estVec)
    }
    theta1 <- theta + MASS::ginv(hessMat) %*% score
    
    diff <- sqrt(sum((theta1-theta)^2))
    theta <- theta1
    ii <- ii+1
  }
  
  cat("======= Number of Iteration =========\n")
  print(ii)
  cat("======= Diff on theta between iterations ======\n")
  print(diff)
  if (ll>=1000) cat("Not converge in thea")
  return(theta1)
}


# application
library(MASS)
# simulate y and x
x <- matrix(rnorm(1000), ncol = 2)
y <- 2 + x %*% c(1,3) + rnorm(500)
x <- cbind(1,x)

# esitimation
initial_value <- runif(ncol(x)+1)
initial_value
theta <- fisher_scoring(y=y,x=x,initial_value)
theta


