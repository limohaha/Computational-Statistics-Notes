innovations.algorithm <- function(acvf,n.max=length(acvf)-1)
  {
    thetas <- vector(mode="list",length=n.max)
    vs <- rep(acvf[1],n.max+1)
    for(n in 1:n.max)
      {
        thetas[[n]] <- rep(0,n)
        thetas[[n]][n] <- acvf[n+1]/vs[1]
        if(n>1)
          {
            for(k in 1:(n-1))
              {
                js <- 0:(k-1)
                thetas[[n]][n-k] <- (acvf[n-k+1] - sum(thetas[[k]][k-js]*thetas[[n]][n-js]*vs[js+1]))/vs[k+1]
              }
          }
        js <- 0:(n-1)
        vs[n+1] <- vs[n+1] - sum(thetas[[n]][n-js]^2*vs[js+1])
      }
    structure(list(vs=vs,thetas=thetas))
  }

### > AR2.acvf <- as.vector(ARMAacf(ar=c(3/4,-1/2),lag=3))*16/9
### > (results <- innovations.algorithm(AR2.acvf))
### $vs
### [1] 1.777778 1.333333 1.000000 1.000000
### 
### $thetas
### $thetas[[1]]
### [1] 0.5
### 
### $thetas[[2]]
### [1]  0.750 -0.125
### 
### $thetas[[3]]
### [1]  0.75000  0.06250 -0.34375
