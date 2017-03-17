#' Create dummy data x
#' 
#' Function used for simulating data for test runs
#' 
#' @param N number of observations
#' @param Nx number of explanatory variables
#' 
#' @export
data_gen <- function(N,Nx=3){
  
  x <- matrix(runif(Nx*N),ncol=Nx)
  y <- data.frame(y1 = runif(N), y2 = factor(sample(letters[1:2],N,replace=TRUE)))

  return(list(x=x,y=y))
}

data_gen_inf <- function(N,Nx=3){
  x <- matrix(runif(Nx*N),ncol=Nx)
  prbs1 <- x[,1]/ (max(x[,1])-min(x[,1]))
  prbs2 <- x[,2]/ (max(x[,2])-min(x[,2]))
  
  yf <- letters[1:3][as.numeric(prbs1 > runif(N,0.4,0.6)) + as.numeric(prbs2 < runif(N,0.4,0.6)) + 1]
  y <- data.frame(y1 = runif(N), y2 = factor(yf)) # factor(letters[1:2][1+(runif(N,0.4,0.6) > prbs)]))
  
  return(list(x=x,y=y))
}

data_gen_tree <- function(N,Ne=0){
  x <- matrix(runif(N*(4+Ne)),ncol=4+Ne)
  yf <- character(N)
  yf[x[,1] < 0.5 & x[,2] >= 0.4] <- 'a'
  yf[x[,1] < 0.5 & x[,2] < 0.4] <- 'b'
  yf[x[,1] >= 0.5 & x[,3] < 0.7] <- 'c'
  yf[x[,1] >= 0.5 & x[,3] >= 0.7 & x[,2] < 0.5] <- 'd'
  yf[x[,1] >= 0.5 & x[,3] >= 0.7 & x[,2] >= 0.5] <- 'e'
  yc <- rnorm(N,as.numeric(x[,4] >= 0.4)*3+1,1)
  y <- data.frame(y1=yc,y2=factor(yf))
  
  return(list(x=x,y=y))
}