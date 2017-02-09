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
  prbs <- x[,1]/ (max(x[,1])-min(x[,1]))
    y <- data.frame(y1 = runif(N), y2 = factor(letters[1:2][1+(runif(N) > prbs)]))
  
  return(list(x=x,y=y))
}