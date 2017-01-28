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