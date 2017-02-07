#' Main function
#' 
#' Function used for fitting mixed random forests.
#' 
#' @param x Data frame containing explanatory variables
#' @param y Data frame containing 2 response variables, one being discrete (factor), the other continuous 
#' @param mtry Number of features examined per tree
#' @param Ntree Number of trees in the forest
#' 
#' @useDynLib JCRForest
#' 
#' @export
jcr_forest <- function(x,y,mtry,Ntree,minsize=5,bla=1:20){
  
  # Data checking and preparation
  if (nrow(x) != nrow(y)) stop("x and y must have the same number of rows\n")
  if (nrow(x) == 0) stop("Data has 0 rows")
  if (ncol(y) != 2 ) stop("y must contain two columns")
  if (!any(c('numeric','integer') %in% sapply(y,class))) stop("y must contain a continuous variable")
  if (!c('factor') %in% sapply(y,class)) stop("y must contain a discrete variable (factor)")
  if (any(is.na(x)) | any(is.na(y))) stop("missing values are not allowed")
  if (mtry < 1 | mtry > ncol(x)) stop("mtry must be between 1 and the number of explanatory variables")
  
  yc <- y[,sapply(y,class) %in% c('numeric','integer')]
  yf <- as.integer(y[,sapply(y,class) %in% c('factor')])
  
  nrnodes <- 10
  nclass <- length(unique(yf))
  cat(Ntree,"\n\n")
  # Forest building
  rfout <- .C("build_jcr_forest",
              x=x,
              yc=yc,
              yf=yf,
              nclass=as.integer(nclass),
              nsample=as.integer(nrow(x)),
              nvar=as.integer(ncol(x)),
              mtry=as.integer(mtry),
              ntree=as.integer(Ntree),
              nrnodes=as.integer(nrnodes),
              ldaughter=matrix(integer(nrnodes * Ntree),ncol=Ntree),
              rdaughter=matrix(integer(nrnodes * Ntree),ncol=Ntree),
              node_status=matrix(integer(nrnodes * Ntree),ncol=Ntree),
              node_var=matrix(integer(nrnodes * Ntree),ncol=Ntree),
              node_xvar=matrix(integer(nrnodes * Ntree),ncol=Ntree),
              dum_vect=as.numeric(bla),
              dum_long=as.integer(length(bla)),
              dum_ind=as.integer((1:length(bla))*2),
              PACKAGE="JCRForest")
  
  return(rfout)
}