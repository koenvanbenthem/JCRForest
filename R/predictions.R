predict.jcr_forest <- function(forest,newdata){
  # tests:
  # 1. number of columns
  # possibly: is x numeric?
  nsample <- nrow(newdata)
  pred.out <- .C("forestpred",
                 x=newdata,
                 yc=numeric(nsample),
                 yf=integer(nsample),
                 nsample=as.integer(nsample),
                 ldaughter=forest$ldaughter,
                 rdaughter=forest$rdaughter,
                 ntree = as.integer(forest$ntree),
                 nrnodes = as.integer(forest$nrnodes),
                 node_var = forest$node_var,
                 node_xvar = forest$node_xvar,
                 PACKAGE="JCRForest"
  )
  return(pred.out)
}
