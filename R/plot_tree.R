arrows_label <- function(x1,y1,x2,y2,lbl,...){
  text(0.5*x1+0.5*x2,0.5*y1+0.5*y2,lbl)
  arrows(x1,y1,x2,y2,...)
}

draw_jcr_tree <- function(forest,n){

  left <- forest$ldaughter[,n] + 1
  right <- forest$rdaughter[,n] + 1
  vars <- forest$node_var[,n] + 1
  thrs <- forest$node_xvar[,n]
  tst <- paste(vars,round(thrs,2),sep = " ? ")
  df <- data.frame(test=tst,left=left,right=right,parent=NA,x=0,y=0)
  #df2 <- df[1:(nrow(df)+1-min(which(rev(df$test) != "0 ? 0"))),]
  df2 <- df[1:max(left),]
  
  for(i in 1:nrow(df2)){
    if(df2$left[i] != 1){
    df2$parent[df2$left[i]] <- i
    df2$parent[df2$right[i]] <- i
    df2$y[df2$left[i]] <- df2$y[i] + 1
    df2$y[df2$right[i]] <- df2$y[i] + 1
    }
  }
  Nlayer <- max(df2$y)
  width <- 2 ^ Nlayer
  
  plot(0,0,type="n",axes = FALSE,xlab="",ylab="",xlim = c(-width/2,width/2),ylim=c(-Nlayer,0))
  for(i in 1:nrow(df2)){

    if(df2$left[i] != 1){    
    offs <- 2^(Nlayer-df2$y[i]-2)
    
    # change + offs and - offs to make the left child also appear on the left side.
    df2$x[df2$left[i]] <- df2$x[i] + offs
    df2$x[df2$right[i]] <- df2$x[i] - offs

    arrows_label(df2$x[i],-df2$y[i],df2$x[df2$left[i]],-df2$y[df2$left[i]],lbl=">",angle=0,col="grey")
    
    arrows_label(df2$x[i],-df2$y[i],df2$x[df2$right[i]],-df2$y[df2$right[i]],lbl="<",angle=0,col="grey")
    
    text(df2$x[i],-df2$y[i],df2$test[i])

    } else if(df2$y[i] != 0) {
      text(df2$x[i],-df2$y[i],paste(forest$yf_pred[i,n],round(forest$yc_mu_pred[i,n],1),sep="\n"),col="green",xpd=TRUE)
    }
  }

}
