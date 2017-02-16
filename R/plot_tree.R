draw_jcr_tree <- function(forest,n){
  left <- forest$ldaughter[,n] + 1
  right <- forest$rdaughter[,n] + 1
  vars <- forest$node_var[,n] + 1
  thrs <- forest$node_xvar[,n]
  tst <- paste(vars,round(thrs,2),sep = " ? ")
  print(tst)
  df <- data.frame(test=tst,left=left,right=right,parent=NA,x=0,y=0)
  df2 <- df[1:(nrow(df)+1-min(which(rev(df$test) != "0 ? 0.0"))),]
  
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
    df2$x[df2$left[i]] <- df2$x[i] - offs
    df2$x[df2$right[i]] <- df2$x[i] + offs

    arrows(df2$x[i],-df2$y[i],df2$x[df2$left[i]],-df2$y[df2$left[i]],angle=0,col="grey")
    arrows(df2$x[i],-df2$y[i],df2$x[df2$right[i]],-df2$y[df2$right[i]],angle=0,col="grey")
    text(df2$x[i],-df2$y[i],df2$test[i])

    } else if(df2$y[i] != 0) {
      text(df2$x[i],-df2$y[i],forest$yf_pred[i],col="green")
    }
  }
  # Determining the coordinates
  
  ##

  #text(0.5,0.5,"oooo")
  return(df2)
}

plot(0,0,xlim=c(-16,16),ylim=c(-5,0),type="n")
for(i in 0:5){
  points((1:2^i)*32/(2^(i)+1)-16,-rep(i,length(1:2^i)))
}    
