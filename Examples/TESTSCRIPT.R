## tree data
dat <- data_gen_tree(2000)
gt <- sample(100,20,replace=TRUE)
output <- jcr_forest(dat$x,dat$y,4,10,bla=gt)


table(apply(predict(output,output$x)$yf,2,which.max),dat$y$y2)
plot(predict(output,output$x)$yc,dat$y$y1)

t(tmp)[!apply(predict(output,output$x)$yf,2,which.max)==as.numeric(dat$y$y2),]
output$x[!apply(predict(output,output$x)$yf,2,which.max)==as.numeric(dat$y$y2),]
output$node_xvar
draw_jcr_tree(output,4)
bla <- rbind(1:nrow(output$yf_pred),output$yf_pred[,2],output$ldaughter[,2]+1,output$rdaughter[,2]+1,1:nrow(output$yf_pred))
t(bla)

plot(dat$y$y1,predict(output,dat$x)$yc)
output$ldaughter
dat2 <- dat$y[1:50,]
aggregate(dat2$y1,by=list(dat2$y2),mean)
aggregate(dat2$y1,by=list(dat2$y2),sd)






# artif cut-off: 0.5
H_c <- function(y){
  pb <- sum(y=='b')/length(y)
  pa <- 1-pb
  return(pb*log(pb)+pa*log(pa))
}
# all
H_c(dat$y$y2)

# x1 > 0.5
for(thr in seq(0.1,0.9,0.05)){

val <- (sum(dat$x[,1]>thr)/nrow(dat$x))*H_c(dat$y$y2[dat$x[,1]>thr]) + (sum(dat$x[,1]<=thr)/nrow(dat$x))*H_c(dat$y$y2[dat$x[,1]<=thr])
cat(thr,"\t",val,"\n")
}

goat <- structure(list(x=7),class="jcr_forest")
predict(goat)


blabla <- c(5,6,1,4)
order(blabla)
blabla[order(blabla)]

blabla[c(1,order(blabla[-1])+1)]

# pb <- sum(dat$y$y2 == 'b')/nrow(dat$y)
# pa <- 1-pb
# pa *log(pa) + pb*log(pb)

#dat$y$y2 <- factor(1:20)
# apply(output$x_bag,1,FUN = function(z) which(apply(dat$x,1 ,FUN = function(p) all(z==p))))
# 
# inds <- sapply(output$yc_bag,FUN = function(z) which(abs(dat$y$y1-z)<1e-6))
# table(output$yf_bag,dat$y$y2[inds])
# output$x_bag - dat$x[inds,]
