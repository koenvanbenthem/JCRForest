dat <- data_gen(20,Nx=3)
dat$x[,1] <- sample(1:20)
dat$x[,2] <- sample(21:40)
dat$x[,3] <- sample(41:60)
#dat$x[,4] <- sample(61:80)

gt <- sample(100,20,replace=TRUE)
output <- jcr_forest(dat$x,dat$y,3,10,bla=gt)

output$dum_ind
output$dum_vect
gt
gt[output$dum_ind/2]

output

gt <- sample(100,20,replace=TRUE)
dat <- data_gen_inf(200,Nx=3)
tmp <- dat$x[,1]
dat$x[,1] <- dat$x[,2]
dat$x[,2] <- tmp
output <- jcr_forest(dat$x,dat$y,3,10,bla=gt)

pred_output <- predict(output,output$x)
table(apply(pred_output$yf,2,which.max),output$yf)
rbind(pred_output$yf,output$yf)

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
