dat <- data_gen(20)
dat$x[,1] <- sample(1:20)
dat$x[,2] <- sample(21:40)
dat$x[,3] <- sample(41:60)


output <- jcr_forest(dat$x,dat$y,3,10)
output
#dat$y$y2 <- factor(1:20)
# apply(output$x_bag,1,FUN = function(z) which(apply(dat$x,1 ,FUN = function(p) all(z==p))))
# 
# inds <- sapply(output$yc_bag,FUN = function(z) which(abs(dat$y$y1-z)<1e-6))
# table(output$yf_bag,dat$y$y2[inds])
# output$x_bag - dat$x[inds,]
