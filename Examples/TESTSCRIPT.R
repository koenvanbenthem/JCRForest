#### REAL DATA
# y <- read.csv("/Users/koen/Dropbox/ML-Scripts/SET3/R249_499.csv")[,-1]
# y$sign <- as.factor(y$sign)
# x <- read.csv("/Users/koen/Dropbox/ML-Scripts/SET3/S150_249.csv")[,-1]
# gt <- sample(100,20,replace=TRUE)
# Ntrain <- 70
# datap <- sample(1:nrow(x))
# training <- datap[1:Ntrain]
# testing <- datap[(Ntrain+1):length(datap)]
# traindat <- x[training,]
# trainresp <- y[training,]
library(JCRForest)
rm(list=ls())
ynobla <- read.csv("/Users/koen/Dropbox/ML-Scripts/SET3/R249_499.csv")[,-1]
ynobla$sign <- factor(ynobla$sign+1,levels=1:2)
#ynobla$sign <- sample(ynobla$sign)
#ynobla$rate <- 10*ynobla$rate
xnobla <- read.csv("/Users/koen/Dropbox/ML-Scripts/SET3/S150_249.csv")[,-1]
# save(ynobla,file="y-nocrash")
# save(xnobla,file="x-nocrash")
# ncol(xnobla)
output <- jcr_forest(xnobla,ynobla,20,Ntree=10,minsize = 3)
# output2 <- jcr_forest(xnobla,ynobla,5,Ntree=500)
# lapply(output,FUN = function(x){
#   if(length(dim(x)==2)){
#     all(apply(x,1,var)==0)
#   }
# })

plot(ynobla$rate,predict(output,xnobla)$yc)
abline(a=0,b=1)
table(ynobla$sign,apply(predict(output,xnobla)$yf,2,which.max))




#output$ldaughter

# output$ldaughter[24:25,]
# output$rdaughter[24:25,]
# 
pdf("tst2.pdf")
xnoblas <- xnobla
xnoblas$silly <- ynobla$rate
output_silly <- jcr_forest(xnoblas,ynobla,23,Ntree=20,minsize=5)
plot(ynobla$rate,predict(output_silly,xnoblas)$yc,xlab="population growth rate",ylab="predicted population growth rate",main="predictions on training data")
abline(a=0,b=1)
#draw_jcr_tree(output_silly,9)
table(ynobla$sign,apply(predict(output_silly,xnoblas)$yf,2,which.max))
table(ynobla$rate <= 1,ynobla$sign)

miss <- ynobla[(ynobla$sign!=apply(predict(output_silly,xnoblas)$yf,2,which.max)),]
points(miss$rate,miss$rate,col="red")
dev.off()
rm(list=ls())
x <- data.frame(x1=as.numeric(rep(1:10,20)),x2=runif(200))
y <- data.frame(y1=rnorm(nrow(x),x$x1,0.01),y2=as.factor(round(10*x$x2)))

lm(y$y1~x)
sil <- jcr_forest(x,y,2,100)

plot(y$y1,predict(sil,x)$yc)
table(y$y2,apply(predict(sil,x)$yf,2,which.max))
abline(b=1,a=0)
mean(y$y1[x$x1==2])
# output2 <- jcr_forest(xnobla,ynobla,5,Ntree=50)
# identical(output,output2)
# for(i in 1:length(output)){
#   cat(i,"_",identical(output[[i]],output2[[i]]),"\n")
# }
# # 
# # 
# # replicate(20,output <- jcr_forest(x,y,4,Ntree=500))
# # output <- jcr_forest(x,y,4,Ntree=500)
# # 
# # summary(output$yf_pred)
# # hist(output$yf_pred)
# # apply(predict(output,traindat)$yf,2,which.max)
# # round(sqrt(ncol(x)))
#
# tree data
dat <- data_gen_tree(2000)
# gt <- sample(100,20,replace=TRUE)
output <- jcr_forest(dat$x,dat$y,4,10)
# output2 <- jcr_forest(dat$x,dat$y,4,500)
# 
# identical(output1,output2)
# lapply(output1,FUN = function(x){
#   if(length(dim(x)==2)){
#     all(apply(x,1,var)==0)
#   }
# })
# output$ldaughter
# 

table(apply(predict(output,output$x)$yf,2,which.max),dat$y$y2)
plot(predict(output,output$x)$yc,dat$y$y1)
# predict(output,matrix(c(0.1,0.9,0.27,0.45),nrow=1))$yf
# draw_jcr_tree(output,4)

# t(tmp)[!apply(predict(output,output$x)$yf,2,which.max)==as.numeric(dat$y$y2),]
# output$x[!apply(predict(output,output$x)$yf,2,which.max)==as.numeric(dat$y$y2),]
# output$node_xvar
# #draw_jcr_tree(output,)
# bla <- rbind(1:nrow(output$yf_pred),output$yf_pred[,2],output$ldaughter[,2]+1,output$rdaughter[,2]+1,1:nrow(output$yf_pred))
# t(bla)
# 
# plot(dat$y$y1,predict(output,dat$x)$yc)
# output$ldaughter
# dat2 <- dat$y[1:50,]
# aggregate(dat2$y1,by=list(dat2$y2),mean)
# aggregate(dat2$y1,by=list(dat2$y2),sd)
# 
# 
# 
# 
# 
# 
# # artif cut-off: 0.5
# H_c <- function(y){
#   pb <- sum(y=='b')/length(y)
#   pa <- 1-pb
#   return(pb*log(pb)+pa*log(pa))
# }
# # all
# H_c(dat$y$y2)
# 
# # x1 > 0.5
# for(thr in seq(0.1,0.9,0.05)){
# 
# val <- (sum(dat$x[,1]>thr)/nrow(dat$x))*H_c(dat$y$y2[dat$x[,1]>thr]) + (sum(dat$x[,1]<=thr)/nrow(dat$x))*H_c(dat$y$y2[dat$x[,1]<=thr])
# cat(thr,"\t",val,"\n")
# }
# 
# goat <- structure(list(x=7),class="jcr_forest")
# predict(goat)
# 
# 
# blabla <- c(5,6,1,4)
# order(blabla)
# blabla[order(blabla)]
# 
# blabla[c(1,order(blabla[-1])+1)]

# pb <- sum(dat$y$y2 == 'b')/nrow(dat$y)
# pa <- 1-pb
# pa *log(pa) + pb*log(pb)

#dat$y$y2 <- factor(1:20)
# apply(output$x_bag,1,FUN = function(z) which(apply(dat$x,1 ,FUN = function(p) all(z==p))))
# 
# inds <- sapply(output$yc_bag,FUN = function(z) which(abs(dat$y$y1-z)<1e-6))
# table(output$yf_bag,dat$y$y2[inds])
# output$x_bag - dat$x[inds,]
