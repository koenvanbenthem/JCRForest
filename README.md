# JCRForest
Here, I will try to implement a joint classification/regression forest in R. 

A very preliminary version now exists. It is entirely based on the randomForest package.

This packages can be installed using:
`devtools::install_github("koenvanbenthem/JCRForest") `

An example use is:
```
library(JCRForest)
dat <- data_gen_tree(200)
gt <- sample(100,20,replace=TRUE)
output <- jcr_forest(dat$x,dat$y,2,10,bla=gt)
draw_jcr_tree(output,2)
```
