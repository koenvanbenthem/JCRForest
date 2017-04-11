# JCRForest
Here, I will try to implement a joint classification/regression forest in R. 

A very preliminary version now exists. It is entirely based on the randomForest package.

This packages can be installed using:
```
devtools::install_github("koenvanbenthem/JCRForest") 
```
or (equivalently)
```
devtools::install_git("git://github.com/koenvanbenthem/JCRForest.git")
```

An example use is:
```
library(JCRForest)
dat <- data_gen_tree(200)
output <- jcr_forest(dat$x,dat$y,2,10)
draw_jcr_tree(output,2)
```
