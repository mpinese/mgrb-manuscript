# install.packages(c("corrplot","forestplot","ggplot2","gplots","knitr","knitrProgressBar","mgcv","openxlsx","plotROC","plyr","pROC","RColorBrewer","reshape2","rmarkdown","viridis","HDInterval","doParallel","viridis"))
setwd("C:/Users/mpinese/Downloads/MGRB_revisions/common_var_burden_expanded")
library(rmarkdown)
knitr::opts_knit$set(dev = "svg")
options(warn = 1, bitmapType = "cairo")
render("common_variants.Rmd")
