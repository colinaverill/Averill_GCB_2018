#Histogram of sites in soil analysis as a function of %forested condition.
#clear environment, load packages.
rm(list=ls())
source('required_products_utilities/master_path_list.r')
library(data.table)

#load the soils data
d <- readRDS(Product_1.path)

#get subset used in soil analysis
dat <-   d[,.(PLT_CN, mat, map, pH_H2O, tot.15, relEM, clay, N.storage, C.storage)]
dat <- dat[complete.cases(dat),]
d <- d[PLT_CN %in% dat$PLT_CN,]
d$forest_proportion <- d$forest_proportion*100

#set save path, size and resolution.
png(filename=Supp_Figure7_all_path,width=4.5,height=4.5,units='in',res=300)

par(mar = c(4.5,4.5,.2,.2))

#make histogram
hist(d$forest_proportion,
     ylab = 'Number of Sites',
     xlab = '% of plot in a forested condition',
     main = NA)

#end plot
dev.off()