#checking for residual spatial autocorrelation of soil C data.
#clearn environemnt, load packages
rm(list=ls())
library (sp)
library(gstat)
library(boot)
library(data.table)
source('required_products_utilities/master_path_list.r')

#load Product_1 for lat/lon
d <- readRDS(Product_1.path)

#load all soils data.
soil <- data.table(read.csv(soil_data.processed.path))

#drop in lat/lon
soil <- merge(soil, d[,.(PLT_CN,latitude,longitude)], by = 'PLT_CN')
soil <- soil[,.(C.storage,latitude,longitude)]

#set coordinates
coordinates(soil) <- ~ longitude + latitude
bubble(soil, zcol = 'C.storage', fill=T, do.sqrt = F, maxsize=3)

#fit variogram
vario.1 <- variogram(log(C.storage)~ 1, data = soil) #get semivariance values
vario.model.1 <- fit.variogram(vario.1, vgm(model='Sph', nugget=T)) #fit a variogram model with a Spherical error structure and a "nugget" (nugget = intercept)

#plot your semivariance and your model.
png(filename=soil_semivariogram_path,width=5,height=5,units='in',res=300)
plot(vario.1, model=vario.model.1, xlab = 'distance (km)')
dev.off()
