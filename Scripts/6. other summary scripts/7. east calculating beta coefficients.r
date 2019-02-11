#calculating standardized (beta) coefficients
#beta coefficient defined as the abosolute value of change in dependent variable 
#in response to an increase in 1SD of a predictor
#relative to the value of the dependent variable at the median value of the predictor
#holding all other predictors constant at their means.
#clear R environment, load packages
rm(list=ls())
library(data.table)
library(boot)

#create an output list, define save path.
output.list <- list()
output.path <- '/fs/data3/caverill/FIA7_Ndep_data.processed/east.data_analysis_output/east_beta_list.rds'

####SOIL C BETA FACTORS####
soil.data  <- readRDS('/fs/data3/caverill/FIA7_Ndep_data.processed/east.data_analysis_output/FIA7_EAST_053017_soilC.filtered.rds')
soil.pred  <- readRDS('/fs/data3/caverill/FIA7_Ndep_data.processed/east.data_analysis_output/FIA7_EAST_053017_soilC.summary.Rdata')
pred.range <- readRDS('/fs/data3/caverill/FIA7_Ndep_data.processed/east.data_analysis_output/FIA7_EAST_053017_soilC.predictor_ranges.rds')

#subset output by predictor
mat       <- soil.pred[grep('y.x1'    , rownames(soil.pred)),];      mat <- exp(mat )
map       <- soil.pred[grep('y.x2'    , rownames(soil.pred)),];      map <- exp(map )
pH        <- soil.pred[grep('y.x3'    , rownames(soil.pred)),];       pH <- exp(pH  )
ndep      <- soil.pred[grep('y.x4'    , rownames(soil.pred)),];     ndep <- exp(ndep)
em        <- soil.pred[grep('y.x5'    , rownames(soil.pred)),];       em <- exp(em  )
clay      <- soil.pred[grep('y.x7'    , rownames(soil.pred)),];     clay <- exp(clay)

#grab standard deviation of each predictor from the data subset for a particular analysis.
sd.C    <- sd(soil.data$C.storage)
sd.mat  <- sd(soil.data$mat)
sd.map  <- sd(soil.data$map)
sd.pH   <- sd(soil.data$pH)
sd.ndep <- sd(soil.data$tot.15)
sd.em   <- sd(soil.data$relEM)
sd.clay <- sd(soil.data$clay)


#C value at center + 1 sd for each predictor, holding others constant
mat.change <-  mat[as.numeric(rownames(pred.range[abs((pred.range[51,1] + sd.mat ) - pred.range$x1.range) == min(abs((pred.range[51,1] + sd.mat ) - pred.range$x1.range)),])),4]
map.change <-  map[as.numeric(rownames(pred.range[abs((pred.range[51,2] + sd.map ) - pred.range$x2.range) == min(abs((pred.range[51,2] + sd.map ) - pred.range$x2.range)),])),4]
pH.change <-   pH[as.numeric(rownames(pred.range[abs((pred.range[51,3] + sd.pH  ) - pred.range$x3.range) == min(abs((pred.range[51,3] + sd.pH  ) - pred.range$x3.range)),])),4]
ndep.change <- ndep[as.numeric(rownames(pred.range[abs((pred.range[51,4] + sd.ndep) - pred.range$x4.range) == min(abs((pred.range[51,4] + sd.ndep) - pred.range$x4.range)),])),4]
em.change <-   em[as.numeric(rownames(pred.range[abs((pred.range[51,5] + sd.em  ) - pred.range$x5.range) == min(abs((pred.range[51,5] + sd.em  ) - pred.range$x5.range)),])),4]
clay.change <- clay[as.numeric(rownames(pred.range[abs((pred.range[51,6] + sd.clay) - pred.range$x7.range) == min(abs((pred.range[51,6] + sd.clay) - pred.range$x7.range)),])),4]

#C value at the center of each predictors range
mat.center <-  mat[51,4]
map.center <-  map[51,4]
pH.center <-   pH[51,4]
ndep.center <- ndep[51,4]
em.center <-   em[51,4]
clay.center <- clay[51,4]

#beta factors
mat.beta <- ( mat.change -  mat.center)#/sd.C
map.beta <- ( map.change -  map.center)#/sd.C
pH.beta <- (  pH.change -   pH.center)#/sd.C
ndep.beta <- (ndep.change - ndep.center)#/sd.C
em.beta <- (  em.change -   em.center)#/sd.C
clay.beta <- (clay.change - clay.center)#/sd.C

#create a data frame of beta factors, add to output list
soil.beta <- data.frame(mat.beta, map.beta, pH.beta, clay.beta, ndep.beta, em.beta)
output.list[[1]] <- soil.beta

####RELATIVE ABUNDANCE BETA FACTORS####
relEM.data  <- readRDS('/fs/data3/caverill/FIA7_Ndep_data.processed/east.data_analysis_output/FIA7_east_053017_rel.abundance.filtered.rds')
relEM.pred  <- readRDS('/fs/data3/caverill/FIA7_Ndep_data.processed/east.data_analysis_output/FIA7_east_053017_rel.abundance.summary.rds')
pred.range  <- readRDS('/fs/data3/caverill/FIA7_Ndep_data.processed/east.data_analysis_output/FIA7_east_053017_rel.abundance.ranges.rds')

#subset output by predictor
mat       <- relEM.pred[grep('x1'    , rownames(relEM.pred)),];      mat <- inv.logit(mat )
map       <- relEM.pred[grep('x2'    , rownames(relEM.pred)),];      map <- inv.logit(map )
cn        <- relEM.pred[grep('x3'    , rownames(relEM.pred)),];       cn <- inv.logit(cn  )
pH        <- relEM.pred[grep('x4'    , rownames(relEM.pred)),];       pH <- inv.logit(pH  )
ndep      <- relEM.pred[grep('x5'    , rownames(relEM.pred)),];     ndep <- inv.logit(ndep)
clay      <- relEM.pred[grep('x6'    , rownames(relEM.pred)),];     clay <- inv.logit(clay)

#get sd
sd.mat  <- sd(relEM.data$mat)
sd.map  <- sd(relEM.data$map)
sd.pH   <- sd(relEM.data$pH)
sd.ndep <- sd(relEM.data$tot.15)
sd.em   <- sd(relEM.data$relEM)
sd.clay <- sd(relEM.data$clay)
sd.cn   <- sd(relEM.data$cn)


#C value at center + 1 sd for each predictor, holding others constant
mat.change  <-  mat[as.numeric(rownames(pred.range[abs((pred.range[51,1] + sd.mat ) - pred.range$mat.range ) == min(abs((pred.range[51,1] + sd.mat ) - pred.range$mat.range )),])),4]
map.change  <-  map[as.numeric(rownames(pred.range[abs((pred.range[51,2] + sd.map ) - pred.range$map.range ) == min(abs((pred.range[51,2] + sd.map ) - pred.range$map.range )),])),4]
cn.change   <-   cn[as.numeric(rownames(pred.range[abs((pred.range[51,3] + sd.cn  ) - pred.range$cn.range  ) == min(abs((pred.range[51,3] + sd.cn  ) - pred.range$cn.range  )),])),4]
pH.change   <-   pH[as.numeric(rownames(pred.range[abs((pred.range[51,4] + sd.pH  ) - pred.range$ph.range  ) == min(abs((pred.range[51,4] + sd.pH  ) - pred.range$ph.range  )),])),4]
ndep.change <- ndep[as.numeric(rownames(pred.range[abs((pred.range[51,5] + sd.ndep) - pred.range$ndep.range) == min(abs((pred.range[51,5] + sd.ndep) - pred.range$ndep.range)),])),4]
clay.change <- clay[as.numeric(rownames(pred.range[abs((pred.range[51,6] + sd.clay) - pred.range$clay.range) == min(abs((pred.range[51,6] + sd.clay) - pred.range$clay.range)),])),4]

#grab center values
mat.center  <-  mat[51,4]
map.center  <-  map[51,4]
pH.center   <-   pH[51,4]
cn.center   <-   cn[51,4]
ndep.center <- ndep[51,4]
clay.center <- clay[51,4]

#beta factors
mat.beta <- ( mat.change -  mat.center)#/sd.C
map.beta <- ( map.change -  map.center)#/sd.C
cn.beta <- (  cn.change -   cn.center)#/sd.C
pH.beta <- (  pH.change -   pH.center)#/sd.C
clay.beta <- (clay.change - clay.center)#/sd.C
ndep.beta <- (ndep.change - ndep.center)#/sd.C

#save beta factors
rel.abundance.beta <- data.frame(mat.beta,map.beta, cn.beta, pH.beta, clay.beta, ndep.beta)
output.list[[2]] <- rel.abundance.beta

####GROWTH BETA FACTORS####
#grab growth data
growth.data  <- readRDS('/fs/data3/caverill/FIA7_Ndep_data.processed/east.data_analysis_output/FIA7_051517_growth.filtered.rds')
growth.pred  <- readRDS('/fs/data3/caverill/FIA7_Ndep_data.processed/east.data_analysis_output/FIA7_051517_growth.summary.Rdata')
pred.range   <- readRDS('/fs/data3/caverill/FIA7_Ndep_data.processed/east.data_analysis_output/FIA7_051517_growth.predictor_ranges.rds')

#subset output by predictor
mat       <- growth.pred[grep('x2'    , rownames(growth.pred)),]
map       <- growth.pred[grep('x3'    , rownames(growth.pred)),]
cn        <- growth.pred[grep('x4'    , rownames(growth.pred)),]
pH        <- growth.pred[grep('x5'    , rownames(growth.pred)),]
clay      <- growth.pred[grep('x8'    , rownames(growth.pred)),]
ndep      <- growth.pred[grep('x7'    , rownames(growth.pred)),]; ndep <- ndep[1:101,]
ndep.am   <- growth.pred[grep('x7.am' , rownames(growth.pred)),]
ndep.em   <- growth.pred[grep('x7.em' , rownames(growth.pred)),]

#get sd
sd.mat  <- sd(relEM.data$mat)
sd.map  <- sd(relEM.data$map)
sd.pH   <- sd(relEM.data$pH)
sd.ndep <- sd(relEM.data$tot.15)
sd.clay <- sd(relEM.data$clay)
sd.cn   <- sd(relEM.data$cn)

#C value at center + 1 sd for each predictor, holding others constant
mat.change     <-     mat[as.numeric(rownames(pred.range[abs((pred.range[51,2] + sd.mat ) - pred.range$x2.range) == min(abs((pred.range[51,2] + sd.mat ) - pred.range$x2.range)),])),4]
map.change     <-     map[as.numeric(rownames(pred.range[abs((pred.range[51,3] + sd.map ) - pred.range$x3.range) == min(abs((pred.range[51,3] + sd.map ) - pred.range$x3.range)),])),4]
cn.change      <-      cn[as.numeric(rownames(pred.range[abs((pred.range[51,4] + sd.cn  ) - pred.range$x4.range) == min(abs((pred.range[51,4] + sd.cn  ) - pred.range$x4.range)),])),4]
pH.change      <-      pH[as.numeric(rownames(pred.range[abs((pred.range[51,5] + sd.pH  ) - pred.range$x5.range) == min(abs((pred.range[51,5] + sd.pH  ) - pred.range$x5.range)),])),4]
clay.change    <-    clay[as.numeric(rownames(pred.range[abs((pred.range[51,8] + sd.clay) - pred.range$x8.range) == min(abs((pred.range[51,8] + sd.clay) - pred.range$x8.range)),])),4]
ndep.change    <-    ndep[as.numeric(rownames(pred.range[abs((pred.range[51,7] + sd.ndep) - pred.range$x7.range) == min(abs((pred.range[51,7] + sd.ndep) - pred.range$x7.range)),])),4]
ndep.am.change <- ndep.am[as.numeric(rownames(pred.range[abs((pred.range[51,7] + sd.ndep) - pred.range$x7.range) == min(abs((pred.range[51,7] + sd.ndep) - pred.range$x7.range)),])),4]
ndep.em.change <- ndep.em[as.numeric(rownames(pred.range[abs((pred.range[51,7] + sd.ndep) - pred.range$x7.range) == min(abs((pred.range[51,7] + sd.ndep) - pred.range$x7.range)),])),4]

#center values
mat.center     <-     mat[51,4]
map.center     <-     map[51,4]
pH.center      <-      pH[51,4]
cn.center      <-      cn[51,4]
clay.center    <-    clay[51,4]
ndep.center    <-    ndep[51,4]
ndep.am.center <- ndep.am[51,4]
ndep.em.center <- ndep.em[51,4]

#beta factors
mat.beta     <- (    mat.change -  mat.center)#/sd.C
map.beta     <- (    map.change -  map.center)#/sd.C
pH.beta      <- (     pH.change -   pH.center)#/sd.C
cn.beta      <- (     cn.change -   cn.center)#/sd.C
clay.beta    <- (   clay.change - clay.center)#/sd.C
ndep.beta    <- (   ndep.change - ndep.center)#/sd.C
ndep.am.beta <- (ndep.am.change - ndep.am.center)#/sd.C
ndep.em.beta <- (ndep.em.change - ndep.em.center)#/sd.C

#save output
growth.beta <- data.frame(mat.beta, map.beta, pH.beta, cn.beta, clay.beta, ndep.beta, ndep.am.beta, ndep.em.beta)
output.list[[3]] <- growth.beta

####AM RECRUITMENT####
recruitAM.data  <- readRDS('/fs/data3/caverill/FIA7_Ndep_data.processed/east.data_analysis_output/FIA7_EAST_052617_AM.adult_EAST_linear_recruit.filtered.rds')
recruitAM.pred  <- readRDS('/fs/data3/caverill/FIA7_Ndep_data.processed/east.data_analysis_output/FIA7_EAST_052617_AM.adult_EAST_linear_recruit.summary.Rdata')
pred.range      <- readRDS('/fs/data3/caverill/FIA7_Ndep_data.processed/east.data_analysis_output/FIA7_EAST_052617_AM.adult_EAST_linear_recruit.predictor_ranges.rds')

#subset output by predictor 
mat       <- recruitAM.pred[grep('mat'   , rownames(recruitAM.pred)),]
map       <- recruitAM.pred[grep('map'   , rownames(recruitAM.pred)),]
cn        <- recruitAM.pred[grep('cn'    , rownames(recruitAM.pred)),]
pH        <- recruitAM.pred[grep('ph'    , rownames(recruitAM.pred)),]
clay      <- recruitAM.pred[grep('clay'  , rownames(recruitAM.pred)),]
ndep      <- recruitAM.pred[grep('ndep'  , rownames(recruitAM.pred)),]

#get sd
sd.mat  <- sd(recruitAM.data$mat)
sd.map  <- sd(recruitAM.data$map)
sd.pH   <- sd(recruitAM.data$pH)
sd.ndep <- sd(recruitAM.data$tot.15)
sd.clay <- sd(recruitAM.data$clay)
sd.cn   <- sd(recruitAM.data$cn)

#Recruitment value at center + 1 sd for each predictor, holding others constant
#these senitivities seem pretty strange. 
mat.change     <-     mat[as.numeric(rownames(pred.range[abs((pred.range[51,3] + sd.mat ) - pred.range$x3.r) == min(abs((pred.range[51,3] + sd.mat ) - pred.range$x3.r)),])),4]
map.change     <-     map[as.numeric(rownames(pred.range[abs((pred.range[51,4] + sd.map ) - pred.range$x4.r) == min(abs((pred.range[51,4] + sd.map ) - pred.range$x4.r)),])),4]
cn.change      <-      cn[as.numeric(rownames(pred.range[abs((pred.range[51,2] + sd.cn  ) - pred.range$x2.r) == min(abs((pred.range[51,2] + sd.cn  ) - pred.range$x2.r)),])),4]
pH.change      <-      pH[as.numeric(rownames(pred.range[abs((pred.range[51,1] + sd.pH  ) - pred.range$x1.r) == min(abs((pred.range[51,1] + sd.pH  ) - pred.range$x1.r)),])),4]
clay.change    <-    clay[as.numeric(rownames(pred.range[abs((pred.range[51,7] + sd.clay) - pred.range$x7.r) == min(abs((pred.range[51,7] + sd.clay) - pred.range$x7.r)),])),4]
ndep.change    <-    ndep[as.numeric(rownames(pred.range[abs((pred.range[51,5] + sd.ndep) - pred.range$x5.r) == min(abs((pred.range[51,5] + sd.ndep) - pred.range$x5.r)),])),4]

#center values
mat.center     <-     mat[51,4]
map.center     <-     map[51,4]
pH.center      <-      pH[51,4]
cn.center      <-      cn[51,4]
clay.center    <-    clay[51,4]
ndep.center    <-    ndep[51,4]

#beta factors
mat.beta     <- (    mat.change -  mat.center)#/sd.C
map.beta     <- (    map.change -  map.center)#/sd.C
pH.beta      <- (     pH.change -   pH.center)#/sd.C
cn.beta      <- (     cn.change -   cn.center)#/sd.C
clay.beta    <- (   clay.change - clay.center)#/sd.C
ndep.beta    <- (   ndep.change - ndep.center)#/sd.C

#save output
am.recruit.beta <- data.frame(mat.beta, map.beta, cn.beta, pH.beta, clay.beta, ndep.beta)
output.list[[4]] <- am.recruit.beta

####EM RECRUITMENT####
recruitEM.data  <- readRDS('/fs/data3/caverill/FIA7_Ndep_data.processed/east.data_analysis_output/FIA7_EAST_052617_EM.adult_EAST_linear_recruit.filtered.rds')
recruitEM.pred  <- readRDS('/fs/data3/caverill/FIA7_Ndep_data.processed/east.data_analysis_output/FIA7_EAST_052617_EM.adult_EAST_linear_recruit.summary.Rdata')
pred.range      <- readRDS('/fs/data3/caverill/FIA7_Ndep_data.processed/east.data_analysis_output/FIA7_EAST_052617_EM.adult_EAST_linear_recruit.predictor_ranges.rds')

#subset output by predictor
mat       <- recruitEM.pred[grep('mat'   , rownames(recruitEM.pred)),]
map       <- recruitEM.pred[grep('map'   , rownames(recruitEM.pred)),]
cn        <- recruitEM.pred[grep('cn'    , rownames(recruitEM.pred)),]
pH        <- recruitEM.pred[grep('ph'    , rownames(recruitEM.pred)),]
clay      <- recruitEM.pred[grep('clay'  , rownames(recruitEM.pred)),]
ndep      <- recruitEM.pred[grep('ndep'  , rownames(recruitEM.pred)),]

#get sd
sd.mat  <- sd(recruitEM.data$mat)
sd.map  <- sd(recruitEM.data$map)
sd.pH   <- sd(recruitEM.data$pH)
sd.ndep <- sd(recruitEM.data$tot.15)
sd.clay <- sd(recruitEM.data$clay)
sd.cn   <- sd(recruitEM.data$cn)

#Recruitment value at center + 1 sd for each predictor, holding others constant
#these senitivities seem pretty strange. 
mat.change     <-     mat[as.numeric(rownames(pred.range[abs((pred.range[51,3] + sd.mat ) - pred.range$x3.r) == min(abs((pred.range[51,3] + sd.mat ) - pred.range$x3.r)),])),4]
map.change     <-     map[as.numeric(rownames(pred.range[abs((pred.range[51,4] + sd.map ) - pred.range$x4.r) == min(abs((pred.range[51,4] + sd.map ) - pred.range$x4.r)),])),4]
cn.change      <-      cn[as.numeric(rownames(pred.range[abs((pred.range[51,2] + sd.cn  ) - pred.range$x2.r) == min(abs((pred.range[51,2] + sd.cn  ) - pred.range$x2.r)),])),4]
pH.change      <-      pH[as.numeric(rownames(pred.range[abs((pred.range[51,1] + sd.pH  ) - pred.range$x1.r) == min(abs((pred.range[51,1] + sd.pH  ) - pred.range$x1.r)),])),4]
clay.change    <-    clay[as.numeric(rownames(pred.range[abs((pred.range[51,7] + sd.clay) - pred.range$x7.r) == min(abs((pred.range[51,7] + sd.clay) - pred.range$x7.r)),])),4]
ndep.change    <-    ndep[as.numeric(rownames(pred.range[abs((pred.range[51,5] + sd.ndep) - pred.range$x5.r) == min(abs((pred.range[51,5] + sd.ndep) - pred.range$x5.r)),])),4]

#center values
mat.center     <-     mat[51,4]
map.center     <-     map[51,4]
pH.center      <-      pH[51,4]
cn.center      <-      cn[51,4]
clay.center    <-    clay[51,4]
ndep.center    <-    ndep[51,4]

#beta factors
mat.beta     <- (    mat.change -  mat.center)#/sd.C
map.beta     <- (    map.change -  map.center)#/sd.C
pH.beta      <- (     pH.change -   pH.center)#/sd.C
cn.beta      <- (     cn.change -   cn.center)#/sd.C
clay.beta    <- (   clay.change - clay.center)#/sd.C
ndep.beta    <- (   ndep.change - ndep.center)#/sd.C

#save output
em.recruit.beta <- data.frame(mat.beta, map.beta, pH.beta, cn.beta, clay.beta, ndep.beta)
output.list[[5]] <- em.recruit.beta

####MORTALITY BETA FACTORS####
mortality.data  <- readRDS('/fs/data3/caverill/FIA7_Ndep_data.processed/east.data_analysis_output/outputFIA7_EAST_052617_runjags_EAST_mortality.filtered.rds')
mortality.pred  <- readRDS('/fs/data3/caverill/FIA7_Ndep_data.processed/east.data_analysis_output/outputFIA7_EAST_052617_runjags_EAST_mortality.summary.Rdata')
pred.range      <- readRDS('/fs/data3/caverill/FIA7_Ndep_data.processed/east.data_analysis_output/outputFIA7_EAST_052617_runjags_EAST_mortality.predictor_ranges.rds')

#subset output by predictor
mat       <- mortality.pred[grep('x2'    , rownames(mortality.pred)),]; mat    <- inv.logit(mat)
map       <- mortality.pred[grep('x3'    , rownames(mortality.pred)),]; map    <- inv.logit(map)
cn        <- mortality.pred[grep('x4'    , rownames(mortality.pred)),];  cn    <- inv.logit(cn)
pH        <- mortality.pred[grep('x5'    , rownames(mortality.pred)),];  pH    <- inv.logit(pH)
clay      <- mortality.pred[grep('x7'    , rownames(mortality.pred)),];clay    <- inv.logit(clay)
ndep      <- mortality.pred[grep('x6'    , rownames(mortality.pred)),];ndep    <- inv.logit(ndep)   ; ndep <- ndep[1:101,]
ndep.am   <- mortality.pred[grep('x6.am' , rownames(mortality.pred)),];ndep.am <- inv.logit(ndep.am)
ndep.em   <- mortality.pred[grep('x6.em' , rownames(mortality.pred)),];ndep.em <- inv.logit(ndep.em)

#get sd
sd.mat  <- sd(mortality.data$mat)
sd.map  <- sd(mortality.data$map)
sd.pH   <- sd(mortality.data$pH)
sd.ndep <- sd(mortality.data$tot.15)
sd.clay <- sd(mortality.data$clay)
sd.cn   <- sd(mortality.data$cn)

#Mortality probability at center + 1 sd for each predictor, holding others constant
mat.change     <-     mat[as.numeric(rownames(pred.range[abs((pred.range[51,2] + sd.mat ) - pred.range$range.x2) == min(abs((pred.range[51,2] + sd.mat ) - pred.range$range.x2)),])),4]
map.change     <-     map[as.numeric(rownames(pred.range[abs((pred.range[51,3] + sd.map ) - pred.range$range.x3) == min(abs((pred.range[51,3] + sd.map ) - pred.range$range.x3)),])),4]
cn.change      <-      cn[as.numeric(rownames(pred.range[abs((pred.range[51,4] + sd.cn  ) - pred.range$range.x4) == min(abs((pred.range[51,4] + sd.cn  ) - pred.range$range.x4)),])),4]
pH.change      <-      pH[as.numeric(rownames(pred.range[abs((pred.range[51,5] + sd.pH  ) - pred.range$range.x5) == min(abs((pred.range[51,5] + sd.pH  ) - pred.range$range.x5)),])),4]
clay.change    <-    clay[as.numeric(rownames(pred.range[abs((pred.range[51,7] + sd.clay) - pred.range$range.x7) == min(abs((pred.range[51,7] + sd.clay) - pred.range$range.x7)),])),4]
ndep.change    <-    ndep[as.numeric(rownames(pred.range[abs((pred.range[51,6] + sd.ndep) - pred.range$range.x6) == min(abs((pred.range[51,6] + sd.ndep) - pred.range$range.x6)),])),4]
ndep.am.change <- ndep.am[as.numeric(rownames(pred.range[abs((pred.range[51,6] + sd.ndep) - pred.range$range.x6) == min(abs((pred.range[51,6] + sd.ndep) - pred.range$range.x6)),])),4]
ndep.em.change <- ndep.em[as.numeric(rownames(pred.range[abs((pred.range[51,6] + sd.ndep) - pred.range$range.x6) == min(abs((pred.range[51,6] + sd.ndep) - pred.range$range.x6)),])),4]

#center values
mat.center     <-     mat[51,4]
map.center     <-     map[51,4]
pH.center      <-      pH[51,4]
cn.center      <-      cn[51,4]
clay.center    <-    clay[51,4]
ndep.center    <-    ndep[51,4]
ndep.em.center <- ndep.am[51,4]
ndep.am.center <- ndep.em[51,4]

#beta factors
mat.beta     <- (    mat.change -     mat.center)#/sd.C
map.beta     <- (    map.change -     map.center)#/sd.C
pH.beta      <- (     pH.change -      pH.center)#/sd.C
cn.beta      <- (     cn.change -      cn.center)#/sd.C
clay.beta    <- (   clay.change -    clay.center)#/sd.C
ndep.beta    <- (   ndep.change -    ndep.center)#/sd.C
ndep.am.beta <- (ndep.am.change - ndep.am.center)#/sd.C
ndep.em.beta <- (ndep.em.change - ndep.em.center)#/sd.C

#save output
mortality.beta <- data.frame(mat.beta, map.beta, cn.beta, pH.beta, clay.beta, ndep.beta, ndep.am.beta, ndep.em.beta)
output.list[[6]] <- mortality.beta

####NAME AND SAVE OUTPUT####
names(output.list) <- c('soil','relative.abundance','growth','recruit.am','recruit.em','mortality')
saveRDS(output.list, output.path)