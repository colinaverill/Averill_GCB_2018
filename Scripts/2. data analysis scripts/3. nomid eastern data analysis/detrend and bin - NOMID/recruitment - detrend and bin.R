#plotting recruitment means.
#detrending binned data to deal with probabilistic component.
rm(list=ls())
library(data.table)
library(runjags)
library(wesanderson) #wes anderson of course. 
library(boot)        #for inv.logit command

     am.fit.path <- '/fs/data3/caverill/FIA7_Ndep_data.processed/nomid_east.data_analysis_output/FIA7_NOMID_EAST_052617_AM.adult_linear_recruit.summary.Rdata'
     em.fit.path <- '/fs/data3/caverill/FIA7_Ndep_data.processed/nomid_east.data_analysis_output/FIA7_NOMID_EAST_052617_EM.adult_linear_recruit.summary.Rdata'
      range.path <- '/fs/data3/caverill/FIA7_Ndep_data.processed/nomid_east.data_analysis_output/FIA7_NOMID_EAST_052617_EM.adult_linear_recruit.predictor_ranges.rds'
am.filtered.path <- '/fs/data3/caverill/FIA7_Ndep_data.processed/nomid_east.data_analysis_output/FIA7_NOMID_EAST_052617_AM.adult_linear_recruit.filtered.rds'
em.filtered.path <- '/fs/data3/caverill/FIA7_Ndep_data.processed/nomid_east.data_analysis_output/FIA7_NOMID_EAST_052617_EM.adult_linear_recruit.filtered.rds'
       bins.path <- '/fs/data3/caverill/FIA7_Ndep_data.processed/nomid_east.data_analysis_output/nomid_recruit_bins.rds'

#load growth fit.
r.out.am    <- readRDS(am.fit.path)
r.out.em    <- readRDS(em.fit.path)
r.ndep.am <- r.out.am[grep('ndep.pred', rownames(r.out.am)),]
r.ndep.em <- r.out.em[grep('ndep.pred', rownames(r.out.em)),]

#we need the Ndep ranges of the original JAGS predictors predictor_ranges
r.range <- readRDS(range.path)
r.range <- r.range$x5.r

#load recruitment multiple regression parameter estimates
r.fit.am <- r.out.am[1:17,]
r.fit.em <- r.out.em[1:17,]

#load recruitment data.
d.am <- data.table(readRDS(am.filtered.path))
d.em <- data.table(readRDS(em.filtered.path))
d.am <- d.am[d.am$tot.15 < 20,]
d.am <- data.table(d.am)
d.em <- d.em[d.em$tot.15 < 20,]
d.em <- data.table(d.em)

#get binned Ndep categories
n.bins <- 6
d.am$cat <- cut(d.am$tot.15,n.bins, labels = F)
d.em$cat <- cut(d.em$tot.15,n.bins, labels = F)
d.am$cat<- as.factor(d.am$cat)
d.em$cat<- as.factor(d.em$cat)

#start binning the raw data.
bins                  <- aggregate(tot.15           ~ cat, data = d.am, FUN = 'mean'  )
bins$recruit.am       <- aggregate(recruit.adult.am ~ cat, data = d.am, FUN = 'mean'  )[,2]
bins$recruit.am.sd    <- aggregate(recruit.adult.am ~ cat, data = d.am, FUN = 'sd'    )[,2]
bins$recruit.am.n     <- aggregate(recruit.adult.am ~ cat, data = d.am, FUN = 'length')[,2]
bins$recruit.em       <- aggregate(recruit.adult.em ~ cat, data = d.em, FUN = 'mean'  )[,2]
bins$recruit.em.sd    <- aggregate(recruit.adult.em ~ cat, data = d.em, FUN = 'sd'    )[,2]
bins$recruit.em.n     <- aggregate(recruit.adult.em ~ cat, data = d.em, FUN = 'length')[,2]
bins$REMPER           <- aggregate(REMPER           ~ cat, data = d.em, FUN = 'mean'  )[,2]
#account for each bins variable rememeasurement period.
bins$recruit.am    <- bins$recruit.am    / bins$REMPER
bins$recruit.em    <- bins$recruit.em    / bins$REMPER
bins$recruit.am.sd <- bins$recruit.am.sd / bins$REMPER
bins$recruit.em.sd <- bins$recruit.em.sd / bins$REMPER

#get binned predictors
bins$mat       <- aggregate(mat     ~ cat, data = d.am, FUN = 'mean')[,2]
bins$map       <- aggregate(map     ~ cat, data = d.am, FUN = 'mean')[,2]
bins$cn        <- aggregate(cn      ~ cat, data = d.am, FUN = 'mean')[,2]
bins$ph        <- aggregate(pH_H2O  ~ cat, data = d.am, FUN = 'mean')[,2]
bins$clay      <- aggregate(clay    ~ cat, data = d.am, FUN = 'mean')[,2]
bins$rel.am    <- aggregate(relAM   ~ cat, data = d.am, FUN = 'mean')[,2]
bins$rel.em    <- aggregate(relEM   ~ cat, data = d.am, FUN = 'mean')[,2]
bins$basal     <- aggregate(log(BASAL)~ cat, data = d.am, FUN = 'mean')[,2]

#detrend the mean values - do the occurence and the poisson components separately.
#important to include intercept and N deposition means in every term because of the transformations. May be relative to other corrections as well.
pois.am   <- (r.fit.am[1,4] + r.fit.am[ 2,4]*bins$ph           + r.fit.am[ 3,4]*bins$cn       + r.fit.am[ 4,4]*bins$mat       + r.fit.am[ 5,4]*bins$map       + r.fit.am[ 6,4]*mean(d.am$tot.15) + r.fit.am[ 7,4]*bins$rel.am      + r.fit.am[ 8,4]*bins$clay      ) #+ r.fit.am[ 8,4]*bins$basal)
bern.am   <- (r.fit.am[9,4] + r.fit.am[10,4]*bins$ph           + r.fit.am[11,4]*bins$cn       + r.fit.am[12,4]*bins$mat       + r.fit.am[13,4]*bins$map       + r.fit.am[14,4]*mean(d.am$tot.15) + r.fit.am[15,4]*bins$rel.am      + r.fit.am[16,4]*bins$clay      ) #+ r.fit.am[16,4]*bins$basal)
pois.am.m <- (r.fit.am[1,4] + r.fit.am[ 2,4]*mean(d.am$pH_H2O) + r.fit.am[ 3,4]*mean(d.am$cn) + r.fit.am[ 4,4]*mean(d.am$mat) + r.fit.am[ 5,4]*mean(d.am$map) + r.fit.am[ 6,4]*mean(d.am$tot.15) + r.fit.am[ 7,4]*mean(d.am$relAM) + r.fit.am[ 8,4]*mean(d.am$clay)) #+ r.fit.am[ 8,4]*mean(log(d.am$BASAL)))
bern.am.m <- (r.fit.am[9,4] + r.fit.am[10,4]*mean(d.am$pH_H2O) + r.fit.am[11,4]*mean(d.am$cn) + r.fit.am[12,4]*mean(d.am$mat) + r.fit.am[13,4]*mean(d.am$map) + r.fit.am[14,4]*mean(d.am$tot.15) + r.fit.am[15,4]*mean(d.am$relAM) + r.fit.am[16,4]*mean(d.am$clay)) #+ r.fit.am[16,4]*mean(log(d.am$BASAL)))
pois.em   <- (r.fit.em[1,4] + r.fit.em[ 2,4]*bins$ph           + r.fit.em[ 3,4]*bins$cn       + r.fit.em[ 4,4]*bins$mat       + r.fit.em[ 5,4]*bins$map       + r.fit.em[ 6,4]*mean(d.em$tot.15) + r.fit.em[ 7,4]*bins$rel.em      + r.fit.em[ 8,4]*bins$clay      ) #+ r.fit.em[ 8,4]*bins$basal)
bern.em   <- (r.fit.em[9,4] + r.fit.em[10,4]*bins$ph           + r.fit.em[11,4]*bins$cn       + r.fit.em[12,4]*bins$mat       + r.fit.em[13,4]*bins$map       + r.fit.em[14,4]*mean(d.em$tot.15) + r.fit.em[15,4]*bins$rel.em      + r.fit.em[16,4]*bins$clay      ) #+ r.fit.em[16,4]*bins$basal)
pois.em.m <- (r.fit.em[1,4] + r.fit.em[ 2,4]*mean(d.em$pH_H2O) + r.fit.em[ 3,4]*mean(d.em$cn) + r.fit.em[ 4,4]*mean(d.em$mat) + r.fit.em[ 5,4]*mean(d.em$map) + r.fit.em[ 6,4]*mean(d.em$tot.15) + r.fit.em[ 7,4]*mean(d.em$relEM) + r.fit.em[ 8,4]*mean(d.em$clay)) #+ r.fit.am[ 8,4]*mean(log(d.am$BASAL)))
bern.em.m <- (r.fit.em[9,4] + r.fit.em[10,4]*mean(d.em$pH_H2O) + r.fit.em[11,4]*mean(d.em$cn) + r.fit.em[12,4]*mean(d.em$mat) + r.fit.em[13,4]*mean(d.em$map) + r.fit.em[14,4]*mean(d.em$tot.15) + r.fit.em[15,4]*mean(d.em$relEM) + r.fit.em[16,4]*mean(d.em$clay)) #+ r.fit.am[16,4]*mean(log(d.am$BASAL)))

#important to include intercept and N deposition means in every term because of the transformations. May be relative to other corrections as well.
pois.am   <- (r.fit.am[1,4] + r.fit.am[ 2,4]*bins$ph           + r.fit.am[ 3,4]*bins$cn       + r.fit.am[ 4,4]*bins$mat       + r.fit.am[ 5,4]*bins$map       + r.fit.am[ 6,4]*(bins$tot.15) + r.fit.am[ 7,4]*bins$rel.am      + r.fit.am[ 8,4]*bins$clay      ) #+ r.fit.am[ 8,4]*bins$basal)
bern.am   <- (r.fit.am[9,4] + r.fit.am[10,4]*bins$ph           + r.fit.am[11,4]*bins$cn       + r.fit.am[12,4]*bins$mat       + r.fit.am[13,4]*bins$map       + r.fit.am[14,4]*(bins$tot.15) + r.fit.am[15,4]*bins$rel.am      + r.fit.am[16,4]*bins$clay      ) #+ r.fit.am[16,4]*bins$basal)
pois.am.m <- (r.fit.am[1,4] + r.fit.am[ 2,4]*mean(d.am$pH_H2O) + r.fit.am[ 3,4]*mean(d.am$cn) + r.fit.am[ 4,4]*mean(d.am$mat) + r.fit.am[ 5,4]*mean(d.am$map) + r.fit.am[ 6,4]*(bins$tot.15) + r.fit.am[ 7,4]*mean(d.am$relAM) + r.fit.am[ 8,4]*mean(d.am$clay)) #+ r.fit.am[ 8,4]*mean(log(d.am$BASAL)))
bern.am.m <- (r.fit.am[9,4] + r.fit.am[10,4]*mean(d.am$pH_H2O) + r.fit.am[11,4]*mean(d.am$cn) + r.fit.am[12,4]*mean(d.am$mat) + r.fit.am[13,4]*mean(d.am$map) + r.fit.am[14,4]*(bins$tot.15) + r.fit.am[15,4]*mean(d.am$relAM) + r.fit.am[16,4]*mean(d.am$clay)) #+ r.fit.am[16,4]*mean(log(d.am$BASAL)))
pois.em   <- (r.fit.em[1,4] + r.fit.em[ 2,4]*bins$ph           + r.fit.em[ 3,4]*bins$cn       + r.fit.em[ 4,4]*bins$mat       + r.fit.em[ 5,4]*bins$map       + r.fit.em[ 6,4]*(bins$tot.15) + r.fit.em[ 7,4]*bins$rel.em      + r.fit.em[ 8,4]*bins$clay      ) #+ r.fit.em[ 8,4]*bins$basal)
bern.em   <- (r.fit.em[9,4] + r.fit.em[10,4]*bins$ph           + r.fit.em[11,4]*bins$cn       + r.fit.em[12,4]*bins$mat       + r.fit.em[13,4]*bins$map       + r.fit.em[14,4]*(bins$tot.15) + r.fit.em[15,4]*bins$rel.em      + r.fit.em[16,4]*bins$clay      ) #+ r.fit.em[16,4]*bins$basal)
pois.em.m <- (r.fit.em[1,4] + r.fit.em[ 2,4]*mean(d.em$pH_H2O) + r.fit.em[ 3,4]*mean(d.em$cn) + r.fit.em[ 4,4]*mean(d.em$mat) + r.fit.em[ 5,4]*mean(d.em$map) + r.fit.em[ 6,4]*(bins$tot.15) + r.fit.em[ 7,4]*mean(d.em$relEM) + r.fit.em[ 8,4]*mean(d.em$clay)) #+ r.fit.am[ 8,4]*mean(log(d.am$BASAL)))
bern.em.m <- (r.fit.em[9,4] + r.fit.em[10,4]*mean(d.em$pH_H2O) + r.fit.em[11,4]*mean(d.em$cn) + r.fit.em[12,4]*mean(d.em$mat) + r.fit.em[13,4]*mean(d.em$map) + r.fit.em[14,4]*(bins$tot.15) + r.fit.em[15,4]*mean(d.em$relEM) + r.fit.em[16,4]*mean(d.em$clay)) #+ r.fit.am[16,4]*mean(log(d.am$BASAL)))

#get corrected values
bins$am.detrend <- bins$recruit.am - exp(pois.am)*inv.logit(bern.am) + exp(pois.am.m)*inv.logit(bern.am.m)
bins$em.detrend <- bins$recruit.em - exp(pois.em)*inv.logit(bern.em) + exp(pois.em.m)*inv.logit(bern.em.m)
bins$am.detrend.upr <- (bins$recruit.am + bins$recruit.am.sd/sqrt(bins$recruit.am.n)) - exp(pois.am)*inv.logit(bern.am) + exp(pois.am.m)*inv.logit(bern.am.m)
bins$am.detrend.lwr <- (bins$recruit.am - bins$recruit.am.sd/sqrt(bins$recruit.am.n)) - exp(pois.am)*inv.logit(bern.am) + exp(pois.am.m)*inv.logit(bern.am.m)
bins$em.detrend.upr <- (bins$recruit.em + bins$recruit.em.sd/sqrt(bins$recruit.em.n)) - exp(pois.em)*inv.logit(bern.em) + exp(pois.em.m)*inv.logit(bern.em.m)
bins$em.detrend.lwr <- (bins$recruit.em - bins$recruit.em.sd/sqrt(bins$recruit.em.n)) - exp(pois.em)*inv.logit(bern.em) + exp(pois.em.m)*inv.logit(bern.em.m)


#save recruitment bins data
saveRDS(bins, bins.path)


#make some plots. 
#some graphic parameters
cols <- wes_palette("Zissou", 5)
trans <- 0.3
par(mfrow=c(1,2))
limy <- c(0,1.1)

#make AM vs EM plot
#AM Plot
plot(bins$am.detrend ~ bins$tot.15, pch = 16, 
     ylim=limy, xlim=c(min(d.am$tot.15),max(d.am$tot.15)))
#error bars
arrows(bins$tot.15, bins$am.detrend.lwr, bins$tot.15, bins$am.detrend.upr, length=0.05, angle=90, code=3)
#Drop fit on top.
lines(smooth.spline(r.ndep.am[,4] ~ r.range), lwd=2, col=cols[5])
polygon(c(r.range, rev(r.range)),c(r.ndep.am[,3], rev(r.ndep.am[,1])), col=adjustcolor(cols[5], trans), lty=0)
lines(smooth.spline(r.ndep.am[,3] ~ r.range), lty=2, lwd=0.5)
lines(smooth.spline(r.ndep.am[,1] ~ r.range), lty=2, lwd=0.5)

#EM Plot
plot(bins$em.detrend ~ bins$tot.15, pch = 16, 
     ylim=limy, xlim=c(min(d.em$tot.15),max(d.em$tot.15)))
#error bars
arrows(bins$tot.15, bins$em.detrend.lwr, bins$tot.15, bins$em.detrend.upr, length=0.05, angle=90, code=3)
#Drop fit on top.
lines(smooth.spline(r.ndep.em[,4] ~ r.range) , lwd=2, col=cols[2])
polygon(c(r.range, rev(r.range)),c(r.ndep.em[,3], rev(r.ndep.em[,1])), col=adjustcolor(cols[2], trans), lty=0)
lines(smooth.spline(r.ndep.em[,3] ~ r.range), lty=2, lwd=0.5)
lines(smooth.spline(r.ndep.em[,1] ~ r.range), lty=2, lwd=0.5)

#mean AM fit R2 = 0.65
am.spline <- smooth.spline(r.ndep.am[,4] ~ r.range)
summary(lm(bins$am.detrend ~ predict(am.spline,bins$tot.15)$y))

#mean EM fit R2 = 0.46
em.spline <- smooth.spline(r.ndep.em[,4] ~ r.range)
summary(lm(bins$em.detrend ~ predict(em.spline,bins$tot.15)$y))