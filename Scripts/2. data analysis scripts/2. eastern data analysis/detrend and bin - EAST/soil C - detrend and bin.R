#plotting soil C by mycorrhizal type and N-deposition, detrending for other factors.
rm(list=ls())
library(data.table)
library(runjags)
library(wesanderson)
library(boot)        #for inv.logit command

#specify paths for the Ndep prediction, the data, and the model parameter estimates. 
fit.path      <- '/fs/data3/caverill/FIA7_Ndep_data.processed/east.data_analysis_output/FIA7_EAST_053017_soilC.summary.Rdata'
filtered.path <- '/fs/data3/caverill/FIA7_Ndep_data.processed/east.data_analysis_output/FIA7_EAST_053017_soilC.filtered.rds'
bins.path     <- '/fs/data3/caverill/FIA7_Ndep_data.processed/east.data_analysis_output/east_soil.myc_bins.rds'
 
#load growth fit and predictions.
s.out    <- readRDS(fit.path)

#subsetresponses to EM abundance at different levels of N-deposition
s.ndep    <-    exp(s.out[grep('y.nmean', rownames(s.out)),])  #mean level of N-dep
s.ndep.lo <-    exp(s.out[grep('y.nlow' , rownames(s.out)),])  #lowest level of N-dep ~ 1.2 kg N / ha / yr
s.ndep.hi <-    exp(s.out[grep('y.nhigh', rownames(s.out)),])  #highest lvel of N-dep ~17.8 kg N / ha / yr


#we need the Ndep ranges of the original JAGS predictors predictor_ranges
em.range <- seq(0,1,by=0.01)

#grab multiple regression parameter estimates
s.fit <- s.out[1:11,1:4]

#load soil data.
d <- data.table(readRDS(filtered.path))
d <- d[d$tot.15 < 20,] #already done in analysis script.
d <- data.table(d)

#get low and high N-dep means used to generate spline fits.
n.lo <- min(d$tot.15)
n.hi <- max(d$tot.15)

#detrend the scores to isolate effects of nitrogen, EM and AM
d$detrend  <- log(d$C.storage) -   (s.fit[1,4]
                                  + s.fit[2,4]*d$mat   
                                  + s.fit[3,4]*d$map      
                                  + s.fit[4,4]*d$pH_H2O       
                                  + s.fit[5,4]*d$tot.15
                                  + s.fit[6,4]*d$relEM
                                  + s.fit[7,4]*log(d$N.storage)
                                  + s.fit[8,4]*log(d$N.storage)*d$relEM
                                  + s.fit[9,4]*d$tot.15*d$relEM
                                  + s.fit[10,4]*d$clay
                                  ) + 
                    (s.fit[1,4]
                   + s.fit[2,4]*mean(d$mat)   
                   + s.fit[3,4]*mean(d$map)     
                   + s.fit[4,4]*mean(d$pH_H2O)       
                   + s.fit[5,4]*mean(d$tot.15)
                   + s.fit[6,4]*d$relEM
                   + s.fit[7,4]*mean(log(d$N.storage))
                   + s.fit[8,4]*mean(log(d$N.storage))*d$relEM
                   + s.fit[9,4]*mean(d$tot.15)*d$relEM
                   + s.fit[10,4]*mean(d$clay)
                   )
d$detrend  <- exp(d$detrend)

#detrend effect of EM at low N-dep
d$detrend.lo  <- log(d$C.storage) -  (s.fit[1,4]
                                    + s.fit[2,4]*d$mat   
                                    + s.fit[3,4]*d$map      
                                    + s.fit[4,4]*d$pH_H2O       
                                    + s.fit[5,4]*d$tot.15
                                    + s.fit[6,4]*d$relEM
                                    + s.fit[7,4]*log(d$N.storage)
                                    + s.fit[8,4]*log(d$N.storage)*d$relEM
                                    + s.fit[9,4]*d$tot.15*d$relEM
                                    + s.fit[10,4]*d$clay
                                    ) + 
              (s.fit[1,4]
             + s.fit[2,4]*mean(d$mat)   
             + s.fit[3,4]*mean(d$map)     
             + s.fit[4,4]*mean(d$pH_H2O)       
             + s.fit[5,4]*n.lo
             + s.fit[6,4]*d$relEM
             + s.fit[7,4]*mean(log(d$N.storage))
             + s.fit[8,4]*mean(log(d$N.storage))*d$relEM
             + s.fit[9,4]*n.lo*d$relEM
             + s.fit[10,4]*mean(d$clay)
             )
d$detrend.lo  <- exp(d$detrend.lo)

#detrend effect of EM at low N-dep
d$detrend.hi  <- log(d$C.storage) -   (s.fit[1,4]
                                       + s.fit[2,4]*d$mat   
                                       + s.fit[3,4]*d$map      
                                       + s.fit[4,4]*d$pH_H2O       
                                       + s.fit[5,4]*d$tot.15
                                       + s.fit[6,4]*d$relEM
                                       + s.fit[7,4]*log(d$N.storage)
                                       + s.fit[8,4]*log(d$N.storage)*d$relEM
                                       + s.fit[9,4]*d$tot.15*d$relEM
                                       + s.fit[10,4]*d$clay
                                       ) + 
              (s.fit[1,4]
               + s.fit[2,4]*mean(d$mat)   
               + s.fit[3,4]*mean(d$map)     
               + s.fit[4,4]*mean(d$pH_H2O)       
               + s.fit[5,4]*n.hi
               + s.fit[6,4]*d$relEM
               + s.fit[7,4]*mean(log(d$N.storage))
               + s.fit[8,4]*mean(log(d$N.storage))*d$relEM
               + s.fit[9,4]*n.hi*d$relEM
               + s.fit[10,4]*mean(d$clay)
               )
d$detrend.hi  <- exp(d$detrend.hi)

#save filtered for downstream plotting, if desired.
saveRDS(d, filtered.path)

#get binned means
n.bins <- 20
d$cat <- cut(d$relEM, n.bins, labels=F)
d$cat <- as.factor(d$cat)


#get response to relEM, at low, high and mean levels of N-deposition.
bins.myc               <- aggregate(relEM      ~ cat, data = d, FUN = 'mean'  )
bins.myc$detrend       <- aggregate(detrend    ~ cat, data = d, FUN = 'mean'  )[,2]
bins.myc$detrend.sd    <- aggregate(detrend    ~ cat, data = d, FUN = 'sd'    )[,2]
bins.myc$detrend.lo    <- aggregate(detrend.lo ~ cat, data = d, FUN = 'mean'  )[,2]
bins.myc$detrend.lo.sd <- aggregate(detrend.lo ~ cat, data = d, FUN = 'sd'    )[,2]
bins.myc$detrend.hi    <- aggregate(detrend.hi ~ cat, data = d, FUN = 'mean'  )[,2]
bins.myc$detrend.hi.sd <- aggregate(detrend.hi ~ cat, data = d, FUN = 'sd'    )[,2]
bins.myc$detrend.n     <- aggregate(detrend    ~ cat, data = d, FUN = 'length')[,2]
bins.myc$detrend.upr <- bins.myc$detrend + bins.myc$detrend.sd/sqrt(bins.myc$detrend.n)
bins.myc$detrend.lwr <- bins.myc$detrend - bins.myc$detrend.sd/sqrt(bins.myc$detrend.n)
bins.myc$detrend.lo.upr <- bins.myc$detrend.lo + bins.myc$detrend.lo.sd/sqrt(bins.myc$detrend.n)
bins.myc$detrend.lo.lwr <- bins.myc$detrend.lo - bins.myc$detrend.lo.sd/sqrt(bins.myc$detrend.n)
bins.myc$detrend.hi.upr <- bins.myc$detrend.hi + bins.myc$detrend.hi.sd/sqrt(bins.myc$detrend.n)
bins.myc$detrend.hi.lwr <- bins.myc$detrend.hi - bins.myc$detrend.hi.sd/sqrt(bins.myc$detrend.n)

#remove any bins that only have one observation - not a problem here.
bins.myc <- bins.myc[!(bins.myc$detrend.n < 2),]

#save growth bins data
saveRDS(bins.myc, bins.path)

#graphic settings
par(mfrow=c(1,2))
cols <- wes_palette("Zissou", 5)
trans <- 0.3
limy <- c(1000,8500)
limx <- c(0,1)

#plot Soil C ~ mycorrhizal abundance at low N-deposition
plot(d$detrend.lo ~ d$relEM, ylim = limy, xlim = limx, col = 'gray', cex = 0.2)
par(new=T)
plot(bins.myc$detrend.lo ~ bins.myc$relEM, pch = 16,
     ylim=limy, xlim=limx)
arrows(bins.myc$relEM, bins.myc$detrend.lo.lwr, bins.myc$relEM, bins.myc$detrend.lo.upr, length=0.05, angle=90, code=3)
#drop fit on top
lines(smooth.spline(s.ndep.lo[,4] ~ em.range) , lwd=2, col=cols[5])
polygon(c(em.range, rev(em.range)),c(s.ndep.lo[,3], rev(s.ndep.lo[,1])), col=adjustcolor(cols[5], trans), lty=0)
lines(smooth.spline(s.ndep.lo[,3] ~ em.range), lty=2, lwd=0.5)
lines(smooth.spline(s.ndep.lo[,1] ~ em.range), lty=2, lwd=0.5)

#plot Soil C ~ mycorrhizal abundance at high N-deposition
plot(d$detrend.hi ~ d$relEM, ylim = limy, xlim = limx, col = 'gray', cex = 0.2)
par(new=T)
plot(bins.myc$detrend.hi ~ bins.myc$relEM, pch = 16,
     ylim=limy, xlim=limx)
arrows(bins.myc$relEM, bins.myc$detrend.hi.lwr, bins.myc$relEM, bins.myc$detrend.hi.upr, length=0.05, angle=90, code=3)
#drop fit on top
lines(smooth.spline(s.ndep.hi[,4] ~ em.range) , lwd=2, col=cols[2])
polygon(c(em.range, rev(em.range)),c(s.ndep.hi[,3], rev(s.ndep.hi[,1])), col=adjustcolor(cols[2], trans), lty=0)
lines(smooth.spline(s.ndep.hi[,3] ~ em.range), lty=2, lwd=0.5)
lines(smooth.spline(s.ndep.hi[,1] ~ em.range), lty=2, lwd=0.5)

#lo-ndep R2 - 0.63
#lo.spline <- smooth.spline(s.ndep.lo[,4] ~ em.range)
#summary(lm(bins.myc$detrend.lo ~ predict(lo.spline, bins.myc$relEM)$y))
#summary(lm(d$detrend.lo ~ predict(lo.spline, d$relEM)$y)) #raw data R2 = 0.08

#hi-ndep R2 - NS - -0.04
#hi.spline <- smooth.spline(s.ndep.hi[,4] ~ em.range)
#summary(lm(bins.myc$detrend.hi ~ predict(hi.spline, bins.myc$relEM)$y))
#summary(lm(d$detrend.hi ~ predict(hi.spline, d$relEM)$y)) #raw data R2 = 0.005 - negative correlation predict vs. obs...