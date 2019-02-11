#plotting growth by mycorrhizal type, detrending for other factors.
rm(list=ls())
library(data.table)
library(runjags)
library(wesanderson)
library(boot)        #for inv.logit command

#specify paths for the Ndep prediction, the data, and the model parameter estimates. 
  fit.path    <- '/fs/data3/caverill/FIA7_Ndep_data.processed/nomid_east.data_analysis_output/FIA7_051517_growth.summary.Rdata'
range.path    <- '/fs/data3/caverill/FIA7_Ndep_data.processed/nomid_east.data_analysis_output/FIA7_051517_growth.predictor_ranges.rds'
filtered.path <- '/fs/data3/caverill/FIA7_Ndep_data.processed/nomid_east.data_analysis_output/FIA7_051517_growth.filtered.rds'
bins.am.path  <- '/fs/data3/caverill/FIA7_Ndep_data.processed/nomid_east.data_analysis_output/nomid_growth.am_bins.rds'
bins.em.path  <- '/fs/data3/caverill/FIA7_Ndep_data.processed/nomid_east.data_analysis_output/nomid_growth.em_bins.rds'

#load growth fit and predictions.
g.out    <- readRDS(fit.path)

#subsetresponses to N. 
g.ndep.am <-    g.out[grep('y.x7.am', rownames(g.out)),]
g.ndep.em <-    g.out[grep('y.x7.em', rownames(g.out)),]
g.ndep    <-    g.out[grep('y.x7'   , rownames(g.out)),]
g.ndep    <-   g.ndep[1:101,]

#we need the Ndep ranges of the original JAGS predictors predictor_ranges
range <- readRDS(range.path)
g.range <- range$x7.range

#grab multiple regression parameter estimates
g.fit <- g.out[1:11,1:4]
#put them back on a log scale.
g.fit <- log(g.fit)

#load growth data.
d <- data.table(readRDS(filtered.path))
d <- d[d$tot.15 < 20,]
d<- data.table(d)

#detrend the scores to isolate effects of nitrogen, EM and AM
d$detrend.am <- log(d$growth) - (g.fit[1,4]
                               + g.fit[2,4]*log(d$basal.m2)       
                               + g.fit[3,4]*d$mat      
                               + g.fit[4,4]*d$map       
                               + g.fit[5,4]*d$cn      
                               + g.fit[6,4]*d$pH_H2O
                               + g.fit[7,4]*d$relEM
                               + g.fit[8,4]*d$tot.15
                               + g.fit[9,4]*d$relEM*d$tot.15
                               + g.fit[10,4]*d$clay) + 
                         (g.fit[1,4]
                        + g.fit[2,4]*mean(log(d$basal.m2)) 
                        + g.fit[3,4]*mean(d$mat) 
                        + g.fit[4,4]*mean(d$map) 
                        + g.fit[5,4]*mean(d$cn) 
                        + g.fit[6,4]*mean(d$pH_H2O)
                        + g.fit[7,4]*0
                        + g.fit[8,4]*d$tot.15
                        + g.fit[9,4]*0*d$tot.15
                        + g.fit[10,4]*mean(d$clay))
d$detrend.am  <- exp(d$detrend.am)

d$detrend.em <- log(d$growth) - (g.fit[1,4]
                                 + g.fit[2,4]*log(d$basal.m2)       
                                 + g.fit[3,4]*d$mat      
                                 + g.fit[4,4]*d$map       
                                 + g.fit[5,4]*d$cn      
                                 + g.fit[6,4]*d$pH_H2O
                                 + g.fit[7,4]*d$relEM
                                 + g.fit[8,4]*d$tot.15
                                 + g.fit[9,4]*d$relEM*d$tot.15
                                 + g.fit[10,4]*d$clay) + 
                (g.fit[1,4]
                 + g.fit[2,4]*mean(log(d$basal.m2)) 
                 + g.fit[3,4]*mean(d$mat) 
                 + g.fit[4,4]*mean(d$map) 
                 + g.fit[5,4]*mean(d$cn) 
                 + g.fit[6,4]*mean(d$pH_H2O)
                 + g.fit[7,4]*1
                 + g.fit[8,4]*d$tot.15
                 + g.fit[9,4]*1*d$tot.15
                 + g.fit[10,4]*mean(d$clay))
d$detrend.em  <- exp(d$detrend.em)

#get binned means
n.bins <- 10
d$cat <- cut(d$tot.15,n.bins, labels = F)
d$cat<- as.factor(d$cat)


#get am vs. em data sets, as not each subset has all of the bins present so mismatches in dataframe size.
bins.em <- aggregate(tot.15 ~ cat, data = d, FUN = 'mean')
bins.em$detrend     <- aggregate(detrend.em    ~ cat, data = d, FUN = 'mean'  )[,2]
bins.em$detrend.sd  <- aggregate(detrend.em    ~ cat, data = d, FUN = 'sd'    )[,2]
bins.em$detrend.n   <- aggregate(detrend.em    ~ cat, data = d, FUN = 'length')[,2]
bins.em$detrend.upr <- bins.em$detrend + bins.em$detrend.sd/sqrt(bins.em$detrend.n)
bins.em$detrend.lwr <- bins.em$detrend - bins.em$detrend.sd/sqrt(bins.em$detrend.n)
bins.am <- aggregate(tot.15 ~ cat, data = d, FUN = 'mean')
bins.am$detrend     <- aggregate(detrend.am    ~ cat, data = d, FUN = 'mean'  )[,2]
bins.am$detrend.sd  <- aggregate(detrend.am    ~ cat, data = d, FUN = 'sd'    )[,2]
bins.am$detrend.n   <- aggregate(detrend.am    ~ cat, data = d, FUN = 'length')[,2]
bins.am$detrend.upr <- bins.am$detrend + bins.am$detrend.sd/sqrt(bins.am$detrend.n)
bins.am$detrend.lwr <- bins.am$detrend - bins.am$detrend.sd/sqrt(bins.am$detrend.n)

#remove any bins that only have one observation
bins.am <- bins.am[!(bins.am$detrend.n < 3),]
bins.em <- bins.em[!(bins.em$detrend.n < 3),]

#save growth bins data
saveRDS(bins.am, bins.am.path)
saveRDS(bins.em, bins.em.path)

#make AM vs EM plot
par(mfrow=c(1,2))
cols <- wes_palette("Zissou", 5)
trans <- 0.3

#AM
plot(bins.am$detrend ~ bins.am$tot.15, pch = 16,
     ylim=c(0,.25), xlim=c(min(d$tot.15),max(d$tot.15)))
arrows(bins.am$tot.15, bins.am$detrend.lwr, bins.am$tot.15, bins.am$detrend.upr, length=0.05, angle=90, code=3)
par(new=T)
#drop fit on top
lines(smooth.spline(g.ndep.am[,4] ~ g.range) , lwd=2, col=cols[5])
polygon(c(g.range, rev(g.range)),c(g.ndep.am[,3], rev(g.ndep.am[,1])), col=adjustcolor(cols[5], trans), lty=0)
lines(smooth.spline(g.ndep.am[,3] ~ g.range), lty=2, lwd=0.5)
lines(smooth.spline(g.ndep.am[,1] ~ g.range), lty=2, lwd=0.5)
#EM
plot(bins.em$detrend ~ bins.em$tot.15, pch = 16,
     ylim=c(0,.25), xlim=c(min(d$tot.15),max(d$tot.15)))
arrows(bins.em$tot.15, bins.em$detrend.lwr, bins.em$tot.15, bins.em$detrend.upr, length=0.05, angle=90, code=3)
#drop fit on top
lines(smooth.spline(g.ndep.em[,4] ~ g.range) , lwd=2, col=cols[2])
polygon(c(g.range, rev(g.range)),c(g.ndep.em[,3], rev(g.ndep.em[,1])), col=adjustcolor(cols[2], trans), lty=0)
lines(smooth.spline(g.ndep.em[,3] ~ g.range), lty=2, lwd=0.5)
lines(smooth.spline(g.ndep.em[,1] ~ g.range), lty=2, lwd=0.5)

#AM R2 - 0.62
am.spline <- smooth.spline(g.ndep.am[,4] ~ g.range)
summary(lm(bins.am$detrend ~ predict(am.spline, bins.am$tot.15)$y))

#EM R2 - 0.45
em.spline <- smooth.spline(g.ndep.em[,4] ~ g.range)
summary(lm(bins.em$detrend ~ predict(em.spline, bins.em$tot.15)$y))
