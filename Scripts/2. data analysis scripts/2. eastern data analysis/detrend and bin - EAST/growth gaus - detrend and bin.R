#NOTE. INSUFFICIENT CONVERGENCE OF THIS MODEL.
#CONSIDER PLOTTING NO-GAUSSIAN, SEE IF GAUS IS NEEDED OR IF DIFFERENT STARTING VALUES.
#Plotting of >90% AM or EM plots doesnt support a gaussian perturbation here.
#Too few all AM or EM plots in the nomid data for the plots to look compelling. 
#plotting gaussian perturbationgrowth by mycorrhizal type, detrending for other factors.
rm(list=ls())
library(data.table)
library(runjags)
library(wesanderson)
library(boot)        #for inv.logit command

#specify paths for the Ndep prediction, the data, and the model parameter estimates. 
     fit.path <- '/fs/data3/caverill/FIA7_Ndep_data.processed/nomid_east.data_analysis_output/FIA7_051517_gaus_growth.model.Rdata'
   range.path <- '/fs/data3/caverill/FIA7_Ndep_data.processed/nomid_east.data_analysis_output/FIA7_051517_gaus_growth.predictor_ranges.rds'
filtered.path <- '/fs/data3/caverill/FIA7_Ndep_data.processed/nomid_east.data_analysis_output/FIA7_051517_gaus_growth.filtered.rds'
 bins.am.path <- '/fs/data3/caverill/FIA7_Ndep_data.processed/nomid_east.data_analysis_output/all_growth.gaus.am_bins.rds'
 bins.em.path <- '/fs/data3/caverill/FIA7_Ndep_data.processed/nomid_east.data_analysis_output/all_growth.gaus.em_bins.rds'
    pred.path <- '/fs/data3/caverill/FIA7_Ndep_data.processed/nomid_east.data_analysis_output/FIA7_051517_gaus_growth.prediction.Rdata'
#load growth fit and predictions.
g.out    <- readRDS(fit.path)

#subset responses to N - no longer the case cause we have to drop a chain
#g.ndep.am <-    g.out[grep('y.x7.am', rownames(g.out)),]
#g.ndep.em <-    g.out[grep('y.x7.em', rownames(g.out)),]
#g.ndep    <-    g.out[grep('y.x7'   , rownames(g.out)),]
#g.ndep    <-   g.ndep[1:101,]

#trim prediction path, only using parameters from the two chains that actually converge.
pred.model <-readRDS(pred.path)
pred.trim <- as.mcmc.list(pred.model)
pred.trim <- mcmc.list(pred.trim[[1]],pred.trim[[3]])
out <- summary(pred.trim)
mean.out <- out$statistics
cred.out <- out$quantiles
pred.out <- data.frame(cred.out[,1],cred.out[,3],cred.out[,5],mean.out[,1])
colnames(pred.out) <- c('2.5%','50%','97.5%','mean')

#grab the mean responses to N in AM and EM conditions
ndep <- pred.out[grep('y.x7',rownames(pred.out)),]
g.ndep.am <- exp(ndep[102:202,])
g.ndep.em <- exp(ndep[203:303,])
g.ndep    <- exp(ndep[1:101,])

#we need the Ndep ranges of the original JAGS predictors predictor_ranges
range <- readRDS(range.path)
g.range <- range$x7.range

#grab multiple regression parameter estimates
#g.fit <- g.out[1:11,1:4]
g.fit <- pred.out[1:11,1:4]
#put them back on a log scale.
#g.fit <- log(g.fit)

#load growth data.
d <- data.table(readRDS(filtered.path))
d <- d[d$tot.15 < 20,]
d<- data.table(d)

#detrend the scores to isolate effects of nitrogen, EM and AM
d$detrend  <- log(d$growth) - (g.fit[2,4]*log(d$basal.m2)       
                               + g.fit[3,4]*d$mat      
                               + g.fit[4,4]*d$map       
                               + g.fit[5,4]*d$cn      
                               + g.fit[6,4]*d$pH_H2O        
                               + g.fit[10,4]*d$clay) + 
  (g.fit[2,4]*mean(log(d$basal.m2)) 
   + g.fit[3,4]*mean(d$mat) 
   + g.fit[4,4]*mean(d$map) 
   + g.fit[5,4]*mean(d$cn) 
   + g.fit[6,4]*mean(d$pH_H2O)  
   + g.fit[10,4]*mean(d$clay))
d$detrend  <- exp(d$detrend)

#get binned means
n.bins <- 9
d$cat <- cut(d$tot.15,n.bins, labels = F)
d$cat<- as.factor(d$cat)

#get plots that are >90% EM or >90% AM for plotting data on interaction
d.am <- d[relEM < 0.1,]
d.em <- d[relEM > 0.9,]

#get am vs. em data sets, as not each subset has all of the bins present so mismatches in dataframe size.
bins.em             <- aggregate(tot.15     ~ cat, data = d.em, FUN = 'mean'  )
bins.em$detrend     <- aggregate(detrend    ~ cat, data = d.em, FUN = 'mean'  )[,2]
bins.em$detrend.sd  <- aggregate(detrend    ~ cat, data = d.em, FUN = 'sd'    )[,2]
bins.em$detrend.n   <- aggregate(detrend    ~ cat, data = d.em, FUN = 'length')[,2]
bins.em$detrend.upr <- bins.em$detrend + bins.em$detrend/sqrt(bins.em$detrend.n)
bins.em$detrend.lwr <- bins.em$detrend - bins.em$detrend/sqrt(bins.em$detrend.n)
bins.am             <- aggregate(tot.15     ~ cat, data = d.am, FUN = 'mean'  )
bins.am$detrend     <- aggregate(detrend    ~ cat, data = d.am, FUN = 'mean'  )[,2]
bins.am$detrend.sd  <- aggregate(detrend    ~ cat, data = d.am, FUN = 'sd'    )[,2]
bins.am$detrend.n   <- aggregate(detrend    ~ cat, data = d.am, FUN = 'length')[,2]
bins.am$detrend.upr <- bins.am$detrend + bins.am$detrend/sqrt(bins.am$detrend.n)
bins.am$detrend.lwr <- bins.am$detrend - bins.am$detrend/sqrt(bins.am$detrend.n)

#remove any bins that only have one observation
bins.am <- bins.am[!(bins.am$detrend.n < 2),]
bins.em <- bins.em[!(bins.em$detrend.n < 2),]

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

#AM R2 - 0.47
#am.spline <- smooth.spline(g.ndep.am[,4] ~ g.range)
#summary(lm(bins.am$detrend ~ predict(am.spline, bins.am$tot.15)$y))

#EM R2 - 0.63
#em.spline <- smooth.spline(g.ndep.em[,4] ~ g.range)
#summary(lm(bins.em$detrend ~ predict(em.spline, bins.em$tot.15)$y))
