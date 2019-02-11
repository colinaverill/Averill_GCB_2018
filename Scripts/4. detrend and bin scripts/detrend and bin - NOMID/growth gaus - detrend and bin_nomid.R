#plotting gaussian perturbationgrowth by mycorrhizal type, detrending for other factors.
rm(list=ls())
library(data.table)
library(runjags)
library(coda)
library(wesanderson)
library(boot)        #for inv.logit command
source('required_products_utilities/master_path_list.r')

#specify paths for the Ndep prediction, the data, and the model parameter estimates.
   model.path <- local_growth.gaus_model_path.nomid
     fit.path <- local_growth.gaus_summary_path.nomid
   range.path <- local_growth.gaus_range_path.nomid
filtered.path <- local_growth.gaus_filtered_path.nomid
    pred.path <- local_growth.gaus_prediction_path.nomid
 bins.am.path <- bins_growth.gaus_nomid_AM_path
 bins.em.path <- bins_growth.gaus_all_EM_path

#load growth fit and predictions.
g.out    <- readRDS(fit.path)

#subset responses to N - no longer the case because we have to drop a chain
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
g.fit <- pred.out[1:14,1:4]
#put them back on a log scale.
#g.fit <- log(g.fit)

#load growth data.
d <- data.table(readRDS(filtered.path))
d <- d[d$tot.15 < 20,]
d <- d[!(PLT_CN == 168931301010661),]
d <- data.table(d)

#detrend the scores to isolate effects of nitrogen, EM and AM
d$detrend.em  <- log(d$growth) -    (g.fit[1,4]
                                   + g.fit[2,4]*log(d$basal.m2)       
                                   + g.fit[3,4]*d$mat      
                                   + g.fit[4,4]*d$map       
                                   + g.fit[5,4]*d$cn      
                                   + g.fit[6,4]*d$pH_H2O
                                   + g.fit[7,4]*d$relEM
                                   + g.fit[8,4]*d$tot.15
                                   + g.fit[9,4]*d$relEM*d$tot.15
                                   + g.fit[10,4]*d$clay
                                   + g.fit[11,4] * exp(- ((d$tot.15 - g.fit[12,4])^2/(2*g.fit[13,4]^2)))) + 
                (g.fit[1,4]
               + g.fit[2,4]*mean(log(d$basal.m2)) 
               + g.fit[3,4]*mean(d$mat) 
               + g.fit[4,4]*mean(d$map) 
               + g.fit[5,4]*mean(d$cn) 
               + g.fit[6,4]*mean(d$pH_H2O)  
               + g.fit[7,4]*1
               + g.fit[8,4]*d$tot.15
               + g.fit[9,4]*1*d$tot.15
               + g.fit[10,4]*mean(d$clay)
               + g.fit[11,4] * exp(- ((d$tot.15 - g.fit[12,4])^2/(2*g.fit[13,4]^2))))
  d$detrend.em  <- exp(d$detrend.em)
  
d$detrend.am  <- log(d$growth) -        (g.fit[1,4]
                                       + g.fit[2,4]*log(d$basal.m2)       
                                       + g.fit[3,4]*d$mat      
                                       + g.fit[4,4]*d$map       
                                       + g.fit[5,4]*d$cn      
                                       + g.fit[6,4]*d$pH_H2O
                                       + g.fit[7,4]*d$relEM
                                       + g.fit[8,4]*d$tot.15
                                       + g.fit[9,4]*d$relEM*d$tot.15
                                       + g.fit[10,4]*d$clay
                                       + g.fit[11,4] * exp(- ((d$tot.15 - g.fit[12,4])^2/(2*g.fit[13,4]^2)))) + 
      (g.fit[1,4]
     + g.fit[2,4]*mean(log(d$basal.m2)) 
     + g.fit[3,4]*mean(d$mat) 
     + g.fit[4,4]*mean(d$map) 
     + g.fit[5,4]*mean(d$cn) 
     + g.fit[6,4]*mean(d$pH_H2O)  
     + g.fit[7,4]*0
     + g.fit[8,4]*d$tot.15
     + g.fit[9,4]*0*d$tot.15
     + g.fit[10,4]*mean(d$clay)
     + g.fit[11,4] * exp(- ((d$tot.15 - g.fit[12,4])^2/(2*g.fit[13,4]^2))))
d$detrend.am  <- exp(d$detrend.am)

#get binned means
n.bins <- 9
d$cat <- cut(d$tot.15,n.bins, labels = F)
d$cat<- as.factor(d$cat)


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
par(new = T)
plot(d$detrend.am ~ d$tot.15, pch = 16, cex = 0.3, col = 'light gray', ylim = c(0,.25),xlim=c(min(d$tot.15),max(d$tot.15)))
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
par(new = T)
plot(d$detrend.em ~ d$tot.15, pch = 16, cex = 0.3, col = 'light gray', ylim = c(0,.25),xlim=c(min(d$tot.15),max(d$tot.15)))
#drop fit on top
lines(smooth.spline(g.ndep.em[,4] ~ g.range) , lwd=2, col=cols[2])
polygon(c(g.range, rev(g.range)),c(g.ndep.em[,3], rev(g.ndep.em[,1])), col=adjustcolor(cols[2], trans), lty=0)
lines(smooth.spline(g.ndep.em[,3] ~ g.range), lty=2, lwd=0.5)
lines(smooth.spline(g.ndep.em[,1] ~ g.range), lty=2, lwd=0.5)

#AM R2 - 0.98
am.spline <- smooth.spline(g.ndep.am[,4] ~ g.range)
summary(lm(bins.am$detrend ~ predict(am.spline, bins.am$tot.15)$y))
am.gaus.r2 <- summary(lm(bins.am$detrend ~ predict(am.spline, bins.am$tot.15)$y))$r.squared

#EM R2 - 0.85
em.spline <- smooth.spline(g.ndep.em[,4] ~ g.range)
summary(lm(bins.em$detrend ~ predict(em.spline, bins.em$tot.15)$y))
em.gaus.r2 <- summary(lm(bins.em$detrend ~ predict(em.spline, bins.em$tot.15)$y))$r.squared

#save R2 values
analyses <- c('growth.gaus.AM','growth.gaus.EM','growth.AM','growth.EM','mortality','beta','recruitment.AM','recruitment.EM','soilC.lo','soilC.hi')
out <- data.frame(analyses,r.sq = NA)

out[analyses == 'growth.gaus.AM',]$r.sq <- am.gaus.r2
out[analyses == 'growth.gaus.EM',]$r.sq <- em.gaus.r2

saveRDS(out, nomid_binned_rsq_summary_path)
