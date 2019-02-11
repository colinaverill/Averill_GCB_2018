#N dep effects on the abundance of EM fungi
rm(list=ls())
library(betareg)
library(data.table)
library(runjags)
library(wesanderson) #wes anderson of course. 
library(boot) #for inv.logit command

     fit.path <- '/fs/data3/caverill/FIA7_Ndep_data.processed/nomid_east.data_analysis_output/FIA7_nomid_east_053017_rel.abundance.summary.rds'
   range.path <- '/fs/data3/caverill/FIA7_Ndep_data.processed/nomid_east.data_analysis_output/FIA7_nomid_east_053017_rel.abundance.ranges.rds'
filtered.path <- '/fs/data3/caverill/FIA7_Ndep_data.processed/nomid_east.data_analysis_output/FIA7_nomid_east_053017_rel.abundance.filtered.rds'
    bins.path <- '/fs/data3/caverill/FIA7_Ndep_data.processed/nomid_east.data_analysis_output/nomid_beta_bins.rds'

#load beta fit.
b.out   <- readRDS(fit.path)
b.ndep  <- b.out[grep('x5'  , rownames(b.out)),]
b.mean  <- inv.logit(b.ndep[,4])
b.upr   <- inv.logit(b.ndep[,3])
b.lwr   <- inv.logit(b.ndep[,1])
b.range <- readRDS(range.path)
b.range <- b.range$ndep.range

#load beta multiple regression parameter estimates
b.fit <- b.out[1:8,]

#load beta data.
d <- data.table(readRDS(filtered.path))
d <- d[d$tot.15 < 20,]
#remove a duplicate column
d <- d[,-4]
d<- data.table(d)

#logit transform rel. abundance scores for detrending with model parameters.
d$score <- logit(d$relEM.beta)

#detrend the scores to isolate nitrogen.
d$detrend <- d$score - (b.fit[2,4]*d$mat       
                        + b.fit[3,4]*d$map       
                        + b.fit[4,4]*d$cn       
                        + b.fit[5,4]*d$pH_H2O       
                        + b.fit[7,4]*d$clay) + 
                       (b.fit[2,4]*mean(d$mat) 
                        + b.fit[3,4]*mean(d$map) 
                        + b.fit[4,4]*mean(d$cn) 
                        + b.fit[5,4]*mean(d$pH_H2O) 
                        + b.fit[7,4]*mean(d$clay))
d$relEM.detrend <- inv.logit(d$detrend)


#get binned means and standard deviation.
d$cat <- cut(d$tot.15,10, labels = F)
d$cat<- as.factor(d$cat)

#start aggregating detrended values
bins                   <- aggregate(tot.15             ~ cat, data = d, FUN = 'mean')
bins$relEM.detrend     <- aggregate(inv.logit(detrend) ~ cat, data = d, FUN = 'mean')[,2]
bins$relEM.detrend.sd  <- aggregate(         (detrend) ~ cat, data = d, FUN = 'sd'  )[,2]
bins$n                 <- aggregate(inv.logit(detrend) ~ cat, data = d, FUN = 'length'  )[,2]
bins$relEM.detrend.upr <- inv.logit(logit(bins$relEM.detrend) + bins$relEM.detrend.sd/sqrt(bins$n))
bins$relEM.detrend.lwr <- inv.logit(logit(bins$relEM.detrend) - bins$relEM.detrend.sd/sqrt(bins$n))
#do the same with the raw values
bins$raw.em  <- aggregate(      relEM.beta  ~ cat, data = d, FUN = 'mean')[,2]
bins$raw.sd  <- aggregate(logit(relEM.beta) ~ cat, data = d, FUN = 'sd'  )[,2]
bins$raw.n   <- aggregate(logit(relEM.beta) ~ cat, data = d, FUN = 'length')[,2]
bins$raw.upr <- inv.logit(logit(bins$raw.em) + bins$raw.sd/sqrt(bins$raw.n))
bins$raw.lwr <- inv.logit(logit(bins$raw.em) - bins$raw.sd/sqrt(bins$raw.n))

#save bins output
saveRDS(bins, bins.path)


#some graphic parameters
cols <- wes_palette("Zissou", 5)
trans <- 0.3
par(mfrow=c(1,1))

#generate a plot of detrended relEM values, means and standard deviations.
#drop continuous fit on top of it.
#raw values.
plot(relEM.detrend ~ tot.15, data = d, cex = 0.3, pch = 16, col = 'gray', ylim=c(0,1), xlim=c(min(d$tot.15),max(d$tot.15)))
#binned means
par(new=T)
plot(bins$relEM.detrend ~ bins$tot.15, pch = 16,
     ylim=c(0,1), xlim=c(min(d$tot.15),max(d$tot.15)))
#error bars of binned means
arrows(bins$tot.15, bins$relEM.detrend.lwr, bins$tot.15, bins$relEM.detrend.upr, length=0.05, angle=90, code=3)

#drop fit on top
lines(smooth.spline(b.mean ~ b.range) , lwd=2, col=cols[3])
polygon(c(b.range, rev(b.range)),c(b.upr, rev(b.lwr)), col=adjustcolor(cols[3], trans), lty=0)
lines(smooth.spline(b.upr ~ b.range), lty=2, lwd=0.5)
lines(smooth.spline(b.lwr ~ b.range), lty=2, lwd=0.5)

#R2 of beta mean fit - 0.85
beta.spline <- smooth.spline(b.mean ~ b.range)
summary(lm(bins$relEM.detrend ~ predict(beta.spline, bins$tot.15)$y))