#plotting mortality means
#detrending using binned means, since we cant correct 0-1 observations.
#clear environemnt, load packages
rm(list=ls())
library(data.table)
library(runjags)
library(wesanderson) #wes anderson of course. 
library(boot)        #for inv.logit command
library(binom)       #for binomial confidence intervals.

#specify paths.
fit.path      <- '/fs/data3/caverill/FIA7_Ndep_data.processed/nomid_east.data_analysis_output/FIA7_NOMID_EAST_052617_runjags_NOMID_EAST_mortality.summary.Rdata'
range.path    <- '/fs/data3/caverill/FIA7_Ndep_data.processed/nomid_east.data_analysis_output/FIA7_NOMID_EAST_052617_runjags_NOMID_EAST_mortality.predictor_ranges.rds'
filtered.path <- '/fs/data3/caverill/FIA7_Ndep_data.processed/nomid_east.data_analysis_output/FIA7_NOMID_EAST_052617_runjags_NOMID_EAST_mortality.filtered.rds'
bins.path     <- '/fs/data3/caverill/FIA7_Ndep_data.processed/nomid_east.data_analysis_output/nomid_mort_bins.rds'

#load mortality fit
m.out    <- readRDS(fit.path)
m.ndep.am <- m.out[grep('x6.am', rownames(m.out)),]
m.ndep.em <- m.out[grep('x6.em', rownames(m.out)),]
m.ndep    <- m.out[grep('x6',    rownames(m.out)),]
m.ndep    <- m.ndep[1:101,]

#we need the Ndep ranges of the original JAGS predictors predictor_ranges
m.range <- readRDS(range.path)
m.range <- m.range$range.x6

#grab mortality multiple regression parameter estimates
m.fit <- m.out[1:11,]

#load mortality data.
d <- data.table(readRDS(filtered.path))
d <- d[tot.15 < 20,]

#get EM vs. AM deaths, total AM indicator
d$am <- 1 - d$em
d$death.em <- d$em * d$death
d$death.am <- d$death - d$death.em

#get binned Ndep categories
d$cat <- cut(d$tot.15, 7, labels = F)
d$cat<- as.factor(d$cat)
d$REMPER <- as.numeric(d$REMPER)

#start binning the raw data.
bins                  <- aggregate(tot.15                  ~ cat, data = d, FUN = 'mean'  )
bins$pre.m.rate       <- aggregate(death                   ~ cat, data = d, FUN = 'mean'  )[,2]
bins$deaths.am        <- aggregate(death.am                ~ cat, data = d, FUN = 'sum'   )[,2]
bins$deaths.em        <- aggregate(death.em                ~ cat, data = d, FUN = 'sum'   )[,2]
bins$n.total          <- aggregate(death                   ~ cat, data = d, FUN = 'length')[,2]
bins$n.em             <- aggregate(em                      ~ cat, data = d, FUN = 'sum'   )[,2]
bins$n.am             <- aggregate(am                      ~ cat, data = d, FUN = 'sum'   )[,2]
bins$REMPER           <- aggregate(REMPER                  ~ cat, data = d, FUN = 'mean'  )[,2]

#get m.rates
bins$pre.m.rate.am <- bins$deaths.am / bins$n.am
bins$pre.m.rate.em <- bins$deaths.em / bins$n.em

#get 95% confidence interval on mortality rate per bin.
bins$pre.m.rate.lwr    <- binom.confint(bins$deaths.am + bins$deaths.em, bins$n.total, conf.level = 0.95, method = 'exact')[,5]
bins$pre.m.rate.upr    <- binom.confint(bins$deaths.am + bins$deaths.em, bins$n.total, conf.level = 0.95, method = 'exact')[,6]
bins$pre.m.rate.am.lwr <- binom.confint(bins$deaths.am, bins$n.am, conf.level = 0.95, method = 'exact')[,5]
bins$pre.m.rate.am.upr <- binom.confint(bins$deaths.am, bins$n.am, conf.level = 0.95, method = 'exact')[,6]
bins$pre.m.rate.em.lwr <- binom.confint(bins$deaths.em, bins$n.em, conf.level = 0.95, method = 'exact')[,5]
bins$pre.m.rate.em.upr <- binom.confint(bins$deaths.em, bins$n.em, conf.level = 0.95, method = 'exact')[,6]

#correct estimated mortality rates using the average REMPER
bins$m.rate        <- 1 - (1-bins$pre.m.rate   )^(1/bins$REMPER)
bins$m.rate.am     <- 1 - (1-bins$pre.m.rate.am)^(1/bins$REMPER)
bins$m.rate.em     <- 1 - (1-bins$pre.m.rate.em)^(1/bins$REMPER)
bins$m.rate.lwr    <- 1 - (1-bins$pre.m.rate.lwr)^(1/bins$REMPER)
bins$m.rate.upr    <- 1 - (1-bins$pre.m.rate.upr)^(1/bins$REMPER)
bins$m.rate.am.lwr <- 1 - (1-bins$pre.m.rate.am.lwr)^(1/bins$REMPER)
bins$m.rate.am.upr <- 1 - (1-bins$pre.m.rate.am.upr)^(1/bins$REMPER)
bins$m.rate.em.lwr <- 1 - (1-bins$pre.m.rate.em.lwr)^(1/bins$REMPER)
bins$m.rate.em.upr <- 1 - (1-bins$pre.m.rate.em.upr)^(1/bins$REMPER)


#bin up associated predictors for detrending.
bins$mat       <- aggregate(mat     ~ cat, data = d, FUN = 'mean')[,2]
bins$map       <- aggregate(map     ~ cat, data = d, FUN = 'mean')[,2]
bins$cn        <- aggregate(cn      ~ cat, data = d, FUN = 'mean')[,2]
bins$ph        <- aggregate(pH_H2O  ~ cat, data = d, FUN = 'mean')[,2]
bins$clay      <- aggregate(clay    ~ cat, data = d, FUN = 'mean')[,2]
bins$DIA       <- aggregate(DIA     ~ cat, data = d, FUN = 'mean')[,2]

#Detrend the means!
resid <- exp(m.fit[1,4] + m.fit[2,4]*bins$DIA)    + exp(-m.fit[3,4]*bins$DIA)    + m.fit[4,4]*bins$mat    + m.fit[5,4]*bins$map    + m.fit[6,4]*bins$ph        + m.fit[7,4]*bins$cn    + m.fit[11,4]*bins$clay
mean  <- exp(m.fit[1,4] + m.fit[2,4]*mean(d$DIA)) + exp(-m.fit[3,4]*mean(d$DIA)) + m.fit[4,4]*mean(d$mat) + m.fit[5,4]*mean(d$map) + m.fit[6,4]*mean(d$pH_H2O) + m.fit[7,4]*mean(d$cn) + m.fit[11,4]*mean(d$clay)
bins$detrend    <- inv.logit(logit(bins$m.rate)    - resid + mean)
bins$detrend.am <- inv.logit(logit(bins$m.rate.am) - resid + mean)
bins$detrend.em <- inv.logit(logit(bins$m.rate.em) - resid + mean)
bins$detrend.lwr    <- inv.logit(log(bins$m.rate.lwr) - resid + mean)
bins$detrend.upr    <- inv.logit(log(bins$m.rate.upr) - resid + mean)
bins$detrend.am.lwr <- inv.logit(logit(bins$m.rate.am.lwr) - resid + mean)
bins$detrend.am.upr <- inv.logit(logit(bins$m.rate.am.upr) - resid + mean)
bins$detrend.em.lwr <- inv.logit(logit(bins$m.rate.em.lwr) - resid + mean)
bins$detrend.em.upr <- inv.logit(logit(bins$m.rate.em.upr) - resid + mean)

#save mortality bins data
saveRDS(bins, bins.path)


#make some plots. 
#some graphic parameters
cols <- wes_palette("Zissou", 5)
trans <- 0.3
par(mfrow=c(1,2))

#AM fit.
plot(bins$detrend.am ~ bins$tot.15, pch = 16, 
     ylim=c(0,0.025), xlim=c(min(d$tot.15),max(d$tot.15)))
arrows(bins$tot.15, bins$detrend.am.lwr, bins$tot.15, bins$detrend.am.upr, length=0.05, angle=90, code=3)
lines(smooth.spline(inv.logit(m.ndep.am[,4]) ~ m.range) , lwd=2, col=cols[5])
polygon(c(m.range, rev(m.range)),c(inv.logit(m.ndep.am[,3]), rev(inv.logit(m.ndep.am[,1]))), col=adjustcolor(cols[5], trans), lty=0)
lines(smooth.spline(inv.logit(m.ndep.am[,3]) ~ m.range), lty=2, lwd=0.5)
lines(smooth.spline(inv.logit(m.ndep.am[,1]) ~ m.range), lty=2, lwd=0.5)

#EM fit.
plot(bins$detrend.em ~ bins$tot.15, pch = 16,
     ylim=c(0,0.025), xlim=c(min(d$tot.15),max(d$tot.15)))
arrows(bins$tot.15, bins$detrend.em.lwr, bins$tot.15, bins$detrend.em.upr, length=0.05, angle=90, code=3)
lines(smooth.spline(inv.logit(m.ndep.em[,4]) ~ m.range) , lwd=2, col=cols[2])
polygon(c(m.range, rev(m.range)),c(inv.logit(m.ndep.em[,3]), rev(inv.logit(m.ndep.em[,1]))), col=adjustcolor(cols[2], trans), lty=0)
lines(smooth.spline(inv.logit(m.ndep.em[,3]) ~ m.range), lty=2, lwd=0.5)
lines(smooth.spline(inv.logit(m.ndep.em[,1]) ~ m.range), lty=2, lwd=0.5)

#plot aggregated job, since these are not actually different between EM and AM trees.
par(mfrow=c(1,1))
plot(bins$detrend ~ bins$tot.15, pch = 16,
     ylim=c(0,0.025), xlim=c(min(d$tot.15),max(d$tot.15)))
arrows(bins$tot.15, bins$detrend.lwr, bins$tot.15, bins$detrend.upr, length=0.05, angle=90, code=3)
lines(smooth.spline(inv.logit(m.ndep[,4]) ~ m.range) , lwd=2, col=cols[2])
polygon(c(m.range, rev(m.range)),c(inv.logit(m.ndep[,3]), rev(inv.logit(m.ndep[,1]))), col=adjustcolor(cols[2], trans), lty=0)
lines(smooth.spline(inv.logit(m.ndep[,3]) ~ m.range), lty=2, lwd=0.5)
lines(smooth.spline(inv.logit(m.ndep[,1]) ~ m.range), lty=2, lwd=0.5)

#This gives you the R2 of the binned means - NA
thing <- smooth.spline(inv.logit(m.ndep[,4]) ~ m.range)
summary(lm(bins$detrend ~ predict(thing, bins$tot.15)$y))