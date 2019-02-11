#N dep effects on the abundance EM vs AM trees - Eastern US only. 
rm(list=ls())
source('required_products_utilities/master_path_list.r')
require(boot) #for inv.logit command

#load output data for growth and mortality. 
#all outputs have been back transformed as required. 
g.out    <- readRDS(local_growth_summary_path.east)
m.out    <- readRDS(local_mortality_summary_path.east)
r.out.am <- readRDS(local_AM_recruitment_summary_path.east)
r.out.em <- readRDS(local_EM_recruitment_summary_path.east)
b.out    <- readRDS(local_beta_summary_path.east)
m.out <- m.out[,1:5]

#inverse logit transform m.out to get annual mortality probability.
m.out <- (inv.logit(m.out))

#transform with inv.logit beta regression output

#subsetresponses to N. 
g.ndep.am <-    g.out[grep('y.x7.am', rownames(g.out)),]
g.ndep.em <-    g.out[grep('y.x7.em', rownames(g.out)),]
m.ndep.am <-    m.out[grep('z.x6.am', rownames(m.out)),]
m.ndep.em <-    m.out[grep('z.x6.em', rownames(m.out)),]
m.ndep    <-    m.out[grep('z.x6'   , rownames(m.out)),]; m.ndep <- m.ndep[1:101,]
r.ndep.am <- r.out.am[grep('ndep.pred', rownames(r.out.am)),]
r.ndep.em <- r.out.em[grep('ndep.pred', rownames(r.out.em)),]
b.ndep    <-    b.out[grep('x5'  , rownames(b.out)),]
b.mean <- inv.logit(b.ndep[,4])
b.upr  <- inv.logit(b.ndep[,3])
b.lwr  <- inv.logit(b.ndep[,1])


#we need the Ndep ranges of the original JAGS predictors predictor_ranges
range <- readRDS(local_growth_range_path.east)
g.range <- range$x7.range
range <- readRDS(local_mortality_range_path.east)
m.range <- range$range.x6
r.range <- readRDS(local_AM_recruitment_range_path.east)
r.range <- r.range$x5.r
b.range <- readRDS(local_beta_range_path.east)
b.range <- b.range$ndep.range


#save dimensions, destination
png(filename=Supp_Figure1_all_path,width=7,height=7,units='in',res=300)

#begin the plots!
require(wesanderson) #wes anderson of course. 
cols <- wes_palette("Zissou", 5)
trans <- 0.3
par(mfrow=c(2,2), oma=c(5,3,1,1), mai=c(.3,.3,.3,.3))
par(mfrow=c(2,2))
#inner/outer label cex (size)
cex.ilab <- 1.0
cex.olab <- 1.2
limx <- c(2.8,20)

#relative abundance EM
limy<- c(0,1)

plot(b.mean ~ b.range, lwd=0, ylim= limy, xlim = limx, cex.axis = cex.olab)
lines(smooth.spline(b.mean ~ b.range) , lwd=2, col=cols[3])
polygon(c(b.range, rev(b.range)),c(b.upr, rev(b.lwr)), col=adjustcolor(cols[3], trans), lty=0)
lines(smooth.spline(b.upr ~ b.range), lty=2, lwd=0.5)
lines(smooth.spline(b.lwr ~ b.range), lty=2, lwd=0.5)

#labels
mtext('a.', adj=0.05,line = -1.7, side=1)
mtext(side=2,expression(paste('EM basal area  (total basal area)'^'-1')), las=0, line = 2)
mtext(expression(paste('kg N ha'^'-2',' yr'^'-1')), side=1, line = 2.5)
mtext('Relative Abundance EM', line = -1.2, cex=cex.ilab)

#mortality plot- all - grouping AM and EM.
limy <- c(0,0.03)
#am line
plot(m.ndep[,4] ~ m.range, lwd=0, ylim= limy, xlim = limx, cex.axis = cex.olab)
lines(smooth.spline(m.ndep[,4] ~ m.range) , lwd=2, col=cols[3], lty=2)
polygon(c(m.range, rev(m.range)),c(m.ndep[,3], rev(m.ndep[,1])), col=adjustcolor(cols[3], trans), lty=0)
lines(smooth.spline(m.ndep[,1] ~ m.range), lty=2, lwd=0.5)
lines(smooth.spline(m.ndep[,3] ~ m.range), lty=2, lwd=0.5)

#labels
mtext('b.', adj=0.05,line = -1.7, side=1)
mtext('Mortality Probability', line = -1.2, cex=cex.ilab)
mtext(side=2,expression(paste('year'^'-1')), las=0, line = 2)
mtext(expression(paste('kg N ha'^'-2',' yr'^'-1')), side=1, line = 2.5)

#Growth plot
limy <- c(0,0.3)
#am line
plot(g.ndep.am[,4] ~ g.range, lwd=0, ylim= limy, xlim = limx, cex.axis = cex.olab)
lines(smooth.spline(g.ndep.am[,4] ~ g.range) , lwd=2, col=cols[1])
polygon(c(g.range, rev(g.range)),c(g.ndep.am[,3], rev(g.ndep.am[,1])), col=adjustcolor(cols[1], trans), lty=0)
lines(smooth.spline(g.ndep.am[,1] ~ g.range) , lty=2, lwd=0.5)
lines(smooth.spline(g.ndep.am[,3] ~ g.range) , lty=2, lwd=0.5)
#em line
lines(smooth.spline(g.ndep.em[,4] ~ g.range) , lwd=2, col=cols[5])
polygon(c(g.range, rev(g.range)),c(g.ndep.em[,3], rev(g.ndep.em[,1])), col=adjustcolor(cols[5], trans), lty=0)
lines(smooth.spline(g.ndep.em[,1] ~ g.range) , lty=2, lwd=0.5)
lines(smooth.spline(g.ndep.em[,3] ~ g.range) , lty=2, lwd=0.5)


#labels
mtext('c.', adj=0.05,line = -1.7, side=1)
#mtext(side=2,'Growth of Surviving Trees',cex.lab = cex.olab, las=0, line = 6)
mtext(side=2,expression(paste('cm'^'2','m'^'-2',' yr'^'-1')), las=0, line = 2)
mtext('Growth of Surviving Trees', line = -1.2, cex=cex.ilab)
mtext(expression(paste('kg N ha'^'-2',' yr'^'-1')), side=1, line = 2.5)

#legend
legend(4, 0.275,c('AM','EM'), lwd=2, col=c(cols[1],cols[5]), bty='n', y.intersp = 1, x.intersp = 0.75, cex=1.2, seg.len=1)


#recruitment plot
limy <- c(0,.65)

#am line
plot(r.ndep.am[,4] ~ r.range, lwd=0, ylim= limy, xlim = limx, cex.axis = cex.olab)
lines(smooth.spline(r.ndep.am[,4] ~ r.range) , lwd=2, col=cols[1])
polygon(c(r.range, rev(r.range)),c(r.ndep.am[,4] + r.ndep.am[,5], rev(r.ndep.am[,4] - r.ndep.am[,5])), col=adjustcolor(cols[1], trans), lty=0)
lines(smooth.spline(r.ndep.am[,4] + r.ndep.am[,5] ~ r.range), lty=2, lwd=0.5)
lines(smooth.spline(r.ndep.am[,4] - r.ndep.am[,5] ~ r.range), lty=2, lwd=0.5)

#em line
lines(smooth.spline(r.ndep.em[,4] ~ r.range) , lwd=2, col=cols[5])
polygon(c(r.range, rev(r.range)),c(r.ndep.em[,4] + r.ndep.em[,5], rev(r.ndep.em[,4] - r.ndep.em[,5])), col=adjustcolor(cols[5], trans), lty=0)
lines(smooth.spline(r.ndep.em[,4] + r.ndep.em[,5] ~ r.range), lty=2, lwd=0.5)
lines(smooth.spline(r.ndep.em[,4] - r.ndep.em[,5] ~ r.range), lty=2, lwd=0.5)

#labels
mtext('d.', adj=0.05,line = -1.7, side=1)
#mtext(side=2,'Recruitment',cex.lab = cex.olab, las=0, line = 6)
mtext(side=2,expression(paste('stems yr'^'-1')), las=0, line = 2)
mtext('Recruitment', line = -1.2, cex=cex.ilab)
mtext(expression(paste('kg N ha'^'-2',' yr'^'-1')), side=1, line = 2.5)


#mortality plot- AM VS EM - effects not significant.
#limy <- c(0,0.07)
#am line
#plot(m.ndep.am[,4] ~ m.range, lwd=0, ylim= limy, cex.axis = cex.olab)
#lines(smooth.spline(m.ndep.am[,4] ~ m.range) , lwd=2, col=cols[1], lty=2)
#polygon(c(m.range, rev(m.range)),c(m.ndep.am[,3], rev(m.ndep.am[,1])), col=adjustcolor(cols[1], trans), lty=0)
#lines(smooth.spline(m.ndep.am[,1] ~ m.range), lty=2, lwd=0.5)
#lines(smooth.spline(m.ndep.am[,3] ~ m.range), lty=2, lwd=0.5)
#em line
#lines(smooth.spline(m.ndep.em[,4] ~ m.range) , lwd=2, col=cols[5], lty=2)
#polygon(c(m.range, rev(m.range)),c(m.ndep.em[,3], rev(m.ndep.em[,1])), col=adjustcolor(cols[5], trans), lty=0)
#lines(smooth.spline(m.ndep.em[,1] ~ m.range), lty=2, lwd=0.5)
#lines(smooth.spline(m.ndep.em[,3] ~ m.range), lty=2, lwd=0.5)

#labels
#mtext('d.', adj=0.05,line = -1.7, side=1)
#mtext(side=2,'Mortality Probability',cex.lab = cex.olab, las=0, line = 6)
#mtext('Mortality Probability', line = -1.2, cex=cex.ilab)
#mtext(expression(paste('kg N ha'^'-2',' yr'^'-1')), side=1, line = 2.5)


#outer labels
mtext('Nitrogen Deposition', side=1, out=T, cex=1.5, line = 3)
mtext('Eastern U.S. only'  , side=3, out=T, cex=1.5, line = -1)

#end plot
dev.off()
