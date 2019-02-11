#Soil C plots with interaction
rm(list=ls())
source('required_products_utilities/master_path_list.r')
raw   <- readRDS(local_soils_filtered_path.nomid)
d     <- readRDS(local_soils_summary_path.nomid)
bins  <- readRDS(bins_soils_nomid_path)
save.path <- Supp_Figure4_all_path
rsq.table <- readRDS(nomid_binned_rsq_summary_path)

#subset output by predictor for MAT, MAP, C:N and pH
ndep.low  <- d[grep('y.nlow'  , rownames(d)),]; ndep.low <- exp(ndep.low)
ndep.mean <- d[grep('y.nmean' , rownames(d)),]; ndep.mean <- exp(ndep.mean)
ndep.high <- d[grep('y.nhigh' , rownames(d)),]; ndep.high <- exp(ndep.high)

em.range <- seq(0,100, by= 1)
raw$relEM <- raw$relEM*100
bins$relEM <- bins$relEM*100

#begin the plots!
require(wesanderson) #wes anderson of course. 
cols <- wes_palette("Moonrise2",4)
trans <- 0.3

#inner/outer label cex (size)
cex.ilab <- 1.0
cex.olab <- 1.0

#ylim values for soil C plots
limy <- c(1000,8000)
limx <- c(0,100)

#save dimensions, destination
png(filename=save.path,width=7,height=4.5,units='in',res=300)

par(mfrow=c(1,2), oma=c(5,5,1,1), mar=c(0,0,0,0))
#1st panel - effect of EM abundance at low N-deposition
#plot(raw$detrend.lo ~ raw$relEM, ylim = limy, xlim = limx, col = 'gray', cex = 0.2)
#par(new=T)
plot(bins$detrend.lo ~ bins$relEM, pch = 16,ylim=limy, xlim=limx)
arrows(bins$relEM, bins$detrend.lo.lwr, bins$relEM, bins$detrend.lo.upr, length=0.05, angle=90, code=3)
#drop fit on top
lines(smooth.spline(ndep.low[,4] ~ em.range) , lwd=2, col=cols[2])
polygon(c(em.range, rev(em.range)),c(ndep.low[,3], rev(ndep.low[,1])), col=adjustcolor(cols[2], trans), lty=0)
lines(smooth.spline(ndep.low[,3] ~ em.range), lty=2, lwd=0.5)
lines(smooth.spline(ndep.low[,1] ~ em.range), lty=2, lwd=0.5)
#axis labels
mtext(expression(paste('g C m'^'-2')),              side=2, line = 2.5, cex=1.2, las=0)
mtext('a.', side = 1, line = -1.3, adj = 0.05, cex = 1.3)
mtext(expression(paste('2.8 kg N ha'^'-1',' yr'^'-1')), cex = 1.3, line = -1.1, adj = 0.95, side = 1)
r.sq <- rsq.table[analyses == 'soilC.lo',]$r.sq
r.sq <- format(round(r.sq,digits = 2), digits = 2)
mylabel = bquote(italic(R)^2 == .(r.sq))
mtext(mylabel, cex = 1.3, line = -2, adj = 0.05)


#2nd panel - effect of EM abundance at high N-deposition
#plot(raw$detrend.hi ~ raw$relEM, ylim = limy, xlim = limx, col = 'gray', cex = 0.2)
#par(new=T)
plot(bins$detrend.hi ~ bins$relEM, pch = 16,ylim=limy, xlim=limx, ylab=NA, xlab=NA, yaxt='n')
arrows(bins$relEM, bins$detrend.hi.lwr, bins$relEM, bins$detrend.hi.upr, length=0.05, angle=90, code=3)
#drop fit on top
lines(smooth.spline(ndep.high[,4] ~ em.range) , lwd=2, col=cols[1])
polygon(c(em.range, rev(em.range)),c(ndep.high[,3], rev(ndep.high[,1])), col=adjustcolor(cols[1], trans), lty=0)
lines(smooth.spline(ndep.high[,3] ~ em.range), lty=2, lwd=0.5)
lines(smooth.spline(ndep.high[,1] ~ em.range), lty=2, lwd=0.5)
#labels
mtext('b.', line = -1.3, side=1, adj = 0.05, cex = 1.3)
mtext(expression(paste('16.9 kg N ha'^'-1',' yr'^'-1')), cex = 1.3, line = -1.1, adj = 0.95, side = 1)
r.sq <- rsq.table[analyses == 'soilC.hi',]$r.sq
r.sq <- format(round(r.sq,digits = 2), digits = 2)
mylabel = bquote(italic(R)^2 == .(r.sq))
mtext(mylabel, cex = 1.3, line = -2, adj = 0.05)

#outer labels
mtext('% EM by basal area', side=1, out=T, cex=1.5, line = 3)

#end plot
dev.off()