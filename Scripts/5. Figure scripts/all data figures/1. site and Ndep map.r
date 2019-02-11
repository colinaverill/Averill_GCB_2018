#mapping all data and demographic data sites using thre ndep15 data.
rm(list=ls())
require(grid)
require(ggplot2)
require(raster)
require(rgdal)
require(RColorBrewer)
source('required_products_utilities/master_path_list.r')
alld <- readRDS(local_beta_filtered_path)
alld$col <- alld$relEM * 100
demo <- readRDS(local_growth_filtered_path)
no.demo <- alld[!(alld$PLT_CN %in% demo$PREV_PLT_CN),]
yes.demo <- alld[(alld$PLT_CN %in% demo$PREV_PLT_CN),]


#load ndep15 (wet+dry 2000-2014) raster 
r.all <- raster(ndep_summary_raster_path)
r.all <- log(r.all)

#put raster in CRS of points
r.proj = projectRaster(r.all, crs=crs('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'))
raster.points <- rasterToPoints(r.proj)
raster.points <- data.frame(raster.points)
colnames(raster.points)<-c('x','y','layer')

#grab points from data set
points <- cbind(alld$longitude, alld$latitude)

#these are colors I like from the wes anderson color pack!
require(wesanderson)
colors.w <- c(wes_palette(5,name="Zissou")[1],wes_palette(5,name="Zissou")[3],wes_palette(5,name="Zissou")[5],'#FF69B4','#FFFFFF')
colors <- wes_palette("Zissou", 21, type='continuous')
colors <- brewer.pal(9,"Blues")
colors.w <- brewer.pal(9,'GnBu')


png(filename=Figure1_all_path,width=10,height=6.5,units='in',res=300)

#map 1 - all sites
mp <- NULL

#latitude <- latitude1
#longitude <- longitude1

#grab world map and choose colors
mapWorld <- borders("usa", colour='white',fill='white', lwd=0.4)
mp <- ggplot(data=raster.points, aes(y=y, x=x))
mp <- mp + mapWorld  
mp <- mp + geom_raster(aes(fill=layer)) 
mp <- mp + scale_fill_gradientn("N deposition",colours =colors.w) 
mp <- mp + theme(legend.position="bottom", legend.box='horizontal')
mp <- mp + geom_point(data= no.demo,aes(x=longitude, y=latitude,color= no.demo$col) , size=1,pch=4) 
mp <- mp + geom_point(data=yes.demo,aes(x=longitude, y=latitude,color=yes.demo$col) , size=1,pch=16)
mp <- mp + scale_colour_gradient2('%EM',low='purple', mid = 'green', high='red', midpoint = 50)

#change up background gridlines, labels and colors.
mp<- mp + theme(axis.text.y=element_blank(),
                axis.text.x=element_blank(),
                axis.title.y=element_blank(),
                axis.title.x=element_blank(),
                axis.ticks=element_blank(), 
                panel.background = element_rect(fill='white'),
                plot.background = element_rect(fill='white'),
                panel.grid.major=element_blank(),
                panel.grid.minor=element_blank())
mp

dev.off()