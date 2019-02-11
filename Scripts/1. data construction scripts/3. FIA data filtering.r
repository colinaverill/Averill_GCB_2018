#This script uses the output of "1. FIA Extract and Filter.r". This output represents all data at the tree-level that meet our filtering criteria, and are plots that are represented in the soils data set. 
#This script uses the output of "2. FIA soil data extraction.r". This output is the profile scale aggregated soil data for each plot that met all filtering criteria. 
#This script builds 3 products
#Product_1: The total basal area of all trees in every plot at the time of soil measurement, and then subsets this basal area by mycorrhizal type and PFT, pairs with soils.
#This is for the soil C storage analysis, and the relative abundance of AM and EM trees analysis.
#Product_2: This is the tree level data for mortality analysis, paired with soils data. Uses the 'future' FIA plot data, generates a death vector, and identifies the mycorrhizal status of each tree.
#Product_3: Plot level growth and recruitment data- plot level basal area increment of all trees that survived the remeasurement interval. I do not subtract death here. Also all new trees that cross 5in DIA threshold counted as new recruits.
rm(list=ls())
require(data.table)
source('required_products_utilities/master_path_list.r')

#load data from Tree and Soil queries
FIA.out    <- readRDS(FIA_extraction_out.path)
FIA.future <- readRDS(FIA_extraction_out_FUTURE.path)
Soils      <- read.csv(soil_data.processed.path)
FIA.states <- data.table(read.csv('required_products_utilities/FIA_state_codes_regions.csv'))

#remove quotes from CN values for both FIA data sets because they mess everything up.
FIA.out$PLT_CN      <- as.numeric(gsub('"', "", FIA.out$PLT_CN     ))
FIA.out$PREV_PLT_CN <- as.numeric(gsub('"', "", FIA.out$PREV_PLT_CN))
FIA.out$TRE_CN      <- as.numeric(gsub('"', "", FIA.out$TRE_CN     ))
FIA.out$PREV_TRE_CN <- as.numeric(gsub('"', "", FIA.out$PREV_TRE_CN))
FIA.future$PLT_CN      <- as.numeric(gsub('"', "", FIA.future$PLT_CN     ))
FIA.future$PREV_PLT_CN <- as.numeric(gsub('"', "", FIA.future$PREV_PLT_CN))
FIA.future$TRE_CN      <- as.numeric(gsub('"', "", FIA.future$TRE_CN     ))
FIA.future$PREV_TRE_CN <- as.numeric(gsub('"', "", FIA.future$PREV_TRE_CN))

#remove any plots that have "clear evidence of artificial regeneration." STDORGCD == 1. 
#294 of 4449 soil sites in intial, 303 of 3221 remeasurement plots.
to.remove        <- unique(   FIA.out[STDORGCD ==1,]$PLT_CN)
to.remove.future <- unique(FIA.future[STDORGCD ==1,]$PLT_CN)
FIA.out <- FIA.out[!(PLT_CN %in% to.remove),]
FIA.future <- FIA.future[!(PLT_CN %in% to.remove.future),]

#remove any plots that have a tree with STATUSCD = 3. These are plots where humans cut down a tree.
#125 of 4155 intiial, and 285/2918 future plots.
to.remove        <- unique(   FIA.out[STATUSCD == 3,]$PLT_CN)
to.remove.future <- unique(FIA.future[STATUSCD == 3,]$PLT_CN)
FIA.out    <-    FIA.out[!(PLT_CN %in% to.remove),]
FIA.future <- FIA.future[!(PLT_CN %in% to.remove.future),]

####Grab sapling data from microplot data for sapling recruitment analysis.
FIA.out.sap        <-    FIA.out[(TPA_UNADJ == 74.965282),]
FIA.out.sap.future <- FIA.future[(TPA_UNADJ == 74.965282),]

####Remove all saplings (DIA < 5inches) based on microplot samplings.
#FIA.out   <-     FIA.out[!(TPA_UNADJ == 74.965282),]
#FIA.future <- FIA.future[!(TPA_UNADJ == 74.965282),]

####Remove one random site that has very strange growth/recruitment numbers.
####Fairly confident this is recovering from a recent clearcut, but was not indicated in other filters.
FIA.out    <-    FIA.out[!(PLT_CN == 65355954010538),]
FIA.future <- FIA.future[!(PREV_PLT_CN == 65355954010538),]

###############################################################
#####Product 1. Basal area of each plot paired with soils######
###############################################################

#Calculate number of species in each plot.
FIA.out[, spp.count := uniqueN(SPCD), by = PLT_CN]

#generate lists of myc types and PFTs.
em.list <- levels(FIA.out$MYCO_ASSO)
pft.list <- levels(FIA.out$PFT)
all.list <- c(em.list,pft.list)

#convert PREVDIA to numeric, because its not for some reason. 
FIA.out[,PREVDIA := as.numeric(PREVDIA)]

#if you are currently dead, your current basal area is assigned NA
FIA.out$BASAL <- ifelse(!is.na(FIA.out$AGENTCD),NA,FIA.out$BASAL)

#current and previous basal area in cm2
FIA.out[,BASAL     := pi*((2.54*DIA    )/2)^2]
FIA.out[,PREVBASAL := pi*((2.54*PREVDIA)/2)^2]

#calculate current basal area of all trees by mycorrhizal type and PFT
for(v in em.list){
  pre <- paste0("FIA.out[,BASAL.",v,":= (MYCO_ASSO =='",v,"') * BASAL]"); eval(parse(text=pre))}
for(v in pft.list){
  pre <- paste0("FIA.out[,BASAL.",v,":= (PFT==       '",v,"') * BASAL]"); eval(parse(text=pre))}


#begin aggregation of tree-level data to plot level. 
scaled<- aggregate(FIA.out$BASAL ~ FIA.out$PLT_CN, FUN='sum', na.rm=T, na.action = na.pass)
colnames(scaled) <- c("PLT_CN", "BASAL")

#aggregate basal area per plot by myctype and PFT
for(v in all.list){ 
  pre <- paste0("scaled$BASAL.",v," <- aggregate(FIA.out$BASAL.",v," ~ FIA.out$PLT_CN, FUN='sum',na.rm=T, na.action=na.pass)[,2]"); eval(parse(text = pre))}

#pop in relevant data from plot table by taking medians
scaled$latitude            <- aggregate(FIA.out$LAT           ~ FIA.out$PLT_CN, FUN='median',na.rm=T,na.action=na.pass)[,2]
scaled$longitude           <- aggregate(FIA.out$LON           ~ FIA.out$PLT_CN, FUN='median',na.rm=T,na.action=na.pass)[,2]
scaled$elevation           <- aggregate(FIA.out$ELEV          ~ FIA.out$PLT_CN, FUN='median',na.rm=T,na.action=na.pass)[,2]
scaled$INVYR               <- aggregate(FIA.out$INVYR         ~ FIA.out$PLT_CN, FUN='median',na.rm=T,na.action=na.pass)[,2]
scaled$STATECD             <- aggregate(FIA.out$STATECD       ~ FIA.out$PLT_CN, FUN='median',na.rm=T,na.action=na.pass)[,2]
scaled$STDAGE              <- aggregate(FIA.out$STDAGE        ~ FIA.out$PLT_CN, FUN='median',na.rm=T,na.action=na.pass)[,2]

#calculate relative abundance of EM trees at time of soil sampling and total EM + AM
scaled$relEM    <-  scaled$BASAL.ECM                    / scaled$BASAL
scaled$relEM.AM <- (scaled$BASAL.ECM + scaled$BASAL.AM) / scaled$BASAL

#merge with soils data. Remember the tree data has both caharacter and numeric PLT_CN values.
Product_1 <- merge(scaled,Soils,by='PLT_CN')
Product_1 <- data.table(Product_1)

#remove sites that have C:N ratios greater than 90. Removes 197 sites.
Product_1[,cn := C.storage/N.storage]
Product_1 <- Product_1[cn < 90,]

#remove sites that are less than 90% AM or EM by basal area. Removes 643 Sites.
Product_1 <- Product_1[!(relEM.AM < 0.9),]

#save this output file. Plot scale tree data at time of soil measurement merged with soils!
saveRDS(Product_1,file=Product_1.path)


############################################################
#####         Product 2. Mortality data               ######
############################################################

#flag whether a tree died or not, for any reason. If it did, there will be a value associated with "AGENTCD" grater than 0.
FIA.future[,death := ifelse(AGENTCD > 0, 1, 0)]

#Trees that are new (cross 5in threshold) during the remeasurement period don't count as surviving the measurement period. give them NA values.
FIA.future[,death := ifelse(PREV_TRE_CN > 0, death, NA)]

#get diameter in centimeters
FIA.future[,DIA.cm := DIA*2.54]

#merge this data with soils. Need to use previous PLT_CN values. Must also merge in relative abundance EM for downstream filtering. 
Product_2 <- merge(FIA.future, Product_1[,.(relEM,relEM.AM,PLT_CN)], by.x = 'PREV_PLT_CN', by.y = 'PLT_CN')
Product_2 <- merge(Product_2, Soils, by.x = 'PREV_PLT_CN', by.y = 'PLT_CN')

#exclude sites in Product 2 that are not in Product 1. This does not exclude anything. 
Product_2 <- Product_2[PREV_PLT_CN %in% Product_1$PLT_CN,]

#save the output! Tree-level mortality data paired with soils!
saveRDS(Product_2,file=Product_2.path)

############################################################
#####      Product 3. Growth and Recruitment data     ######
############################################################

#In this section I calculate:
#1. Plot area from TPA_UNADJ column (expansion factors)
#2. calculate current basal area total and by myc/pft
#3. Basal area growth of surviving trees by mycorrhizal type and PFT
#4. Calculate total recruitment and by mycorrhizal type. 
#5. Aggregate individual tree data to the plot scale. Add site level data. 

###########################################################
#1. Plot area from TPA_UNADJ column (expansion factors)
###########################################################
#Note- all microplot observations of saplings have already been removed. 
#assuming they sample everything within the plot at the biggest resolution, the correct TPA_UNADJ is the smallest non-zero, non-NA TPA_UNADJ number.
FIA.future[, TPA.fixed := ifelse(TPA_UNADJ ==  0, NA, TPA_UNADJ)]
FIA.future[, TPA.fixed := ifelse(TPA_UNADJ == "", NA, TPA.fixed)]
FIA.future[, TPA.fixed := min(TPA.fixed, na.rm=T), by = PLT_CN]

#17 sites do not have a TPA_UNADJ values that is not NA. remove them. 
length(unique(FIA.future$PLT_CN))
length(unique(FIA.future[is.na(TPA.fixed)]$PLT_CN))
FIA.future <- FIA.future[!(is.na(TPA.fixed)),]

#double check that everything has a single value of TPA.fixed
FIA.future[,n.expansion  := sum(!is.na(unique(TPA.fixed))), by=PLT_CN]
max(FIA.future$n.expansion) == 1

#convert back to numeric.
FIA.future[,TPA.fixed := as.numeric(TPA.fixed)]

##calculate plot area as 1/TPA_UNADJ, which returns the area in acres. Convert to m2 by multiplying by 4046.86
#Note- some entries will be NA. These represent individual trees that died during the recensus interval. 
FIA.future[,area.m2 := 1/TPA.fixed * 4046.86]

#########################################################################
#2. calculate current diameter and basal area by mycorrhizal type and PFT
#########################################################################

#generate lists of myc types and PFTs.
em.list  <- levels(FIA.future$MYCO_ASSO)
pft.list <- levels(FIA.future$PFT)
all.list <- c(em.list,pft.list)

#convert PREVDIA to numeric because its not.
FIA.future[,PREVDIA := as.numeric(PREVDIA)]

#current and previous basal area in cm2
FIA.future[,BASAL     := pi*((2.54*DIA    )/2)^2]
FIA.future[,PREVBASAL := pi*((2.54*PREVDIA)/2)^2]

#if you are currently dead, your current basal area is assigned NA
FIA.future$BASAL <- ifelse(FIA.future$AGENTCD > 0,NA,FIA.future$BASAL)

#calculate current basal area of all trees by mycorrhizal type and PFT
for(v in em.list){
  pre <- paste0("FIA.future[,BASAL.",v,":= (MYCO_ASSO =='",v,"') * BASAL]"); eval(parse(text=pre))}
for(v in pft.list){
  pre <- paste0("FIA.future[,BASAL.",v,":= (PFT==       '",v,"') * BASAL]"); eval(parse(text=pre))}


###########################################################
#3. calculate basal area increment- total and by myc/PFT
###########################################################
#convert REMPER to numeric because its not
FIA.future[,REMPER:= as.numeric(REMPER)]

#calculate basal area increment for each tree. If else statement required, as trees w/o a previous basal crossed 5in threshold during census interval. 
FIA.future[,GROWTH.total := ifelse(!is.na(PREVBASAL), (BASAL-PREVBASAL)/REMPER, BASAL/REMPER)]
#mark trees that died as NA - this should already be addressed by previous code. 
#FIA.future$GROWTH.total <- ifelse(!is.na(FIA.future$AGENTCD),NA,FIA.future$GROWTH.total)

#Repeat process for every myc type and pft
#loop to get growth for each myc type. 
for(v in em.list){
  pre <- paste0("FIA.future[,GROWTH.",v," := ifelse(MYCO_ASSO == '",v,"',ifelse(!is.na(PREVBASAL), (BASAL-PREVBASAL)/REMPER, BASAL/REMPER), NA)]")
  eval(parse(text = pre))
}


###########################################################
#4. calculate recruitment, by MYC
###########################################################

#get a recruitment vector. Anything that does not have a PREV_TRE_CN value.
#all recruits, regardless of adult or sapling recruit.
FIA.future[, recruit := ifelse(is.na(PREV_TRE_CN), 1, 0)]
#adult recruits are not from microplots
FIA.future[, recruit.adult := ifelse(is.na(PREV_TRE_CN),
                                     ifelse(!(TPA_UNADJ==74.965282),
                                            1,0), 
                                     0)
           ]
#sapling recruits are from microplots.
#FIA.future[, recruit.sapling := ifelse(is.na(PREV_TRE_CN),
#                                     ifelse((TPA_UNADJ==74.965282),
#                                            1,0), 
#                                     0)
#           ]

#get EM and AM recruitment vectors. 
FIA.future[,recruit.em         := ifelse(MYCO_ASSO == 'ECM',recruit        , 0)]
FIA.future[,recruit.am         := ifelse(MYCO_ASSO ==  'AM',recruit        , 0)]
FIA.future[,recruit.adult.em   := ifelse(MYCO_ASSO == 'ECM',recruit.adult  , 0)]
FIA.future[,recruit.adult.am   := ifelse(MYCO_ASSO ==  'AM',recruit.adult  , 0)]
#FIA.future[,recruit.sapling.em := ifelse(MYCO_ASSO == 'ECM',recruit.sapling, 0)]
#FIA.future[,recruit.sapling.am := ifelse(MYCO_ASSO ==  'AM',recruit.sapling, 0)]

########################################################################
#5. Begin aggregating tree-level data by PFT and myc type to plot-level. 
########################################################################

scaled2 <- aggregate(FIA.future$BASAL ~ FIA.future$PLT_CN, FUN='sum', na.rm=T, na.action = na.pass)
colnames(scaled2) <- c("PLT_CN", "BASAL")

#grab total growth and mortality
scaled2$GROWTH.total <- aggregate(FIA.future$GROWTH.total ~ FIA.future$PLT_CN, FUN='sum', na.rm=T, na.action = na.pass)[,2]
scaled2$PREVBASAL    <- aggregate(FIA.future$PREVBASAL    ~ FIA.future$PLT_CN, FUN='sum', na.rm=T, na.action = na.pass)[,2]

#grab total recruitment
scaled2$recruit              <- aggregate(FIA.future$recruit              ~ FIA.future$PLT_CN, FUN='sum', na.rm=T, na.action = na.pass)[,2]
scaled2$recruit.adult        <- aggregate(FIA.future$recruit.adult        ~ FIA.future$PLT_CN, FUN='sum', na.rm=T, na.action = na.pass)[,2]
#scaled2$recruit.sapling      <- aggregate(FIA.future$recruit.sapling      ~ FIA.future$PLT_CN, FUN='sum', na.rm=T, na.action = na.pass)[,2]
scaled2$recruit.em           <- aggregate(FIA.future$recruit.em           ~ FIA.future$PLT_CN, FUN='sum', na.rm=T, na.action = na.pass)[,2]
scaled2$recruit.am           <- aggregate(FIA.future$recruit.am           ~ FIA.future$PLT_CN, FUN='sum', na.rm=T, na.action = na.pass)[,2]
scaled2$recruit.adult.am     <- aggregate(FIA.future$recruit.adult.am     ~ FIA.future$PLT_CN, FUN='sum', na.rm=T, na.action = na.pass)[,2]
scaled2$recruit.adult.em     <- aggregate(FIA.future$recruit.adult.em     ~ FIA.future$PLT_CN, FUN='sum', na.rm=T, na.action = na.pass)[,2]
#scaled2$recruit.sapling.am   <- aggregate(FIA.future$recruit.sapling.am   ~ FIA.future$PLT_CN, FUN='sum', na.rm=T, na.action = na.pass)[,2]
#scaled2$recruit.sapling.em   <- aggregate(FIA.future$recruit.sapling.em   ~ FIA.future$PLT_CN, FUN='sum', na.rm=T, na.action = na.pass)[,2]

#account for the fact that plots with no previous measurement will now have 0s for Growth.total and PREVBASAL
#go through and changes these 0s to NAs
scaled2$GROWTH.total <- ifelse(scaled2$PREVBASAL == 0, NA, scaled2$GROWTH.total)
scaled2$PREVBASAL    <- ifelse(scaled2$PREVBASAL == 0, NA, scaled2$PREVBASAL)

#basal area per plot by myctype and PFT
for(v in all.list){ 
  pre <- paste0("scaled2$BASAL.",v," <- aggregate(FIA.future$BASAL.",v," ~ FIA.future$PLT_CN, FUN='sum',na.rm=T, na.action=na.pass)[,2]"); eval(parse(text = pre))}

#pop in other details by taking the median
scaled2$area.m2             <- aggregate(FIA.future$area.m2     ~ FIA.future$PLT_CN, FUN='median',na.rm=T,na.action=na.pass)[,2]
scaled2$REMPER              <- aggregate(FIA.future$REMPER      ~ FIA.future$PLT_CN, FUN='median',na.rm=T,na.action=na.pass)[,2]
#scaled2$TPA.fixed           <- aggregate(FIA.future$TPA_UNADJ   ~ FIA.future$PLT_CN, FUN='median',na.rm=T,na.action=na.pass)[,2]
scaled2$latitude            <- aggregate(FIA.future$LAT         ~ FIA.future$PLT_CN, FUN='median',na.rm=T,na.action=na.pass)[,2]
scaled2$longitude           <- aggregate(FIA.future$LON         ~ FIA.future$PLT_CN, FUN='median',na.rm=T,na.action=na.pass)[,2]
scaled2$elevation           <- aggregate(FIA.future$ELEV        ~ FIA.future$PLT_CN, FUN='median',na.rm=T,na.action=na.pass)[,2]
scaled2$INVYR               <- aggregate(FIA.future$INVYR       ~ FIA.future$PLT_CN, FUN='median',na.rm=T,na.action=na.pass)[,2]
scaled2$PREV_PLT_CN         <- aggregate(FIA.future$PREV_PLT_CN ~ FIA.future$PLT_CN, FUN='median',na.rm=T,na.action=na.pass)[,2]
scaled2$STATECD             <- aggregate(FIA.future$STATECD     ~ FIA.future$PLT_CN, FUN='median',na.rm=T,na.action=na.pass)[,2]

#account for the fact that plots with no REMPER observation will now have 0s from Growth.total, because it was summed w/ a na.rm=T flag
scaled2$GROWTH.total <- ifelse(is.na(scaled2$REMPER), NA, scaled2$GROWTH.total)

#merge with soils- note, soils match with the "future" scaled 2 plots 'PREV_PLT_CN'. Also grab initial relative abundance of AM+EM basal area from Product_1 for filtering
#note, merging with previous Product 1 sites drops 830 sites. 
#but all of them have soils data. Must be prior filtering criteria above, require > 90% AM.EM
scaled2   <- merge(scaled2, Product_1[,.(PLT_CN,relEM.AM)]    ,by.x='PREV_PLT_CN', by.y ="PLT_CN")
Product_3 <- merge(scaled2,Soils                              ,by.x='PREV_PLT_CN', by.y ="PLT_CN")

#double check that Product_3 only includes sites present in Product_1
Product_3 <- Product_3[Product_3$PREV_PLT_CN %in% Product_1$PLT_CN,]
Product_3 <- data.table(Product_3)

#save this output file. Gross basal increment of surviving trees at the plot scale, paired with soils!
saveRDS(Product_3,file=Product_3.path)


#####Make Eastern US subset using FIA STATECD column. 
Product_1.E <- Product_1[STATECD %in% FIA.states[east==1,STATECD],]
Product_2.E <- Product_2[STATECD %in% FIA.states[east==1,STATECD],]
Product_3.E <- Product_3[STATECD %in% FIA.states[east==1,STATECD],]

#Make Eastern US subset using FIA STATECD column, also excluding the midwest.
Product_1.E_nomid <- Product_1[STATECD %in% FIA.states[east==1 & midwest == 0,STATECD],]
Product_2.E_nomid <- Product_2[STATECD %in% FIA.states[east==1 & midwest == 0,STATECD],]
Product_3.E_nomid <- Product_3[STATECD %in% FIA.states[east==1 & midwest == 0,STATECD],]

#save output.
saveRDS(Product_1.E, Product_1.E.path)
saveRDS(Product_2.E, Product_2.E.path)
saveRDS(Product_3.E, Product_3.E.path)

saveRDS(Product_1.E_nomid, Product_1.E_nomid.path)
saveRDS(Product_2.E_nomid, Product_2.E_nomid.path)
saveRDS(Product_3.E_nomid, Product_3.E_nomid.path)

##end script.