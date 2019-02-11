#which trees were most abundant at low and high N-deposition?
#clear environment, load packakges
rm(list=ls())
library(data.table)

#load all rel abundance data that went into the fit.
d <- data.table(readRDS('/fs/data3/caverill/FIA7_Ndep_data.processed/all.data_analysis_output/rel.abundance.filtered_FIA7_051517.rds'))
#load species codes to FIA codes - in FIA traits folder.
codes <- data.table(read.csv('/home/caverill/FIA_trait_analyses/raw_data/FIA_codes_all.csv'))

#grab low and high ndep plots - 200 lowest and highest N-dep plots, respectively
d <- d[order(d$tot.15),]
lo <- d[1:200,]
hi <- d[(nrow(d) - 199):nrow(d),]

#now we have to load up the thing used to calculate basal area in the first place.
tree <- readRDS('/fs/data3/caverill/FIA7_Ndep_data.processed/soilC.FIA.out.rds')
tree$PLT_CN      <- as.numeric(gsub('"', "", tree$PLT_CN     ))
tree$PREV_PLT_CN <- as.numeric(gsub('"', "", tree$PREV_PLT_CN))
tree$TRE_CN      <- as.numeric(gsub('"', "", tree$TRE_CN     ))
tree$PREV_TRE_CN <- as.numeric(gsub('"', "", tree$PREV_TRE_CN))

#convert diameter in incehs to basal area per meter squared.
tree$basal.cm2 <- pi*((tree$DIA*2.54)/2)^2

#check all your plots are in here. They are.
sum(lo$PLT_CN %in% tree$PLT_CN)
sum(hi$PLT_CN %in% tree$PLT_CN)

#subset tree table by lo and hi observations
tree.lo <- tree[PLT_CN %in% lo$PLT_CN,]
tree.hi <- tree[PLT_CN %in% hi$PLT_CN,]

#aggregate lo and hi tables by species. Get the sum of basal area.
ag.lo <- aggregate(basal.cm2 ~ SPCD, data = tree.lo, FUN='sum')
ag.hi <- aggregate(basal.cm2 ~ SPCD, data = tree.hi, FUN='sum')

#drop in species codes
ag.lo <- merge(ag.lo, codes[,.(FIA_code,genus,species)], by.x = 'SPCD', by.y = 'FIA_code', all.x=T)
ag.hi <- merge(ag.hi, codes[,.(FIA_code,genus,species)], by.x = 'SPCD', by.y = 'FIA_code', all.x=T)

#sort by basal area represented.
ag.lo <- ag.lo[order(-ag.lo$basal.cm2),]
ag.hi <- ag.hi[order(-ag.hi$basal.cm2),]

#do it for all data as well
ag.all <- aggregate(basal.cm2 ~ SPCD, data = tree, FUN='sum')
ag.all <- merge(ag.all, codes[,.(FIA_code,genus,species)], by.x = 'SPCD',by.y = 'FIA_code', all.x=T)
ag.all <- ag.all[order(-ag.all$basal.cm2),]