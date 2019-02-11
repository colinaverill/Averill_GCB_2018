#calculating number of plots and trees used for each analysis and data set.
#clear environment, load packages
rm(list=ls())
library(data.table)

#relative abundance analysis data
abundance <- readRDS('/fs/data3/caverill/FIA7_Ndep_data.processed/all.data_analysis_output/rel.abundance.filtered_FIA7_051517.rds')
#growth data
growth <- readRDS('/fs/data3/caverill/FIA7_Ndep_data.processed/all.data_analysis_output/FIA7_051517_growth.filtered.rds')
#recruitment
recruit <- readRDS('/fs/data3/caverill/FIA7_Ndep_data.processed/all.data_analysis_output/FIA7_051517_AM.adult_linear_recruit.filtered.rds')
#mortality
mortality <- readRDS('/fs/data3/caverill/FIA7_Ndep_data.processed/all.data_analysis_output/FIA7_051517_mortality.filtered.rds')

#load FIA.out that has the full number of trees by plot - and make sure to remove saplings. 
FIA.out <- readRDS('/fs/data3/caverill/FIA7_Ndep_data.processed/soilC.FIA.out.rds')
FIA.future <- readRDS('/fs/data3/caverill/FIA7_Ndep_data.processed/soilC.FIA.out.FUTURE.rds')
#remove quotes from CN values for both FIA data sets because they mess everything up.
FIA.out$PLT_CN      <- as.numeric(gsub('"', "", FIA.out$PLT_CN     ))
FIA.out$PREV_PLT_CN <- as.numeric(gsub('"', "", FIA.out$PREV_PLT_CN))
FIA.out$TRE_CN      <- as.numeric(gsub('"', "", FIA.out$TRE_CN     ))
FIA.out$PREV_TRE_CN <- as.numeric(gsub('"', "", FIA.out$PREV_TRE_CN))
FIA.future$PLT_CN      <- as.numeric(gsub('"', "", FIA.future$PLT_CN     ))
FIA.future$PREV_PLT_CN <- as.numeric(gsub('"', "", FIA.future$PREV_PLT_CN))
FIA.future$TRE_CN      <- as.numeric(gsub('"', "", FIA.future$TRE_CN     ))
FIA.future$PREV_TRE_CN <- as.numeric(gsub('"', "", FIA.future$PREV_TRE_CN))
FIA.future$REMPER      <- as.numeric(gsub('"', "", FIA.future$REMPER     ))

#kill saplings
FIA.out    <-    FIA.out[!(TPA_UNADJ == 74.965282),]
FIA.future <- FIA.future[!(TPA_UNADJ == 74.965282),]

#list of things for number of sites and trees
names <- c('remeasurement','abundance','growth','recruitment','mortality')
plots <- c(length(unique(FIA.future$PLT_CN)),
           length(unique( abundance$PLT_CN)),
           length(unique(    growth$PLT_CN)),
           length(unique(   recruit$PLT_CN)),
           length(unique( mortality$PLT_CN))
           )
#grab number of trees per site
trees <- c(nrow(FIA.future),
           nrow(   FIA.out[PLT_CN %in% abundance$PLT_CN,]),
           nrow(FIA.future[PLT_CN %in%    growth$PLT_CN,]),
           nrow(FIA.future[PLT_CN %in%   recruit$PLT_CN,]),
           nrow(FIA.future[PLT_CN %in% mortality$PLT_CN,])
          )

#number of plots and trees by subset
data.frame(names,plots,trees)

#date range of initial and final visits
c(min(   FIA.out$INVYR), max(   FIA.out[INVYR<9000,]$INVYR))
c(min(FIA.future$INVYR), max(FIA.future$INVYR))

#average remeasurement interval
mean(FIA.future$REMPER, na.rm=T)

#number of plots lost due to missing expansion factor data
17