#Querying the FIA database.
#Colin downloaded the most recent PLOT, COND, TREE_GRM_EST, and TREE tables from the FIA on May 12, 2017
#This uses RSQlite to query the respective tables, rather than pSQL as done previously by R. Kelly.
#This script takes a while to run (>20 minutes).
#clear environment, load packages.
rm(list=ls())
library(data.table)
library(RSQLite)
source('required_products_utilities/master_path_list.r')

#Some functions to make life easier
sumNA  = function(x)  sum(x,na.rm=T)
meanNA = function(x) mean(x,na.rm=T)
maxNA  = function(x)  max(x,na.rm=T)
tic = function() assign("timer", Sys.time(), envir=.GlobalEnv)
toc = function() print(Sys.time()-timer)

#Connect to FIA7 database.
con <- dbConnect(SQLite(), dbname = FIAdb.path)

#File paths to data
file.pft = "required_products_utilities/gcbPFT.csv"
file.myc = "required_products_utilities/mycorrhizal_SPCD_data.csv"
file.soil = read.csv(soil_data.processed.path)
file.soil$plt_cn <- file.soil$PLT_CN

#Where to save filtered outputs. 
file.out        = FIA_extraction_out.path
file.out.future = FIA_extraction_out_FUTURE.path


###---Query PLOT table
cat("Query PLOT...\n")
tic()
PLOT = dbGetQuery(con, "select 
                        CN, STATECD, PREV_PLT_CN, REMPER, LAT, LON, ELEV 
                        from 
                        PLOT 
                        where 
                        STATECD <= 56")
PLOT = data.table(PLOT)
toc()
setnames(PLOT,"CN","PLT_CN")
states  = sort(unique(PLOT$STATECD))
n.state = length(states)

# Remove this one miscellaneous plot, per Trevor Andrews
PLOT = PLOT[ PLT_CN!= 134680578010854 ]


###---Query COND table
cat("Query COND...\n")
tic()
COND = dbGetQuery(con, "select 
                        PLT_CN, CONDID, STDORGCD, STDAGE, CONDPROP_UNADJ, AFFORESTATION_CD
                        from 
                        COND")
COND = data.table(COND)
toc()

###---Query SUBP_COND table
cat("Query SUBP_COND...\n")
tic()
SUBP_COND = dbGetQuery(con, "select 
                       CN, PLT_CN, CONDID, SUBP, SUBPCOND_PROP
                       from 
                       SUBP_COND")
SUBP_COND = data.table(SUBP_COND)
toc()

#Calculate forested condition proportion.
forest_prop <- SUBP_COND[CONDID == 1,]
forest_prop <- data.table(aggregate(SUBPCOND_PROP ~ PLT_CN, data = forest_prop, FUN = 'sum'))
forest_prop[,SUBPCOND_PROP := SUBPCOND_PROP / 4]
colnames(forest_prop)[2] <- 'forest_proportion'
#merge into COND table
COND <- merge(COND,forest_prop, all.x = T)


#need to have at least one plot with condition =1 (forested)
COND <- COND[CONDID == 1,]

#Mod for GCB revision 1 for reviewier 1. Only include sites that are completely condition = 1 (forested).
#Exclude any plots that are not condition 1, or are multiple conditions (i.e. forested and not-forested).
#to.remove <- COND[!(CONDID == 1),]
#COND <- COND[!(PLT_CN %in% to.remove$PLT_CN),]

###---merge PLOT and COND tables
PC = merge(COND, PLOT, by="PLT_CN")

#remove PLOT, COND and SUBP_COND tables from memory.
rm(PLOT,COND, SUBP_COND)


#---Query GRM table
#This is to remove sites that have evidence of harvesting.
# --- Query
query = paste("select PLT_CN, INVYR, TPAGROW_UNADJ, DIA_BEGIN, DIA_END, COMPONENT, TRE_CN, REMPER, STATECD from TREE_GRM_ESTN where ESTN_TYPE='\"AL\"'")

cat("Query TREE_GRM_ESTN...\n")
tic()
GRM = dbGetQuery(con, query)
GRM = data.table(GRM)
toc()
nrow(GRM)


# --- Filtering GRM table
cat("Filtering TREE_GRM_ESTN...\n")

# By plot/cond criteria- 378,658 unique sites 
GRM = GRM[ PLT_CN %in% PC$PLT_CN ]

# Assign GRM$START + GRM$CUT and restrict to cut==0, start>0
GRM[, START      := INVYR - REMPER                                  ]
GRM[, REMPER     := NULL                                            ]
GRM[, CUT1TPA    := (COMPONENT=="CUT1") * TPAGROW_UNADJ             ]
GRM[, CUT2TPA    := (COMPONENT=="CUT2") * TPAGROW_UNADJ             ]
GRM[, CUT        := sumNA(CUT2TPA + CUT1TPA), by=PLT_CN             ]

# Assign Reversion/Diversion, and exclude plots with either.
GRM[, DIVERSION1TPA  := (COMPONENT=="DIVERSION1") * TPAGROW_UNADJ   ]
GRM[, DIVERSION2TPA  := (COMPONENT=="DIVERSION2") * TPAGROW_UNADJ   ]
GRM[, REVERSION1TPA  := (COMPONENT=="REVERSION1") * TPAGROW_UNADJ   ]
GRM[, REVERSION2TPA  := (COMPONENT=="REVERSION2") * TPAGROW_UNADJ   ]
GRM[, REDIV          := sumNA(REVERSION2TPA+REVERSION1TPA+DIVERSION2TPA+DIVERSION1TPA), by=PLT_CN]

cat("Calculating TPA and Diameter...\n")
# Compute TPA
GRM[, INGROWTHTPA    := (COMPONENT=="INGROWTH") * TPAGROW_UNADJ     ]
GRM[, MORTALITY1TPA  := (COMPONENT=="MORTALITY1") * TPAGROW_UNADJ   ]
GRM[, MORTALITY2TPA  := (COMPONENT=="MORTALITY2") * TPAGROW_UNADJ   ]
GRM[, MORTALITYTPA   := MORTALITY1TPA + MORTALITY2TPA               ]


#remove plots from PC table with evidence of cutting, reversion or diversion. 
#this retains 4,333/5,257 of the soil profiles.
toBeRemoved<- rbind(GRM[ REDIV>0, ],GRM[CUT>0,])
PC = PC[!(PC$PLT_CN %in% toBeRemoved$PLT_CN),]


#---Query TREE table
#Colin is modidfying this query to only grab sites with soils data. 
#Colin did this because the full query started crashing pecan2. boring.
# 1. colin has removed p2a_grm_flg!=\'N\' from the query. This killed the west coast.
# 2. statuscd=1 restricts query to live trees. We take anything that is alive in the current or previous measurement period.
cat("Query TREE...\n")


#first, grab all PLT_CN values where soil was measured, as well as remeasurements
     PC$PLT_CN_filter <- as.numeric(gsub('"', "", PC$PLT_CN))
PC$PREV_PLT_CN_filter <- as.numeric(gsub('"', "", PC$PREV_PLT_CN))
a <- PC[     PLT_CN_filter %in% file.soil$PLT_CN,     PLT_CN]
b <- PC[PREV_PLT_CN_filter %in% file.soil$PLT_CN,     PLT_CN]
of_interest <- data.frame(c(a,b))
colnames(of_interest)<- c('test')
initial <- data.frame(a)
revisit <- data.frame(b)

#build an empty data.table to store output. 
out <- data.frame(matrix(NA, nrow = 0, ncol = 18))
colnames(out) <- c('cn','prev_tre_cn','plt_cn','invyr','condid','dia','tpa_unadj','carbon_ag','carbon_bg','spcd','stocking','statuscd','prevdia','prev_status_cd','p2a_grm_flg','reconcilecd','agentcd','tpamort_unadj')
out <- data.table(out)
setnames(out, toupper(names(out)))

#RI is 44, CT is 9

#write a for loop because everything is the worst.
#really I should just query ones that match the file.soil PLT_CN vector. But. SQL queries hate me. So I'm doing this.
#querying based on the 'of_interest' PLT_CN values would probably speed this up a ton.
#could also query in parallel using the doParallel package and a foreach loop.
tic()
for(i in 1:length(states)){
  query = paste('select CN, PREV_TRE_CN, PLT_CN, INVYR, CONDID, DIA, TPA_UNADJ, CARBON_AG, CARBON_BG,
              SPCD, STOCKING, STATUSCD, PREVDIA, PREV_STATUS_CD, P2A_GRM_FLG, RECONCILECD, AGENTCD, TPAMORT_UNADJ 
                from TREE 
                WHERE (PREVDIA>5 OR DIA>5) AND (STATUSCD=1 OR PREV_STATUS_CD=1) AND 
                STATECD IN (', paste(states[i],collapse=','), ')')
  pre.tree = as.data.table(dbGetQuery(con, query))
  #only save data that has soil information
  pre.tree <- pre.tree[PLT_CN %in% of_interest$test,]
  out <- rbind(out,pre.tree)
  cat(paste0(i,' states queried.'))
}
toc()
TREE <-out

#how many plots and trees?
length(unique(TREE$PLT_CN))
nrow(TREE)


# --- Filter TREE
cat("Filter TREE ...\n")
# By plot/cond criteria
TREE = TREE[ PLT_CN %in% PC$PLT_CN ]

# CONDID ("Remove edge effects" --TA)
TREE[, CONmax := maxNA(CONDID), by=PLT_CN]

# STATUSCD
# *** RK: Next line looks wrong. It's a sum, not max, despite the name. I did rewrite the line but this is equivalent to what Travis had so keeping for now.
TREE[, STATUSCDmax := sumNA(3*as.integer(STATUSCD==3)), by=PLT_CN]

# RECONCILECD
TREE[is.na(RECONCILECD), RECONCILECD :=0] # Set NA values to 0 (unused)

# Filter
#TREE = TREE[ CONmax==1 & STATUSCDmax!=3 & STATUSCD!=0 & RECONCILECD<=4 ]
TREE = TREE[STATUSCDmax!=3 & STATUSCD!=0 & RECONCILECD<=4 ]

#CALCULATE number of trees in a plot. 
TREE[,n.trees  := length(TPA_UNADJ), by=PLT_CN]

# --- Merge in PFTs and mycorrhizal associations
cat("Merge in PFTs and mycorrhizal associations...\n")
MCDPFT     = as.data.table(read.csv(file.pft, header = T)) 
CA_myctype = as.data.table(read.csv(file.myc, header = T)) #colin loads in mycorrhizal associations
CA_myctype = CA_myctype[,c("SPCD","MYCO_ASSO"),with=F]     #colin loads in mycorrhizal associations
TREE = merge(TREE, MCDPFT    , all.x=T, by = "SPCD")
TREE = merge(TREE, CA_myctype, all.x=T, by = "SPCD")


###subset to only include TREE sites with soil profiles, merge with complementary PC keys. 
setnames(TREE, 'CN','TRE_CN')
PC.soil   <-   PC[PLT_CN %in% initial$a,]
TREE.soil <- TREE[PLT_CN %in% initial$a,]

#grab the set of future remeasurements of all soil FIA plots. 
PC.soil.future <-     PC[PLT_CN %in% revisit$b,]
TREE.soil.future <- TREE[PLT_CN %in% revisit$b,]


#merge current soils-trees-plot data
cat("Final merge...\n")
ALL = merge(TREE.soil, PC.soil, by='PLT_CN')
#kill observations that are actually remeasurements. 
#this happens because some sites have had soils measured multiple times. ~186 sites. 
ALL = ALL[!PLT_CN %in% PREV_PLT_CN,]

#merge future soils-trees-plot data
ALL.future = merge(TREE.soil.future, PC.soil.future, by='PLT_CN')

#save outputs
cat("Save.../n")
tic()
saveRDS(ALL       , file = file.out       )
saveRDS(ALL.future, file = file.out.future)
toc()
###end script.