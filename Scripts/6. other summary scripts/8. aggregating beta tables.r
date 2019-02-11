#aggregate beta factors for publication. 
#clear environment, load packages and data
rm(list=ls())

       all <- readRDS('/fs/data3/caverill/FIA7_Ndep_data.processed/all.data_analysis_output/beta_list.rds')
      east <- readRDS('/fs/data3/caverill/FIA7_Ndep_data.processed/east.data_analysis_output/east_beta_list.rds')
nomid_east <- readRDS('/fs/data3/caverill/FIA7_Ndep_data.processed/nomid_east.data_analysis_output/nomid_east_beta_list.rds')

#Soil beta table
soil.table <- t(rbind(all$soil, east$soil, nomid_east$soil))
colnames(soil.table) <- c('all','east','nomid_east')

#relative abundance beta table
rel.table <- t(rbind(all$relative.abundance, east$relative.abundance, nomid_east$relative.abundance))
colnames(rel.table) <- c('all','east','nomid_east')

#growth beta table
growth.table <- t(rbind(all$growth, east$growth, nomid_east$growth))
colnames(growth.table) <- c('all','east','nomid_east')

#AM recruitment beta table
am.recruit.table <- t(rbind(all$recruit.am, east$recruit.am, nomid_east$recruit.am))
colnames(am.recruit.table) <- c('all','east','nomid_east')

#EM recruitment beta table
em.recruit.table <- t(rbind(all$recruit.em, east$recruit.em, nomid_east$recruit.em))
colnames(em.recruit.table) <- c('all','east','nomid_east')

#mortality beta table
mortality.table <- t(rbind(all$mortality, east$mortality, nomid_east$mortality))
colnames(mortality.table) <- c('all','east','nomid_east')
