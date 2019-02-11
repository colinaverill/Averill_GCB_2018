#summary parameter values and p < 0.05 'significance' from JAGS models.
#For ALL data subset.
#in future summaries will be save with predictor or interaction names.
#clear environment, load packages.
rm(list=ls())
library(runjags)
library(data.table)
source('required_products_utilities/master_path_list.r')
p.val_fun <- function(x) {ifelse(x$Lower95*x$Upper95 < 0, 'N.S.','p < 0.05')}

#paths
      beta <- readRDS(local_beta_summary_path.nomid)
    growth <- readRDS(local_growth.gaus_summary_path.nomid)
recruit.AM <- readRDS(local_AM_recruitment_summary_path.nomid)
recruit.EM <- readRDS(local_EM_recruitment_summary_path.nomid)
      mort <- readRDS(local_mortality_summary_path.nomid)
      soil <- readRDS(local_soils_summary_path.nomid)
  out.path <- nomid_model.table_path
outCI.path <- nomid_model.table_CI_path
  
  #beta model summary.
  beta <- data.frame(beta)
  beta <- beta[grep('a',rownames(beta)),]
  beta$p.val <- p.val_fun(beta)
  beta <- beta[1:nrow(beta) -1,]
  #add names
  beta$predictor <- c('intercept','MAT','MAP','C:N','pH','N depositon','clay')
  beta <- data.table(beta)
  betaCI <- beta[,.(predictor,Mean,Lower95,Upper95)]
  beta <- beta[,.(predictor,Mean,SD,p.val)]
  
  #growth model summary
  growth <- data.table(growth)
  growth <- growth[1:13,]
  growth$p.val <- p.val_fun(growth)
  growth$predictor <- c('intercept','log(basal area)','MAT','MAP','C:N','pH','relEM','N deposition',
                        'relEM : N deposition','clay',
                        'gaus1 - N deposition',
                        'gaus2 - N deposition',
                        'gaus3 - N deposition')
  growthCI <- growth[,.(predictor,Mean,Lower95,Upper95)]
  growth <- growth[,.(predictor,Mean,SD,p.val)]
  
  #recruitment AM model summary
  recruit.AM <- data.table(recruit.AM)
  recruit.AM <- recruit.AM[1:16,]
  recruit.AM$p.val <- p.val_fun(recruit.AM)
  recruit.AM$predictor <- c('intercept_pois','pH_pois','C:N_pois','MAT_pois','MAP_pois','N deposition_pois','relAM_pois','clay_pois',
                            'intercept_bern','pH_bern','C:N_bern','MAT_bern','MAP_bern','N deposition_bern','relAM_bern','clay_bern')
  recruit.AMCI <- recruit.AM[,.(predictor,Mean,Lower95,Upper95)]
  recruit.AM <- recruit.AM[,.(predictor,Mean,SD,p.val)]
  
  #recruitment EM model summary
  recruit.EM <- data.table(recruit.EM)
  recruit.EM <- recruit.EM[1:16,]
  recruit.EM$p.val <- p.val_fun(recruit.EM)
  recruit.EM$predictor <- c('intercept_pois','pH_pois','C:N_pois','MAT_pois','MAP_pois','N deposition_pois','relEM_pois','clay_pois',
                            'intercept_bern','pH_bern','C:N_bern','MAT_bern','MAP_bern','N deposition_bern','relEM_bern','clay_bern')
  recruit.EMCI <- recruit.EM[,.(predictor,Mean,Lower95,Upper95)]
  recruit.EM <- recruit.EM[,.(predictor,Mean,SD,p.val)]
  
  #mortality model summary
  mort <- data.table(mort)
  mort <- mort[1:11,]
  mort$p.val <- p.val_fun(mort)
  mort$predictor <- c('intercept','exp1_diameter','exp2_diameter','MAT','MAP','pH','C:N','N deposition','EM','N deposition : EM','clay')
  mortCI <- mort[,.(predictor,Mean,Lower95,Upper95)]
  mort <- mort[,.(predictor,Mean,SD,p.val)]
  
  #soils model summary
  soil <- data.table(soil)
  soil <- soil[1:10,]
  soil$p.val <- p.val_fun(soil)
  soil$predictor <- c('intercept','MAT','MAP','pH','N deposition','relEM','log(N.storage)','relEM : log(N.storage)','relEM : N deposition','clay')
  soilCI <- soil[,.(predictor,Mean,Lower95,Upper95)]
  soil <- soil[,.(predictor,Mean,SD,p.val)]
  
  #combine summaries as list, name and save.
  out <- list(beta,growth,recruit.AM,recruit.EM,mort,soil)
  outCI <- list(betaCI,growthCI,recruit.AMCI,recruit.EMCI,mortCI,soilCI)
  names(out) <- c('beta','growth','recruit.AM','recruit.EM','mortality','soil')
  names(outCI) <- c('beta','growth','recruit.AM','recruit.EM','mortality','soil')
  saveRDS(out,out.path)
  saveRDS(outCI,outCI.path)