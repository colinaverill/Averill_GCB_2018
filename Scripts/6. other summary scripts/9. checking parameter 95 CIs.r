#clear environment, load paths.
rm(list=ls())
source('required_products_utilities/master_path_list.r')
library(runjags)
check_significance <- function(model){
  mod <- summary(model)
  sig_vector <- mod[,1]*mod[,3]
  sig_vector <- ifelse(sig_vector > 0,'yes','no')
  mod <- as.data.frame(mod)
  mod$p_.05 <- sig_vector
  return(mod)
}

#this is just form looking at the models. 
#I should write code in the actual models that stores which parameters are tied to which variables.
beta.names <- c('intercept','mat','map','cn','pH_H2O','tot.15','clay','NA')
grow.names <- c('intercept','basal.m2','mat','map','cn','pH_H2O','relEM','tot.15','relEM:tot.15','clay','g1','g2','g3','NA')
r.AM.names <- c(c('intercept','pH_H2O','cn','mat','map','tot.15','rel.AM','clay'),c('intercept','pH_H2O','cn','mat','map','tot.15','rel.AM','clay','NA'))
r.EM.names <- c(c('intercept','pH_H2O','cn','mat','map','tot.15','rel.EM','clay'),c('intercept','pH_H2O','cn','mat','map','tot.15','rel.EM','clay','NA'))
mort.names <- c('intercept','DIA1','DIA2','mat','map','pH_H2O','cn','tot.15','em','tot.15:em','clay')
soil.names <- c('intercept','mat','map','pH_H2O','tot.15','relEM','N.storage','relEM:N.storage','tot.15:relEM','clay','NA')


#All data models for checking significance
soil <- check_significance(readRDS(local_soils_model_path))
beta <- check_significance(readRDS(local_beta_model_path))
grow <- check_significance(readRDS(local_growth.gaus_model_path))
r.am <- check_significance(readRDS(local_AM_recruitment_model_path))
r.em <- check_significance(readRDS(local_EM_recruitment_model_path))
mort <- check_significance(readRDS(local_mortality_model_path))
soil$names <- soil.names
beta$names <- beta.names
grow$names <- grow.names
r.am$names <- r.AM.names
r.em$names <- r.EM.names
mort$names <- mort.names
all.list <- list(soil,beta,grow,r.am,r.em,mort)
names(all.list) <- c('soil','beta','growth','AM.recruitment','EM.recruitment','mortality')

#EAST data models for checking significance
soil <- check_significance(readRDS(local_soils_model_path.east))
beta <- check_significance(readRDS(local_beta_model_path.east))
grow <- check_significance(readRDS(local_growth.gaus_model_path.east))
r.am <- check_significance(readRDS(local_AM_recruitment_model_path.east))
r.em <- check_significance(readRDS(local_EM_recruitment_model_path.east))
mort <- check_significance(readRDS(local_mortality_model_path.east))
soil$names <- soil.names
beta$names <- beta.names
grow$names <- grow.names
r.am$names <- r.AM.names
r.em$names <- r.EM.names
mort$names <- mort.names
east.list <- list(soil,beta,grow,r.am,r.em,mort)
names(east.list) <- c('soil','beta','growth','AM.recruitment','EM.recruitment','mortality')


#NOMID EAST data models for checking significance
soil <- check_significance(readRDS(local_soils_model_path.nomid))
beta <- check_significance(readRDS(local_beta_model_path.nomid))
grow <- check_significance(readRDS(local_growth.gaus_model_path.nomid))
r.am <- check_significance(readRDS(local_AM_recruitment_model_path.nomid))
r.em <- check_significance(readRDS(local_EM_recruitment_model_path.nomid))
mort <- check_significance(readRDS(local_mortality_model_path.nomid))
soil$names <- soil.names
beta$names <- beta.names
grow$names <- grow.names
r.am$names <- r.AM.names
r.em$names <- r.EM.names
mort$names <- mort.names
nomid.list <- list(soil,beta,grow,r.am,r.em,mort)
names(nomid.list) <- c('soil','beta','growth','AM.recruitment','EM.recruitment','mortality')

#wrap up and save
out <- list(all.list, east.list, nomid.list)
names(out) <- c('all','east','nomid')
saveRDS(out, significance_table_path)