#build model parameter, sd and p value tables.
#clear environment.
rm(list=ls())
source('required_products_utilities/master_path_list.r')

#load all, east, nomid.
a <- readRDS  (all_model.table_path)
e <- readRDS (east_model.table_path)
n <- readRDS(nomid_model.table_path)

#model tables
beta.table <- cbind(  a$beta, e$beta  [,2:4], n$beta  [,2:4])
grow.table <- cbind(a$growth, e$growth[,2:4], n$growth[,2:4])
am.r.table <- cbind(a$recruit.AM, e$recruit.AM[,2:4], n$recruit.AM[,2:4])
em.r.table <- cbind(a$recruit.EM, e$recruit.EM[,2:4], n$recruit.EM[,2:4])
mort.table <- cbind(a$mortality, e$mortality[,2:4], n$mortality[,2:4])
soil.table <- cbind(a$soil, e$soil[,2:4], n$soil[,2:4])

#model CI tables
a <- readRDS  (all_model.table_CI_path)
e <- readRDS (east_model.table_CI_path)
n <- readRDS(nomid_model.table_CI_path)

beta.table <- cbind(  a$beta, e$beta  [,2:4], n$beta  [,2:4])
grow.table <- cbind(a$growth, e$growth[,2:4], n$growth[,2:4])
am.r.table <- cbind(a$recruit.AM, e$recruit.AM[,2:4], n$recruit.AM[,2:4])
em.r.table <- cbind(a$recruit.EM, e$recruit.EM[,2:4], n$recruit.EM[,2:4])
mort.table <- cbind(a$mortality, e$mortality[,2:4], n$mortality[,2:4])
soil.table <- cbind(a$soil, e$soil[,2:4], n$soil[,2:4])