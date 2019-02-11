#Paths. Giving this a shot.
#This is the main project data directory. This is where all data products are saved.
#You probably still want to run this on a server, as you're going to be storing a lot of data here (>15GB).
main <- '/fs/data3/caverill/Averill_2018_GCB/'
#where R project is home based.
wd_home <- '/home/caverill/Averill_FIA_GCB_2018/'

#This is the remote project directory for running analyses in parallel on the cluster.
#Colin's server isn't a cluster, which is why he has local and remote directories, as well as transfer scripts.
#You must create this folder on the remote computer, as well as the sub folders "input_data" and "output_data"
remote <- '/project/talbot-lab-data/caverill/Averill_GCB_2018/'

#The below lines will setup the project directory. ONLY DO THIS THE FIRST TIME.
#Otherwise you may delete anything already saved in these directories.
#system(paste0('mkdir ',main))
jags.directory <- paste0(main,'JAGS_output/')
detrend.directory <- paste0(main,'detrend_data')
summary.directory <- paste0(main,'summary_data')
#system(paste0('mkdir ',jags.directory))
#system(paste0('mkdir ',detrend.directory))
#system(paste0('mkdir ',summary.directory))

#Colin would like to write scripts that download these particular products, put in a project folder on local.
#PRISM paths. necessary for spatial climate extraction.
map.PRISM.path <- '/fs/data3/caverill/PRISM/PRISM_ppt_30yr_normal_800mM2_annual_bil/PRISM_ppt_30yr_normal_800mM2_annual_bil.bil'
mat.PRISM.path <- '/fs/data3/caverill/PRISM/PRISM_tmean_30yr_normal_800mM2_annual_bil/PRISM_tmean_30yr_normal_800mM2_annual_bil.bil'

#NACP soil texture products.
NACP.sand.raster.path <- '/fs/data3/caverill/NACP_soil_map/NACP_MSTMIP_UNIFIED_NA_SOIL_MA_1242/data/Unified_NA_Soil_Map_Topsoil_Sand_Fraction.tif'
NACP.silt.raster.path <- '/fs/data3/caverill/NACP_soil_map/NACP_MSTMIP_UNIFIED_NA_SOIL_MA_1242/data/Unified_NA_Soil_Map_Topsoil_Silt_Fraction.tif'
NACP.clay.raster.path <- '/fs/data3/caverill/NACP_soil_map/NACP_MSTMIP_UNIFIED_NA_SOIL_MA_1242/data/Unified_NA_Soil_Map_Topsoil_Clay_Fraction.tif'

#Summary Ndep raster for plotting
#Colin would like to write a script that generates this layer for the user on local.
ndep_summary_raster_path <- '/fs/data3/caverill/CASTNET_Ndep/total_15year_Ndep.gri'

#FIA paths. Requires downloading the FIA database tables and building a sqlite database.
#There is a script in data construction that creates these directories, downloads the FIA data, and builds the sqlite database.
#Also require separately downloading the FIA soils data.
     FIA7.dir.path <- '/fs/data3/caverill/FIA7/'
        FIAdb.path <- paste0(FIA7.dir.path,'FIA7.sqlite')
FIA.soil.chem.path <- paste0(FIA7.dir.path,'soils/SOILS_LAB.csv')
 FIA.soil.loc.path <- paste0(FIA7.dir.path,'soils/SOILS_SAMPLE_LOC.csv')
 
##########################################################################################
 ##################Everything from here on is internal directory paths.###################
 ##################You *should* not have to change any of these. #########################
 #########################################################################################
 
#Data product outputs
      #soils
      soil_data.processed.path <- paste0(main,'FIA7soil_output.csv')
      #initial FIA data extractions, visit 1 and visit 2.
       FIA_extraction_out.path <- paste0(main,'soilC.FIA.out.rds')
FIA_extraction_out_FUTURE.path <- paste0(main,'soilC.FIA.out.FUTURE.rds')

#Products used for final downstream analyses of soil, relative abundance, growth, recruitment and mortality.
#broken out by regional subsets (All U.S., East of Missippi, or East of Missippi w/o midwest)
Product_1.path         <- paste0(main,"Product_1.rds")
Product_2.path         <- paste0(main,"Product_2.rds")
Product_3.path         <- paste0(main,"Product_3.rds")
Product_1.E.path       <- paste0(main,"Product_1.E.rds")
Product_2.E.path       <- paste0(main,"Product_2.E.rds")
Product_3.E.path       <- paste0(main,"Product_3.E.rds")
Product_1.E_nomid.path <- paste0(main,"Product_1.E_nomid.rds")
Product_2.E_nomid.path <- paste0(main,"Product_2.E_nomid.rds")
Product_3.E_nomid.path <- paste0(main,"Product_3.E_nomid.rds")

#Remote data product paths
remote_Product_1.path         <- paste0(remote,"input_data/Product_1.rds")
remote_Product_2.path         <- paste0(remote,"input_data/Product_2.rds")
remote_Product_3.path         <- paste0(remote,"input_data/Product_3.rds")
remote_Product_1.E.path       <- paste0(remote,"input_data/Product_1.E.rds")
remote_Product_2.E.path       <- paste0(remote,"input_data/Product_2.E.rds")
remote_Product_3.E.path       <- paste0(remote,"input_data/Product_3.E.rds")
remote_Product_1.E_nomid.path <- paste0(remote,"input_data/Product_1.E_nomid.rds")
remote_Product_2.E_nomid.path <- paste0(remote,"input_data/Product_2.E_nomid.rds")
remote_Product_3.E_nomid.path <- paste0(remote,"input_data/Product_3.E_nomid.rds")

####################################
###remote output paths. ALL data.###
####################################

#significance table path
significance_table_path <- paste0(main,'significance_table.rds')

output.dir <- paste0(remote,'output_data/')
  #model paths
          remote_beta_model_path <- paste0(output.dir,'beta_model_FIA_all.rds')
        remote_growth_model_path <- paste0(output.dir,'growth_model_FIA_all.rds')
   remote_growth.gaus_model_path <- paste0(output.dir,'growth.gaus_model_FIA_all.rds')
remote_AM_recruitment_model_path <- paste0(output.dir,'AM_recruitment_model_FIA_all.rds')
remote_EM_recruitment_model_path <- paste0(output.dir,'EM_recruitment_model_FIA_all.rds')
     remote_mortality_model_path <- paste0(output.dir,'mortality_model_FIA_all.rds')
         remote_soils_model_path <- paste0(output.dir,'soils_model_FIA_all.rds')
      remote_soils.ff_model_path <- paste0(output.dir,'soils.ff_model_FIA_all.rds')
         
  #filtered data paths
          remote_beta_filtered_path <- paste0(output.dir,'beta_filtered_FIA_all.rds')
        remote_growth_filtered_path <- paste0(output.dir,'growth_filtered_FIA_all.rds')
   remote_growth.gaus_filtered_path <- paste0(output.dir,'growth.gaus_filtered_FIA_all.rds')
remote_AM_recruitment_filtered_path <- paste0(output.dir,'AM_recruitment_filtered_FIA_all.rds')
remote_EM_recruitment_filtered_path <- paste0(output.dir,'EM_recruitment_filtered_FIA_all.rds')  
     remote_mortality_filtered_path <- paste0(output.dir,'mortality_filtered_FIA_all.rds')
         remote_soils_filtered_path <- paste0(output.dir,'soils_filtered_FIA_all.rds')
      remote_soils.ff_filtered_path <- paste0(output.dir,'soils.ff_filtered_FIA_all.rds')
         
  #data predictor range paths
          remote_beta_range_path <- paste0(output.dir,'beta_range_FIA_all.rds')
        remote_growth_range_path <- paste0(output.dir,'growth_range_FIA_all.rds')
   remote_growth.gaus_range_path <- paste0(output.dir,'growth.gaus_range_FIA_all.rds')
remote_AM_recruitment_range_path <- paste0(output.dir,'AM_recruitment_range_FIA_all.rds')
remote_EM_recruitment_range_path <- paste0(output.dir,'EM_recruitment_range_FIA_all.rds')
     remote_mortality_range_path <- paste0(output.dir,'mortality_range_FIA_all.rds')
         remote_soils_range_path <- paste0(output.dir,'soils_range_FIA_all.rds')
      remote_soils.ff_range_path <- paste0(output.dir,'soils.ff_range_FIA_all.rds')
         
  #model summary paths
          remote_beta_summary_path <- paste0(output.dir,'beta_summary_FIA_all.rds')
        remote_growth_summary_path <- paste0(output.dir,'growth_summary_FIA_all.rds')
   remote_growth.gaus_summary_path <- paste0(output.dir,'growth.gaus_summary_FIA_all.rds')
remote_AM_recruitment_summary_path <- paste0(output.dir,'AM_recruitment_summary_FIA_all.rds')
remote_EM_recruitment_summary_path <- paste0(output.dir,'EM_recruitment_summary_FIA_all.rds')
     remote_mortality_summary_path <- paste0(output.dir,'mortality_summary_FIA_all.rds')
         remote_soils_summary_path <- paste0(output.dir,'soils_summary_FIA_all.rds')  
      remote_soils.ff_summary_path <- paste0(output.dir,'soils.ff_summary_FIA_all.rds')  
         
  #model prediction paths
          remote_beta_prediction_path <- paste0(output.dir,'beta_prediction_FIA_all.rds')
        remote_growth_prediction_path <- paste0(output.dir,'growth_prediction_FIA_all.rds')
   remote_growth.gaus_prediction_path <- paste0(output.dir,'growth.gaus_prediction_FIA_all.rds')
remote_AM_recruitment_prediction_path <- paste0(output.dir,'AM_recruitment_prediction_FIA_all.rds')
remote_EM_recruitment_prediction_path <- paste0(output.dir,'EM_recruitment_prediction_FIA_all.rds')
     remote_mortality_prediction_path <- paste0(output.dir,'mortality_prediction_FIA_all.rds')
         remote_soils_prediction_path <- paste0(output.dir,'soils_prediction_FIA_all.rds')
      remote_soils.ff_prediction_path <- paste0(output.dir,'soils.ff_prediction_FIA_all.rds')
         
      
#repeat all these for the local paths.
output.dir <- paste0(main,'JAGS_output/')
  #model paths
          local_beta_model_path <- paste0(output.dir,'beta_model_FIA_all.rds')
        local_growth_model_path <- paste0(output.dir,'growth_model_FIA_all.rds')
   local_growth.gaus_model_path <- paste0(output.dir,'growth.gaus_model_FIA_all.rds')
local_AM_recruitment_model_path <- paste0(output.dir,'AM_recruitment_model_FIA_all.rds')
local_EM_recruitment_model_path <- paste0(output.dir,'EM_recruitment_model_FIA_all.rds')
     local_mortality_model_path <- paste0(output.dir,'mortality_model_FIA_all.rds')
         local_soils_model_path <- paste0(output.dir,'soils_model_FIA_all.rds')
      local_soils.ff_model_path <- paste0(output.dir,'soils.ff_model_FIA_all.rds')
         #filtered data paths
          local_beta_filtered_path <- paste0(output.dir,'beta_filtered_FIA_all.rds')
        local_growth_filtered_path <- paste0(output.dir,'growth_filtered_FIA_all.rds')
   local_growth.gaus_filtered_path <- paste0(output.dir,'growth.gaus_filtered_FIA_all.rds')
local_AM_recruitment_filtered_path <- paste0(output.dir,'AM_recruitment_filtered_FIA_all.rds')
local_EM_recruitment_filtered_path <- paste0(output.dir,'EM_recruitment_filtered_FIA_all.rds')  
     local_mortality_filtered_path <- paste0(output.dir,'mortality_filtered_FIA_all.rds')
         local_soils_filtered_path <- paste0(output.dir,'soils_filtered_FIA_all.rds')
      local_soils.ff_filtered_path <- paste0(output.dir,'soils.ff_filtered_FIA_all.rds')
      #data predictor range paths
          local_beta_range_path <- paste0(output.dir,'beta_range_FIA_all.rds')
        local_growth_range_path <- paste0(output.dir,'growth_range_FIA_all.rds')
   local_growth.gaus_range_path <- paste0(output.dir,'growth.gaus_range_FIA_all.rds')
local_AM_recruitment_range_path <- paste0(output.dir,'AM_recruitment_range_FIA_all.rds')
local_EM_recruitment_range_path <- paste0(output.dir,'EM_recruitment_range_FIA_all.rds')
     local_mortality_range_path <- paste0(output.dir,'mortality_range_FIA_all.rds')
         local_soils_range_path <- paste0(output.dir,'soils_range_FIA_all.rds')
      local_soils.ff_range_path <- paste0(output.dir,'soils.ff_range_FIA_all.rds')
         #model summary paths
          local_beta_summary_path <- paste0(output.dir,'beta_summary_FIA_all.rds')
        local_growth_summary_path <- paste0(output.dir,'growth_summary_FIA_all.rds')
   local_growth.gaus_summary_path <- paste0(output.dir,'growth.gaus_summary_FIA_all.rds')
local_AM_recruitment_summary_path <- paste0(output.dir,'AM_recruitment_summary_FIA_all.rds')
local_EM_recruitment_summary_path <- paste0(output.dir,'EM_recruitment_summary_FIA_all.rds')
     local_mortality_summary_path <- paste0(output.dir,'mortality_summary_FIA_all.rds')
         local_soils_summary_path <- paste0(output.dir,'soils_summary_FIA_all.rds')  
      local_soils.ff_summary_path <- paste0(output.dir,'soils.ff_summary_FIA_all.rds')  
         #model prediction paths
          local_beta_prediction_path <- paste0(output.dir,'beta_prediction_FIA_all.rds')
        local_growth_prediction_path <- paste0(output.dir,'growth_prediction_FIA_all.rds')
   local_growth.gaus_prediction_path <- paste0(output.dir,'growth.gaus_prediction_FIA_all.rds')
local_AM_recruitment_prediction_path <- paste0(output.dir,'AM_recruitment_prediction_FIA_all.rds')
local_EM_recruitment_prediction_path <- paste0(output.dir,'EM_recruitment_prediction_FIA_all.rds')
     local_mortality_prediction_path <- paste0(output.dir,'mortality_prediction_FIA_all.rds')
         local_soils_prediction_path <- paste0(output.dir,'soils_prediction_FIA_all.rds')
      local_soils.ff_prediction_path <- paste0(output.dir,'soils.ff_prediction_FIA_all.rds')
         
#####################################
###remote output paths. EAST data.###
#####################################
output.dir <- paste0(remote,'output_data/')
  #model paths east
          remote_beta_model_path.east <- paste0(output.dir,'beta_model_FIA_east.rds')
        remote_growth_model_path.east <- paste0(output.dir,'growth_model_FIA_east.rds')
   remote_growth.gaus_model_path.east <- paste0(output.dir,'growth.gaus_model_FIA_east.rds')
remote_AM_recruitment_model_path.east <- paste0(output.dir,'AM_recruitment_model_FIA_east.rds')
remote_EM_recruitment_model_path.east <- paste0(output.dir,'EM_recruitment_model_FIA_east.rds')
     remote_mortality_model_path.east <- paste0(output.dir,'mortality_model_FIA_east.rds')
         remote_soils_model_path.east <- paste0(output.dir,'soils_model_FIA_east.rds')
      remote_soils.ff_model_path.east <- paste0(output.dir,'soils.ff_model_FIA_east.rds')
         
         #filtered data path.easts
          remote_beta_filtered_path.east <- paste0(output.dir,'beta_filtered_FIA_east.rds')
        remote_growth_filtered_path.east <- paste0(output.dir,'growth_filtered_FIA_east.rds')
   remote_growth.gaus_filtered_path.east <- paste0(output.dir,'growth.gaus_filtered_FIA_east.rds')
remote_AM_recruitment_filtered_path.east <- paste0(output.dir,'AM_recruitment_filtered_FIA_east.rds')
remote_EM_recruitment_filtered_path.east <- paste0(output.dir,'EM_recruitment_filtered_FIA_east.rds')  
     remote_mortality_filtered_path.east <- paste0(output.dir,'mortality_filtered_FIA_east.rds')
         remote_soils_filtered_path.east <- paste0(output.dir,'soils_filtered_FIA_east.rds')
      remote_soils.ff_filtered_path.east <- paste0(output.dir,'soils.ff_filtered_FIA_east.rds')
         
         #data predictor range path.easts
          remote_beta_range_path.east <- paste0(output.dir,'beta_range_FIA_east.rds')
        remote_growth_range_path.east <- paste0(output.dir,'growth_range_FIA_east.rds')
   remote_growth.gaus_range_path.east <- paste0(output.dir,'growth.gaus_range_FIA_east.rds')
remote_AM_recruitment_range_path.east <- paste0(output.dir,'AM_recruitment_range_FIA_east.rds')
remote_EM_recruitment_range_path.east <- paste0(output.dir,'EM_recruitment_range_FIA_east.rds')
     remote_mortality_range_path.east <- paste0(output.dir,'mortality_range_FIA_east.rds')
         remote_soils_range_path.east <- paste0(output.dir,'soils_range_FIA_east.rds')
      remote_soils.ff_range_path.east <- paste0(output.dir,'soils.ff_range_FIA_east.rds')
         
         #model summary path.easts
          remote_beta_summary_path.east <- paste0(output.dir,'beta_summary_FIA_east.rds')
        remote_growth_summary_path.east <- paste0(output.dir,'growth_summary_FIA_east.rds')
   remote_growth.gaus_summary_path.east <- paste0(output.dir,'growth.gaus_summary_FIA_east.rds')
remote_AM_recruitment_summary_path.east <- paste0(output.dir,'AM_recruitment_summary_FIA_east.rds')
remote_EM_recruitment_summary_path.east <- paste0(output.dir,'EM_recruitment_summary_FIA_east.rds')
     remote_mortality_summary_path.east <- paste0(output.dir,'mortality_summary_FIA_east.rds')
         remote_soils_summary_path.east <- paste0(output.dir,'soils_summary_FIA_east.rds')  
      remote_soils.ff_summary_path.east <- paste0(output.dir,'soils.ff_summary_FIA_east.rds')  
         
         #model prediction path.easts
          remote_beta_prediction_path.east <- paste0(output.dir,'beta_prediction_FIA_east.rds')
        remote_growth_prediction_path.east <- paste0(output.dir,'growth_prediction_FIA_east.rds')
   remote_growth.gaus_prediction_path.east <- paste0(output.dir,'growth.gaus_prediction_FIA_east.rds')
remote_AM_recruitment_prediction_path.east <- paste0(output.dir,'AM_recruitment_prediction_FIA_east.rds')
remote_EM_recruitment_prediction_path.east <- paste0(output.dir,'EM_recruitment_prediction_FIA_east.rds')
     remote_mortality_prediction_path.east <- paste0(output.dir,'mortality_prediction_FIA_east.rds')
         remote_soils_prediction_path.east <- paste0(output.dir,'soils_prediction_FIA_east.rds')
      remote_soils.ff_prediction_path.east <- paste0(output.dir,'soils.ff_prediction_FIA_east.rds')

    #repeat east these for the local path.easts.
output.dir <- paste0(main,'JAGS_output/')
  #model path.easts
          local_beta_model_path.east <- paste0(output.dir,'beta_model_FIA_east.rds')
        local_growth_model_path.east <- paste0(output.dir,'growth_model_FIA_east.rds')
   local_growth.gaus_model_path.east <- paste0(output.dir,'growth.gaus_model_FIA_east.rds')
local_AM_recruitment_model_path.east <- paste0(output.dir,'AM_recruitment_model_FIA_east.rds')
local_EM_recruitment_model_path.east <- paste0(output.dir,'EM_recruitment_model_FIA_east.rds')
     local_mortality_model_path.east <- paste0(output.dir,'mortality_model_FIA_east.rds')
         local_soils_model_path.east <- paste0(output.dir,'soils_model_FIA_east.rds')
      local_soils.ff_model_path.east <- paste0(output.dir,'soils.ff_model_FIA_east.rds')
         
  #filtered data path.easts
          local_beta_filtered_path.east <- paste0(output.dir,'beta_filtered_FIA_east.rds')
        local_growth_filtered_path.east <- paste0(output.dir,'growth_filtered_FIA_east.rds')
   local_growth.gaus_filtered_path.east <- paste0(output.dir,'growth.gaus_filtered_FIA_east.rds')
local_AM_recruitment_filtered_path.east <- paste0(output.dir,'AM_recruitment_filtered_FIA_east.rds')
local_EM_recruitment_filtered_path.east <- paste0(output.dir,'EM_recruitment_filtered_FIA_east.rds')  
     local_mortality_filtered_path.east <- paste0(output.dir,'mortality_filtered_FIA_east.rds')
         local_soils_filtered_path.east <- paste0(output.dir,'soils_filtered_FIA_east.rds')
      local_soils.ff_filtered_path.east <- paste0(output.dir,'soils.ff_filtered_FIA_east.rds')
         
  #data predictor range path.easts
          local_beta_range_path.east <- paste0(output.dir,'beta_range_FIA_east.rds')
        local_growth_range_path.east <- paste0(output.dir,'growth_range_FIA_east.rds')
   local_growth.gaus_range_path.east <- paste0(output.dir,'growth.gaus_range_FIA_east.rds')
local_AM_recruitment_range_path.east <- paste0(output.dir,'AM_recruitment_range_FIA_east.rds')
local_EM_recruitment_range_path.east <- paste0(output.dir,'EM_recruitment_range_FIA_east.rds')
     local_mortality_range_path.east <- paste0(output.dir,'mortality_range_FIA_east.rds')
         local_soils_range_path.east <- paste0(output.dir,'soils_range_FIA_east.rds')
      local_soils.ff_range_path.east <- paste0(output.dir,'soils.ff_range_FIA_east.rds')
         
  #model summary path.easts
          local_beta_summary_path.east <- paste0(output.dir,'beta_summary_FIA_east.rds')
        local_growth_summary_path.east <- paste0(output.dir,'growth_summary_FIA_east.rds')
   local_growth.gaus_summary_path.east <- paste0(output.dir,'growth.gaus_summary_FIA_east.rds')
local_AM_recruitment_summary_path.east <- paste0(output.dir,'AM_recruitment_summary_FIA_east.rds')
local_EM_recruitment_summary_path.east <- paste0(output.dir,'EM_recruitment_summary_FIA_east.rds')
     local_mortality_summary_path.east <- paste0(output.dir,'mortality_summary_FIA_east.rds')
         local_soils_summary_path.east <- paste0(output.dir,'soils_summary_FIA_east.rds')  
      local_soils.ff_summary_path.east <- paste0(output.dir,'soils.ff_summary_FIA_east.rds')  
         
  #model prediction path.easts
          local_beta_prediction_path.east <- paste0(output.dir,'beta_prediction_FIA_east.rds')
        local_growth_prediction_path.east <- paste0(output.dir,'growth_prediction_FIA_east.rds')
   local_growth.gaus_prediction_path.east <- paste0(output.dir,'growth.gaus_prediction_FIA_east.rds')
local_AM_recruitment_prediction_path.east <- paste0(output.dir,'AM_recruitment_prediction_FIA_east.rds')
local_EM_recruitment_prediction_path.east <- paste0(output.dir,'EM_recruitment_prediction_FIA_east.rds')
     local_mortality_prediction_path.east <- paste0(output.dir,'mortality_prediction_FIA_east.rds')
         local_soils_prediction_path.east <- paste0(output.dir,'soils_prediction_FIA_east.rds')
      local_soils.ff_prediction_path.east <- paste0(output.dir,'soils.ff_prediction_FIA_east.rds')
         
######################################
###remote output paths. NOMID data.###
######################################
output.dir <- paste0(remote,'output_data/')
         
  #model paths nomid
          remote_beta_model_path.nomid <- paste0(output.dir,'beta_model_FIA_nomid.rds')
        remote_growth_model_path.nomid <- paste0(output.dir,'growth_model_FIA_nomid.rds')
   remote_growth.gaus_model_path.nomid <- paste0(output.dir,'growth.gaus_model_FIA_nomid.rds')
remote_AM_recruitment_model_path.nomid <- paste0(output.dir,'AM_recruitment_model_FIA_nomid.rds')
remote_EM_recruitment_model_path.nomid <- paste0(output.dir,'EM_recruitment_model_FIA_nomid.rds')
     remote_mortality_model_path.nomid <- paste0(output.dir,'mortality_model_FIA_nomid.rds')
         remote_soils_model_path.nomid <- paste0(output.dir,'soils_model_FIA_nomid.rds')
      remote_soils.ff_model_path.nomid <- paste0(output.dir,'soils.ff_model_FIA_nomid.rds')
         
  #filtered data path.nomids
          remote_beta_filtered_path.nomid <- paste0(output.dir,'beta_filtered_FIA_nomid.rds')
        remote_growth_filtered_path.nomid <- paste0(output.dir,'growth_filtered_FIA_nomid.rds')
   remote_growth.gaus_filtered_path.nomid <- paste0(output.dir,'growth.gaus_filtered_FIA_nomid.rds')
remote_AM_recruitment_filtered_path.nomid <- paste0(output.dir,'AM_recruitment_filtered_FIA_nomid.rds')
remote_EM_recruitment_filtered_path.nomid <- paste0(output.dir,'EM_recruitment_filtered_FIA_nomid.rds')  
     remote_mortality_filtered_path.nomid <- paste0(output.dir,'mortality_filtered_FIA_nomid.rds')
         remote_soils_filtered_path.nomid <- paste0(output.dir,'soils_filtered_FIA_nomid.rds')
      remote_soils.ff_filtered_path.nomid <- paste0(output.dir,'soils.ff_filtered_FIA_nomid.rds')
         
  #data predictor range path.nomids
          remote_beta_range_path.nomid <- paste0(output.dir,'beta_range_FIA_nomid.rds')
        remote_growth_range_path.nomid <- paste0(output.dir,'growth_range_FIA_nomid.rds')
   remote_growth.gaus_range_path.nomid <- paste0(output.dir,'growth.gaus_range_FIA_nomid.rds')
remote_AM_recruitment_range_path.nomid <- paste0(output.dir,'AM_recruitment_range_FIA_nomid.rds')
remote_EM_recruitment_range_path.nomid <- paste0(output.dir,'EM_recruitment_range_FIA_nomid.rds')
     remote_mortality_range_path.nomid <- paste0(output.dir,'mortality_range_FIA_nomid.rds')
         remote_soils_range_path.nomid <- paste0(output.dir,'soils_range_FIA_nomid.rds')
      remote_soils.ff_range_path.nomid <- paste0(output.dir,'soils.ff_range_FIA_nomid.rds')
         
  #model summary path.nomids
          remote_beta_summary_path.nomid <- paste0(output.dir,'beta_summary_FIA_nomid.rds')
        remote_growth_summary_path.nomid <- paste0(output.dir,'growth_summary_FIA_nomid.rds')
   remote_growth.gaus_summary_path.nomid <- paste0(output.dir,'growth.gaus_summary_FIA_nomid.rds')
remote_AM_recruitment_summary_path.nomid <- paste0(output.dir,'AM_recruitment_summary_FIA_nomid.rds')
remote_EM_recruitment_summary_path.nomid <- paste0(output.dir,'EM_recruitment_summary_FIA_nomid.rds')
     remote_mortality_summary_path.nomid <- paste0(output.dir,'mortality_summary_FIA_nomid.rds')
         remote_soils_summary_path.nomid <- paste0(output.dir,'soils_summary_FIA_nomid.rds')  
      remote_soils.ff_summary_path.nomid <- paste0(output.dir,'soils.ff_summary_FIA_nomid.rds')  
         
  #model prediction path.nomids
          remote_beta_prediction_path.nomid <- paste0(output.dir,'beta_prediction_FIA_nomid.rds')
        remote_growth_prediction_path.nomid <- paste0(output.dir,'growth_prediction_FIA_nomid.rds')
   remote_growth.gaus_prediction_path.nomid <- paste0(output.dir,'growth.gaus_prediction_FIA_nomid.rds')
remote_AM_recruitment_prediction_path.nomid <- paste0(output.dir,'AM_recruitment_prediction_FIA_nomid.rds')
remote_EM_recruitment_prediction_path.nomid <- paste0(output.dir,'EM_recruitment_prediction_FIA_nomid.rds')
     remote_mortality_prediction_path.nomid <- paste0(output.dir,'mortality_prediction_FIA_nomid.rds')
         remote_soils_prediction_path.nomid <- paste0(output.dir,'soils_prediction_FIA_nomid.rds')
      remote_soils.ff_prediction_path.nomid <- paste0(output.dir,'soils.ff_prediction_FIA_nomid.rds')
         
         
#repeat nomid these for the local path.nomids.
output.dir <- paste0(main,'JAGS_output/')
  #model path.nomids
          local_beta_model_path.nomid <- paste0(output.dir,'beta_model_FIA_nomid.rds')
        local_growth_model_path.nomid <- paste0(output.dir,'growth_model_FIA_nomid.rds')
   local_growth.gaus_model_path.nomid <- paste0(output.dir,'growth.gaus_model_FIA_nomid.rds')
local_AM_recruitment_model_path.nomid <- paste0(output.dir,'AM_recruitment_model_FIA_nomid.rds')
local_EM_recruitment_model_path.nomid <- paste0(output.dir,'EM_recruitment_model_FIA_nomid.rds')
     local_mortality_model_path.nomid <- paste0(output.dir,'mortality_model_FIA_nomid.rds')
         local_soils_model_path.nomid <- paste0(output.dir,'soils_model_FIA_nomid.rds')
      local_soils.ff_model_path.nomid <- paste0(output.dir,'soils.ff_model_FIA_nomid.rds')
         
  #filtered data path.nomids
          local_beta_filtered_path.nomid <- paste0(output.dir,'beta_filtered_FIA_nomid.rds')
        local_growth_filtered_path.nomid <- paste0(output.dir,'growth_filtered_FIA_nomid.rds')
   local_growth.gaus_filtered_path.nomid <- paste0(output.dir,'growth.gaus_filtered_FIA_nomid.rds')
local_AM_recruitment_filtered_path.nomid <- paste0(output.dir,'AM_recruitment_filtered_FIA_nomid.rds')
local_EM_recruitment_filtered_path.nomid <- paste0(output.dir,'EM_recruitment_filtered_FIA_nomid.rds')  
     local_mortality_filtered_path.nomid <- paste0(output.dir,'mortality_filtered_FIA_nomid.rds')
         local_soils_filtered_path.nomid <- paste0(output.dir,'soils_filtered_FIA_nomid.rds')
      local_soils.ff_filtered_path.nomid <- paste0(output.dir,'soils.ff_filtered_FIA_nomid.rds')
         
  #data predictor range path.nomids
          local_beta_range_path.nomid <- paste0(output.dir,'beta_range_FIA_nomid.rds')
        local_growth_range_path.nomid <- paste0(output.dir,'growth_range_FIA_nomid.rds')
   local_growth.gaus_range_path.nomid <- paste0(output.dir,'growth.gaus_range_FIA_nomid.rds')
local_AM_recruitment_range_path.nomid <- paste0(output.dir,'AM_recruitment_range_FIA_nomid.rds')
local_EM_recruitment_range_path.nomid <- paste0(output.dir,'EM_recruitment_range_FIA_nomid.rds')
     local_mortality_range_path.nomid <- paste0(output.dir,'mortality_range_FIA_nomid.rds')
         local_soils_range_path.nomid <- paste0(output.dir,'soils_range_FIA_nomid.rds')
      local_soils.ff_range_path.nomid <- paste0(output.dir,'soils.ff_range_FIA_nomid.rds')
         
  #model summary path.nomids
          local_beta_summary_path.nomid <- paste0(output.dir,'beta_summary_FIA_nomid.rds')
        local_growth_summary_path.nomid <- paste0(output.dir,'growth_summary_FIA_nomid.rds')
   local_growth.gaus_summary_path.nomid <- paste0(output.dir,'growth.gaus_summary_FIA_nomid.rds')
local_AM_recruitment_summary_path.nomid <- paste0(output.dir,'AM_recruitment_summary_FIA_nomid.rds')
local_EM_recruitment_summary_path.nomid <- paste0(output.dir,'EM_recruitment_summary_FIA_nomid.rds')
     local_mortality_summary_path.nomid <- paste0(output.dir,'mortality_summary_FIA_nomid.rds')
         local_soils_summary_path.nomid <- paste0(output.dir,'soils_summary_FIA_nomid.rds')  
      local_soils.ff_summary_path.nomid <- paste0(output.dir,'soils.ff_summary_FIA_nomid.rds')  
         
  #model prediction path.nomids
          local_beta_prediction_path.nomid <- paste0(output.dir,'beta_prediction_FIA_nomid.rds')
        local_growth_prediction_path.nomid <- paste0(output.dir,'growth_prediction_FIA_nomid.rds')
   local_growth.gaus_prediction_path.nomid <- paste0(output.dir,'growth.gaus_prediction_FIA_nomid.rds')
local_AM_recruitment_prediction_path.nomid <- paste0(output.dir,'AM_recruitment_prediction_FIA_nomid.rds')
local_EM_recruitment_prediction_path.nomid <- paste0(output.dir,'EM_recruitment_prediction_FIA_nomid.rds')
     local_mortality_prediction_path.nomid <- paste0(output.dir,'mortality_prediction_FIA_nomid.rds')
         local_soils_prediction_path.nomid <- paste0(output.dir,'soils_prediction_FIA_nomid.rds')
      local_soils.ff_prediction_path.nomid <- paste0(output.dir,'soils.ff_prediction_FIA_nomid.rds')
         
#Detrend and bins paths
  #ALL data.
  bins_growth.gaus_all_AM_path <- paste0(detrend.directory,'/bins_growth.gaus_AM_all.rds')
  bins_growth.gaus_all_EM_path <- paste0(detrend.directory,'/bins_growth.gaus_EM_all.rds')
     bins_recruitment_all_path <- paste0(detrend.directory,'/bins_recruitment_all.rds')
       bins_growth_all_AM_path <- paste0(detrend.directory,'/bins_growth_AM_all.rds')
       bins_growth_all_EM_path <- paste0(detrend.directory,'/bins_growth_EM_all.rds')
            bins_beta_all_path <- paste0(detrend.directory,'/bins_beta_all.rds')
       bins_mortality_all_path <- paste0(detrend.directory,'/bins_mortality_all.rds')
           bins_soils_all_path <- paste0(detrend.directory,'/bins_soils_all.rds')
   all_binned_rsq_summary_path <- paste0(detrend.directory,'/bins_rsq_all.rds')
   
   #EAST data.
   bins_growth.gaus_east_AM_path <- paste0(detrend.directory,'/bins_growth.gaus_AM_east.rds')
   bins_growth.gaus_east_EM_path <- paste0(detrend.directory,'/bins_growth.gaus_EM_east.rds')
      bins_recruitment_east_path <- paste0(detrend.directory,'/bins_recruitment_east.rds')
        bins_growth_east_AM_path <- paste0(detrend.directory,'/bins_growth_AM_east.rds')
        bins_growth_east_EM_path <- paste0(detrend.directory,'/bins_growth_EM_east.rds')
             bins_beta_east_path <- paste0(detrend.directory,'/bins_beta_east.rds')
        bins_mortality_east_path <- paste0(detrend.directory,'/bins_mortality_east.rds')
            bins_soils_east_path <- paste0(detrend.directory,'/bins_soils_east.rds')
    east_binned_rsq_summary_path <- paste0(detrend.directory,'/bins_rsq_east.rds')
    
    #NOMID data.
    bins_growth.gaus_nomid_AM_path <- paste0(detrend.directory,'/bins_growth.gaus_AM_nomid.rds')
    bins_growth.gaus_nomid_EM_path <- paste0(detrend.directory,'/bins_growth.gaus_EM_nomid.rds')
       bins_recruitment_nomid_path <- paste0(detrend.directory,'/bins_recruitment_nomid.rds')
         bins_growth_nomid_AM_path <- paste0(detrend.directory,'/bins_growth_AM_nomid.rds')
         bins_growth_nomid_EM_path <- paste0(detrend.directory,'/bins_growth_EM_nomid.rds')
              bins_beta_nomid_path <- paste0(detrend.directory,'/bins_beta_nomid.rds')
         bins_mortality_nomid_path <- paste0(detrend.directory,'/bins_mortality_nomid.rds')
             bins_soils_nomid_path <- paste0(detrend.directory,'/bins_soils_nomid.rds')
     nomid_binned_rsq_summary_path <- paste0(detrend.directory,'/bins_rsq_nomid.rds')
     
#Figure paths
           Figure1_all_path <- 'Figures/all_data_figures/Fig_1._Site_Map.png'
           Figure2_all_path <- 'Figures/all_data_figures/Fig_2._binned_mean_fits.png'
           Figure3_all_path <- 'Figures/all_data_figures/Fig_3._soilC_EMxN_interaction.png'
      Supp_Figure1_all_path <- 'Figures/all_data_figures/Supp_Fig_1.east_forest_4panel.png'
      Supp_Figure2_all_path <- 'Figures/all_data_figures/Supp_Fig_2.nomid_forest_4panel.png'
      Supp_Figure3_all_path <- 'Figures/all_data_figures/Supp_Fig_3.east_soilC_ndep.png'
      Supp_Figure4_all_path <- 'Figures/all_data_figures/Supp_Fig_4.east_soilC_ndep.png'
      Supp_Figure5_all_path <- 'Figures/all_data_figures/Supp_Fig_5._Soils_4-panel.png'
      Supp_Figure6_all_path <- 'Figures/all_data_figures/Supp_Fig_6._20-panel_summary.png'
      Supp_Figure7_all_path <- 'Figures/all_data_figures/Supp_Fig_7._soil_forest_proportion_histogram.png'
    soil_semivariogram_path <- 'Figures/soil_semivariogram.png'
      
#paths for beta factor, parameter value and p-value tables.
            all_beta_path <- paste0(summary.directory,'/all_beta_list.rds')
           east_beta_path <- paste0(summary.directory,'/east_beta_list.rds')
          nomid_beta_path <- paste0(summary.directory,'/nomid_beta_list.rds')
     all_model.table_path <- paste0(summary.directory,'/all_model.table_list.rds')
    east_model.table_path <- paste0(summary.directory,'/east_model.table_list.rds')
   nomid_model.table_path <- paste0(summary.directory,'/nomid_model.table_list.rds')
  all_model.table_CI_path <- paste0(summary.directory,'/all_model.table_CI_list.rds')
 east_model.table_CI_path <- paste0(summary.directory,'/east_model.table_CI_list.rds')
nomid_model.table_CI_path <- paste0(summary.directory,'/nomid_model.table_CI_list.rds')
           