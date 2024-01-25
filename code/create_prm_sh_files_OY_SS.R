# 12/21/2023
# Alberto Rovellini
# Create a series of harvest.prm files where mFC for one species changes while mFC for the other species is fixed to 1/4 FOFL (light fishing)
# NB: this uses the calibrated value of mFC that results in 1/4 FOFL catch. See run 768 and all subsequent runs that inherit this (including 1328)
# This is for the single-species F testing
# Do it on Tier 3 species to start
# For now, we only have one fleet, so only the first value of the mFC vector is of interest
# Indexing for the harvest.prm and run.sh files is in increasing integers to be used by the foreach loop in doAzureParallel

# This needs to now incorporate a burn-in of 30 years
# we need to use the parameters:

# flagchangeF 1

# flagFchange_XXX 1

# XXX_mFC_changes	33
# 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ## because it is BG fishing that changes

# mFCchange_start_XXX	1
# 10950 ## this is 30 * 365
 
# mFCchange_period_XXX
# 1 ## do it over 1 day

# mFCchange_mult_XXX
# this is the tricky part

# We need to set burn-in F to 1/4 FOFL (calibrated)
# So start from 1473 harvest.prm
# The F ranges to test are multipliers of FOFL in this single-species experiment, that have the purpose of identifying F35
# The range of multipliers on FOFL is important - if we push it to 4 FOFL we lack the resolution to get F35, but if we only do up to 2 * FOFL we miss what happens at high F
# Options are:
# 1. Do a lot of runs (unweildy and unnecessary)
# 2. Do uneven spacing of the multipliers: 0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 2.0, 2.5, 3.0, 4.0
# this would be 16 * 12 = 192 runs - seems reasonable
# The multipliers will need to be:
# 1 for all non-target species
# 4*0, 4*0.2, 4*0.4, ... 4*4.0 for the target species. This brings calibration F back to FOFL, and then scales it

library(readxl)
library(tidyverse)
library(here)

nfleet <- 33 # how many fleets in the harvest.prm file?
select <- dplyr::select

# Create a folder from the Base model that contains all the harvest.prm and run.sh files
# This folder needs to be cloned on each node of the foreach loop
# if(!file.exists(here('goa_runs','AtlantisGOA_OY_SS'))){
#   system(paste('sudo cp -r', here('goa_runs','AtlantisGOA_Base'), here('goa_runs','AtlantisGOA_F_test_4')))
#   system(paste('sudo chmod -R a+rwx ', here('goa_runs','AtlantisGOA_F_test_4'))) # add permission
#   # remove harvest.prm and run.sh files as we will create new ones, also remove any output file that may have been copied
#   system(paste('sudo rm -r ', 
#                here('goa_runs','AtlantisGOA_F_test_4','outputFolder/'),
#                here('goa_runs','AtlantisGOA_F_test_4','out14'),
#                here('goa_runs','AtlantisGOA_F_test_4','GOA_harvest_background.prm'),
#                here('goa_runs','AtlantisGOA_F_test_4','RunAtlantis.sh')))
#   system(paste('sudo rm -rf ', here('goa_runs','AtlantisGOA_F_test_4','.git'))) # get rid of the git tracking from the base model folder
# }

# T3 stocks from FMP and Halibut
stocks <- c('POL','COD','SBF','FFS','FFD','REX','ATF','FHS','POP','RFS','RFP','HAL')

# make mFC vectors based on a range of F:
# start from the mFC values in the harvest.prm
# we do not need to modify these values directly
# we act on the multiplier
prm_file <- here('goa_runs','AtlantisGOA_OY_SS','GOA_harvest_background.prm')
prm_vals <- readLines(prm_file)

# get set of multipliers
mult <- c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 2.0, 2.5, 3.0, 4.0)

# set the path to the Atlantis folder
f_path <- here('goa_runs','AtlantisGOA_OY_SS')

# Produce a look-up table for which species and F each file corresponds to
# do for each species
idx <- 0
f_lookup_ls_SS <- list()

for(sp in stocks){
  
  for(m in 1:length(mult)){
    
    idx <- idx+1
    
    # create file
    newfile <-  paste0(f_path, '/GOA_harvest_', idx, '.prm')
    
    file.create(newfile)
    
    this_sp <- sp
    this_mult <- mult[m] * 4 # because the default is 1/4 FOFL (calibrated)
    
    # find the lines that have parameters for the species of interest
    # flagchangeF 1
    l1 <- grep('flagchangeF', prm_vals)
    new_l1 <- 'flagchangeF	1	 ## F mortality rate forcing changes [1] Yes or [0] remains constant'
    # flagFchange_XXX 1
    l2 <- grep(paste0('flagFchange_',this_sp), prm_vals) + 1
    new_l2 <- paste(as.character(c(1, rep(0, nfleet-1))), collapse = ' ') # turns the species-specific flag to 1
    # XXX_mFC_changes 
    l3 <- grep(paste0(this_sp,'_mFC_changes'), prm_vals) + 1
    new_l3 <- paste(as.character(c(1, rep(0, nfleet-1))), collapse = ' ')
    # mFCchange_start_XXX
    l4 <- grep(paste0('mFCchange_start_',this_sp), prm_vals) + 1
    new_l4 <- as.character(10950) # 30 * 365 # because the burn-in is 30 years
    # mFCchange_period_XXX
    l5 <- grep(paste0('mFCchange_period_',this_sp), prm_vals) + 1
    new_l5 <- as.character(1) # do the change over 1 day (instantly)
    # mFCchange_mult_XXX
    l6 <- grep(paste0('mFCchange_mult_',this_sp), prm_vals) + 1
    new_l6 <- as.character(this_mult) 
    
    # make a new harvest.prm object and modify it
    prm_new <- prm_vals
    prm_new[l1] <- new_l1
    prm_new[l2] <- new_l2
    prm_new[l3] <- new_l3
    prm_new[l4] <- new_l4
    prm_new[l5] <- new_l5
    prm_new[l6] <- new_l6
    
    # write to file
    writeLines(prm_new, con = newfile)
    
    # for later lookup, get mFC for this group and its multiplier, and convert to F
    mfc_line <- grep(paste0('mFC_', this_sp), prm_vals) + 1
    base_mfc <- prm_new[mfc_line] %>% strsplit(" ") %>% unlist() %>% head(1) %>% as.numeric()
    mult_mfc <- base_mfc * this_mult
    # Formula for mFC is mfc = 1-exp(-F / 365)
    # F = -365 * log(1-mFC)
    base_f <- -365 * log(1 - base_mfc)
    mult_f <- -365 * log(1 - mult_mfc)
    
    # add to index f_lookup
    f_lookup_ls_SS[[idx]] <- data.frame('idx'=idx, 'species'=this_sp, 'base_mfc'= base_mfc, 'mult_mfc' = mult_mfc, 'base_f' = base_f, 'mult_f'= mult_f, 'mult_idx' = m, 'fofl_mult' = mult[m])
    
  }
  
}

f_lookup <- bind_rows(f_lookup_ls_SS) # save this 
write.csv(f_lookup, here('NOAA_Azure','data','f_lookup_OY_SS.csv'),row.names = F)

# now make the RunAtlantis.sh scripts pointing to the correct path and creating the correct output
for (idx in f_lookup$idx) {
  
  # create Atlantis sh script
  filenm <- paste0(f_path, "/RunAtlantis_",idx,".sh")
  fileConn<-file(filenm,open="w")
  cat(paste0("atlantisMerged -i GOA_cb_summer.nc  0 -o output_", 
             idx,
             ".nc -r GOA_run.prm -f GOA_force_base.prm -p GOA_physics.prm -b GOAbioparam_test_OY_SS.prm -h GOA_harvest_",
             idx,
             ".prm -m GOAMigrations.csv -s GOA_Groups.csv -q GOA_fisheries.csv -d output\n"),file = fileConn,append=T)
  close(fileConn)
  system(paste0("chmod 775 ",filenm))
}
