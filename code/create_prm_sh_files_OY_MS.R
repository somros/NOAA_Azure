# 12/21/2023
# Alberto Rovellini
# Create a series of harvest.prm files where mFC for all T3 stocks is dialed up or down with a scalar
# NB: this uses the F35 (or proxy thereof) from the sinle-species experiment
# As a first stab we use the old vector of F35 (stocks had different steepness, runtime was different, this will need to change)
# We can still learn something, in particular if we do not compare too much with the SS runs but rather between MS scenarios (until the F35 vector being permuted here comes from the new SS runs)

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
# the multipliers here cannot be the same as in the SS experiment (those were applied to SS FOFL, there are applied to MS F35 or proxy)
# The meaning of FOFL and F35 in SS and MS context, respectively, are not the same
# However, the multipliers here should still give us a good spread of F to compare between SS and MS runs
# The SS multipliers (on FOFL) are: mult <- c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 2.0, 2.5, 3.0, 4.0)
# Let's use the same (13) for dialing MS F35

# So start from f35_vector_PROXY.csv for now
# multiply that vector by the scalar
# write a new harvest.prm for each permutation (13)

# prepare one set of harvest.prm files where F for ATF is fixed (another 13)

# each of these has to be ran with and without CC (13*4 = 52 runs)
library(tidyverse)
select <- dplyr::select

# read groups 
grp <- read.csv('NOAA_Azure/data/GOA_Groups.csv') # functional groups

# Read F35 vector
f_vec <- read.csv("NOAA_Azure/data/f35_vector_PROXY_OY_SS.csv")

# convert to mFC
# Formula for mFC is mfc = 1-exp(-F / 365)

mfc_vec <- f_vec %>%
  left_join(grp %>% select(Code, LongName)) %>%
  mutate(mfc = 1-exp(-atlantis_fmsy / 365)) %>%
  select(Code, LongName, mfc)

# permute
mult <- c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 2.0, 2.5, 3.0, 4.0)

perm_frame <- data.frame(matrix(nrow = length(mfc_vec$Code), ncol = length(mult)+1))
colnames(perm_frame) <- c('Code',1:length(mult)) # p6 = F35

for(i in 1:length(mfc_vec$Code)){
  
  this_sp <- mfc_vec$Code[i]
  this_mfc <- mfc_vec %>% filter(Code == this_sp) %>% pull(mfc)
  mfc_perm <- this_mfc * mult
  perm_frame[i,1] <- this_sp
  perm_frame[i,2:ncol(perm_frame)] <- mfc_perm
  
}

# Write mFC matrix --------------------------------------------------------
nfleet <- 33 # how many fleets in the harvest.prm file?

# we need the burn-in period here
# burn in should occur under old background F (calibrated 1/4 FOFL)
# because of the Atlantis setup, we need to relate the burn-in mFC (base aka background aka 1/4 FOFL) to the MS F35 value
# Step 1: pull background mFC (calibration value) for the species
# Step 2: relate background mFC to mfC in the permutation frame created above, as m1 = mFC@perm / mFC@base
# These steps need to occur for all species simultaneously

# T3 stocks from FMP and Halibut
stocks <- c('POL','COD','SBF','FFS','FFD','REX','ATF','FHS','POP','RFS','RFP','HAL')

# get base values from 1517 harvest.prm
prm_file <- here('goa_runs','AtlantisGOA_OY_SS','GOA_harvest_background.prm')
prm_vals <- readLines(prm_file)

# make a path for multispecies f
f_path <- here('goa_runs','AtlantisGOA_OY_MS')

idx <- 0

# for each permutation
for(m in 1:length(mult)){
  
  idx <- idx+1
  
  # create file
  newfile <-  paste0(f_path, '/GOA_harvest_', idx, '.prm')
  
  file.create(newfile)
  
  this_mult <- mult[m]
  
  prm_new <- prm_vals
  
  # one stock at a time 
  for(j in 1:length(stocks)){
    
    this_sp <- stocks[j]
    
    # get background mFC for this stock
    bg_mfc_line <- grep(paste0('mFC_', this_sp), prm_new) + 1 # get the line where the background mFC value is stored
    bg_mfc <- prm_new[bg_mfc_line] %>% strsplit(" ") %>% unlist() %>% head(1) %>% as.numeric()
    
    # get new mFC in this permutation for this stock
    new_mfc <- perm_frame[perm_frame$Code == this_sp,(m+1)]
    
    # get scaling factor to go from bg to permutation
    scl_f <- new_mfc / bg_mfc
    
    # now manipulate the mFCchange parameters for this species
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
    new_l6 <- as.character(scl_f) 
    
    # make a new harvest.prm object and modify it
    prm_new[l1] <- new_l1
    prm_new[l2] <- new_l2
    prm_new[l3] <- new_l3
    prm_new[l4] <- new_l4
    prm_new[l5] <- new_l5
    prm_new[l6] <- new_l6
    
  }
  
  # write to file
  writeLines(prm_new, con = newfile)
  
}

# now create a copy of these files, but with ATF mFC fixed to background value
# this entails changing the parameter mFCchange_mult_ATF to 1, so that the burn-in value is kept

harvest_files <- list.files("goa_runs/AtlantisGOA_OY_MS/", pattern = "GOA_harvest_[0-9]*.prm", full.names = T)

# reorder these based on the number in the filename
num_idx <- as.numeric(sub(".*_([0-9]+)\\.prm", "\\1", harvest_files))
harvest_files <- harvest_files[order(num_idx)]

for(i in 1:length(harvest_files)){
  
  # pick file
  this_file <- harvest_files[i]
  this_file_label <- gsub(".prm", "", gsub("goa_runs/AtlantisGOA_OY_MS//", "", this_file))
  
  # copy file
  new_file <- paste0("goa_runs/AtlantisGOA_OY_MS//",
                    this_file_label,
                    "_ATF.prm")
  # create new file
  file.copy(this_file, new_file)
  
  # open new file
  # get base values from 1517 harvest.prm
  new_prm_vals <- readLines(new_file)
  
  # update line for ATF mortality
  # mFCchange_mult_XXX
  line_to_change <- grep(paste0('mFCchange_mult_ATF'), new_prm_vals) + 1
  new_prm_vals[line_to_change] <- as.character(1)
  
  writeLines(new_prm_vals, con = new_file)
  
}

# now write sh files so that they account for correct forcing and ATF regime
ms_path <- "goa_runs/AtlantisGOA_OY_MS/"

# forcing files
# copy forcings.prm file from AtlantisGOA_OY_SS
# file modified manually. There is a burn-in on 1999 conditions
# Arguably the base conditions should be a climatology and not 1999, we've had problems before
# Now the SS runs have been done with 1999. Changing this means that the reference points will change
# There is an argument for keeping 1999 because of the burn-in but it is not very convincing
# you should do the burn-in at 1999 (base) but then get reference points and stock status under climatologies for the SS runs

# create sh files
# here is where the numbering becomes important and we should make a key with:
# Harvest.prm and corresponding multiplier
# Force.prm
# relist harves files from new folder
harvest_files_atf <- list.files(ms_path, pattern = "GOA_harvest.*ATF", full.names = T)
harvest_files_atf <- harvest_files_atf[order(num_idx)]

key_ms <- data.frame("run" = c(rep("base",length(mult)),
                               rep("atf", length(mult)),
                               rep("climate", length(mult)),
                               rep("atf_climate", length(mult))),
                     "mult" = rep(mult, 4),
                     "harvest" = rep(c(harvest_files, harvest_files_atf), 2),
                     "force" = c(rep("GOA_force_base.prm", 26), rep("GOA_force_warm.prm", 26)),
                     "idx" = 1:(4*length(mult)))

# save
# write.table(key_ms, "NOAA_Azure/data/oy_key.csv", row.names = F)

# now make the RunAtlantis.sh scripts pointing to the correct path and creating the correct output
for (i in key_ms$idx) {
  
  # get the harvest.prm and force.prm that we need for each run
  this_harvest <- key_ms %>% filter(idx == i) %>% pull(harvest)
  this_harvest <- gsub("goa_runs/AtlantisGOA_OY_MS//", "", this_harvest)
  this_force <- key_ms %>% filter(idx == i) %>% pull(force)
  
  # create Atlantis sh script
  filenm <- paste0(ms_path, "/RunAtlantis_",i,".sh")
  fileConn<-file(filenm,open="w")
  cat(paste0("atlantisMerged -i GOA_cb_summer.nc  0 -o output_", 
             i,
             ".nc -r GOA_run.prm -f ",
             this_force,
             " -p GOA_physics.prm -b GOAbioparam_test_OY_SS.prm -h ",
             this_harvest,
             " -m GOAMigrations.csv -s GOA_Groups.csv -q GOA_fisheries.csv -d output\n"),file = fileConn,append=T)
  close(fileConn)
  system(paste0("chmod 775 ",filenm))
}
