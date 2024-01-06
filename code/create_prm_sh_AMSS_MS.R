# Alberto Rovellini
# 1/5/2024
# This code creates a folder with Atlantis model inputs for multispecies fishing experiments in batch
# Current set of experiments is:
# 1. Base MS runs with a set of F35 permutations
# 2. Same as 1 but with ATF F fixed to 1/4 FOFL
# 3. Same as 1 but with warm forcings and decreased productivity
# 4. Same as 2 but with warm forcings and decreased productivity
# All these runs can be performed with the right combination of input files, namely harvest and climate forcing files
# The harvest files have already been created in make_f35_permutations.R
# This code will:
# 1. Copy an existing Atlantis folder to use as template (ensure we are using the correct PRM) - this has to be consistent with the SS experiment if you want to compare
# In this case, consistency with the SS experiment exposes us to incorrect FSPB settings
# FOR AMSS ONLY: there is an argument to drift off from SS biological parameters 
# NOTE: this means no comparisons with SS, which is fine for AMSS
# Just upload 1517 for the purpose of AMSS and the correct biological parameters 
# 2. Copy over the harvest.prm from the code/f35_perm folders
# 3. Create a set of (2) forcing files with and without climate. In case of no burn-in, just copy over the file from 1331 (this allows burn-in for climate but not for F)
# 4. Create a set of 44 run.sh files that point to the correct combinations of ATF and forcings


library(tidyverse)
ms_path <- "goa_runs/AtlantisGOA_AMSS/"

# copy F files
files_4 <- list.files("NOAA_Azure/code/f35_perms/4", full.names = T)
files_4_atf <- list.files("NOAA_Azure/code/f35_perms/4-ATF-low/", full.names = T)

for(i in 1:length(files_4)){
  copyFile(files_4[i], ms_path)
  copyFile(files_4_atf[i], ms_path)
}

# 3. Copied force.prm from 1331

# create sh files
# here is where the numbering becomes important and we should make a key with:
# Harvest.prm and corresponding multiplier
# Force.prm
mult <- seq(0,4,length.out=11)
# relist harves files from new folder
harvest_files_atf <- list.files(ms_path, pattern = "GOA_harvest.*ATF")
harvest_files <- setdiff(list.files(ms_path, pattern = "GOA_harvest"), harvest_files_atf)

# reorder these based on the number in the filename
num_idx <- as.numeric(sub(".*_([0-9]+)\\.prm", "\\1", harvest_files))
harvest_files_atf <- harvest_files_atf[order(num_idx)]
harvest_files <- harvest_files[order(num_idx)]

key_ms <- data.frame("run" = c(rep("base",length(mult)),
                                   rep("atf", length(mult)),
                                   rep("climate", length(mult)),
                                   rep("atf_climate", length(mult))),
                     "mult" = rep(mult, 4),
                     "harvest" = rep(c(harvest_files, harvest_files_atf), 2),
                     "force" = c(rep("GOA_force.prm", 22), rep("GOA_force_warm.prm", 22)),
                     "idx" = 1:(4*length(mult)))

# now make the RunAtlantis.sh scripts pointing to the correct path and creating the correct output
for (i in key_ms$idx) {
  
  # get the harvest.prm and force.prm that we need for each run
  this_harvest <- key_ms %>% filter(idx == i) %>% pull(harvest)
  this_force <- key_ms %>% filter(idx == i) %>% pull(force)
  
  # create Atlantis sh script
  filenm <- paste0(ms_path, "/RunAtlantis_",i,".sh")
  fileConn<-file(filenm,open="w")
  cat(paste0("atlantisMerged -i GOA_cb_summer.nc  0 -o output_", 
             i,
             ".nc -r GOA_run.prm -f ",
             this_force,
             " -p GOA_physics.prm -b GOAbioparam_test.prm -h ",
             this_harvest,
             " -m GOAMigrations.csv -s GOA_Groups.csv -q GOA_fisheries.csv -d output\n"),file = fileConn,append=T)
  close(fileConn)
  system(paste0("chmod 775 ",filenm))
}

