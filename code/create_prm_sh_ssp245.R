# Alberto Rovellini
# 04/03/2024
# this code adds the sh files and the force.prm files to be added to the Atlantis_GOA_OY_MS folder, to run the ssp245 runs
# To re-do this clone the Atlantis_GOA_OY_MS repo:
# git clone https://github.com/somros/Atlantis_GOA_OY_MS.git
library(tidyverse)

# what's the path to the folder
ms_path <- "C:/Users/Alberto Rovellini/Documents/GOA/Atlantis_GOA_OY_MS/"

# force.prm files
force_file <- paste0(ms_path, "/GOA_force_warm.prm")
force_warm <- readLines(force_file)

# change the paths
force_warm <- gsub("hydro/hydro/", "hydro/hydro_ssp245/", force_warm) # hydro
force_warm <- gsub("goa_roms_temp_2075_2085", "goa_roms_temp_2075_2085_ssp245", force_warm) # temp
force_warm <- gsub("goa_roms_salt_2075_2085", "goa_roms_salt_2075_2085_ssp245", force_warm) # salt

# write out
writeLines(force_warm, "C:/Users/Alberto Rovellini/Documents/GOA/Atlantis_GOA_OY_MS/GOA_force_warm_ssp245.prm")

# now create the sh files. Let's identify which runs are using the climate forcings (ssp585 originally)
# this "key" file also needs updating
key <- read.csv("data/oy_key.csv")
clim_idx <- key %>% filter(run %in% c("climate", "atf_climate")) %>% pull(idx)

# add rows
key_245 <- key %>% rbind(
  key %>% 
    filter(run %in% c("climate", "atf_climate")) %>%
    mutate(run = ifelse(run == "climate", "climate_245", "atf_climate_245"),
           force = "GOA_force_warm_ssp245.prm",
           idx = 53:78)
)

# write out
write.csv(key_245, "data/oy_key_ssp245.csv", row.names = F)

# now make the RunAtlantis.sh scripts pointing to the correct path and creating the correct output
key_ms <- key_245

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
  #system(paste0("chmod 775 ",filenm))
}
