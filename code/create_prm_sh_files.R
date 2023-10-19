# 10/16/2023
# Alberto Rovellini
# Create a series of harvest.prm files where mFC for one species changes
# This is for the single-species F testing
# Do it on Tier 3 species to start
# For now, we only have one fleet, so only the first value of the mFC vector is of interest
# Indexing for the harvest.prm and run.sh files is in increasing integers to be used by the foreach loop in doAzureParallel
library(readxl)
library(tidyverse)
library(here)

nfleet <- 33 # how many fleets in the harvest.prm file?
select <- dplyr::select

# Create a folder from the Base model that contains all the harvest.prm and run.sh files
# This folder needs to be cloned on each node of the foreach loop
if(!file.exists(here('goa_runs','AtlantisGOA_F_test'))){
  system(paste('sudo cp -r', here('goa_runs','AtlantisGOA_Base'), here('goa_runs','AtlantisGOA_F_test')))
  system(paste('sudo chmod -R a+rwx ', here('goa_runs','AtlantisGOA_F_test'))) # add permission
  # remove harvest.prm and run.sh files as we will create new ones, also remove any output file that may have been copied
  system(paste('sudo rm -r ', 
               here('goa_runs','AtlantisGOA_F_test','outputFolder/'),
               here('goa_runs','AtlantisGOA_F_test','out14'),
               here('goa_runs','AtlantisGOA_F_test','GOA_harvest_background.prm'),
               here('goa_runs','AtlantisGOA_F_test','RunAtlantis.sh')))
  system(paste('sudo rm -rf ', here('goa_runs','AtlantisGOA_F_test','.git'))) # get rid of the git tracking from the base model folder
}

# make mFC vectors based on a range of F:
# start from 0 and go up to 2*FOFL - test 10 values (closer together the closer to 0 we are?)
# what to do with imposed and realized f?

# Formula for mFC is mfc = 1-exp(-F / 365)

# read FOFL table
# read in MSY information (from FMP)
# setting the cap for the range of F values to explore to 2FOFL and testing 8 equally-spaced values between 0 and this cap
f_cap <- read_xlsx(here('NOAA_Azure','data','msy.xlsx'), sheet = 1, range = 'A3:J19') %>%
  select(Stock, FOFL) %>%
  mutate(Code = c('POL','COD','SBF','FFS','FFS','FFS','FFS','FFD',
                  'REX','REX','ATF','FHS','POP','RFS','RFS','RFP')) %>%
  group_by(Code) %>%
  summarise(FOFL = mean(FOFL)) %>%
  ungroup() %>%
  mutate(cap = 2*FOFL) %>%
  select(-FOFL)

# add halibut (2M)
f_cap <- rbind(f_cap, data.frame('Code'='HAL', cap = 0.4)) # 0.2*2
species_to_test <- unique(f_cap$Code)

# make range of values, stored in a matrix
n_f <- 8
mfc_tab <- f_tab <- matrix(NA, nrow = length(species_to_test), ncol = n_f)

# make table of F and mFC values for each species
for(i in 1:length(species_to_test)){
  sp <- species_to_test[i]
  f <- f_cap %>% filter(Code == sp) %>% pull(cap)
  fvec <- seq(0,f,length.out=n_f)
  mfcvec <- 1-exp(-fvec / 365)
  mfcrange <- matrix(mfcvec, nrow = 1)
  
  # store
  mfc_tab[i,] <- mfcrange
  
  # also make df with F values, we need to store it for plotting
  frange <- matrix(fvec, nrow = 1)
  
  # store
  f_tab[i,] <- frange
}

# set names
mfc_tab <- mfc_tab %>% 
  data.frame() %>% 
  mutate(Code = species_to_test) %>%
  set_names(c(1:8, 'Code'))

f_tab <- f_tab %>% 
  data.frame() %>% 
  mutate(Code = species_to_test) %>%
  set_names(c(1:8, 'Code'))

prm_file <- here('NOAA_Azure','data','GOA_harvest_background.prm')
prm_vals <- readLines(prm_file)

# set the path to the Atlantis folder
f_path <- here('goa_runs','AtlantisGOA_F_test')

# Produce a look-up table for which species and F each file corresponds to
# do for each species
idx <- 0
f_lookup_ls <- list()

for(sp in species_to_test){
  
  for(f in 1:n_f){
    
    idx <- idx+1
    
    # create file
    newfile <-  paste0(f_path, '/GOA_harvest_', idx, '.prm')
    
    file.create(newfile)
    
    this_sp <- sp
    this_mfc <- mfc_tab %>% filter(Code == sp) %>% pull(as.character(f))
    this_f <- f_tab %>% filter(Code == sp) %>% pull(as.character(f)) # for f_lookup and bookkeeping
    
    # find the line that has the mFC vector for the species of interest
    this_line <- grep(paste0('mFC_',this_sp), prm_vals) + 1
    
    new_mFC_vector <- paste(as.character(c(this_mfc, rep(0, nfleet-1))), collapse = ' ')
    
    # make a new harvest.prm object and modify it
    prm_new <- prm_vals
    prm_new[this_line] <- new_mFC_vector
    
    # write to file
    writeLines(prm_new, con = newfile)
    
    # add to index f_lookup
    f_lookup_ls[[idx]] <- data.frame('idx'=idx, 'species'=this_sp, 'mfc'= this_mfc, 'f' = this_f, 'fidx' = f)
    
  }
  
}

f_lookup <- bind_rows(f_lookup_ls) # save this 
write.csv(f_lookup, here('NOAA_Azure','data','f_lookup.csv'),row.names = F)

# now make the RunAtlantis.sh scripts pointing to the correct path and creating the correct output
for (idx in f_lookup$idx) {

  # create Atlantis sh script
  filenm <- paste0(f_path, "/RunAtlantis_",idx,".sh")
  fileConn<-file(filenm,open="w")
  cat(paste0("atlantisMerged -i GOA_cb_summer.nc  0 -o output_", 
             idx,
             ".nc -r GOA_run.prm -f GOA_force.prm -p GOA_physics.prm -b GOAbioparam_test.prm -h GOA_harvest_",
             idx,
             ".prm -m GOAMigrations.csv -s GOA_Groups.csv -q GOA_fisheries.csv -d output\n"),file = fileConn,append=T)
  close(fileConn)
  system(paste0("chmod 775 ",filenm))
}
