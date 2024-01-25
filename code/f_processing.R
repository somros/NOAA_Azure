# Alberto Rovellini
# 11/7/2022
# this code pulls terminal biomass and catch from each run, and it calculates realized F based on catch in the last 5 years
# these variables are stored in a table with one row and written to a csv file, to be processed later
# NOTE: this version is to pull results from a folder where they have been downloaded from the blob storage
# If we get the merging to wor, it will be a different scoping possibly

# here is an example to check out with azcopy:
# sudo azcopy cp "https://atlantisgoastorage.blob.core.windows.net/job20231228020254?sp=rl&st=2024-01-05T20:02:16Z&se=2024-01-06T04:02:16Z&spr=https&sv=2022-11-02&sr=c&sig=7U2qIPMQTSGr8HpbDsrnRW6b5DvGeEwi3T%2Bfw8yhFUg%3D" "." --recursive --include-pattern="*result.rds"
library(dplyr)
library(tidyr)
library(here)

this_job <- "job20240109052342"
outdir <- paste0('NOAA_Azure/results/post-processing/',this_job)
dir.create(outdir)

burnin <- 30 # years of burn-in

folder_path <- here("NOAA_Azure","results","pre-processing",this_job,"results")
# list the results files
results_list <- list.files(folder_path, full.names = T)

# list the data
grp_path <- here('NOAA_Azure/data/GOA_Groups.csv') # functional groups
mat_path <- here('NOAA_Azure/data/age_at_mat.csv') # maturity at 50%
fspb_path <- here('NOAA_Azure/data/fspb.csv') # proportion mature at age
selex_path <- here('NOAA_Azure/data/age_at_selex_new.csv') # selectivity pattern, after adjusting startage
lookup_path <- here('NOAA_Azure/data/f_lookup_OY_SS.csv')

atlantis_fg <- read.csv(grp_path)
mat <- read.csv(mat_path, header = T) # age at maturity
fspb <- read.csv(fspb_path, header = T) # proportion mature
# reshape fspb
fspb <- fspb %>%
  pivot_longer(-Code, names_to = "Age", values_to = "fspb") %>%
  mutate(Age = gsub("X","",Age))

selex <- read.csv(selex_path, header = T) # age at selectivity
f_lookup <- read.csv(lookup_path) # lookup for species and F per run

for(i in 1:length(results_list)){
  
  print(paste("Doing",i))
  
  # load results
  this_result <- readRDS(results_list[i])
  
  # this idx - CAREFUL! This is not the same as i because of alphabetical sorting of the reuslts files
  this_idx <- as.numeric(names(this_result))

  # first identify the species and the level of fishing for this run. These are unrelated from runname
  sp <- f_lookup %>% filter(idx == this_idx) %>% pull(species)
  fidx <- f_lookup %>% filter(idx == this_idx) %>% pull(mult_idx)
  
  # rename the result object in the list to avoid problems with indexing
  names(this_result) <- "res"
  
  # extract tables from results
  biomage <- this_result$res[[paste0("biomage_",this_idx)]]
  catch <- this_result$res[[paste0("catch_",this_idx)]]
  
  # # check catch soemthis in wrong
  # catch %>%
  #   select(Time:BIV) %>%
  #   pivot_longer(-Time, names_to = "Code", values_to = "mt") %>%
  #   filter(Time > 0) %>%
  #   ggplot()+
  #   geom_line(aes(x = Time, y = mt))+
  #   geom_vline(xintercept = 30*365, color = "red")+
  #   facet_wrap(~Code, scales = "free")
  
  # now extract data
  #biodat_age_tmp <- read.table('../Parametrization/output_files/data/out_1342/outputGOA01342_testAgeBiomIndx.txt', sep = ' ', header = T)
  # SSB to plot and report in tables
  spawning_biomass <- biomage %>% 
    slice_tail(n = 5) %>% # use last xxx years
    summarise(across(-"Time", ~ mean(.x, na.rm = TRUE))) %>%
    ungroup() %>%
    pivot_longer(everything(), names_to = 'Code.Age', values_to = 'biomass_mt') %>%
    separate_wider_delim(Code.Age, delim = '.', names = c('Code', 'Age')) %>%
    filter(Code == sp) %>%
    #left_join(mat, by = 'Code') %>%
    # mutate(idx = as.numeric(Age) - as.numeric(age_class_mat)) %>%
    # filter(is.na(idx) | idx >= 0) %>%
    left_join(fspb, by = c('Code','Age')) %>%
    mutate(biomass_mt = biomass_mt * fspb) %>%
    group_by(Code) %>%
    summarise(biomass_mt = sum(biomass_mt)) %>%
    ungroup() %>%
    pull(biomass_mt)
  
  # total catch
  # taking mean of the last 5 years
  catch_val <- catch %>% 
    slice_tail(n = 5) %>% 
    pull(sp) %>%
    mean()
  
  # # calculate realized F after 1 year of data
  # For runs with a burn-in, this has to be the biomass at the end of the burn-in, when we start fishing with the new scalar
  # # get initial biomass for the selected age classes
  biom_age_t1 <- biomage %>% 
    filter(Time == 365 * burnin) %>%# this is the burn-in years
    pivot_longer(-Time, names_to = 'Code.Age', values_to = 'biomass') %>%
    separate_wider_delim(Code.Age, delim = '.', names = c('Code', 'Age')) %>%
    left_join(selex, by = 'Code') %>%
    mutate(idx = as.numeric(Age) - as.numeric(age_class_selex)) %>%
    filter(is.na(idx) | idx >= 0) %>%
    group_by(Code) %>%
    summarise(biomass = sum(biomass)) %>%
    ungroup() %>% 
    filter(Code == sp)
  # 
  # # catch (one time step after biomass: how much did we catch in this time?)
  catch_t1 <- catch %>% 
    select(Time, all_of(sp)) %>% 
    filter(Time == 365 * (burnin + 1)) %>% # careful - there is a small transition phase
    summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) %>%
    pivot_longer(-Time, names_to = 'Code', values_to = 'catch') %>%
    select(-Time)
  # 
  # # calc realized f
  f_t1 <- biom_age_t1 %>% left_join(catch_t1) %>%
    mutate(exp_rate = catch/biomass,
           f = -log(1-exp_rate),
           fidx = fidx) %>% # need this for joining later on
    select(Code, f, fidx) 
  
  # # bind all
  f_frame <- f_t1 %>%
    mutate(biomass = spawning_biomass,
           catch = catch_val)
  
  # write out to be then brought together with all other runs
  write.csv(f_frame, paste(outdir,paste(sp,fidx,'f.csv',sep='_'), sep = "/"), row.names = F)
  
}
