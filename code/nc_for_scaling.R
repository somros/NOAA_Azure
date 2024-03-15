# code to copy in the appropriate pipeline
# Open out.nc file
# get WAA and NAA of vertebrates by age class by box by time step
# from that, get biomass
# Make two frames from this:
# 1. One with biomass by age class and resulting proportions for AK and BC to scale the biomage.txt output
# 2. One with aggregate biomass of selected age classes and resulting proportions for AK and BC to scale the catch.txt output which is not by age class
# These will then used to scale the respective CSV files every time that thse are called for analyses
# Why does this work? Because mFC extracts catch proportionally to biomass and there are no further limitations to fishing
# this will need to be done on each netcdf file, i.e. there will be a scaling factor for each run
# they can all be combined into one df with an index that is then fileterd and applied for each set of CSV outputs


library(tidyverse)
library(here)
library(tidyr)
library(readxl)
library(viridis)
library(tidync)
library(ncdf4)

# get example netcdf from base
this_nc_file <- "../Parametrization/output_files/data/out_1517/outputGOA01517_test.nc"
this_tidync <- tidync(this_nc_file)
this_nc <- ncdf4::nc_open(this_nc_file)

grps <- read.csv('data/GOA_Groups.csv', header = T) # read in groups
selex <- read.csv("data/age_at_selex_new.csv", header = T) # age at selectivity

# stocks of interest
t3 <- c("POL","COD","ATF","POP","SBF","FFS","FHS","FFD","REX","RFS","RFP","HAL")

all_fg <- grps %>% filter(Code %in% t3) %>% pull(Name) # only focus on vertebrates

# make an empty list to populate with biomass for each group
catch_prop_ls <- list()

for(n in 1:length(all_fg)){
  
  # args for the function below
  fg <- all_fg[n] # this needs to use the "Name" to pull from the NC file
  out <- this_tidync
  this.nc <- this_nc

  #Extract from the output .nc file the appropriate reserve N time series variables
  resN_vars <- hyper_vars(out) %>% # all variables in the .nc file active grid
    filter(grepl("_ResN",name)) %>% # filter for reserve N
    filter(grepl(fg,name)) # filter for specific functional group
  
  #Extract from the output .nc file the appropriate structural N time series variables
  strucN_vars <- hyper_vars(out) %>% # all variables in the .nc file active grid
    filter(grepl("_StructN",name)) %>% # filter for structural N
    filter(grepl(fg,name)) # filter for specific functional group
  
  # Get numbers by box
  abun_vars <- hyper_vars(out) %>% # all variables in the .nc file active grid
    filter(grepl("_Nums",name)) %>% # filter for abundance variables
    filter(grepl(fg,name)) # filter for specific functional group
  
  if(nrow(resN_vars)==0) {return("no data.")}
  else {
    # # Actually pull the data from the .nc
    # here we can collapse the depth layers but need to keep the boxes
    resN <- purrr::map(resN_vars$name,ncdf4::ncvar_get,nc=this.nc) 
    strucN <- purrr::map(strucN_vars$name,ncdf4::ncvar_get,nc=this.nc)
    nums <-purrr::map(abun_vars$name,ncdf4::ncvar_get,nc=this.nc) #numbers by age group,box,layer,time
    
    # Formula for biomass in each cell will be (resN+structN)*nums
    # this step maintains the spatial and temporal structure of the nc arrays
    biom <- list()
    for(i in 1:length(resN)){
      biom[[i]] <- (resN[[i]] + strucN[[i]]) * nums[[i]] * 20 * 5.7 / 1000000000 # go from mgN to tons
    }
    
    # now collapse over depth
    biom_box <- lapply(biom, function(x) apply(x, c(2, 3), sum))
    
    # now turn to data frame
    # create empty list
    biom_box_ls <- list()
    
    # Loop over each matrix in the collapsed list to fill the data frame (i.e. over each age class)
    for(i in 1:length(biom_box)) {
      # Get the current matrix
      mat <- biom_box[[i]]
      
      # Convert matrix to a data frame
      mat_df <- as.data.frame((mat))
      
      # add box_id
      mat_df <- mat_df %>% mutate(box_id = 0:108)
      
      # reshape
      mat_df_long <- mat_df %>%
        pivot_longer(-box_id, names_to = "ts", values_to = "mt")
      
      # turn ts column to integer
      mat_df_long <- mat_df_long %>%
        mutate(ts = gsub("V","",ts)) %>%
        mutate(ts = as.numeric(ts)) %>%
        mutate(ts = ts - 1) # start numbering ts from 0
      
      # add age
      mat_df_long <- mat_df_long %>%
        mutate(age = i-1) # number from 0 for consistency with age at selex and age mat
      
      biom_box_ls[[i]] <- mat_df_long
      
    }
    
    # turn list to data frame
    biom_box_df <- bind_rows(biom_box_ls)
    
    # to avoid immense tables when we do it for all species (and once per run), work out the proportions here
    # do it by age class first
    # this will be used to scale the biomass of each age class
    # I am not even convinced we need this - if we are keeping everything to total population level except for the global yield plot
    # 92 is the first box in BC
    # biom_ak_prop <- biom_box_df %>%
    #   mutate(area = ifelse(box_id < 92, "ak", "bc")) %>%
    #   group_by(ts, age) %>% # get total biomass by time step by age across the model domain
    #   mutate(mt_tot = sum(mt)) %>%
    #   group_by(ts, age, area) %>%
    #   mutate(mt_area = sum(mt)) %>%
    #   ungroup() %>%
    #   mutate(prop = mt_area / mt_tot) %>%
    #   select(ts, age, area, prop) %>%
    #   distinct()
    
    # view:
    # biom_ak_prop %>%
    #   filter(area == "ak") %>%
    #   ggplot(aes(x = ts, y = prop, color = factor(age))) +
    #   geom_line()
    
    # now bring in selex, filter at or above it, sum, and get prop for total
    catch_prop_from_ak <- biom_box_df %>%
      mutate(area = ifelse(box_id < 92, "ak", "bc")) %>% # define areas based on boxes
      mutate(Name = fg) %>% # add name
      left_join(grps %>% dplyr::select(Code, Name), by = "Name") %>% 
      left_join(selex, by = "Code") %>% # bring in age at selectivity (which counts from 0)
      mutate(idx = age - age_class_selex) %>% # is age >= age at selectivity?
      filter(idx >= 0) %>% # keep only age classes >= age at selex
      group_by(ts) %>% # get total biomass by time step across the model domain (aggregate age classes)
      mutate(mt_tot = sum(mt)) %>%
      group_by(ts, area) %>%
      mutate(mt_area = sum(mt)) %>%
      ungroup() %>%
      mutate(prop = mt_area / mt_tot) %>%
      select(ts, Name, area, prop) %>%
      distinct()
    
    # handle time:
    # drop t0 for easier filtering
    # subsample at every 5 time steps to have annual values
    # keep the last 5 of the series only
    # average
    # produce one value per run per species
    
    catch_prop_from_ak <- catch_prop_from_ak %>%
      filter(area == "ak") %>% # keep proportion in AK only
      filter(ts > 0) %>% # drop t0
      filter(ts %in% seq(5,250,5)) %>% # resample to have annual time steps instead of 73 days
      mutate(ts = ts / 5) %>% # reindex the time step accordingly 
      slice_tail(n = 5) %>% # keep end of the run for consistency with the catch sampling
      group_by(Name) %>%
      summarize(ak_prop = mean(prop))
    
  }
  
  catch_prop_ls[[n]] <- catch_prop_from_ak
  
}

# now bind all species
catch_prop_df <- bind_rows(catch_prop_ls)

# output this data frame from each run (how do we index it?)






