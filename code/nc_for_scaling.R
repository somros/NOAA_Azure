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
library(ggh4x)
library(viridis)
library(tidync)
library(ncdf4)

# get example netcdf from base
this_nc_file <- "../Parametrization/output_files/data/out_1517/outputGOA01517_test.nc"
this_tidync <- tidync(this_nc_file)
this_nc <- ncdf4::nc_open(this_nc_file)

# read in groups
grps <- read.csv('data/GOA_Groups.csv', header = T)

# set up a functional group types table
vertebrate_groups <- grps %>% filter(GroupType%in%c("FISH","SHARK","BIRD","MAMMAL")) %>% mutate(BiomassType="vertebrate")
plankton_groups <- grps %>% filter(GroupType %in% c("PWN",'CEP','LG_ZOO','MED_ZOO','SM_ZOO','LG_PHY','SM_PHY')) %>%
  mutate(BiomassType="plankton")
bottom_groups <- grps %>% filter(GroupType %in% c("MOB_EP_OTHER",'SED_EP_FF','SED_EP_OTHER','PHYTOBEN')) %>%
  mutate(BiomassType="2D")
other_groups <- grps %>% filter(GroupType %in% c("LG_INF","MICROPHTYBENTHOS","SED_BACT","PL_BACT","SM_INF","CARRION","LAB_DET","REF_DET"))%>%
  mutate(BiomassType="other")
biomass_groups <- bind_rows(vertebrate_groups,plankton_groups,bottom_groups,other_groups)

# add to grps df
grps <- grps %>% left_join(biomass_groups)

selex <- read.csv("data/age_at_selex_new.csv", header = T) # age at selectivity

all_fg <- grps %>% filter(GroupType %in% c("FISH", "SHARK", "MAMMAL", "BIRD")) %>% pull(Name)


biom_spatial_ls <- list()

for(n in 1:length(all_fg)){
  
  # args for the function below
  fg <- all_fg[n] # this needs to use the "Name" to pull from the NC file
  out <- this_tidync
  this.nc <- this_nc
  run <- 1
  
  fg_atts <- grps %>% filter(Name==fg)
  
  if(fg_atts$BiomassType!="vertebrate") stop("weight at age only for vertebrates.")
  
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
    biom <- list()
    for(i in 1:length(resN)){
      biom[[i]] <- (resN[[i]] + strucN[[i]]) * nums[[i]] * 20 * 5.7 / 1000000000 # go from mgN to tons
    }
    
    # now collapse over depth
    biom_box <- lapply(biom, function(x) apply(x, c(2, 3), sum))
    
    # now turn to data frame
    biom_box_ls <- list()
    
    # Loop over each matrix in the collapsed list to fill the data frame
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
    
    biom_box_df <- bind_rows(biom_box_ls)
    
    # write this out for comparison with catch
    biom_box_df <- biom_box_df %>% 
      mutate(Name = fg)
    
    biom_spatial_ls[[n]] <- biom_box_df
    
    # to avoid immense tables when we do it for all species (and once per run), work out the proportions here
    # do it by age class first
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
    
    # now add name of the species
    
  }
  
}

biom_spatial_df <- bind_rows(biom_spatial_ls)


# now compare spatial dist of biomass to spatial distribution of the catch
# use the TOTCATCH.nc file for the comparison
# script check_catch_outputs.R demonstrates that the 3 outputs for catch are rather comparable

# nc tot file
catch_nc_tot_file <- "../Parametrization/output_files/data/out_1517/outputGOA01517_testTOTCATCH.nc"
this_tidync <- tidync(catch_nc_tot_file)
this_nc <- ncdf4::nc_open(catch_nc_tot_file)

all_fg <- grps %>% filter(GroupType %in% c("FISH", "SHARK", "MAMMAL", "BIRD")) %>% pull(Code)

catch_TOTCATCH <- list()

for(i in 1:length(all_fg)){
  fg <- all_fg[i] # this needs to use the "Code" to pull from the NC file
  out <- this_tidync
  this.nc <- this_nc
  run <- 1
  
  fg_atts <- grps %>% filter(Name==fg)
  
  #Extract from the output .nc file the appropriate reserve N time series variables
  catch_vars <- out %>%
    activate("D1,D0") %>%
    hyper_vars() %>% # all variables in the .nc file active grid
    filter(grepl("^Tot_.*_Catch$",name)) %>% # filter for reserve N
    filter(grepl(fg,name)) # filter for specific functional group
  
  # # Actually pull the data from the .nc
  # here we can collapse the depth layers but need to keep the boxes
  catch <- purrr::map(catch_vars$name,ncdf4::ncvar_get,nc=this.nc)[[1]] # only one variable per group becasue no ages
  
  # now turn to data frame
  # Convert matrix to a data frame
  mat_df <- as.data.frame((catch))
  
  # add box_id
  mat_df <- mat_df %>% mutate(box_id = 0:108)
  
  # reshape
  mat_df_long <- mat_df %>%
    pivot_longer(-box_id, names_to = "ts", values_to = "mt")
  
  # turn ts column to integer
  catch_box_df <- mat_df_long %>%
    mutate(ts = gsub("V","",ts)) %>%
    mutate(ts = as.numeric(ts)) %>%
    mutate(ts = ts - 1) %>% # start numbering ts from 0
    mutate(Code = fg)
  
  catch_TOTCATCH[[i]] <- catch_box_df
  
}

catch_TOTCATCH <- bind_rows(catch_TOTCATCH)

# bring in selectivity, because juvenile biomass will have the distribution of the adults mostly
# flatten biomass over age classes
biom_spatial_df_noage <- biom_spatial_df %>%
  left_join(grps %>% select(Code, Name), by = "Name") %>%
  left_join(selex, by = "Code") %>%
  mutate(idx = age - age_class_selex) %>%
  filter(idx >= 0) %>% # keep only selected age classes
  group_by(ts, box_id, Code) %>%
  summarise(mt = sum(mt)) %>%
  filter(ts > 0) %>%
  filter(ts %in% seq(5,250,5)) %>% # keep every 5th time step
  mutate(ts = ts / 5)

# now get proportions, which is what we need to compare
biom_props <- biom_spatial_df_noage %>%
  group_by(ts, Code) %>%
  mutate(tot_mt = sum(mt)) %>%
  ungroup() %>%
  mutate(prop_biom = mt / tot_mt)

catch_props <- catch_TOTCATCH %>%
  filter(ts > 0) %>%
  group_by(ts, Code) %>%
  mutate(tot_mt = sum(mt)) %>%
  ungroup() %>%
  mutate(prop_catch = mt / tot_mt)

# only do t3
comp <- biom_props %>%
  left_join(catch_props, by = c("ts","box_id","Code")) %>%
  filter(Code %in% c("POL","COD","ATF","POP","SBF","FFS","FHS","FFD","REX","RFS","RFP","HAL")) %>%
  mutate(comp = prop_biom / prop_catch)

hist(comp$comp)
# this looks pretty good for the stocks we are interested in - props of biomass trach with props of catch

# now calculate prop in AK for both
comp_ak <- comp %>%
  filter(Code %in% c("POL","COD","ATF","POP","SBF","FFS","FHS","FFD","REX","RFS","RFP","HAL")) %>%
  filter(box_id < 92) %>%
  group_by(ts, Code) %>%
  summarise(prop_biom_ak = sum(prop_biom, na.rm = T),
            prop_catch_ak = sum(prop_catch, na.rm = T)) %>%
  pivot_longer(-c(ts,Code), names_to = "type", values_to = "prop")

# view
comp_ak %>%
  ggplot(aes(x = ts, y = prop, color = type))+
  geom_line()+
  facet_wrap(~Code, scales = "free")



