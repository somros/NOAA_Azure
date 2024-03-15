# this script checks that 3 output files for catch contain the same information:
# Catch.txt, CATCH.nc, and TOTCATCH.nc

catch_flat_file <- "../Parametrization/output_files/data/out_1517/outputGOA01517_testCatch.txt"
catch_nc_file <- "../Parametrization/output_files/data/out_1517/outputGOA01517_testCATCH.nc"
catch_nc_tot_file <- "../Parametrization/output_files/data/out_1517/outputGOA01517_testTOTCATCH.nc"

# read in files and make tables
# flat file
catch_flat <- read.delim(catch_flat_file, sep = " ")

# nc tot file
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

# compare flat text to TOTCATCH.nc
catch_flat_long <- catch_flat %>%
  pivot_longer(-Time, names_to = "Code", values_to = "mt") %>%
  filter(Code %in% all_fg) %>%
  mutate(Time = Time / 365) %>%
  rename(ts = Time) %>%
  mutate(file = "flat")

# collapse boxes
catch_TOTCATCH_flattened <- catch_TOTCATCH %>%
  group_by(ts, Code) %>%
  summarise(mt = sum(mt)) %>%
  mutate(file = "tot_nc")

comp <- catch_flat_long %>%
  left_join(catch_TOTCATCH_flattened, by = c("ts", "Code")) %>%
  mutate(comp = mt.x/mt.y)

min(comp$comp, na.rm = T)
max(comp$comp, na.rm = T)
mean(comp$comp, na.rm = T)

# flat and TOTCATCH.nc are the same
# now let's bring in the long nc files. We are going to need the bio.nc files here too, because catches are reported in numbers
# so, you need to: read in the WAA, sum them up, collapse by depth (weighting by numbers...), sample every 5th on them, multiply by the catches in numbers
# then we need to compare this, in space, to the TOTCATCH.
# Finally, we can take TOTCATCH.nc and use it to compare to biomass distributions in the biomass nc file

# nc files
bio_nc_file <- "../Parametrization/output_files/data/out_1517/outputGOA01517_test.nc"

# catch
this_tidync <- tidync(catch_nc_file)
this_nc <- ncdf4::nc_open(catch_nc_file)
#bio
this_tidync_bio <- tidync(bio_nc_file)
this_nc_bio <- ncdf4::nc_open(bio_nc_file)

all_fg <- grps %>% filter(GroupType %in% c("FISH", "SHARK", "MAMMAL", "BIRD")) %>% pull(Name)

catch_nc_ls <- list()

for(n in 1:length(all_fg)){
  fg <- all_fg[n] # this needs to use the "Name" to pull from the NC file
  out <- this_tidync
  this.nc <- this_nc
  out_bio <- this_tidync_bio
  this.nc_bio <- this_nc_bio
  
  fg_atts <- grps %>% filter(Name==fg)
  
  if(fg_atts$BiomassType!="vertebrate") stop("weight at age only for vertebrates.")
  
  #Extract from the output .nc file the appropriate catch time series variables
  catch_vars <- out %>%
    activate("D1,D0") %>%
    hyper_vars() %>% # all variables in the .nc file active grid
    filter(grepl("_Catch",name)) %>% # filter for reserve N
    filter(grepl(fg,name)) # filter for specific functional group
  
  #Extract from the output .nc file the appropriate reserve N time series variables
  resN_vars <- hyper_vars(out_bio) %>% # all variables in the .nc file active grid
    filter(grepl("_ResN",name)) %>% # filter for reserve N
    filter(grepl(fg,name)) # filter for specific functional group
  
  #Extract from the output .nc file the appropriate structural N time series variables
  strucN_vars <- hyper_vars(out_bio) %>% # all variables in the .nc file active grid
    filter(grepl("_StructN",name)) %>% # filter for structural N
    filter(grepl(fg,name)) # filter for specific functional group
  
  # Get numbers by box
  # going to need these for weighting over depth layers (which we need to average over because catch is by box and not by depth)
  abun_vars <- hyper_vars(out_bio) %>% # all variables in the .nc file active grid
    filter(grepl("_Nums",name)) %>% # filter for abundance variables
    filter(grepl(fg,name)) # filter for specific functional group
  
  
  if(nrow(resN_vars)==0) {return("no data.")}
  else {
    # # Actually pull the data from the .nc
    # start from bio information
    resN <- purrr::map(resN_vars$name,ncdf4::ncvar_get,nc=this.nc_bio) 
    strucN <- purrr::map(strucN_vars$name,ncdf4::ncvar_get,nc=this.nc_bio)
    nums <-purrr::map(abun_vars$name,ncdf4::ncvar_get,nc=this.nc_bio) #numbers by age group,box,layer,time
    
    # sum res and struct and then take averages over the water column, ending up with a value per box per time step
    weighted_avg_list <- vector("list", length = length(resN)) # Initialize the list to store the weighted averages
    
    # Perform the operation (A+B) and then calculate the weighted average using C
    weighted_mean_func <- function(x, w) {
      sum(x * w, na.rm = TRUE) / sum(w, na.rm = TRUE)
    }
    
    weighted_waa <- list()
    for(i in 1:length(resN)) {
      # Calculate the sum of A and B
      sum_ab <- (resN[[i]] + strucN[[i]]) * 20 * 5.7 / 1000000000 # from from mgN to tons
      this_nums <- nums[[i]] # numbers at this age over the spatial structure
      
      # Initialize an array to store the results
      this_weighted_waa <- array(dim = c(109, 251))
      
      # Calculate weighted mean across the first dimension
      for(j in 1:109) {
        for(k in 1:251) {
          # Extract the slice for the current i, j across all 7 values
          sum_ab_slice <- sum_ab[, j, k]
          nums_slice <- this_nums[, j, k]
          
          # Calculate the weighted mean for this slice
          this_weighted_waa[j, k] <- weighted_mean_func(sum_ab_slice, nums_slice)
        }
      }
      
      # turn to df and pivot
      this_weighted_waa <- this_weighted_waa %>%
        as.data.frame() %>%
        mutate(box_id = 0:108) %>%
        pivot_longer(-box_id, values_to = "mt", names_to = "ts") %>%
        mutate(ts = gsub("V","",ts)) %>%
        mutate(ts = as.numeric(ts)) %>%
        mutate(ts = ts - 1) %>%
        mutate(age = i - 1) # count from 0
      
      weighted_waa[[i]] <- this_weighted_waa
    }
    
    # turn it ito a df
    weighted_waa_df <- bind_rows(weighted_waa)
    
    # drop t 0 and then take every 5th 
    weighted_waa_df <- weighted_waa_df %>%
      filter(ts > 0) %>%
      filter(ts %in% seq(5,250,5)) %>%
      mutate(ts = ts / 5)
    
    # now let's extract the catch for this fg
    # here we can collapse the depth layers but need to keep the boxes
    catch <- purrr::map(catch_vars$name,ncdf4::ncvar_get,nc=this.nc) 
    
    # now turn to data frame
    catch_box_ls <- list()
    
    # Loop over each matrix in the collapsed list to fill the data frame
    for(i in 1:length(catch)) {
      # Get the current matrix
      mat <- catch[[i]]
      
      # Convert matrix to a data frame
      mat_df <- as.data.frame((mat))
      
      # add box_id
      mat_df <- mat_df %>% mutate(box_id = 0:108)
      
      # reshape
      mat_df_long <- mat_df %>%
        pivot_longer(-box_id, names_to = "ts", values_to = "nums")
      
      # turn ts column to integer
      mat_df_long <- mat_df_long %>%
        mutate(ts = gsub("V","",ts)) %>%
        mutate(ts = as.numeric(ts)) %>%
        mutate(ts = ts - 1) # start numbering ts from 0
      
      # add age
      mat_df_long <- mat_df_long %>%
        mutate(age = i-1) # number from 0 for consistency with age at selex and age mat
      
      catch_box_ls[[i]] <- mat_df_long
      
    }
    
    catch_box_df <- bind_rows(catch_box_ls) # this is in numbers... so we'd need to multiply this by WAA, pulling that in
    
    # drop ts = 0
    catch_box_df <- catch_box_df %>%
      filter(ts > 0)
    
    # now join and multiply to get mt per box
    catch_box_df <- catch_box_df %>%
      left_join(weighted_waa_df) %>%
      mutate(mt_tot = nums * mt) %>%
      mutate(Name = fg)
    
  }
  
  catch_nc_ls[[n]] <- catch_box_df
}

catch_nc <- bind_rows(catch_nc_ls)

# turn NaN to NA
catch_nc$mt[is.nan(catch_nc$mt)] <- NA
catch_nc$mt_tot[is.nan(catch_nc$mt_tot)] <- NA

# collapse ages
catch_nc_noage <- catch_nc %>%
  group_by(ts, box_id, Name) %>%
  summarise(mt_tot = sum(mt_tot, na.rm = T)) %>%
  left_join(grps %>% select(Code, Name))

# compare to totcatch
comp_catch_nc <- catch_TOTCATCH %>%
  filter(ts > 0) %>%
  left_join(catch_nc_noage, by = c("ts", "box_id", "Code")) %>%
  filter(Code %in% setdiff(unique(catch_TOTCATCH$Code), c("SCH","SCO","SPI","SCM","SSO","HAK"))) %>% # drop migrating species, reporting gets wonky
  filter(mt > 0) %>% # keep only non-zero catches from tot nc
  mutate(comp = mt / mt_tot)

max(comp_catch_nc$comp)
min(comp_catch_nc$comp)
mean(comp_catch_nc$comp)

hist(comp_catch_nc$comp)

# the following step will be to check that one of the spatially-explicit ones (TOTCATCH.nc) has the same spatial distributions as biomass in out.nc
