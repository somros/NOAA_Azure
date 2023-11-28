# Alberto Rovellini
# 11/13/2023
# This code takes output of the F35 permutation runs (AgeBiomass and Catch) and:
# 1. Plots total yield per value of the permutation, with lines for FMP OY and internal "Atlantis" OY
# 2. Plots biomass and yield curves next to the old ones
# 3. Plots a number of diagnostics (biomass of unselected age classes, biomass of predators and forage fish, yield functions, etc.)

library(tidyverse)
library(readxl)
library(ggh4x)
library(viridis)

# Set up env and read data ------------------------------------------------

# identify which data we want to work on
batch <- "maxF_4f35"
this_job <- "job20231103012611" # and which job to use for the single-species stuff
runs <- 1448:1458
maxmult <- 4

# set the clock to date plots
t <- format(Sys.time(),'%Y-%m-%d %H-%M-%S')

# read in Groups.csv file
grp_path <- here('NOAA_Azure/data/GOA_Groups.csv') # functional groups
grps <- read.csv(grp_path)

# read maturity and selectivity information
mat <- read.csv("NOAA_Azure/data/age_at_mat.csv", header = T) # age at maturity
selex <- read.csv("NOAA_Azure/data/age_at_selex_new.csv", header = T) # age at selectivity

# read in lookup tables and f35 proxy values
f_lookup <- read.csv("NOAA_Azure/data/f_lookup.csv")
f35_vector <- read.csv("NOAA_Azure/data/f35_vector_PROXY.csv")

# list Tier 3 stocks 
t3_fg <- f_lookup %>% pull(species) %>% unique()

# get the multispecies files
f35_path <- paste0("NOAA_Azure/results/f35/",batch)
f35_biomass_files <- list.files(f35_path, pattern = "AgeBiom", full.names = T)
f35_catch_files <- list.files(f35_path, pattern = "Catch", full.names = T)

# Calculate global yield from F35 permutation runs ------------------------

# make a lookup table (should find a better naming system for the runs)
f35_key <- data.frame("run" = runs, mult = seq(0,maxmult,length.out=11))

# list
catch_list <- list()
for(i in 1:nrow(f35_key)){
  
  this_catch_file <- f35_catch_files[i]
  this_run <- as.numeric(substr((gsub(paste0(f35_path,"/outputGOA0"), "", this_catch_file)), 1, 4))
  this_mult <- f35_key %>% filter(run == this_run) %>% pull(mult)
  
  this_catch <- read.table(this_catch_file, sep=" ", header = T)
  this_catch <- this_catch %>%
    slice_tail(n = 5) %>%
    summarise(across(all_of(t3_fg), ~mean(.x, na.rm = T))) %>%
    mutate(mult = this_mult)
  
  catch_list[[i]] <- this_catch
}

catch_df <- bind_rows(catch_list)

# reshape and calculate total
catch_df_long <- catch_df %>%
  pivot_longer(-mult, names_to = "Code", values_to = "mt") %>%
  group_by(mult) %>%
  mutate(total_yield = sum(mt),
         prop = mt / total_yield) %>%
  ungroup() %>%
  left_join(grps %>% select(Code, LongName), by = "Code")

# make Atlantis-derived OY: sum single-species MSY estimates and discount by 8%
# atlantis_oy <- f_df %>%
#   filter(type == "Catch") %>%
#   group_by(Code) %>%
#   slice_max(mt) %>%
#   ungroup() %>%
#   pull(mt) %>%
#   sum()

# discount it by 8%
# atlantis_oy <- atlantis_oy - 0.08*atlantis_oy # this is REALLY close to OY
  
# plot
catch_df_long %>%
  ggplot(aes(x = mult, y = total_yield))+
  geom_point()+
  geom_line()+
  theme_bw()

# at 1.2 F35 vector we have maximum yield

global_yield_plot <- catch_df_long %>%
  ggplot(aes(x = mult, y = mt / 1000, fill = LongName))+
  geom_bar(stat = "identity", position = "stack")+
  scale_fill_viridis_d()+
  scale_x_continuous(breaks = seq(0,maxmult,length.out=11))+
  geom_hline(yintercept = 116, color = "red", linetype = "dashed")+
  geom_hline(yintercept = 800, color = "red", linetype = "dashed")+
  #geom_hline(yintercept = atlantis_oy / 1000, color = "blue", linetype = "dashed")+
  theme_bw()+
  labs(x = "F35 multiplier", y = "Catch (1000 mt)", fill = "Stock")
global_yield_plot

ggsave(paste0("NOAA_Azure/results/figures/", batch, "/global_yield.png"), global_yield_plot, width = 8, height = 6)
  
# proportions remain fairly constant over time. This is not a good sign, trophic interactions are really weak
# negligible amounts of POP (which is bad in the GOA, it is a large fishery)
# shooting past the real-world cap at f35 (never happened).
# need to calculate the Atlantis cap

# get single species files
# now we need to compare the yield curves obtained with this method to the single-species curves
# tricy part is what goes on the x-axis
# need to calculate F from catch and spawning stock biomass like before. Similar process but need to do it for all species at the same time
# TODO: it can become a function with the f processing code

ms_yield_list <- list()

for(i in 1:nrow(f35_key)){
  
  this_catch_file <- f35_catch_files[i]
  this_run <- as.numeric(substr((gsub(paste0(f35_path,"/outputGOA0"), "", this_catch_file)), 1, 4))
  this_mult <- f35_key %>% filter(run == this_run) %>% pull(mult)
  
  # extract tables from results
  biomage <- read.table(f35_biomass_files[i], sep = " ", header = T)
  catch <- read.table(this_catch_file , sep = " ", header = T)
  
  # now extract data
  # SSB to plot and report in tables
  spawning_biomass <- biomage %>% 
    slice_tail(n = 5) %>% # use last xxx years
    summarise(across(-"Time", ~ mean(.x, na.rm = TRUE))) %>%
    ungroup() %>%
    pivot_longer(everything(), names_to = 'Code.Age', values_to = 'biomass_mt') %>%
    separate_wider_delim(Code.Age, delim = '.', names = c('Code', 'Age')) %>%
    filter(Code %in% t3_fg) %>%
    left_join(mat, by = 'Code') %>%
    mutate(idx = as.numeric(Age) - as.numeric(age_class_mat)) %>%
    filter(is.na(idx) | idx >= 0) %>%
    group_by(Code) %>%
    summarise(biomass_mt = sum(biomass_mt)) %>%
    ungroup()
  
  # total catch
  # taking mean of the last 5 years
  catch_vals <- catch %>% 
    slice_tail(n = 5) %>%
    summarise(across(all_of(t3_fg), ~mean(.x, na.rm = T))) %>%
    pivot_longer(cols = everything(), names_to = "Code", values_to = "catch_mt")
  
  # # calculate realized F after 1 year of data
  # # get initial biomass for the selected age classes
  biom_age_t1 <- biomage %>% 
    filter(Time == 0) %>% 
    pivot_longer(-Time, names_to = 'Code.Age', values_to = 'biomass') %>%
    separate_wider_delim(Code.Age, delim = '.', names = c('Code', 'Age')) %>%
    left_join(selex, by = 'Code') %>%
    mutate(idx = as.numeric(Age) - as.numeric(age_class_selex)) %>%
    filter(is.na(idx) | idx >= 0) %>%
    group_by(Code) %>%
    summarise(biomass = sum(biomass)) %>%
    ungroup() %>% 
    filter(Code %in% t3_fg)
  # 
  # # catch (one time step after biomass: how much did we catch in this time?)
  catch_t1 <- catch %>% 
    select(Time, all_of(t3_fg)) %>% 
    filter(Time == 365) %>%
    summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) %>%
    pivot_longer(-Time, names_to = 'Code', values_to = 'catch') %>%
    select(-Time)
  # 
  # # calc realized f
  f_t1 <- biom_age_t1 %>% left_join(catch_t1, by = "Code") %>%
    mutate(exp_rate = catch/biomass,
           f = -log(1-exp_rate)) %>%#,
           #fidx = fidx) %>% # need this for joining later on
    select(Code, f)#, fidx) 
  
  # # bind all
  f_frame <- f_t1 %>%
    left_join(spawning_biomass) %>%
    left_join(catch_vals) %>%
    mutate(mult = this_mult)
  
  # add to multispecies yield list
  ms_yield_list[[i]] <- f_frame
}

ms_yield_df <- bind_rows(ms_yield_list)

# add species long name and reshape catch and biomass and add a label that this is the ms approach
ms_yield_long <- df_mult <- ms_yield_df %>% # one of these is for later plots that need mult
  left_join(grps %>% select(Code, LongName), by = "Code") %>%
  rename(Biomass = biomass_mt, Catch = catch_mt) %>%
  pivot_longer(cols = -c(Code, LongName, f, mult), names_to = "type", values_to = "mt") %>%
  mutate(experiment = "ms") %>%
  select(Code, LongName, experiment, f, mult, type, mt)
  
# drop mult
ms_yield_long <- ms_yield_long %>%
  dplyr::select(-mult)


# Bring in the single-species experiment (96 runs) ------------------------

# now create a similar object from the single-species f experiments

# list all csv files we need to read
f_files <- list.files(paste0("NOAA_Azure/results/post-processing/",this_job), full.names = T)

# create empty list to fill with data frame for the yield curve
f_df_ls <- list()

for(i in 1:length(f_files)){
  
  this_f_files <- f_files[[i]]
  
  # read all csv files
  f_ls <- list()
  for(j in 1:length(this_f_files)){
    this_file <- this_f_files[j]
    f_ls[[j]] <- read.csv(this_file)
  }
  
  # bind into a data frame
  f_df <- f_ls %>% bind_rows() %>% rename(Biomass = biomass, Catch = catch)
  
  # clean up and format
  f_df <- f_df %>%
    pivot_longer(-c(Code, f, fidx), values_to = 'mt', names_to = 'type') %>%
    left_join(grps %>% select(Code, LongName), by = 'Code')
  
  f_df_ls[[i]] <- f_df
  
}

f_df <- bind_rows(f_df_ls) 

ss_yield_long <- f_df %>%
  mutate(experiment = "ss") %>%
  select(Code, LongName, f, experiment, type, mt)

# produce a dataset of 35% B0, to be used to plot horizontal lines that will intersect the yield curve
# but B35% will now be different between runs with different steepness? Not between runs with smaller selex in theory
b35 <- f_df %>%
  filter(f == 0) %>% # producing it off of v1 (v2 should be the same, v3 should be similar)
  rowwise() %>%
  mutate(b35 = ifelse(type== 'Biomass', mt * 0.35, NA)) %>%
  ungroup() %>%
  select(LongName, Code, type, b35)

# read in MSY information (from FMP)
tier3 <- read_xlsx('NOAA_Azure/data/msy.xlsx', sheet = 1, range = 'A3:J19') %>%
  select(Stock, FOFL) %>%
  set_names(c('Stock', 'FMSY'))

tier4_5 <- read_xlsx('NOAA_Azure/data/msy.xlsx', sheet = 2, range = 'A3:I10') %>%
  select(`Stock/Stock complex`, `M or FMSY`)%>%
  set_names(c('Stock', 'FMSY'))

tier_3_4_5 <- rbind(tier3, tier4_5)

# make key
tier_3_4_5 <- tier_3_4_5 %>%
  mutate(Code = c('POL','COD','SBF','FFS','FFS','FFS','FFS','FFD',
                  'REX','REX','ATF','FHS','POP','RFS','RFS','RFP',
                  'FFS','RFD','RFD','RFD','RFD','THO','DOG')) %>%
  group_by(Code) %>%
  summarise(FMSY = mean(FMSY))

# as soon as we introduce species that are not in the FMP, including forage species, we will need estimates of M from the parameters file
all_f <- tier_3_4_5 # placeholder for now

# find groups to plot
to_plot <- unique(f_df$Code)

# bind FMSY information
fmsy <- data.frame('Code' = to_plot) %>%
  left_join(all_f) %>%
  left_join(grps %>% select(Code, LongName))

# add halibut (M from IPHC assessment)
fmsy[fmsy$Code=='HAL',]$FMSY <- 0.2 # this is M

# get f that returned the highest yield, and level of depletion for that F
sp <- unique(f_df$LongName)

atlantis_fmsy_ls <- list()

for(i in 1:length(sp)){
  this_f_df <- f_df %>% filter(LongName == sp[i])
  
  atlantis_fmsy <- this_f_df %>% filter(type == 'Catch') %>%
    slice_max(mt) %>%
    pull(f)
  
  b0 <- this_f_df %>%
    filter(f == 0, type == 'Biomass') %>%
    pull(mt)
  
  b_fmsy <- this_f_df %>%
    filter(f == atlantis_fmsy, type == 'Biomass') %>%
    pull(mt)
  
  depletion_fmsy <- b_fmsy / b0
  
  fidx_fmsy <- this_f_df %>%
    filter(f == atlantis_fmsy) %>%
    pull(fidx) %>%
    unique()
  
  atlantis_fmsy_ls[[i]] <- data.frame('LongName' = sp[i], 
                                      'atlantis_fmsy' = atlantis_fmsy,
                                      'b_fmsy' = b_fmsy,
                                      'depletion' = depletion_fmsy,
                                      'fidx' = fidx_fmsy)
}

atlantis_fmsy <- bind_rows(atlantis_fmsy_ls)

# save this for future calculations
# write.csv(atlantis_fmsy, "NOAA_Azure/data/f35_vector_PROXY.csv", row.names = F)

# annotations for the plots (atlantis depletion)
annotations <- atlantis_fmsy %>% 
  mutate(depletion=round(depletion,digits=2), atlantis_fmsy=round(atlantis_fmsy,digits = 2))

# plot
to_plot <- ms_yield_long %>%
  rbind(ss_yield_long) %>%
  mutate(experiment = ifelse(experiment == "ms", "Multispecies", "Single-species"))

# spaces
to_plot$LongNamePlot <- gsub(" ", "\n", to_plot$LongName)
fmsy$LongNamePlot <- gsub(" ", "\n", fmsy$LongName)
atlantis_fmsy$LongNamePlot <- gsub(" ", "\n", atlantis_fmsy$LongName)
b35$LongNamePlot <- gsub(" ", "\n", b35$LongName)

f_plot <- to_plot %>%
  ggplot(aes(x = f, y = mt/1000, color = experiment))+
  geom_line()+
  geom_point(size = 2)+
  scale_color_viridis_d(begin = 0.2, end = 0.8)+
  geom_vline(data = fmsy, aes(xintercept = FMSY, group = LongNamePlot), linetype = 'dashed', color = 'orange')+
  geom_vline(data = atlantis_fmsy, aes(xintercept = atlantis_fmsy, group = LongNamePlot), linetype = 'dashed', color = 'blue')+
  geom_hline(data = b35 , aes(yintercept = b35/1000, group = LongNamePlot), linetype = 'dashed', color = 'red')+
  theme_bw()+
  scale_y_continuous(limits = c(0, NA))+
  labs(x = 'F as perceived by the model', y = '1000\'s of tons')+
  facet_grid2(LongNamePlot~type, scales = 'free', independent = 'all')+
  theme(strip.text.y = element_text(angle=0))
f_plot

# make a figure
ggsave(paste0(paste0('NOAA_Azure/results/figures/',batch,'/yield_curves'),t,'_MS.png'), f_plot, width = 8, height = 16)

# make figures for slides (break into two columns)
# plot
grp1 <- unique(to_plot$LongNamePlot)[1:6]
f_plot1 <- to_plot %>%
  filter(LongNamePlot %in% grp1) %>%
  ggplot(aes(x = f, y = mt/1000, color = experiment))+
  geom_line(linewidth = 1)+
  #geom_point(size = 1.6)+
  scale_color_viridis_d(begin = 0.2, end = 0.8)+
  geom_vline(data = fmsy %>% filter(LongNamePlot %in% grp1), aes(xintercept = FMSY, group = LongNamePlot), linetype = 'dashed', color = 'orange')+
  geom_vline(data = atlantis_fmsy %>% filter(LongNamePlot %in% grp1), aes(xintercept = atlantis_fmsy, group = LongNamePlot), linetype = 'dashed', color = 'blue')+
  geom_hline(data = b35 %>% filter(LongNamePlot %in% grp1), aes(yintercept = b35/1000, group = LongNamePlot), linetype = 'dashed', color = 'red')+
  theme_bw()+
  scale_y_continuous(limits = c(0, NA))+
  labs(x = 'F (as perceived by the model)', y = '1000\'s of tons')+
  facet_grid2(LongNamePlot~type, scales = 'free', independent = 'all')+
  theme(strip.text.y = element_text(angle=0))

grp2 <- unique(to_plot$LongNamePlot)[7:12]
f_plot2 <- to_plot %>%
  filter(LongNamePlot %in% grp2) %>%
  ggplot(aes(x = f, y = mt/1000, color = experiment))+
  geom_line(linewidth = 1)+
  #geom_point(size = 1.6)+
  scale_color_viridis_d(begin = 0.2, end = 0.8)+
  geom_vline(data = fmsy %>% filter(LongNamePlot %in% grp2), aes(xintercept = FMSY, group = LongNamePlot), linetype = 'dashed', color = 'orange')+
  geom_vline(data = atlantis_fmsy %>% filter(LongNamePlot %in% grp2), aes(xintercept = atlantis_fmsy, group = LongNamePlot), linetype = 'dashed', color = 'blue')+
  geom_hline(data = b35 %>% filter(LongNamePlot %in% grp2), aes(yintercept = b35/1000, group = LongNamePlot), linetype = 'dashed', color = 'red')+
  theme_bw()+
  scale_y_continuous(limits = c(0, NA))+
  labs(x = 'F (as perceived by the model)', y = '1000\'s of tons')+
  facet_grid2(LongNamePlot~type, scales = 'free', independent = 'all')+
  theme(strip.text.y = element_text(angle=0))

# make a figure
ggsave(paste0('NOAA_Azure/results/figures/',batch,'/yield_curves',t,'_MS_1.png'), f_plot1, width = 6, height = 6)
ggsave(paste0('NOAA_Azure/results/figures/',batch,'/yield_curves',t,'_MS_2.png'), f_plot2, width = 6, height = 6)

# save the file as a csv to compare it to other runs
write.csv(to_plot %>% mutate(batch = batch), 
          paste0('NOAA_Azure/results/for_comp/',batch,'_ms_vs_ss.csv'), row.names = F)
##############################################################################################

# Multispecies yield is always higher than single species yield. 
# Most species do better when other species are also fished. 
# To get a full picture we need to look also at:
# Biomass of unselected age classes
# Biomass of forage fish and predators
# Overall, this is yet another symptom of low consumption across the model
# Maximum priority to mapping consumption from PROD.nc files

# plot the above as a function of the multiplier for the permutations
to_plot <- df_mult 

# spaces
to_plot$LongNamePlot <- gsub(" ", "\n", to_plot$LongName)
b35$LongNamePlot <- gsub(" ", "\n", b35$LongName)

# plot
grp1 <- unique(to_plot$LongNamePlot)[1:6]
mult_plot1 <- to_plot %>%
  filter(LongNamePlot %in% grp1) %>%
  ggplot(aes(x = mult, y = mt / 1000))+
  geom_line(linewidth = 1)+
  geom_point(size = 1.6)+
  geom_vline(xintercept = 1, color = "blue", linetype = "dashed")+
  geom_hline(data = b35 %>% filter(LongNamePlot %in% grp1), aes(yintercept = b35/1000, group = LongNamePlot), linetype = 'dashed', color = 'red')+
  theme_bw()+
  scale_y_continuous(limits = c(0, NA))+
  labs(x = 'Multiplier of F35', y = '1000\'s of tons')+
  facet_grid2(LongNamePlot~type, scales = 'free', independent = 'all')+
  theme(strip.text.y = element_text(angle=0))

grp2 <- unique(to_plot$LongNamePlot)[7:12]
mult_plot2 <- to_plot %>%
  filter(LongNamePlot %in% grp2) %>%
  ggplot(aes(x = mult, y = mt / 1000))+
  geom_line(linewidth = 1)+
  geom_point(size = 1.6)+
  geom_vline(xintercept = 1, color = "blue", linetype = "dashed")+
  geom_hline(data = b35 %>% filter(LongNamePlot %in% grp2), aes(yintercept = b35/1000, group = LongNamePlot), linetype = 'dashed', color = 'red')+
  theme_bw()+
  scale_y_continuous(limits = c(0, NA))+
  labs(x = 'Multiplier of F35', y = '1000\'s of tons')+
  facet_grid2(LongNamePlot~type, scales = 'free', independent = 'all')+
  theme(strip.text.y = element_text(angle=0))

ggsave(paste0('NOAA_Azure/results/figures/',batch,'/mult_curves',t,'_MS_1.png'), mult_plot1, width = 6, height = 6)
ggsave(paste0('NOAA_Azure/results/figures/',batch,'/mult_curves',t,'_MS_2.png'), mult_plot2, width = 6, height = 6)

# Diagnostics: biomass of unselected age classes --------------------------

# look at unselected biomass in both ms and ss
ms_unselected_list <- list()

for(i in 1:nrow(f35_key)){
  
  this_catch_file <- f35_catch_files[i]
  this_run <- as.numeric(substr((gsub(paste0("NOAA_Azure/results/f35/",batch,"/outputGOA0"), "", this_catch_file)), 1, 4))
  this_mult <- f35_key %>% filter(run == this_run) %>% pull(mult)
  
  # extract tables from results
  biomage <- read.table(f35_biomass_files[i], sep = " ", header = T)
  catch <- read.table(this_catch_file , sep = " ", header = T)
  
  # now extract data
  # SSB to plot and report in tables
  unselected_biomass <- biomage %>% 
    slice_tail(n = 5) %>% # use last xxx years
    summarise(across(-"Time", ~ mean(.x, na.rm = TRUE))) %>%
    ungroup() %>%
    pivot_longer(everything(), names_to = 'Code.Age', values_to = 'biomass_mt') %>%
    separate_wider_delim(Code.Age, delim = '.', names = c('Code', 'Age')) %>%
    filter(Code %in% t3_fg) %>%
    left_join(selex, by = 'Code') %>%
    mutate(idx = as.numeric(Age) - as.numeric(age_class_selex)) %>%
    filter(is.na(idx) | idx < 0) %>%
    group_by(Code) %>%
    summarise(biomass_mt = sum(biomass_mt)) %>%
    ungroup()
  
  # # calculate realized F after 1 year of data
  # # get initial biomass for the selected age classes
  biom_age_t1 <- biomage %>% 
    filter(Time == 0) %>% 
    pivot_longer(-Time, names_to = 'Code.Age', values_to = 'biomass') %>%
    separate_wider_delim(Code.Age, delim = '.', names = c('Code', 'Age')) %>%
    left_join(selex, by = 'Code') %>%
    mutate(idx = as.numeric(Age) - as.numeric(age_class_selex)) %>%
    filter(is.na(idx) | idx >= 0) %>%
    group_by(Code) %>%
    summarise(biomass = sum(biomass)) %>%
    ungroup() %>% 
    filter(Code %in% t3_fg)
  # 
  # # catch (one time step after biomass: how much did we catch in this time?)
  catch_t1 <- catch %>% 
    select(Time, all_of(t3_fg)) %>% 
    filter(Time == 365) %>%
    summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) %>%
    pivot_longer(-Time, names_to = 'Code', values_to = 'catch') %>%
    select(-Time)
  # 
  # # calc realized f
  f_t1 <- biom_age_t1 %>% left_join(catch_t1, by = "Code") %>%
    mutate(exp_rate = catch/biomass,
           f = -log(1-exp_rate)) %>%#,
    #fidx = fidx) %>% # need this for joining later on
    select(Code, f)#, fidx) 
  
  # # bind all
  f_frame <- f_t1 %>%
    left_join(unselected_biomass) %>%
    #left_join(catch_vals) %>%
    mutate(mult = this_mult)
  
  # add to multispecies yield list
  ms_unselected_list[[i]] <- f_frame
}

ms_unselected_df <- bind_rows(ms_unselected_list)

ms_unselected_df <- ms_unselected_df %>%
  left_join(grps %>% select(Code, LongName), by = "Code") %>%
  mutate(experiment = "ms") %>%
  select(Code, LongName, experiment, f, biomass_mt)
  
# create empty list to fill with data frame for the yield curve
folder_path <- here("NOAA_Azure","results","pre-processing",this_job,"results")
# list the results files
results_list <- list.files(folder_path, full.names = T)

ss_unselected_list <- list()
for(i in 1:length(results_list)){ # TODO: turn all of these into a function with a few args (e.g., fished, unselected, top preds, ms, ss, etc.)
  
  print(paste("Doing",i))
  
  # load results
  this_result <- readRDS(results_list[i])
  
  # this idx - CAREFUL! This is not the same as i because of alphabetical sorting of the reuslts files
  this_idx <- as.numeric(names(this_result))
  
  # first identify the species and the level of fishing for this run. These are unrelated from runname
  sp <- f_lookup %>% filter(idx == this_idx) %>% pull(species)
  fidx <- f_lookup %>% filter(idx == this_idx) %>% pull(fidx)
  
  # rename the result object in the list to avoid problems with indexing
  names(this_result) <- "res"
  
  # extract tables from results
  biomage <- this_result$res[[paste0("biomage_",this_idx)]]
  catch <- this_result$res[[paste0("catch_",this_idx)]]
  
  # now extract data
  # SSB to plot and report in tables
  unselected_biomass <- biomage %>% 
    slice_tail(n = 5) %>% # use last xxx years
    summarise(across(-"Time", ~ mean(.x, na.rm = TRUE))) %>%
    ungroup() %>%
    pivot_longer(everything(), names_to = 'Code.Age', values_to = 'biomass_mt') %>%
    separate_wider_delim(Code.Age, delim = '.', names = c('Code', 'Age')) %>%
    filter(Code == sp) %>%
    left_join(selex, by = 'Code') %>%
    mutate(idx = as.numeric(Age) - as.numeric(age_class_selex)) %>%
    filter(is.na(idx) | idx < 0) %>%
    group_by(Code) %>%
    summarise(biomass_mt = sum(biomass_mt)) %>%
    ungroup() %>%
    pull(biomass_mt)
  
  # # calculate realized F after 1 year of data
  # # get initial biomass for the selected age classes
  biom_age_t1 <- biomage %>% 
    filter(Time == 0) %>% 
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
    filter(Time == 365) %>%
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
    mutate(biomass = unselected_biomass)#,
           #catch = catch_val)
  
  ss_unselected_list[[i]] <- f_frame
  
  # # write out to be then brought together with all other runs
  # write.csv(f_frame, paste(outdir,paste(sp,fidx,'f.csv',sep='_'), sep = "/"), row.names = F)
  
}

ss_unselected_df <- bind_rows(ss_unselected_list) %>%
  left_join(grps %>% select(Code, LongName), by= "Code") %>%
  mutate(experiment = "ss") %>%
  select(Code, LongName, experiment, f, biomass) %>%
  rename(biomass_mt = biomass)

# plot
unselected_plot <- ms_unselected_df %>%
  rbind(ss_unselected_df) %>%
  mutate(experiment = ifelse(experiment == "ms", "Multispecies", "Single-species")) %>%
  ggplot(aes(x = f, y = biomass_mt/1000, color = experiment))+
  geom_line()+
  geom_point(size = 2)+
  scale_color_viridis_d(begin = 0.2, end = 0.8)+
  geom_vline(data = fmsy, aes(xintercept = FMSY, group = LongName), linetype = 'dashed', color = 'orange')+
  geom_vline(data = atlantis_fmsy, aes(xintercept = atlantis_fmsy, group = LongName), linetype = 'dashed', color = 'blue')+
  theme_bw()+
  scale_y_continuous(limits = c(0, NA))+
  labs(x = 'F as perceived by the model', y = '1000\'s of tons')+
  facet_wrap(~LongName, scales = 'free', ncol = 3)+
  theme(strip.text.y = element_text(angle=0))
unselected_plot

ggsave(paste0("NOAA_Azure/results/figures/",batch,"/unselected.png"), unselected_plot, width = 8, height = 6)

# save the file as a csv to compare it to other runs
write.csv(ms_unselected_df %>% mutate(batch = batch), 
          paste0('NOAA_Azure/results/for_comp/',batch,'_unselected.csv'), row.names = F)

# In the MS scenario, the biomass of the unselected age classes is in general higher across the board. 
# One possible explanation is that their adult predators are getting fished out
# This does not really hold if we look at their biomass, which is not lower than in the ss runs
# So, in the MS runs we have higher biomass and higher catches across all T3 species.
# In other words, when we fish a stock at FMSY, if we fish all stocks at FMSY we get more of all of them
# Since they can't all be winners, who is losing?
# look at top predators and forage fish across the ms runs

# Diagnostics: top predators and forage fish ------------------------------

top_preds <- grps %>% filter(GroupType %in% c("MAMMAL","BIRD","SHARK")) %>% pull(Code)
forage <- c("CAP","SAN","HER","EUL","FOS")
other_fg <- c(top_preds, forage)
  
ms_other_list <- list()

for(i in 1:nrow(f35_key)){
  
  # get the multiplier
  this_biomage_file <- f35_biomass_files[i]
  this_run <- as.numeric(substr((gsub(paste0("NOAA_Azure/results/f35/",batch,"/outputGOA0"), "", this_biomage_file)), 1, 4))
  this_mult <- f35_key %>% filter(run == this_run) %>% pull(mult)
  
  # extract tables from results
  biomage <- read.table(this_biomage_file, sep = " ", header = T)

  # now extract data
  # SSB to plot and report in tables
  other_biomass <- biomage %>% 
    slice_tail(n = 5) %>% # use last xxx years
    summarise(across(-"Time", ~ mean(.x, na.rm = TRUE))) %>%
    ungroup() %>%
    pivot_longer(everything(), names_to = 'Code.Age', values_to = 'biomass_mt') %>%
    separate_wider_delim(Code.Age, delim = '.', names = c('Code', 'Age')) %>%
    filter(Code %in% other_fg) %>%
    group_by(Code) %>%
    summarise(biomass_mt = sum(biomass_mt)) %>%
    ungroup() %>%
    mutate(mult = this_mult)
  
  # add to multispecies yield list
  ms_other_list[[i]] <- other_biomass
}

ms_other_df <- bind_rows(ms_other_list) %>%
  left_join(grps %>% select(Code, LongName), by = "Code")

# get b0
b0_other <- ms_other_df %>% filter(mult == 0) %>% dplyr::select(LongName, biomass_mt) %>% rename(b0 = biomass_mt)

ms_other_df <- ms_other_df %>%
  left_join(b0_other, by = "LongName") %>%
  mutate(biomchange = biomass_mt / b0)

# plot (separate top preds and forage)
other_plot_top <- ms_other_df %>%
  filter(Code %in% c("KWT","KWR","WHT","WHH","WHB","WHG","DOL","SSL","PIN","BDF","BDI","BSF","BSI","SHD","SHP")) %>%
  ggplot(aes(x = mult, y = biomchange))+
  geom_line()+
  geom_point()+
  geom_hline(yintercept = 1, color = "red", linetype = "dashed")+
  theme_bw()+
  labs(x = "F35 permutation", y = "Change in biomass from unfished")+
  facet_wrap(~LongName)

other_plot_forage <- ms_other_df %>%
  filter(Code %in% c("CAP","SAN","EUL","HER","FOS")) %>%
  ggplot(aes(x = mult, y = biomchange))+
  geom_line()+
  geom_point()+
  geom_hline(yintercept = 1, color = "red", linetype = "dashed")+
  theme_bw()+
  labs(x = "F35 permutation", y = "Change in biomass from unfished")+
  facet_wrap(~LongName)

ggsave(paste0("NOAA_Azure/results/figures/",batch,"/other_top.png"), other_plot_top, width = 8, height = 6)
ggsave(paste0("NOAA_Azure/results/figures/",batch,"/other_forage.png"), other_plot_forage, width = 8, height = 6)

# save the file as a csv to compare it to other runs
write.csv(ms_other_df %>% mutate(batch = batch), 
          paste0('NOAA_Azure/results/for_comp/',batch,'_top_pred.csv'), row.names = F)
write.csv(ms_other_df %>% mutate(batch = batch), 
          paste0('NOAA_Azure/results/for_comp/',batch,'_forage.csv'), row.names = F)

# overall pretty minimal differences, but in general prey increases (a little) and predators decrease (a little) under higher fishing
# perhaps more of a case for KWR (halved), SSL (-25% or so), and sharks (small declines), the rest does fairly well under fishing

# Diagnostics: M -----------------------------------------------------------------------

f35_mort_files <- list.files(f35_path, pattern = "Mort", full.names = T)

# for each group, calculate the relationship between F and M as per mort.txt file
f_to_m <- data.frame()
for(i in 1:nrow(f35_key)){
  
  this_m_file <- f35_mort_files[i]
  this_run <- as.numeric(substr(gsub(paste0("NOAA_Azure/results/f35/",batch,"/outputGOA0"), "", this_m_file), 1, 4))
  this_mult <- f35_key %>% filter(run == this_run) %>% pull(mult)
  
  # extract tables from results
  m <- read.table(f35_mort_files[i], sep = " ", header = T)
  
  # the values of M and F are only meaningful relative to one another
  m_long <- m %>%
    slice_tail(n = 5) %>% # use last xxx years
    summarise(across(-"Time", ~ mean(.x, na.rm = TRUE))) %>%
    ungroup() %>%
    pivot_longer(everything(), names_to = 'Code.MType', values_to = 'rawM') %>%
    separate_wider_delim(Code.MType, delim = '.', names = c('Code', 'MType'))
  
  # now get M and F on a line each and get the proportion
  
  this_f_to_m <- m_long %>%
    pivot_wider(id_cols = Code, names_from = MType, values_from = rawM) %>%
    mutate(f_mult = M / `F`,
           mult = this_mult) %>%
    dplyr::select(Code, mult, f_mult)

  # add to multispecies yield list
  f_to_m <- rbind(f_to_m, this_f_to_m)
} 

# join to data frame with f values from the multispecies experiment
m_df <- df_mult %>%
  left_join(f_to_m, by = c("Code","mult")) %>%
  mutate(m = f * f_mult)

# plot
mp <- m_df %>%
  filter(mult > 0) %>%
  ggplot(aes(x = mult, y = m))+
  geom_line()+
  theme_bw()+
  facet_wrap(~LongName, scales = "free")
mp

# M values are all over the place. Also we don't get to plot the values for f = 0, because the multiplier is not calculated
# if we need a baseline, apply very low mFC
# start from mult = 0.2
# apart from the scale that is all over the place, M decreases with increasing F for ATF, FFD, FHS, COD, RFS, and POL - likely all ATF's prey items
# More surprisingly, M increases with F for HAL, POP, REX, RFP, SBF. Need to investigate, but one possible explanation is that their (very few) predators (like SSL and KWR) switch from some of the depleted groundfish to these other species

# to convey the same message without worrying about the scale, do this relative to mult = 2.0
m_df_mult02 <- m_df %>% filter(mult == 0.2) %>% dplyr::select(LongName, m) %>% rename(m02 = m) %>% distinct()
m_df_rel <- m_df %>%
  left_join(m_df_mult02, by = "LongName") %>%
  mutate(rel_m = m / m02)

# plot
rel_m_plot <- m_df_rel %>%
  filter(mult > 0) %>%
  ggplot(aes(x = mult, y = rel_m))+
  geom_line()+
  geom_hline(yintercept = 1, color = "red", linetype = "dashed")+
  theme_bw()+
  labs(x = "F35 permutation", y = "Relative change in M (from mult = 0.2)")+
  facet_wrap(~LongName)
rel_m_plot

ggsave(paste0("NOAA_Azure/results/figures/", batch, "/relative_m.png"), rel_m_plot, width = 8, height = 4)

write.csv(m_df_rel %>% mutate(batch = batch), 
          paste0('NOAA_Azure/results/for_comp/',batch,'_rel_m.csv'), row.names = F)

# Diagnostics: yield functions --------------------------------------------

# get two data frames: one for b0 and one for maximum yield
# b0
b0_ms <- ms_yield_long %>% filter(f == 0, type == "Biomass") %>% dplyr::select(LongName, experiment, mt) %>% rename(b0 = mt)
b0_ss <- ss_yield_long %>% filter(f == 0, type == "Biomass") %>% dplyr::select(LongName, experiment, mt) %>% rename(b0 = mt)
b0 <- rbind(b0_ms, b0_ss)

# max yield
ymax_ms <- ms_yield_long %>% 
  filter(type == "Catch") %>% 
  group_by(LongName, experiment) %>%
  slice_max(mt) %>%
  ungroup() %>%
  dplyr::select(LongName, experiment, mt) %>% 
  rename(ymax = mt) 

ymax_ss <- ss_yield_long %>% 
  filter(type == "Catch") %>% 
  group_by(LongName, experiment) %>%
  slice_max(mt) %>%
  ungroup() %>%
  dplyr::select(LongName, experiment, mt) %>% 
  rename(ymax = mt) 

ymax <- rbind(ymax_ms, ymax_ss)

# we are plotting yield fraction against depletion
yield_func <- ms_yield_long %>%
  rbind(ss_yield_long) %>%
  dplyr::select(LongName, experiment, type, f, mt) %>%
  pivot_wider(id_cols = c(LongName, experiment, f), names_from = type, values_from = mt) %>%
  left_join(b0, by = c("LongName", "experiment")) %>%
  mutate(depletion = Biomass / b0) %>%
  left_join(ymax, by = c("LongName", "experiment")) %>%
  mutate(yfrac = Catch / ymax) %>%
  #mutate(experiment = ifelse(experiment == "ms", "Multispecies", "Single-species")) %>%
  dplyr::select(LongName, experiment, yfrac, depletion)

# prepare data frames to write the following quantities on the plot:
# final depletion, final yield fraction, biomass corresponding to final depletion, biomass corresponding to final yield fraction
yfun_terminal <- yield_func %>%
  arrange(LongName, experiment) %>%
  group_by(LongName, experiment) %>%
  slice_tail(n = 1) %>%
  mutate(depletion = depletion * 100,
         yfrac = yfrac * 100) # turn depletion and yield to percentages for easier interpretation

annotations <- b0 %>%
  left_join(ymax) %>%
  left_join(yfun_terminal) %>%
  mutate(catch = ymax / 100 * yfrac / 1000,
         ssb = b0 / 100 * depletion / 1000) %>%
  filter(experiment == "ms")

# plot
yield_func_plot <- yield_func %>%
  filter(experiment == "ms") %>%
  ggplot(aes(x = depletion, y = yfrac))+
  geom_line()+
  scale_x_reverse()+
  geom_text(data = annotations,
            aes(x = 0.5, y = 0.5, hjust=0.5, vjust=1,
                label=paste0('Depletion(%)=', round(depletion,2),
                             '\n',
                             'Yield fraction(%)=', round(yfrac, 2),
                             '\n',
                             'SSB(1000mt)=', round(ssb, 2),
                             '\n',
                             'Catch(1000mt)=', round(catch, 2))), 
            size = 3)+
  theme_bw()+
  labs(x = "Depletion", y = "Yield fraction")+
  facet_wrap(~LongName)
yield_func_plot

ggsave(paste0("NOAA_Azure/results/figures/", batch, "/yield_functions.png"), yield_func_plot, width = 8, height = 5)

# The yield curves illustrate fundamental issues with either model productivity or my calculations here
# Stocks are way depleted by the catch won't go down
# to the extent that POL SSB is 200,000 mt, and catch is 300,000 mt.
# same applies to all other stocks
# For POL, startage < age_mat, so there is no selectivity refuge here

# Diagnostics: equilibrium ------------------------------------------------

# one thing to check is - are we at equilibrium here or not?
eq_catch_list <- list()
for(i in 1:nrow(f35_key)){
  
  this_catch_file <- f35_catch_files[i]
  this_run <- as.numeric(substr((gsub(paste0("NOAA_Azure/results/f35/",batch,"/outputGOA0"), "", this_catch_file)), 1, 4))
  this_mult <- f35_key %>% filter(run == this_run) %>% pull(mult)
  
  this_catch <- read.table(this_catch_file, sep=" ", header = T)
  this_catch <- this_catch %>%
    slice_tail(n = 20) %>% # look at the  last 20 years without averageing
    dplyr::select(Time, all_of(t3_fg)) %>%
    pivot_longer(-Time, names_to = "Code", values_to = "mt") %>%
    mutate(mult = this_mult)

  eq_catch_list[[i]] <- this_catch
}

eq_catch_df <- bind_rows(eq_catch_list) %>%
  mutate(yr = Time / 365)

# plot
eq_catch_df %>%
  ggplot(aes(x = yr, y = mt, color = factor(mult)))+
  geom_line()+
  theme_bw()+
  facet_wrap(~Code, scales = "free_y")

# catches in slight decline over time but not to extent of suggesting a collapse immediately after the end of the run (50 years)
# now biomass
eq_biom_list <- list()
for(i in 1:nrow(f35_key)){
  
  this_biom_file <- f35_biomass_files[i]
  this_run <- as.numeric(substr((gsub(paste0("NOAA_Azure/results/f35/",batch,"/outputGOA0"), "", this_biom_file)), 1, 4))
  this_mult <- f35_key %>% filter(run == this_run) %>% pull(mult)
  
  this_biom <- read.table(this_biom_file, sep=" ", header = T)
  this_biom <- this_biom %>%
    slice_tail(n = 20) %>% # look at the  last 20 years without averageing
    pivot_longer(-Time, names_to = 'Code.Age', values_to = 'biomass_mt') %>%
    separate_wider_delim(Code.Age, delim = '.', names = c('Code', 'Age')) %>%
    filter(Code %in% t3_fg) %>%
    left_join(mat, by = 'Code') %>%
    mutate(idx = as.numeric(Age) - as.numeric(age_class_mat)) %>%
    filter(is.na(idx) | idx >= 0) %>%
    group_by(Time, Code) %>%
    summarize(mt = sum(biomass_mt))%>%
    mutate(mult = this_mult)
    
  
  eq_biom_list[[i]] <- this_biom
}

eq_biom_df <- bind_rows(eq_biom_list) %>%
  mutate(yr = Time / 365)

# plot
eq_biom_df %>%
  ggplot(aes(x = yr, y = mt, color = factor(mult)))+
  geom_line()+
  theme_bw()+
  facet_wrap(~Code, scales = "free_y")

# similarly to catch, we are fairly close to equilibrium for SSB and there are no signs of imminent collapse
# This points to excessive productivity, possibly mediated by excessive recruitment
# It may be worth repeating the set with (much) lower steepness. Reducing steepness did not help when stocks had mature age classes below the knife edge, but it may help in cases where that is not the case.

# Plot global yield but from single-species (it won't be the same as the other one)


# Plot changes in weight at age
# this one will be as little more complex but we have all the components (take from code we use)
