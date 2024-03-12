# Alberto Rovellini
# 01/16/2024
# This code takes output of the F35 permutation runs from 4 scenarios and plots:
# Production functions
# Catch and biomass curves
# Numbers at age
# biomass of forage fish and predators

# Update 12/14/2023
# Dropping age at maturity entirely and using FSPB to get SSB
# TODO: change single-species stuff too, but first see how they compare

library(tidyverse)
library(here)
library(tidyr)
library(readxl)
library(ggh4x)
library(viridis)
library(tidync)
library(ncdf4)

# Set up env and read data ------------------------------------------------

burnin <- 30 # years of burn-in

# identify which data we want to work on
# there are a variety of paths here because many runs failed in the batch job and were (after much trying and tribulation) done manually on the fast eScience machines
batch_res1 <- "job20240210013612/results" # ms_oy - these are the biomage, catch, and mort files
batch_res2 <- "job20240210013612-part2/results" # ms_oy - these are the biomage, catch, and mort files

batch_nc1 <- "job20240210013612/oy-ms-3" # these are the full out.nc files
batch_nc2 <- "job20240210013612-part2/oy-ms-3-part2" # these are the full out.nc files

ss_job <- "job20240109052342" # this is SS runs
maxmult <- 4 # this is the full range of explored F

# set the clock to date plots
t <- format(Sys.time(),'%Y-%m-%d %H-%M-%S')

# read in Groups.csv file
grp_path <- here('NOAA_Azure/data/GOA_Groups.csv') # functional groups
grps <- read.csv(grp_path)

# read maturity and selectivity information
mat <- read.csv("NOAA_Azure/data/age_at_mat.csv", header = T) # age at maturity
selex <- read.csv("NOAA_Azure/data/age_at_selex_new.csv", header = T) # age at selectivity
fspb <- read.csv("NOAA_Azure/data/fspb.csv", header = T) # proportion of spawning biomass per age class
# reshape fspb
fspb <- fspb %>%
  pivot_longer(-Code, names_to = "Age", values_to = "fspb") %>%
  mutate(Age = gsub("X","",Age))

# read in lookup tables and f35 proxy values
oy_key <- read.csv("NOAA_Azure/data/oy_key.csv")
f_lookup <- read.csv("NOAA_Azure/data/f_lookup_OY_SS.csv")
f35_vector <- read.csv("NOAA_Azure/data/f35_vector_PROXY.csv")

# list Tier 3 stocks 
t3_fg <- f_lookup %>% pull(species) %>% unique() %>% sort()

# get the results files from the multispecies runs
f35_path1 <- paste0("NOAA_Azure/results/f35/",batch_res1)
f35_path2 <- paste0("NOAA_Azure/results/f35/",batch_res2)
# list the rds files
f35_results <- c(list.files(f35_path1, pattern = ".rds", full.names = T),
                 list.files(f35_path2, pattern = ".rds", full.names = T))

# order them correctly
# reorder these based on the number in the filename
num_idx <- as.numeric(gsub("([0-9]+)-result\\.rds", "\\1", 
                           c(list.files(f35_path1, pattern = ".rds", full.names = F), 
                             list.files(f35_path2, pattern = ".rds", full.names = F))))
f35_results <- f35_results[order(num_idx)]

# get the nc files
f35_nc_path1 <- paste0("NOAA_Azure/results/f35/",batch_nc1)
f35_nc_path2 <- paste0("NOAA_Azure/results/f35/",batch_nc2)
# list the nc files
f35_nc <- c(list.files(f35_nc_path1, pattern = ".nc", full.names = T),
            list.files(f35_nc_path2, pattern = ".nc", full.names = T))
# order them correctly
# reorder these based on the number in the filename
num_idx <- as.numeric(gsub("output_([0-9]+)\\.nc", "\\1", 
                           c(list.files(f35_nc_path1, pattern = ".nc", full.names = F),
                             list.files(f35_nc_path2, pattern = ".nc", full.names = F))))
f35_nc <- f35_nc[order(num_idx)]

# extract biomass and catch from the MS runs
ms_yield_list <- list()

for(i in 1:length(f35_results)){
  
  print(paste("Doing", f35_results[i]))
  
  # grab the index from the file name
  this_idx <- as.numeric(gsub("-result.rds", "", gsub("NOAA.*2/results/", "", f35_results[i])))
  
  # run information based on the index
  this_run <- oy_key %>% filter(idx == this_idx) %>% pull(run)
  this_mult <- oy_key %>% filter(idx == this_idx) %>% pull(mult)
  
  # extract tables from results
  this_result <- readRDS(f35_results[i])
  # the packaging of the RDS object was different between the eScience runs and the batch (doAzureParallel) runs
  if(length(this_result)==1) {
    this_result <- this_result[[1]]
  }
  
  biomage <- this_result[[2]]
  catch <- this_result[[3]]
  mort <- this_result[[4]]
  
  # now extract data
  # SSB to plot and report in tables
  spawning_biomass <- biomage %>% 
    slice_tail(n = 5) %>% # use last xxx years
    summarise(across(-"Time", ~ mean(.x, na.rm = TRUE))) %>%
    ungroup() %>%
    pivot_longer(everything(), names_to = 'Code.Age', values_to = 'biomass_mt') %>%
    # separate_wider_delim(Code.Age, delim = '.', names = c('Code', 'Age')) %>% # not working on eScience
    separate(Code.Age, into = c('Code', 'Age'), sep = '\\.') %>%
    filter(Code %in% t3_fg) %>%
    #left_join(mat, by = 'Code') %>%
    # mutate(idx = as.numeric(Age) - as.numeric(age_class_mat)) %>%
    # filter(is.na(idx) | idx >= 0) %>%
    left_join(fspb, by = c('Code','Age')) %>%
    mutate(biomass_mt = biomass_mt * fspb) %>%
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
  # For runs with a burn-in, this has to be the biomass at the end of the burn-in, when we start fishing with the new scalar
  # # get initial biomass for the selected age classes
  biom_age_t1 <- biomage %>% 
    filter(Time == 365 * burnin) %>%# this is the burn-in years
    pivot_longer(-Time, names_to = 'Code.Age', values_to = 'biomass') %>%
    # separate_wider_delim(Code.Age, delim = '.', names = c('Code', 'Age')) %>%
    separate(Code.Age, into = c('Code', 'Age'), sep = '\\.') %>%
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
    filter(Time == 365 * (burnin + 1)) %>% # careful - there is a small transition phase
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
    mutate(run = this_run,
           mult = this_mult)
  
  # add to multispecies yield list
  ms_yield_list[[i]] <- f_frame
}

ms_yield_df <- bind_rows(ms_yield_list) # Issues with FHS at high F - catch is higher than SSB at T0, so F from exploitation rate is infinite

# add species long name and reshape catch and biomass and add a label that this is the ms approach
ms_yield_long <- df_mult <- ms_yield_df %>% # one of these is for later plots that need mult
  left_join(grps %>% select(Code, LongName), by = "Code") %>%
  rename(Biomass = biomass_mt, Catch = catch_mt) %>%
  pivot_longer(cols = -c(Code, LongName, f, mult, run), names_to = "type", values_to = "mt") %>%
  select(Code, LongName, run, f, mult, type, mt)

# Reference points
# get two data frames: one for b0 and one for maximum yield
# b0
# For AMSS and OY climate scenarios:
# What is B0 under different climate and ATF fishing regimes?
# b0 <- ms_yield_long %>% filter(mult == 0, type == "Biomass") %>% dplyr::select(LongName, run, mt) %>% rename(b0 = mt)
# or leave it fixed to base conditions
b0 <- ms_yield_long %>% filter(mult == 0, type == "Biomass", run == "base") %>% dplyr::select(LongName, mt) %>% rename(b0 = mt)
# b0_ss <- ss_yield_long %>% filter(f == 0, type == "Biomass") %>% dplyr::select(LongName, experiment, mt) %>% rename(b0 = mt)
# b0 <- rbind(b0_ms, b0_ss)

# max yield
ymax_ms <- ms_yield_long %>% 
  filter(type == "Catch") %>% 
  group_by(LongName, run) %>%
  slice_max(mt) %>%
  ungroup() %>%
  dplyr::select(LongName, run, mt, f) %>% 
  rename(ymax = mt) 

ymax <- ymax_ms

# handle the NaN's from FHS
ms_yield_long <- as.data.frame(ms_yield_long)
ms_yield_long$f[is.nan(ms_yield_long$f)] <- NA

# Plot yield functions ----------------------------------------------------
# These are unclear to me when thinking of the different scenarios - what is b0?
# Seems like there is no sure answer (see work from Cody)
# Using b0 from the base run means that depletion of, for instance, Pcod would be really low at the lowest F as well
# but using b0 from the respective climate scenario seems equally problematic for other reasons, one being it will hide obvious depletion 
# the information we may get from such a case is how the shape of the curve changes: under different scenarios, how does fishing intensity affect the stock and its yield?
# If the curves from relative reference points are identical, F affects the stock in the smame way across scenarios
# If the curve changes (i.e. the peak is in a different position), the maximum yield is attained at different exploitation rates in different scenarios

# we are plotting yield fraction against depletion
yield_func <- ms_yield_long %>%
  drop_na() %>%
  #rbind(ss_yield_long) %>%
  dplyr::select(LongName, run, type, f, mt) %>%
  pivot_wider(id_cols = c(LongName, run, f), names_from = type, values_from = mt) %>%
  left_join(b0, by = c("LongName")) %>% # if you are keep static reference point
  # left_join(b0, by = c("LongName", "run")) %>%
  mutate(depletion = Biomass / b0) %>%
  left_join(ymax %>% select(-f), by = c("LongName", "run")) %>%
  mutate(yfrac = Catch / ymax) %>%
  #mutate(experiment = ifelse(experiment == "ms", "Multispecies", "Single-species")) %>%
  dplyr::select(LongName, run, yfrac, depletion, f)

# prepare data frames to write the following quantities on the plot:
# final depletion, final yield fraction, biomass corresponding to final depletion, biomass corresponding to final yield fraction
yfun_terminal <- yield_func %>%
  arrange(LongName, run) %>%
  group_by(LongName, run) %>%
  #slice_tail(n = 1) %>%
  mutate(depletion = depletion * 100,
         yfrac = yfrac * 100) # turn depletion and yield to percentages for easier interpretation

annotations <- b0 %>%
  left_join(ymax) %>%
  left_join(yfun_terminal) %>%
  mutate(catch = ymax / 100 * yfrac / 1000,
         ssb = b0 / 100 * depletion / 1000) 

# prepare for visualization as Cool vs Warm and as ATF varying vs fixed
yield_func <- yield_func %>%
  mutate(`F on\narrowtooth` = ifelse(run %in% c("atf","atf_climate"), "Fixed (1/4 FMSY)", "Varying"),
         Climate = ifelse(run %in% c("climate","atf_climate"), "ssp585 (2075-2085)", "Historical (1999)"))

# reorder ATF F
yield_func$`F on\narrowtooth` <- factor(yield_func$`F on\narrowtooth`, levels = c("Varying", "Fixed (1/4 FMSY)"))

# plot
yield_func_plot <- yield_func %>%
  #filter(LongName %in% c("Walleye pollock", "Pacific cod", "Arrowtooth flounder", "Pacific halibut")) %>%
  ggplot(aes(x = depletion, y = yfrac, color = Climate, linetype = `F on\narrowtooth`))+
  geom_point()+
  geom_line()+
  scale_color_manual(values = c("blue3", "red3"))+
  scale_x_reverse()+
  # geom_text(data = annotations,
  #           aes(x = 0.5, y = 0.5, hjust=0.5, vjust=1,
  #               label=paste0('Depletion(%)=', round(depletion,2),
  #                            '\n',
  #                            'Yield fraction(%)=', round(yfrac, 2),
  #                            '\n',
  #                            'SSB(1000mt)=', round(ssb, 2),
  #                            '\n',
  #                            'Catch(1000mt)=', round(catch, 2))), 
  #           size = 3)+
  theme_bw()+
  labs(x = "Depletion", y = "Yield fraction")+
  facet_wrap(~LongName)
yield_func_plot

# is there a way to make this plot useful?
# ggsave(paste0("NOAA_Azure/results/figures/oy/yield_functions.png"), yield_func_plot, width = 11, height = 6.5)

# Biomass and catch curves ------------------------------------------------
# plot catch and biomass curves
to_plot <- ms_yield_long

# spaces
to_plot$LongNamePlot <- gsub(" ", "\n", to_plot$LongName)
ymax$LongNamePlot <- gsub(" ", "\n", ymax$LongName)
ymax$type <- "Catch"

# rename scenarios and order them
# key_scenarios <- data.frame("Scenario" = c("Base", "Warm", "Fixed ATF", "Fixed ATF + Warm"),
#                             "run" = c("base", "climate", "atf", "atf_climate"))
# to_plot <- to_plot %>%
#   left_join(key_scenarios)
# 
# to_plot$Scenario <- factor(to_plot$Scenario, levels = c("Base", "Warm", "Fixed ATF", "Fixed ATF + Warm"))

# prepare for visualization as Cool vs Warm and as ATF varying vs fixed
# to_plot <- to_plot %>%
#   mutate(`F on\narrowtooth` = ifelse(run %in% c("atf","atf_climate"), "Fixed (1/4 FMSY)", "Varying"),
#          Climate = ifelse(run %in% c("climate","atf_climate"), "ssp585 (2075-2085)", "Historical (1999)"))
to_plot <- to_plot %>%
  mutate(Fishing = ifelse(run %in% c("atf","atf_climate"),
                                     "1/4 FMSY on arrowtooth flounder,\nMFMSY varying for all other stocks",
                                     "MFMSY varying for all stocks"),
         Climate = ifelse(run %in% c("climate","atf_climate"), "ssp585 (2075-2085)", "Historical (1999)"))

# ymax <- ymax %>%
#   mutate(`F on\narrowtooth` = ifelse(run %in% c("atf","atf_climate"), "Fixed (1/4 FMSY)", "Varying"),
#          Climate = ifelse(run %in% c("climate","atf_climate"), "ssp585 (2075-2085)", "Historical (1999)"))
ymax <- ymax %>%
  mutate(Fishing = ifelse(run %in% c("atf","atf_climate"),
                          "1/4 FMSY on arrowtooth flounder,\nMFMSY varying for all other stocks",
                          "MFMSY varying for all stocks"),
         Climate = ifelse(run %in% c("climate","atf_climate"), "ssp585 (2075-2085)", "Historical (1999)"))

# reorder ATF F
# to_plot$`F on\narrowtooth` <- factor(to_plot$`F on\narrowtooth`,
#                           levels = c("Varying",
#                                      "Fixed (1/4 FMSY)"))
# ymax$`F on\narrowtooth` <- factor(ymax$`F on\narrowtooth`,
#                                   levels = c("Varying",
#                                              "Fixed (1/4 FMSY)"))
to_plot$Fishing <- factor(to_plot$Fishing,
                          levels = c("MFMSY varying for all stocks",
                                     "1/4 FMSY on arrowtooth flounder,\nMFMSY varying for all other stocks"))
ymax$Fishing <- factor(ymax$Fishing,
                          levels = c("MFMSY varying for all stocks",
                                     "1/4 FMSY on arrowtooth flounder,\nMFMSY varying for all other stocks"))

# make figures for slides (break into two columns)
# plot
grp1 <- unique(to_plot$LongNamePlot)[1:6]
f_plot1 <- to_plot %>%
  filter(LongNamePlot %in% grp1) %>%
  filter(type == "Catch") %>%
  ggplot(aes(x = f, y = mt/1000, color = Climate, linetype = Fishing))+
  geom_line(linewidth = 1)+
  scale_color_viridis_d(begin = 0.2, end = 0.8)+
  geom_vline(data = ymax %>% filter(LongNamePlot %in% grp1), aes(xintercept = f, color = Climate, linetype = Fishing))+
  theme_bw()+
  scale_y_continuous(limits = c(0, NA))+
  labs(x = 'Fishing mortality (F)', y = '1000\'s of tons')+
  facet_grid2(LongNamePlot~type, scales = 'free', independent = 'all')+
  theme(strip.text.y = element_text(angle=0))
f_plot1

grp2 <- unique(to_plot$LongNamePlot)[7:12]
f_plot2 <- to_plot %>%
  filter(LongNamePlot %in% grp2) %>%
  filter(type == "Catch") %>%
  ggplot(aes(x = f, y = mt/1000, color = Climate, linetype = Fishing))+
  geom_line(linewidth = 1)+
  scale_color_viridis_d(begin = 0.2, end = 0.8)+
  geom_vline(data = ymax %>% filter(LongNamePlot %in% grp2), aes(xintercept = f, color = Climate, linetype = Fishing))+
  theme_bw()+
  scale_y_continuous(limits = c(0, NA))+
  labs(x = 'Fishing mortality (F)', y = '1000\'s of tons')+
  facet_grid2(LongNamePlot~type, scales = 'free', independent = 'all')+
  theme(strip.text.y = element_text(angle=0))
f_plot2

# make a figure
ggsave(paste0('NOAA_Azure/results/figures/oy_paper/biomass_catch',t,'_MS_1.png'), f_plot1, width = 7.5, height = 7)
ggsave(paste0('NOAA_Azure/results/figures/oy_paper/biomass_catch',t,'_MS_2.png'), f_plot2, width = 7.5, height = 7)

# look at a couple of key species only
# for plots
key_grps <- grps %>% filter(Code %in% c("POL", "COD", "ATF", "HAL", "SBF", "POP", "FFS")) %>% pull(LongName)
f_plot_ms <- to_plot %>%
  filter(LongName %in% key_grps, type == "Catch") %>%
  ggplot(aes(x = f, y = mt/1000, color = Climate, linetype = Fishing))+
  geom_line(linewidth = 1)+
  # geom_point(size = 1.6)+
  scale_color_viridis_d(begin = 0.2, end = 0.8)+
  geom_vline(data = ymax %>% 
               filter(LongName %in% key_grps) %>% 
               filter(!(LongNamePlot == "Arrowtooth\nflounder" & Fishing == "Low F on arrowtooth flounder,\nMFMSY varying for all other stocks")), 
             aes(xintercept = f, color = Climate, linetype = Fishing))+
  theme_bw()+
  scale_y_continuous(limits = c(0, NA))+
  labs(x = 'Fishing mortality (F)', y = 'Catch (1000 mt)')+
  facet_grid2(LongNamePlot~type, scales = 'free', independent = 'all')+
  #facet_wrap(~ LongName, scales = "free", ncol = 1)+
  theme(strip.text.y = element_text(angle=0))
f_plot_ms

ggsave(paste0('NOAA_Azure/results/figures/oy_paper/catch',t,'_MS_ms.png'), f_plot_ms, width = 6, height = 6)

# Diagnostics: top predators and forage fish ------------------------------

top_preds <- grps %>% filter(GroupType %in% c("MAMMAL","BIRD","SHARK")) %>% pull(Code) %>% as.character()
forage <- c("CAP","SAN","HER","EUL","FOS")
other_fg <- c(top_preds, forage)

ms_other_list <- list()

for(i in 1:length(f35_results)){
  
  print(paste("Doing", f35_results[i]))
  
  # grab the index from the file name
  this_idx <- as.numeric(gsub("-result.rds", "", gsub("NOAA.*2/results/", "", f35_results[i])))
  
  # run information based on the index
  this_run <- oy_key %>% filter(idx == this_idx) %>% pull(run)
  this_mult <- oy_key %>% filter(idx == this_idx) %>% pull(mult)
  
  # extract tables from results
  this_result <- readRDS(f35_results[i])
  # the packaging of the RDS object was different between the eScience runs and the batch (doAzureParallel) runs
  if(length(this_result)==1) {
    this_result <- this_result[[1]]
  }
  
  biomage <- this_result[[2]]
  
  # now extract data
  # SSB to plot and report in tables
  other_biomass <- biomage %>% 
    slice_tail(n = 5) %>% # use last xxx years
    summarise(across(-"Time", ~ mean(.x, na.rm = TRUE))) %>%
    ungroup() %>%
    pivot_longer(everything(), names_to = 'Code.Age', values_to = 'biomass_mt') %>%
    # separate_wider_delim(Code.Age, delim = '.', names = c('Code', 'Age')) %>%
    separate(Code.Age, into = c('Code', 'Age'), sep = '\\.') %>%
    filter(Code %in% other_fg) %>%
    group_by(Code) %>%
    summarise(biomass_mt = sum(biomass_mt)) %>%
    ungroup() %>%
    mutate(run = this_run,
           mult = this_mult)
  
  # add to multispecies yield list
  ms_other_list[[i]] <- other_biomass
}

ms_other_df <- bind_rows(ms_other_list) %>%
  left_join(grps %>% select(Code, LongName), by = "Code")

# get b0
# note that 
b0_other <- ms_other_df %>% filter(mult == 0) %>% dplyr::select(LongName, run, biomass_mt) %>% rename(b0 = biomass_mt)

ms_other_df <- ms_other_df %>%
  left_join(b0_other, by = c("LongName", "run")) %>%
  mutate(biomchange = biomass_mt / b0)

# add factors for plot
# ms_other_df <- ms_other_df %>%
#   mutate(`F on\narrowtooth` = ifelse(run %in% c("atf","atf_climate"), "Fixed (1/4 FMSY)", "Varying"),
#          Climate = ifelse(run %in% c("climate","atf_climate"), "ssp585 (2075-2085)", "Historical (1999)"))
# 
# # reorder ATF F
# ms_other_df$`F on\narrowtooth` <- factor(ms_other_df$`F on\narrowtooth`, levels = c("Varying", "Fixed (1/4 FMSY)"))

ms_other_df <- ms_other_df %>%
  mutate(Fishing = ifelse(run %in% c("atf","atf_climate"),
                          "Low F on arrowtooth flounder,\nMFMSY varying for all other stocks",
                          "MFMSY varying for all stocks"),
         Climate = ifelse(run %in% c("climate","atf_climate"), "ssp585 (2075-2085)", "Historical (1999)"))

# reorder ATF F
ms_other_df$Fishing <- factor(ms_other_df$Fishing,
                          levels = c("MFMSY varying for all stocks",
                                     "Low F on arrowtooth flounder,\nMFMSY varying for all other stocks"))

# plot (separate top preds and forage)
other_plot_top <- ms_other_df %>%
  filter(Code %in% c("KWT","KWR","DOL","SSL","PIN","BDF","BDI","BSF","BSI")) %>%
  ggplot(aes(x = mult, y = biomchange, color = Climate, linetype = Fishing))+
  geom_line()+
  geom_point()+
  #scale_color_manual(values = c("blue3", "red3"))+
  scale_color_viridis_d(begin = 0.2, end = 0.8)+
  geom_hline(yintercept = 1, color = "grey", linetype = "dotted", linewidth = 1)+
  geom_vline(xintercept = 1, color = 'black', linetype = "dotted", linewidth = 1)+
  theme_bw()+
  labs(x = expression(MF[MSY] ~ "multiplier"), y = "Change in biomass from unfished")+
  facet_wrap(~LongName)

other_plot_forage <- ms_other_df %>%
  filter(Code %in% c("CAP","SAN","EUL","HER")) %>%
  ggplot(aes(x = mult, y = biomchange, color = Climate, linetype = Fishing))+
  geom_line()+
  geom_point()+
  #scale_color_manual(values = c("blue3", "red3"))+
  scale_color_viridis_d(begin = 0.2, end = 0.8)+
  geom_hline(yintercept = 1, color = "grey", linetype = "dotted", linewidth = 1)+
  geom_vline(xintercept = 1, color = 'black', linetype = "dotted", linewidth = 1)+
  theme_bw()+
  labs(x = expression(MF[MSY] ~ "multiplier"), y = "Change in biomass from unfished")+
  facet_wrap(~LongName)

ggsave(paste0("NOAA_Azure/results/figures/oy_wfc/other_top.png"), other_plot_top, width = 8, height = 6)
ggsave(paste0("NOAA_Azure/results/figures/oy_wfc/other_forage.png"), other_plot_forage, width = 8, height = 4)


# Global yield ------------------------------------------------------------

# First, get the global yield from single-species runs (156)

# list all csv files we need to read
f_files <- list.files(paste0("NOAA_Azure/results/post-processing/",ss_job), full.names = T)

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

ss_yield_long <- ss_yield_long_fidx <- f_df %>%
  mutate(experiment = "ss") %>%
  filter(LongName != "Pacific halibut") %>% # drop halibut as it is not relevant to OY
  select(Code, LongName, f, fidx, experiment, type, mt)

# calcuate the Atlantis internal OY
# this is the procedure:
# for the single-species runs, get the highest yield at equilibrum and interpret that as MSY for a stock
# MSY for different stocks may occur at different modifiers of FOFL from FMP
# Then, add them up and discount them by 8%
atlantis_oy <- ss_yield_long_fidx %>%
  filter(type == "Catch") %>%
  group_by(LongName) %>%
  slice_max(mt) %>%
  pull(mt) %>%
  sum()

# discount by 8%
atlantis_oy <- atlantis_oy - (atlantis_oy * 0.08)

# the result is 597469.1 mt, much lower than 800,000
# This means that MSY calculated in Atlantis are lower than MSY from the assessments in 1987
# This is fine. It means that stocks are not as productive in Atlantis
# A caveat is that F on the nontarget species is at 1/4 FMSY, so this OY accounts for other stocks at a very specific state
# This state may be different from the state that other stocks had in 1987, and also their F is different
# So, the different result is not surprising
# More importantly, we should put this into the context of the multispecies MFSMY runs

ss_yield_long_fidx$LongNamePlot <- gsub(" - ", "\n", ss_yield_long_fidx$LongName)

global_yield_ss <- ss_yield_long_fidx %>%
  filter(type == "Catch") %>%
  left_join(f_lookup, by = c("Code"="species","fidx"="mult_idx")) %>%
  ggplot(aes(x = fofl_mult, y = mt / 1000, fill = LongNamePlot))+
  geom_bar(stat = "identity", position = "stack")+
  scale_fill_viridis_d()+
  scale_x_continuous(breaks = seq(0,maxmult,length.out=11))+
  scale_y_continuous(limits = c(0,900))+
  #geom_hline(yintercept = 116, color = "red", linetype = "dashed")+
  geom_hline(yintercept = 800, color = "red", linetype = "dashed")+
  geom_hline(yintercept = atlantis_oy / 1000, color = "blue", linetype = "dashed")+
  theme_bw()+
  labs(x = expression("F"[OFL]~" (from single-species assessment) multiplier"), y = "Catch (1000 mt)", fill = "Stock")
global_yield_ss

ggsave(paste0("NOAA_Azure/results/figures/oy_paper/global_yield_ss.png"), global_yield_ss, width = 6, height = 4)

# now for the multispecies runs
# list
catch_list <- list()
for(i in 1:length(f35_results)){
  
  print(paste("Doing", f35_results[i]))
  
  # grab the index from the file name
  this_idx <- as.numeric(gsub("-result.rds", "", gsub("NOAA.*2/results/", "", f35_results[i])))
  
  # run information based on the index
  this_run <- oy_key %>% filter(idx == this_idx) %>% pull(run)
  this_mult <- oy_key %>% filter(idx == this_idx) %>% pull(mult)
  
  # extract tables from results
  this_result <- readRDS(f35_results[i])
  # the packaging of the RDS object was different between the eScience runs and the batch (doAzureParallel) runs
  if(length(this_result)==1) {
    this_result <- this_result[[1]]
  }
  
  this_catch <- this_result[[3]]
  this_catch <- this_catch %>%
    slice_tail(n = 5) %>%
    summarise(across(all_of(t3_fg), ~mean(.x, na.rm = T))) %>%
    mutate(mult = this_mult,
           run = this_run)
  
  catch_list[[i]] <- this_catch
  
}

catch_df <- bind_rows(catch_list)

# reshape and calculate total
catch_df_long <- catch_df %>%
  pivot_longer(-c(run, mult), names_to = "Code", values_to = "mt") %>%
  filter(Code != "HAL") %>%
  group_by(run, mult) %>%
  mutate(total_yield = sum(mt),
         prop = mt / total_yield) %>%
  ungroup() %>%
  left_join(grps %>% select(Code, LongName), by = "Code")

# add scenario information
# pretty verbose for WFC now
catch_df_long <- catch_df_long %>%
  mutate(`F on\narrowtooth` = ifelse(run %in% c("atf","atf_climate"), 
                                     "1/4 FMSY on arrowtooth flounder,\nMFMSY varying for all other stocks", 
                                     "MFMSY varying for all stocks"),
         Climate = ifelse(run %in% c("climate","atf_climate"), "ssp585 (2075-2085)", "Historical (1999)"))

# reorder ATF F
catch_df_long$`F on\narrowtooth` <- factor(catch_df_long$`F on\narrowtooth`, 
                                           levels = c("MFMSY varying for all stocks",
                                                      "1/4 FMSY on arrowtooth flounder,\nMFMSY varying for all other stocks"))

# spaces
catch_df_long$LongNamePlot <- gsub(" - ", "\n", catch_df_long$LongName)

global_yield_ms <- catch_df_long %>%
  #filter(run %in% c("base", "climate")) %>%
  ggplot(aes(x = mult, y = mt / 1000, fill = LongNamePlot))+
  geom_bar(stat = "identity", position = "stack")+
  scale_fill_viridis_d()+
  #scale_x_continuous(breaks = seq(0,maxmult,length.out=11))+
  scale_y_continuous(limits = c(0,900))+
  # geom_hline(yintercept = 116, color = "red", linetype = "dashed")+
  geom_hline(yintercept = 800, color = "red", linetype = "dashed")+
  # geom_vline(xintercept = 1, color = "black", linetype = "dotted")+
  # geom_hline(yintercept = atlantis_oy / 1000, color = "blue", linetype = "dashed")+
  theme_bw()+
  labs(x = expression(MF[MSY] ~ "multiplier"), y = "Catch (1000 mt)", fill = "Stock") +
  theme(legend.position="bottom",
        legend.spacing.x = unit(0.1, 'cm'))+
  guides(fill = guide_legend(nrow = 4))+
  facet_grid(Climate~`F on\narrowtooth`)
global_yield_ms

ggsave(paste0("NOAA_Azure/results/figures/oy_paper/global_yield_ms.png"), global_yield_ms, width = 6.1, height = 7)

# test <- catch_df_long %>% filter(run == "base")
# test %>% group_by(mult) %>% summarize(mt = sum(mt))

# make a similar plot but with catch against F as areas or the likes
# global_yield_ms_area <- to_plot %>%
#   filter(type == "Catch") %>%
#   ggplot(aes(x = f, y = mt / 1000, fill = LongName))+
#   geom_area(stat = "align", position = "stack")+
#   scale_fill_viridis_d()+
#   #scale_x_continuous(breaks = seq(0,maxmult,length.out=11))+
#   scale_y_continuous(limits = c(0,1000))+
#   geom_hline(yintercept = 116, color = "red", linetype = "dashed")+
#   geom_hline(yintercept = 800, color = "red", linetype = "dashed")+
#   geom_vline(xintercept = 1, color = "black", linetype = "dotted")+
#   # geom_hline(yintercept = atlantis_oy / 1000, color = "blue", linetype = "dashed")+
#   theme_bw()+
#   labs(x = expression(MF[MSY] ~ "multiplier"), y = "Catch (1000 mt)", fill = "Stock") +
#   facet_grid(Climate~Fishing)
# global_yield_ms_area
# this shows the abrupt end of the range we are testing

# this is fairly interesting I think:
# varying F, cold: we get past the Atlantis OY when fishing multispecies
# varying F, warm: we are very close to the Atlantis OY
# fixed F, cold (kinda status quo): fairly close to Atlantis OY
# fixed F, warm (possible future): well below Atlantis OY
# Take home: under projected warming and low pressure on predators, we are way below Atlantis OY
# Part of the reason is stocks tanking under warming (cod), part is predation (pollock)

# ggsave(paste0("NOAA_Azure/results/figures/oy_wfc/global_yield_wfc_blank.png"), global_yield_ms, width = 8, height = 4)

# make a table with max catch per scenario for the report
max_catch <- catch_df_long %>%
  select(run, total_yield) %>%
  distinct() %>%
  group_by(run) %>%
  slice_max(total_yield)

# Single-species vs multispecies fishing ----------------------------------
ms_vs_ss <- ss_yield_long %>%
  mutate(run = "ss") %>%
  select(LongName, run, f, type, mt) %>%
  rbind(ms_yield_long %>% select(LongName, run, f, type, mt)) %>%
  filter(run %in% c("base", "atf", "ss"))

# change long names for plots
ms_vs_ss$LongNamePlot <- gsub(" ", "\n", ms_vs_ss$LongName)

# clarify labels
ms_vs_ss <- ms_vs_ss %>%
  mutate(Fishing = ifelse(run == "base", "MS - varying for all", ifelse(run == "ss", "SS", "MS - fixed for ATF")))

# reorder factors for legend
ms_vs_ss$Fishing <- factor(ms_vs_ss$Fishing, levels = c("SS", "MS - varying for all", "MS - fixed for ATF"))

# make some vertical lines for FMSY
ymax_comp <- ms_vs_ss %>%
  filter(type == "Catch") %>%
  group_by(run, LongNamePlot) %>%
  slice_max(mt)

# view
ms_vs_ss_plot_1 <- ms_vs_ss %>%
  filter(LongNamePlot %in% grp1, type == "Catch") %>%
  ggplot(aes(x = f, y = mt/1000, color = Fishing))+
  geom_line(linewidth = 1)+
  geom_vline(data = ymax_comp %>% 
               filter(LongNamePlot %in% grp1) %>% 
               filter(!(LongNamePlot == "Arrowtooth\nflounder" & Fishing == "MS - fixed for ATF")), 
             aes(xintercept = f, color =  Fishing))+
  #scale_color_manual(values = c("blue3", "red3", "orange"))+
  scale_color_viridis_d(begin = 0.1, end = 0.9)+
  theme_bw()+
  scale_y_continuous(limits = c(0, NA))+
  labs(x = 'Fishing mortality (F)', y = '1000\'s of tons')+
  facet_grid2(LongNamePlot~type, scales = 'free', independent = 'all')+
  theme(strip.text.y = element_text(angle=0))

ms_vs_ss_plot_2 <- ms_vs_ss %>%
  filter(LongNamePlot %in% grp2, type == "Catch") %>%
  ggplot(aes(x = f, y = mt/1000, color = Fishing))+
  geom_line(linewidth = 1)+
  geom_vline(data = ymax_comp %>% 
               filter(LongNamePlot %in% grp2) %>% 
               filter(!(LongNamePlot == "Arrowtooth\nflounder" & Fishing == "MS - fixed for ATF")), 
             aes(xintercept = f, color =  Fishing))+
  #scale_color_manual(values = c("blue3", "red3", "orange"))+
  scale_color_viridis_d(begin = 0.1, end = 0.9)+
  theme_bw()+
  scale_y_continuous(limits = c(0, NA))+
  labs(x = 'Fishing mortality (F)', y = '1000\'s of tons')+
  facet_grid2(LongNamePlot~type, scales = 'free', independent = 'all')+
  theme(strip.text.y = element_text(angle=0))

# make a figure
ggsave(paste0('NOAA_Azure/results/figures/oy_wfc/catch',t,'_comp_1.png'), ms_vs_ss_plot_1, width = 5.5, height = 7)
ggsave(paste0('NOAA_Azure/results/figures/oy_wfc/catch',t,'_comp_2.png'), ms_vs_ss_plot_2, width = 5.5, height = 7)


# Numbers at age from nc files --------------------------------------------

# expected to decline and be fairly close to 0 for older age classes when SSB is near 0
# should that not be the case, there is an iddues with how we count biomass

# function sum over depth layers in each array slice
collapse_array <- function(mat){
  mat2 <- apply(mat, 3, colSums)
  mat3 <- data.frame(t(mat2))
  colnames(mat3) <- 0:108
  mat3
}

fl <- 'NOAA_Azure/data/GOA_WGS84_V4_final.bgm'
bgm <- rbgm::read_bgm(fl)
goa_sf <- rbgm::box_sf(bgm)
boundary_boxes <- goa_sf %>% sf::st_set_geometry(NULL) %>% filter(boundary == TRUE) %>% pull(box_id) # get boundary boxes
# function to set values in the boundary boxes to NA
setNA <- function(mat) {
  mat2 <- mat
  if(length(dim(mat2))==3) mat2[,(boundary_boxes+1),]<-NA
  if(length(dim(mat2))==2) mat2[(boundary_boxes+1),] <- NA
  mat2
}

# get t3 names, make sure you maintain the same order as t3_fg
t3_names <- grps %>% 
  filter(Code %in% t3_fg) %>% 
  mutate(Code = factor(Code, levels = t3_fg)) %>%
  arrange(Code) %>%
  pull(Name) #%>%
  #sort()

extract_naa <- function(ncfile){
  
  # get run number and corresponding multiplier for F35
  this_idx <- as.numeric(gsub(".*output_", "", gsub(".nc","", ncfile)))
  this_mult <- oy_key %>% filter(idx == this_idx) %>% pull(mult)
  this_run <- oy_key %>% filter(idx == this_idx) %>% pull(run)
  
  this_ncfile <- tidync(ncfile)
  this_ncdata <- nc_open(ncfile)
  
  ts <- ncdf4::ncvar_get(this_ncdata,varid = "t") %>% as.numeric
  tyrs <- ts/(60*60*24*365)
  
  # do one fg at a time, then bring them back together
  naa_frame <- data.frame()
  for (i in 1:length(t3_names)){
    
    fg <- t3_names[i]
    
    # Get numbers by box
    abun_vars <- hyper_vars(this_ncfile) %>% # all variables in the .nc file active grid
      filter(grepl("_Nums",name)) %>% # filter for abundance variables
      filter(grepl(fg,name)) # filter for specific functional group
    
    abun1 <- purrr::map(abun_vars$name,ncdf4::ncvar_get,nc=this_ncdata) %>% 
      lapply(setNA) %>%
      purrr::map(apply,MARGIN=3,FUN=sum,na.rm=T) %>% 
      bind_cols() %>% 
      suppressMessages() %>% 
      set_names(abun_vars$name) %>% 
      mutate(t=tyrs)
    
    abun2 <- abun1 %>%
      pivot_longer(cols = -t,names_to = 'age_group',values_to = 'abun') %>%
      mutate(age=parse_number(age_group)) %>%
      mutate(year = ceiling(t)) %>%
      group_by(year, age_group, age) %>%
      summarise(abun = mean(abun)) %>%
      ungroup() %>%
      mutate(Name = t3_names[i]) %>%
      dplyr::select(year, Name, age, abun)
    
    # get end of the time series (last 5 years average)
    abun3 <- abun2 %>%
      slice_max(year, n = 5) %>%
      group_by(Name, age) %>%
      summarize(abun = mean(abun)) %>%
      ungroup()
    
    naa_frame <- rbind(naa_frame, abun3)
    
  }
  
  # add multiplier for the run
  naa_frame <- naa_frame %>%
    mutate(mult = this_mult,
           run = this_run)
  
  return(naa_frame)
  
}

# apply function to the nc files
naa <- bind_rows(lapply(f35_nc, extract_naa)) # this is slow with 40

# bring in long names
naa <- naa %>%
  left_join(grps %>% select(Name, LongName), by = "Name")

# spaces
naa$LongNamePlot <- gsub(" ", "\n", naa$LongName)

# add scenario information
naa <- naa %>%
  mutate(`F on\narrowtooth` = ifelse(run %in% c("atf","atf_climate"), "Fixed (1/4 FMSY)", "Varying"),
         Climate = ifelse(run %in% c("climate","atf_climate"), "ssp585 (2075-2085)", "Historical (1999)"))

# reorder ATF F
naa$`F on\narrowtooth` <- factor(naa$`F on\narrowtooth`, levels = c("Varying", "Fixed (1/4 FMSY)"))

# add column for age as factor
naa$`Age class` <- factor(naa$age)

# plot
grp1 <- unique(naa$LongNamePlot)[1:6]
naa_plot1 <- naa %>%
  filter(LongNamePlot %in% grp1) %>%
  ggplot(aes(x = mult, y = abun/1000000, color = `Age class`, linetype = `F on\narrowtooth`))+
  geom_line()+
  geom_vline(xintercept = 1, color = "black", linetype = "dotted")+
  scale_color_viridis_d()+
  theme_bw()+
  labs(x = expression(MF[MSY] ~ "multiplier"), y = 'Individuals (millions)')+
  facet_grid2(LongNamePlot~Climate, scales = 'free')+
  theme(strip.text.y = element_text(angle=0))
naa_plot1

grp2 <- unique(naa$LongNamePlot)[7:12]
naa_plot2 <- naa %>%
  filter(LongNamePlot %in% grp2) %>%
  ggplot(aes(x = mult, y = abun/1000000, color = `Age class`, linetype = `F on\narrowtooth`))+
  geom_line()+
  geom_vline(xintercept = 1, color = "black", linetype = "dotted")+
  scale_color_viridis_d()+
  theme_bw()+
  labs(x = expression(MF[MSY] ~ "multiplier"), y = 'Individuals (millions)')+
  facet_grid2(LongNamePlot~Climate, scales = 'free')+
  theme(strip.text.y = element_text(angle=0))
naa_plot2

# make a figure
# ggsave(paste0('NOAA_Azure/results/figures/oy/naa',t,'_1.png'), naa_plot1, width = 7, height = 7)
# ggsave(paste0('NOAA_Azure/results/figures/oy/naa',t,'_2.png'), naa_plot2, width = 7, height = 7)

grp3 <- c("Walleye\npollock", "Pacific\ncod", "Pacific\nhalibut")
naa_plot3 <- naa %>%
  filter(LongNamePlot %in% grp3) %>%
  ggplot(aes(x = mult, y = abun/1000000, color = `Age class`, linetype = `F on\narrowtooth`))+
  geom_line()+
  geom_vline(xintercept = 1, color = "black", linetype = "dotted")+
  scale_color_viridis_d()+
  theme_bw()+
  labs(x = expression(F[MSY] ~ "multiplier"), y = 'Individuals (millions)')+
  facet_grid2(LongNamePlot~Climate, scales = 'free')+
  theme(strip.text.y = element_text(angle=0))
naa_plot3

# ggsave(paste0("NOAA_Azure/results/figures/oy/NAA.png"), naa_plot, width = 8, height = 5)


# Weight at age from nc files -----------------------------------------------------------

extract_waa <- function(ncfile){
  
  # get run number and corresponding multiplier for F35
  this_idx <- as.numeric(gsub(".*output_", "", gsub(".nc","", ncfile)))
  this_mult <- oy_key %>% filter(idx == this_idx) %>% pull(mult)
  this_run <- oy_key %>% filter(idx == this_idx) %>% pull(run)
  
  this_ncfile <- tidync(ncfile)
  this_ncdata <- nc_open(ncfile)
  
  ts <- ncdf4::ncvar_get(this_ncdata,varid = "t") %>% as.numeric
  tyrs <- ts/(60*60*24*365)
  
  # do one fg at a time, then bring them back together
  waa_frame <- data.frame()
  for (i in 1:length(t3_names)){
    
    fg <- t3_names[i]
    
    #Extract from the output .nc file the appropriate reserve N time series variables
    resN_vars <- hyper_vars(this_ncfile) %>% # all variables in the .nc file active grid
      filter(grepl("_ResN",name)) %>% # filter for reserve N
      filter(grepl(fg,name)) # filter for specific functional group
    
    #Extract from the output .nc file the appropriate structural N time series variables
    strucN_vars <- hyper_vars(this_ncfile) %>% # all variables in the .nc file active grid
      filter(grepl("_StructN",name)) %>% # filter for structural N
      filter(grepl(fg,name)) # filter for specific functional group
    
    # Get numbers by box
    abun_vars <- hyper_vars(this_ncfile) %>% # all variables in the .nc file active grid
      filter(grepl("_Nums",name)) %>% # filter for abundance variables
      filter(grepl(fg,name)) # filter for specific functional group
    
    resN <- purrr::map(resN_vars$name,ncdf4::ncvar_get,nc=this_ncdata) 
    strucN <- purrr::map(strucN_vars$name,ncdf4::ncvar_get,nc=this_ncdata)
    nums <-purrr::map(abun_vars$name,ncdf4::ncvar_get,nc=this_ncdata) #numbers by age group,box,layer,time
    totnums <-nums %>% purrr::map(apply,MARGIN=3,FUN=sum) # total numbers by age group, time
    relnums <- purrr::map2(nums,totnums,sweep,MARGIN=3,FUN=`/`) # divide nums by totnums along the time axis to get relative annual nums per age group/box/layer
    
    # add the two matrices to get total nitrogen weight
    rnsn <- purrr::map2(resN,strucN,`+`)
    
    # multiply and sum to get abundance-weighted mean weight at age
    rnsn_summ <- purrr::map2(rnsn,relnums,`*`) %>% 
      purrr::map(apply,MARGIN=3,FUN=sum) %>% # mean total N by time
      bind_cols() %>% # bind age groups elements together
      suppressMessages() %>% 
      set_names(resN_vars$name) %>% 
      mutate(t=tyrs) %>%
      # pivot to long form
      pivot_longer(cols = -t,names_to = 'age_group',values_to = 'totN') %>%
      mutate(age=parse_number(age_group)) %>% 
      mutate(weight=totN*20*5.7/1000000) %>%   # convert totN to weight/individual in kg
      dplyr::filter(t>0) %>%
      mutate(year = ceiling(t)) %>%
      group_by(year, age_group, age) %>%
      summarise(weight = mean(weight)) %>%
      ungroup() %>%
      mutate(Name = t3_names[i]) %>%
      dplyr::select(year, Name, age, weight)
    
    # get end of the time series (last 5 years average)
    rnsn_summ <- rnsn_summ %>%
      slice_max(year, n = 5) %>%
      group_by(Name, age) %>%
      summarize(weight = mean(weight)) %>%
      ungroup()
    
    waa_frame <- rbind(waa_frame, rnsn_summ)
    
  }
  
  # add multiplier for the run
  waa_frame <- waa_frame %>%
    mutate(mult = this_mult,
           run = this_run)
  
  return(waa_frame)
  
}

# apply function to the nc files
waa <- bind_rows(lapply(f35_nc, extract_waa))

# bring in long names
waa <- waa %>%
  left_join(grps %>% select(Name, LongName), by = "Name")

# spaces
waa$LongNamePlot <- gsub(" ", "\n", waa$LongName)

# add scenario information
waa <- waa %>%
  mutate(`F on\narrowtooth` = ifelse(run %in% c("atf","atf_climate"), "Fixed (1/4 FMSY)", "Varying"),
         Climate = ifelse(run %in% c("climate","atf_climate"), "ssp585 (2075-2085)", "Historical (1999)"))

# reorder ATF F
waa$`F on\narrowtooth` <- factor(waa$`F on\narrowtooth`, levels = c("Varying", "Fixed (1/4 FMSY)"))

# add column for age as factor
waa$`Age class` <- factor(waa$age)

# plot
grp1 <- unique(waa$LongNamePlot)[1:6]
waa_plot1 <- waa %>%
  filter(LongNamePlot %in% grp1) %>%
  ggplot(aes(x = mult, y = weight, color = `Age class`, linetype = `F on\narrowtooth`))+
  geom_line()+
  geom_vline(xintercept = 1, color = "black", linetype = "dotted")+
  scale_color_viridis_d()+
  theme_bw()+
  labs(x = expression(MF[MSY] ~ "multiplier"), y = 'Weight (kg)')+
  facet_grid2(LongNamePlot~Climate, scales = 'free')+
  theme(strip.text.y = element_text(angle=0))
waa_plot1

grp2 <- unique(naa$LongNamePlot)[7:12]
waa_plot2 <- waa %>%
  filter(LongNamePlot %in% grp2) %>%
  ggplot(aes(x = mult, y = weight, color = `Age class`, linetype = `F on\narrowtooth`))+
  geom_line()+
  geom_vline(xintercept = 1, color = "black", linetype = "dotted")+
  scale_color_viridis_d()+
  theme_bw()+
  labs(x = expression(MF[MSY] ~ "multiplier"), y = 'Weight (kg)')+
  facet_grid2(LongNamePlot~Climate, scales = 'free')+
  theme(strip.text.y = element_text(angle=0))
waa_plot2

# make a figure
# ggsave(paste0('NOAA_Azure/results/figures/oy/waa',t,'_1.png'), waa_plot1, width = 7, height = 7)
# ggsave(paste0('NOAA_Azure/results/figures/oy/waa',t,'_2.png'), waa_plot2, width = 7, height = 7)

# Below the target of 35% --------------------------------------------------------------

# How many stocks are below 35% B0 for each scenario?
# Use static B0 from Base Scenario for this (but why?)
below_target <- ms_yield_long %>%
  filter(type == "Biomass") %>%
  filter(!(run %in% c("atf", "atf_climate") & LongName == "Arrowtooth flounder")) %>%
  left_join(b0, by = "LongName") %>%
  mutate(depletion = mt / b0) %>%
  mutate(below_target = ifelse(depletion < 0.35, 1, 0)) %>%
  group_by(run, mult) %>%
  mutate(n_below_target = sum(below_target)) %>%
  mutate(prop_below_target = n_below_target / length(unique(LongName))) %>%
  ungroup() %>%
  select(run, mult, prop_below_target) %>%
  distinct()

# add scenario information
below_target <- below_target %>%
  mutate(`F on\narrowtooth` = ifelse(run %in% c("atf","atf_climate"), "Fixed (1/4 FMSY)", "Varying"),
         Climate = ifelse(run %in% c("climate","atf_climate"), "ssp585 (2075-2085)", "Historical (1999)"))

# reorder ATF F
below_target$`F on\narrowtooth` <- factor(below_target$`F on\narrowtooth`, levels = c("Varying", "Fixed (1/4 FMSY)"))

p_below_target <- below_target %>%
  ggplot(aes(x = mult, y = prop_below_target*100))+
  geom_col()+
  theme_bw()+
  geom_vline(xintercept = 1, color = "black", linetype = "dotted")+
  labs(x = expression(MF[MSY] ~ "multiplier"), y = "SSB < B35% (%)") +
  facet_grid(Climate~`F on\narrowtooth`)
p_below_target

# ggsave(paste0("NOAA_Azure/results/figures/oy/below_target_nprb.png"), p_below_target, width = 5, height = 4)








# # Diagnostics: biomass of unselected age classes --------------------------
# 
# # look at unselected biomass in both ms and ss
# ms_unselected_list <- list()
# 
# for(i in 1:nrow(f35_key)){
#   
#   this_catch_file <- f35_catch_files[i]
#   this_run <- as.numeric(substr((gsub(paste0("NOAA_Azure/results/f35/",batch,"/outputGOA0"), "", this_catch_file)), 1, 4))
#   this_mult <- f35_key %>% filter(run == this_run) %>% pull(mult)
#   
#   # extract tables from results
#   biomage <- read.table(f35_biomass_files[i], sep = " ", header = T)
#   catch <- read.table(this_catch_file , sep = " ", header = T)
#   
#   # now extract data
#   # SSB to plot and report in tables
#   unselected_biomass <- biomage %>% 
#     slice_tail(n = 5) %>% # use last xxx years
#     summarise(across(-"Time", ~ mean(.x, na.rm = TRUE))) %>%
#     ungroup() %>%
#     pivot_longer(everything(), names_to = 'Code.Age', values_to = 'biomass_mt') %>%
#     separate_wider_delim(Code.Age, delim = '.', names = c('Code', 'Age')) %>%
#     filter(Code %in% t3_fg) %>%
#     left_join(selex, by = 'Code') %>%
#     mutate(idx = as.numeric(Age) - as.numeric(age_class_selex)) %>%
#     filter(is.na(idx) | idx < 0) %>%
#     group_by(Code) %>%
#     summarise(biomass_mt = sum(biomass_mt)) %>%
#     ungroup()
#   
#   # # calculate realized F after 1 year of data
#   # # get initial biomass for the selected age classes
#   biom_age_t1 <- biomage %>% 
#     filter(Time == 0) %>% 
#     pivot_longer(-Time, names_to = 'Code.Age', values_to = 'biomass') %>%
#     separate_wider_delim(Code.Age, delim = '.', names = c('Code', 'Age')) %>%
#     left_join(selex, by = 'Code') %>%
#     mutate(idx = as.numeric(Age) - as.numeric(age_class_selex)) %>%
#     filter(is.na(idx) | idx >= 0) %>%
#     group_by(Code) %>%
#     summarise(biomass = sum(biomass)) %>%
#     ungroup() %>% 
#     filter(Code %in% t3_fg)
#   # 
#   # # catch (one time step after biomass: how much did we catch in this time?)
#   catch_t1 <- catch %>% 
#     select(Time, all_of(t3_fg)) %>% 
#     filter(Time == 365) %>%
#     summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) %>%
#     pivot_longer(-Time, names_to = 'Code', values_to = 'catch') %>%
#     select(-Time)
#   # 
#   # # calc realized f
#   f_t1 <- biom_age_t1 %>% left_join(catch_t1, by = "Code") %>%
#     mutate(exp_rate = catch/biomass,
#            f = -log(1-exp_rate)) %>%#,
#     #fidx = fidx) %>% # need this for joining later on
#     select(Code, f)#, fidx) 
#   
#   # # bind all
#   f_frame <- f_t1 %>%
#     left_join(unselected_biomass) %>%
#     #left_join(catch_vals) %>%
#     mutate(mult = this_mult)
#   
#   # add to multispecies yield list
#   ms_unselected_list[[i]] <- f_frame
# }
# 
# ms_unselected_df <- bind_rows(ms_unselected_list)
# 
# ms_unselected_df <- ms_unselected_df %>%
#   left_join(grps %>% select(Code, LongName), by = "Code") %>%
#   mutate(experiment = "ms") %>%
#   select(Code, LongName, experiment, f, biomass_mt)
# 
# # create empty list to fill with data frame for the yield curve
# folder_path <- here("NOAA_Azure","results","pre-processing",this_job,"results")
# # list the results files
# results_list <- list.files(folder_path, full.names = T)
# 
# ss_unselected_list <- list()
# for(i in 1:length(results_list)){ # TODO: turn all of these into a function with a few args (e.g., fished, unselected, top preds, ms, ss, etc.)
#   
#   print(paste("Doing",i))
#   
#   # load results
#   this_result <- readRDS(results_list[i])
#   
#   # this idx - CAREFUL! This is not the same as i because of alphabetical sorting of the reuslts files
#   this_idx <- as.numeric(names(this_result))
#   
#   # first identify the species and the level of fishing for this run. These are unrelated from runname
#   sp <- f_lookup %>% filter(idx == this_idx) %>% pull(species)
#   fidx <- f_lookup %>% filter(idx == this_idx) %>% pull(fidx)
#   
#   # rename the result object in the list to avoid problems with indexing
#   names(this_result) <- "res"
#   
#   # extract tables from results
#   biomage <- this_result$res[[paste0("biomage_",this_idx)]]
#   catch <- this_result$res[[paste0("catch_",this_idx)]]
#   
#   # now extract data
#   # SSB to plot and report in tables
#   unselected_biomass <- biomage %>% 
#     slice_tail(n = 5) %>% # use last xxx years
#     summarise(across(-"Time", ~ mean(.x, na.rm = TRUE))) %>%
#     ungroup() %>%
#     pivot_longer(everything(), names_to = 'Code.Age', values_to = 'biomass_mt') %>%
#     separate_wider_delim(Code.Age, delim = '.', names = c('Code', 'Age')) %>%
#     filter(Code == sp) %>%
#     left_join(selex, by = 'Code') %>%
#     mutate(idx = as.numeric(Age) - as.numeric(age_class_selex)) %>%
#     filter(is.na(idx) | idx < 0) %>%
#     group_by(Code) %>%
#     summarise(biomass_mt = sum(biomass_mt)) %>%
#     ungroup() %>%
#     pull(biomass_mt)
#   
#   # # calculate realized F after 1 year of data
#   # # get initial biomass for the selected age classes
#   biom_age_t1 <- biomage %>% 
#     filter(Time == 0) %>% 
#     pivot_longer(-Time, names_to = 'Code.Age', values_to = 'biomass') %>%
#     separate_wider_delim(Code.Age, delim = '.', names = c('Code', 'Age')) %>%
#     left_join(selex, by = 'Code') %>%
#     mutate(idx = as.numeric(Age) - as.numeric(age_class_selex)) %>%
#     filter(is.na(idx) | idx >= 0) %>%
#     group_by(Code) %>%
#     summarise(biomass = sum(biomass)) %>%
#     ungroup() %>% 
#     filter(Code == sp)
#   # 
#   # # catch (one time step after biomass: how much did we catch in this time?)
#   catch_t1 <- catch %>% 
#     select(Time, all_of(sp)) %>% 
#     filter(Time == 365) %>%
#     summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) %>%
#     pivot_longer(-Time, names_to = 'Code', values_to = 'catch') %>%
#     select(-Time)
#   # 
#   # # calc realized f
#   f_t1 <- biom_age_t1 %>% left_join(catch_t1) %>%
#     mutate(exp_rate = catch/biomass,
#            f = -log(1-exp_rate),
#            fidx = fidx) %>% # need this for joining later on
#     select(Code, f, fidx) 
#   
#   # # bind all
#   f_frame <- f_t1 %>%
#     mutate(biomass = unselected_biomass)#,
#   #catch = catch_val)
#   
#   ss_unselected_list[[i]] <- f_frame
#   
#   # # write out to be then brought together with all other runs
#   # write.csv(f_frame, paste(outdir,paste(sp,fidx,'f.csv',sep='_'), sep = "/"), row.names = F)
#   
# }
# 
# ss_unselected_df <- bind_rows(ss_unselected_list) %>%
#   left_join(grps %>% select(Code, LongName), by= "Code") %>%
#   mutate(experiment = "ss") %>%
#   select(Code, LongName, experiment, f, biomass) %>%
#   rename(biomass_mt = biomass)
# 
# # plot
# unselected_plot <- ms_unselected_df %>%
#   rbind(ss_unselected_df) %>%
#   mutate(experiment = ifelse(experiment == "ms", "Multispecies", "Single-species")) %>%
#   ggplot(aes(x = f, y = biomass_mt/1000, color = experiment))+
#   geom_line()+
#   geom_point(size = 2)+
#   scale_color_viridis_d(begin = 0.2, end = 0.8)+
#   geom_vline(data = fmsy, aes(xintercept = FMSY, group = LongName), linetype = 'dashed', color = 'orange')+
#   geom_vline(data = atlantis_fmsy, aes(xintercept = atlantis_fmsy, group = LongName), linetype = 'dashed', color = 'blue')+
#   theme_bw()+
#   scale_y_continuous(limits = c(0, NA))+
#   labs(x = 'F as perceived by the model', y = '1000\'s of tons')+
#   facet_wrap(~LongName, scales = 'free', ncol = 3)+
#   theme(strip.text.y = element_text(angle=0))
# unselected_plot
# 
# ggsave(paste0("NOAA_Azure/results/figures/",batch,"/unselected.png"), unselected_plot, width = 8, height = 6)
# 
# # save the file as a csv to compare it to other runs
# write.csv(ms_unselected_df %>% mutate(batch = batch), 
#           paste0('NOAA_Azure/results/for_comp/',batch,'_unselected.csv'), row.names = F)
# 
# # In the MS scenario, the biomass of the unselected age classes is in general higher across the board. 
# # One possible explanation is that their adult predators are getting fished out
# # This does not really hold if we look at their biomass, which is not lower than in the ss runs
# # So, in the MS runs we have higher biomass and higher catches across all T3 species.
# # In other words, when we fish a stock at FMSY, if we fish all stocks at FMSY we get more of all of them
# # Since they can't all be winners, who is losing?
# # look at top predators and forage fish across the ms runs
# 



# # Diagnostics: M -----------------------------------------------------------------------
# 
# f35_mort_files <- list.files(f35_path, pattern = "Mort", full.names = T)
# 
# # for each group, calculate the relationship between F and M as per mort.txt file
# f_to_m <- data.frame()
# for(i in 1:nrow(f35_key)){
#   
#   this_m_file <- f35_mort_files[i]
#   this_run <- as.numeric(substr(gsub(paste0("NOAA_Azure/results/f35/",batch,"/outputGOA0"), "", this_m_file), 1, 4))
#   this_mult <- f35_key %>% filter(run == this_run) %>% pull(mult)
#   
#   # extract tables from results
#   m <- read.table(f35_mort_files[i], sep = " ", header = T)
#   
#   # the values of M and F are only meaningful relative to one another
#   m_long <- m %>%
#     slice_tail(n = 5) %>% # use last xxx years
#     summarise(across(-"Time", ~ mean(.x, na.rm = TRUE))) %>%
#     ungroup() %>%
#     pivot_longer(everything(), names_to = 'Code.MType', values_to = 'rawM') %>%
#     separate_wider_delim(Code.MType, delim = '.', names = c('Code', 'MType'))
#   
#   # now get M and F on a line each and get the proportion
#   
#   this_f_to_m <- m_long %>%
#     pivot_wider(id_cols = Code, names_from = MType, values_from = rawM) %>%
#     mutate(f_mult = M / `F`,
#            mult = this_mult) %>%
#     dplyr::select(Code, mult, f_mult)
#   
#   # add to multispecies yield list
#   f_to_m <- rbind(f_to_m, this_f_to_m)
# } 
# 
# # join to data frame with f values from the multispecies experiment
# m_df <- df_mult %>%
#   left_join(f_to_m, by = c("Code","mult")) %>%
#   mutate(m = f * f_mult)
# 
# # plot
# mp <- m_df %>%
#   filter(mult > 0) %>%
#   ggplot(aes(x = mult, y = m))+
#   geom_line()+
#   theme_bw()+
#   facet_wrap(~LongName, scales = "free")
# mp
# 
# # M values are all over the place. Also we don't get to plot the values for f = 0, because the multiplier is not calculated
# # if we need a baseline, apply very low mFC
# # start from mult = 0.2
# # apart from the scale that is all over the place, M decreases with increasing F for ATF, FFD, FHS, COD, RFS, and POL - likely all ATF's prey items
# # More surprisingly, M increases with F for HAL, POP, REX, RFP, SBF. Need to investigate, but one possible explanation is that their (very few) predators (like SSL and KWR) switch from some of the depleted groundfish to these other species
# 
# # to convey the same message without worrying about the scale, do this relative to mult = 2.0 (or 0.4)
# if(maxmult == 2) {minmult <- 0.2} else {minmult <- 0.4}
# m_df_mult02 <- m_df %>% filter(mult == minmult) %>% dplyr::select(LongName, m) %>% rename(m02 = m) %>% distinct()
# m_df_rel <- m_df %>%
#   left_join(m_df_mult02, by = "LongName") %>%
#   mutate(rel_m = m / m02)
# 
# # plot
# rel_m_plot <- m_df_rel %>%
#   filter(mult > 0) %>%
#   ggplot(aes(x = mult, y = rel_m))+
#   geom_line()+
#   geom_hline(yintercept = 1, color = "red", linetype = "dashed")+
#   theme_bw()+
#   labs(x = "F35 permutation", y = "Relative change in M (from mult = 0.2)")+
#   facet_wrap(~LongName)
# rel_m_plot
# 
# ggsave(paste0("NOAA_Azure/results/figures/", batch, "/relative_m.png"), rel_m_plot, width = 8, height = 4)
# 
# write.csv(m_df_rel %>% mutate(batch = batch), 
#           paste0('NOAA_Azure/results/for_comp/',batch,'_rel_m.csv'), row.names = F)
# 

# # Diagnostics: equilibrium ------------------------------------------------
# 
# # one thing to check is - are we at equilibrium here or not?
# eq_catch_list <- list()
# for(i in 1:nrow(f35_key)){
#   
#   this_catch_file <- f35_catch_files[i]
#   this_run <- as.numeric(substr((gsub(paste0("NOAA_Azure/results/f35/",batch,"/outputGOA0"), "", this_catch_file)), 1, 4))
#   this_mult <- f35_key %>% filter(run == this_run) %>% pull(mult)
#   
#   this_catch <- read.table(this_catch_file, sep=" ", header = T)
#   this_catch <- this_catch %>%
#     slice_tail(n = 20) %>% # look at the  last 20 years without averageing
#     dplyr::select(Time, all_of(t3_fg)) %>%
#     pivot_longer(-Time, names_to = "Code", values_to = "mt") %>%
#     mutate(mult = this_mult)
#   
#   eq_catch_list[[i]] <- this_catch
# }
# 
# eq_catch_df <- bind_rows(eq_catch_list) %>%
#   mutate(yr = Time / 365)
# 
# # plot
# eq_catch_df %>%
#   ggplot(aes(x = yr, y = mt, color = factor(mult)))+
#   geom_line()+
#   theme_bw()+
#   facet_wrap(~Code, scales = "free_y")
# 
# # catches in slight decline over time but not to extent of suggesting a collapse immediately after the end of the run (50 years)
# # now biomass
# eq_biom_list <- list()
# for(i in 1:nrow(f35_key)){
#   
#   this_biom_file <- f35_biomass_files[i]
#   this_run <- as.numeric(substr((gsub(paste0("NOAA_Azure/results/f35/",batch,"/outputGOA0"), "", this_biom_file)), 1, 4))
#   this_mult <- f35_key %>% filter(run == this_run) %>% pull(mult)
#   
#   this_biom <- read.table(this_biom_file, sep=" ", header = T)
#   this_biom <- this_biom %>%
#     slice_tail(n = 20) %>% # look at the  last 20 years without averageing
#     pivot_longer(-Time, names_to = 'Code.Age', values_to = 'biomass_mt') %>%
#     separate_wider_delim(Code.Age, delim = '.', names = c('Code', 'Age')) %>%
#     filter(Code %in% t3_fg) %>%
#     #left_join(mat, by = 'Code') %>%
#     # mutate(idx = as.numeric(Age) - as.numeric(age_class_mat)) %>%
#     # filter(is.na(idx) | idx >= 0) %>%
#     left_join(fspb, by = c('Code','Age')) %>%
#     mutate(biomass_mt = biomass_mt * fspb) %>%
#     group_by(Time, Code) %>%
#     summarize(mt = sum(biomass_mt))%>%
#     mutate(mult = this_mult)
#   
#   
#   eq_biom_list[[i]] <- this_biom
# }
# 
# eq_biom_df <- bind_rows(eq_biom_list) %>%
#   mutate(yr = Time / 365)
# 
# # plot
# eq_biom_df %>%
#   ggplot(aes(x = yr, y = mt, color = factor(mult)))+
#   geom_line()+
#   theme_bw()+
#   facet_wrap(~Code, scales = "free_y")
# 
# # similarly to catch, we are fairly close to equilibrium for SSB and there are no signs of imminent collapse
# # This points to excessive productivity, possibly mediated by excessive recruitment
# # It may be worth repeating the set with (much) lower steepness. Reducing steepness did not help when stocks had mature age classes below the knife edge, but it may help in cases where that is not the case.
# 
# # Plot global yield but from single-species
# # to make sense, this will be the sum of the single-species yield at FOFL
# # The original range explored was in 8 increments between 0 and 2*FOFL
# # If we bring in idx, we get to explore these rough permutations of F for each stock:
# # 0.0000000 0.2857143 0.5714286 0.8571429 1.1428571 1.4285714 1.7142857 2.0000000
# 
# idx_to_mult <- data.frame("fidx" = 1:8, "mult" = seq(0,maxmult,length.out=8))
# 
# global_yield_ss <- ss_yield_long_fidx %>%
#   filter(type == "Catch") %>%
#   left_join(f_lookup, by = c("Code"="species","fidx")) %>%
#   left_join(idx_to_mult, by = "fidx") %>%
#   ggplot(aes(x = mult, y = mt / 1000, fill = LongName))+
#   geom_bar(stat = "identity", position = "stack")+
#   scale_fill_viridis_d()+
#   scale_x_continuous(breaks = seq(0,maxmult,length.out=11))+
#   scale_y_continuous(limits = c(0,1100))+
#   geom_hline(yintercept = 116, color = "red", linetype = "dashed")+
#   geom_hline(yintercept = 800, color = "red", linetype = "dashed")+
#   #geom_hline(yintercept = atlantis_oy / 1000, color = "blue", linetype = "dashed")+
#   theme_bw()+
#   labs(x = "FOFL (SS) multiplier", y = "Catch (1000 mt)", fill = "Stock")
# global_yield_ss
# 
# ggsave(paste0("NOAA_Azure/results/figures/", batch, "/global_yield_SS.png"), global_yield_ss, width = 8, height = 6)
# 
# 
# # Recruitment -------------------------------------------------------------
# 
# # Above we have biomass of the unselected age classes, but to have more of a pulse on what happens with recruitment under different
# # configurations we need to extract numbers at age of the first age class from the nc files
# 
# # list nc files
# nc_files <- list.files(paste0("NOAA_Azure/results/f35/",batch), pattern = ".nc", full.names = T)
# 
# # get t3 names, make sure you maintain the same order as t3_fg
# t3_names <- grps %>% 
#   filter(Code %in% t3_fg) %>% 
#   mutate(Code = factor(Code, levels = t3_fg)) %>%
#   arrange(Code) %>%
#   pull(Name) %>%
#   sort()
# 
# # function to pull recruits for the T3 stocks
# extract_recruits <- function(ncfile){
#   
#   # get run number and corresponding multiplier for F35
#   this_run <- as.numeric(gsub("_test.nc", "", gsub(".*outputGOA0","", ncfile)))
#   this_mult <- f35_key %>% filter(run == this_run) %>% pull(mult)
#   
#   this_ncfile <- tidync(ncfile)
#   this_ncdata <- nc_open(ncfile)
#   
#   # get variable names for recruitment (numbers of age 0)
#   t3_pattern <- paste(t3_names, collapse = "|")
#   recruit_vars <- this_ncfile %>%
#     hyper_vars() %>%
#     filter(grepl("1_Nums", name)) %>%
#     filter(grepl(t3_pattern, name)) %>%
#     arrange(name)
#   
#   # pull data
#   recs <-purrr::map(recruit_vars$name,ncdf4::ncvar_get,nc=this_ncdata) #numbers by fg,box,layer,time
#   
#   # for each group, collapse all spatial dimensions and get 251 values, then cut to the last (5x5) time steps and average for an average of the last 5 years
#   terminal_recruits_list <- lapply(recs, function(arr) {
#     # dims are (7,109,251), i.e. (z,b,t)
#     # Sum over the 'z' dimension (1st dimension)
#     sum_over_z <- apply(arr, c(2, 3), sum)
#     
#     # Sum over the 'b' dimension (now the 1st dimension)
#     sum_over_b <- apply(sum_over_z, 2, sum)
#     
#     # now take avg of last 25 vals (for last 5 years, 5 data points a year)
#     termrec <- mean(tail(sum_over_b), 25)
#     
#     return(termrec)
#   })
#   
#   terminal_recruits <- unlist(terminal_recruits_list)
#   
#   # make a data frame
#   df_rec <- data.frame("Name" = t3_names, "Code" = t3_fg, "R" = terminal_recruits)
#   
#   # tie in the respective mult for this run
#   df_rec <- df_rec %>% mutate(mult = this_mult)
#   
#   return(df_rec)
#   
# }
# 
# # apply function to the nc files
# rec_by_mult <- bind_rows(lapply(nc_files, extract_recruits))
# 
# # tie to SSB data
# rec_by_mult_df <- rec_by_mult %>%
#   left_join(ms_yield_df, by = c("Code", "mult"))
# 
# # make R0 df
# r0 <- rec_by_mult_df %>%
#   filter(mult == 0) %>%
#   select(Code, R) %>%
#   rename(r0 = R)
# 
# # join to rec df and get R as proportion of R0
# rec_by_mult_df <- rec_by_mult_df %>%
#   left_join(r0, by = "Code") %>%
#   mutate(rprop = R / r0)
# 
# # join b0 df and get SSB as proportion of B0 (depletion)
# rec_by_mult_df <- rec_by_mult_df %>%
#   left_join((b0 %>% 
#                filter(experiment == "ms") %>%
#                left_join(grps %>% 
#                            select(LongName, Code), by = "LongName")) %>%
#               select(Code, b0), 
#             by ="Code") %>%
#   mutate(depletion = biomass_mt / b0)
# 
# # plot R vs depletion
# recruitment_plot <- rec_by_mult_df %>%
#   ggplot(aes(x = depletion, y = rprop))+
#   geom_line()+
#   geom_point()+
#   geom_hline(yintercept = 0.5, linetype = "dashed")+
#   geom_vline(xintercept = 0.125, color = "red", linetype = "dashed")+
#   theme_bw()+
#   scale_x_continuous(limits = c(0,1))+
#   scale_y_continuous(limits = c(0,1))+
#   labs(x = "SSB/S0", y = "R/R0")+
#   facet_wrap(~Name)
# 
# ggsave(paste0("NOAA_Azure/results/figures/", batch, "/SR_curve.png"), recruitment_plot, width = 8, height = 6)
# 
# # write out file for comparison plots
# write.csv(rec_by_mult_df %>% mutate(batch = batch), 
#           paste0('NOAA_Azure/results/for_comp/',batch,'_sr.csv'), row.names = F)
# 
