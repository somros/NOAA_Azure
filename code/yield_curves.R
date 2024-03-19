# Alberto Rovellini
# 9/25/2023
# This script produces yield curves from the individual tables produced in the parallel runs
# This script is not run in parallel, just once after the batch that extracted catch and biomass and F etc.
# This is to create panels with three curves each: 
# 1. original
# 2. reduced age at first selex for POP, HAL, RFP, RFS
# 3. reduced steepness for POP, HAL, SBF, RFP, RFS

# loop over batches of output csv files depending on label ('', 'v2', 'v3')

library(tidyverse)
library(readxl)
library(ggh4x)
library(viridis)

# read in Groups.csv file
grp_path <- here('NOAA_Azure/data/GOA_Groups.csv') # functional groups
grps <- read.csv(grp_path)

# list all csv files we need to read
this_job <- "job20240109052342"
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

# produce a dataset of 35% B0, to be used to plot horizontal lines that will intersect the yield curve
# but B35% will now be different between runs with different steepness? Not between runs with smaller selex in theory
b35 <- f_df %>%
  group_by(LongName) %>%
  slice_min(f) %>% # producing it off of v1 (v2 should be the same, v3 should be similar)
  ungroup() %>%
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
    filter(type == 'Biomass') %>%
    slice_min(f) %>%
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
# write.csv(atlantis_fmsy, "NOAA_Azure/data/f35_vector_PROXY_OY_SS.csv", row.names = F)

# annotations for the plots (atlantis depletion)
annotations <- atlantis_fmsy %>% 
  mutate(depletion=round(depletion,digits=2), atlantis_fmsy=round(atlantis_fmsy,digits = 2))

# plot
# only focus on POP, HAL, SBF, RFP, RFS
f_plot <- f_df %>%
  #filter(Code %in% c('POP','HAL','RFP','RFS')) %>%
  ggplot(aes(x = f, y = mt/1000))+
  geom_line()+
  geom_point(size = 2)+
  geom_vline(data = fmsy, aes(xintercept = FMSY, group = LongName), linetype = 'dashed', color = 'orange')+
  geom_vline(data = atlantis_fmsy, aes(xintercept = atlantis_fmsy, group = LongName), linetype = 'dashed', color = 'blue')+
  geom_hline(data = b35 , aes(yintercept = b35/1000, group = LongName), linetype = 'dashed', color = 'red')+
  geom_hline(data = atlantis_fmsy %>% mutate(type = 'Biomass'), 
             aes(yintercept = b_fmsy/1000, group = LongName), linetype = 'dashed', color = 'blue')+
  geom_text(data = annotations %>% mutate(type = 'Biomass'), 
            aes(x=Inf,y=Inf,hjust=1,vjust=1,label=paste0('Depletion=',depletion)), color = 'blue')+
  theme_bw()+
  scale_y_continuous(limits = c(0, NA))+
  labs(x = 'Fishing mortality (F)', y = '1000\'s of tons')+
  facet_grid2(LongName~type, scales = 'free', independent = 'all')+
  theme(strip.text.y = element_text(angle=0))
f_plot

# make a figure
t <- format(Sys.time(),'%Y-%m-%d %H-%M-%S')
ggsave(paste0('NOAA_Azure/results/figures/oy_paper/yield_curves',t,'.png'), f_plot, width = 7, height = 9)

# make figures for OY paper
# plot
f_df_ms <- f_df

f_df_ms$LongNamePlot <- gsub(" ", "\n", f_df_ms$LongName)
fmsy$LongNamePlot <- gsub(" ", "\n", fmsy$LongName)
atlantis_fmsy$LongNamePlot <- gsub(" ", "\n", atlantis_fmsy$LongName)
b35$LongNamePlot <- gsub(" ", "\n", b35$LongName)
annotations$LongNamePlot <- gsub(" ", "\n", annotations$LongName)

grp1 <- unique(f_df_ms$LongNamePlot)[1:6]
f_plot1 <- f_df_ms %>%
  filter(LongNamePlot %in% grp1) %>%
  ggplot(aes(x = f, y = mt/1000))+
  geom_line()+
  geom_point(size = 2)+
  geom_vline(data = fmsy %>% filter(LongNamePlot %in% grp1), aes(xintercept = FMSY, group = LongNamePlot), linetype = 'dashed', color = 'orange')+
  geom_vline(data = atlantis_fmsy %>% filter(LongNamePlot %in% grp1), aes(xintercept = atlantis_fmsy, group = LongNamePlot), linetype = 'dashed', color = 'blue')+
  #geom_hline(data = b35 %>% filter(LongNamePlot %in% grp1), aes(yintercept = b35/1000, group = LongNamePlot), linetype = 'dashed', color = 'red')+
  geom_hline(data = atlantis_fmsy %>% filter(LongNamePlot %in% grp1) %>% mutate(type = 'Biomass'),
             aes(yintercept = b_fmsy/1000, group = LongNamePlot), linetype = 'dashed', color = 'blue')+
  geom_text(data = annotations %>% filter(LongNamePlot %in% grp1) %>% mutate(type = 'Biomass'),
            aes(x=Inf,y=Inf,hjust=1,vjust=1,label=paste0('Depletion=',depletion)), color = 'blue')+
  theme_bw()+
  scale_y_continuous(limits = c(0, NA))+
  labs(x = 'Fishing mortality (F)', y = '1000\'s of tons')+
  facet_grid2(LongNamePlot~type, scales = 'free', independent = 'all')+
  theme(strip.text.y = element_text(angle=0))

grp2 <- unique(f_df_ms$LongNamePlot)[7:12]
f_plot2 <- f_df_ms %>%
  filter(LongNamePlot %in% grp2) %>%
  ggplot(aes(x = f, y = mt/1000))+
  geom_line()+
  geom_point(size = 2)+
  geom_vline(data = fmsy %>% filter(LongNamePlot %in% grp2), aes(xintercept = FMSY, group = LongNamePlot), linetype = 'dashed', color = 'orange')+
  geom_vline(data = atlantis_fmsy %>% filter(LongNamePlot %in% grp2), aes(xintercept = atlantis_fmsy, group = LongNamePlot), linetype = 'dashed', color = 'blue')+
  #geom_hline(data = b35 %>% filter(LongNamePlot %in% grp2), aes(yintercept = b35/1000, group = LongNamePlot), linetype = 'dashed', color = 'red')+
  geom_hline(data = atlantis_fmsy %>% filter(LongNamePlot %in% grp2) %>% mutate(type = 'Biomass'),
             aes(yintercept = b_fmsy/1000, group = LongNamePlot), linetype = 'dashed', color = 'blue')+
  geom_text(data = annotations %>% filter(LongNamePlot %in% grp2) %>% mutate(type = 'Biomass'),
            aes(x=Inf,y=Inf,hjust=1,vjust=1,label=paste0('Depletion=',depletion)), color = 'blue')+
  theme_bw()+
  scale_y_continuous(limits = c(0, NA))+
  labs(x = 'Fishing mortality (F)', y = '1000\'s of tons')+
  facet_grid2(LongNamePlot~type, scales = 'free', independent = 'all')+
  theme(strip.text.y = element_text(angle=0))

# make a figure
# t <- format(Sys.time(),'%Y-%m-%d %H-%M-%S')
ggsave(paste0('NOAA_Azure/results/figures/oy_paper/yield_curves',t,'_OY_1.png'), f_plot1, width = 7, height = 7)
ggsave(paste0('NOAA_Azure/results/figures/oy_paper/yield_curves',t,'_OY_2.png'), f_plot2, width = 7, height = 7)

# # make a figure for the OY paper
# focus on key stocks: POL, COD, ATF, HAL, SBF, POP, FFS
key_grps <- grps %>% filter(Code %in% c("POL", "COD", "ATF", "HAL", "SBF", "POP", "FFS")) %>% pull(LongName)
p_ms <- f_df_ms %>%
  filter(LongName %in% key_grps) %>%
  ggplot(aes(x = f, y = mt/1000))+
  geom_line()+
  geom_point(size = 1.5)+
  geom_vline(data = fmsy %>% filter(LongName %in% key_grps), aes(xintercept = FMSY, group = LongNamePlot), linetype = 'dashed', color = 'orange')+
  geom_vline(data = atlantis_fmsy %>% filter(LongName %in% key_grps), aes(xintercept = atlantis_fmsy, group = LongNamePlot), linetype = 'dashed', color = 'blue')+
  geom_hline(data = atlantis_fmsy %>% filter(LongName %in% key_grps) %>% mutate(type = 'Biomass'),
             aes(yintercept = b_fmsy/1000, group = LongNamePlot), linetype = 'dashed', color = 'blue')+
  geom_text(data = annotations %>% filter(LongName %in% key_grps) %>% mutate(type = 'Biomass'),
            aes(x=Inf,y=Inf,hjust=1,vjust=1,label=paste0('Depletion=',depletion)), color = 'blue')+
  theme_bw()+
  scale_y_continuous(limits = c(0, NA))+
  labs(x = 'Fishing mortality (F)', y = '1000\'s of tons')+
  facet_grid2(LongNamePlot~type, scales = 'free', independent = 'all')+
  theme(strip.text.y = element_text(angle=0))
ggsave(paste0('NOAA_Azure/results/figures/oy_paper/biom_catch_key_stocks.png'), p_ms, width = 7, height = 7)

# make simple catch curves for methods figure
forplot <- "Pacific halibut"
tt <- f_df_ms %>%
  filter(LongName == forplot, type == "Catch") %>%
  ggplot(aes(x = f, y = mt/1000))+
  geom_line(linewidth = 1.5)+
  geom_vline(data = atlantis_fmsy %>% filter(LongName == forplot), aes(xintercept = atlantis_fmsy, group = LongNamePlot), 
             linetype = 'dashed', color = 'blue', linewidth = 1.5)+
  theme_bw()+
  scale_y_continuous(limits = c(0, NA))+
  labs(x = 'Fishing mortality (F)', y = '1000\'s of tons')
ggsave(paste0("NOAA_Azure/results/figures/oy_paper/",forplot,"_methods.png"),tt,width=3, height = 2)


# Yield functions ---------------------------------------------------------
ss_yield_long <- f_df %>%
  select(Code, LongName, f, fidx, type, mt)

# get b0
b0 <- ss_yield_long %>% 
  group_by(Code) %>%
  slice_min(f) %>% 
  ungroup() %>%
  filter(type == "Biomass") %>% 
  dplyr::select(LongName, mt) %>% 
  rename(b0 = mt)

# get max yield
ymax <- ss_yield_long %>% 
  filter(type == "Catch") %>% 
  group_by(LongName) %>%
  slice_max(mt) %>%
  ungroup() %>%
  dplyr::select(LongName, mt, f) %>% 
  rename(ymax = mt) 

# we are plotting yield fraction against depletion
yield_func <- ss_yield_long %>%
  drop_na() %>%
  dplyr::select(LongName, type, f, mt) %>%
  pivot_wider(id_cols = c(LongName, f), names_from = type, values_from = mt) %>%
  left_join(b0, by = c("LongName")) %>% # if you are keep static reference point
  mutate(depletion = Biomass / b0) %>%
  left_join(ymax %>% select(-f), by = c("LongName")) %>%
  mutate(yfrac = Catch / ymax) %>%
  #mutate(experiment = ifelse(experiment == "ms", "Multispecies", "Single-species")) %>%
  dplyr::select(LongName, yfrac, depletion, f)

# prepare data frames to write the following quantities on the plot:
# final depletion, final yield fraction, biomass corresponding to final depletion, biomass corresponding to final yield fraction
yfun_terminal <- yield_func %>%
  group_by(LongName) %>%
  slice_max(f) %>% # for each group, get highest f
  ungroup() %>%
  mutate(depletion = depletion * 100,
         yfrac = yfrac * 100) # turn depletion and yield to percentages for easier interpretation

annotations <- b0 %>%
  left_join(ymax, by = "LongName") %>%
  left_join(yfun_terminal, by = "LongName") %>%
  mutate(catch = ymax / 100 * yfrac / 1000,
         ssb = b0 / 100 * depletion / 1000) 

# plot
yield_func_plot <- yield_func %>%
  #filter(LongName %in% c("Walleye pollock", "Pacific cod", "Arrowtooth flounder", "Pacific halibut")) %>%
  ggplot(aes(x = depletion, y = yfrac))+
  geom_point()+
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

ggsave(paste0("NOAA_Azure/results/figures/oy_paper/yield_functions_ss.png"), yield_func_plot, width = 8, height = 6.5)

# only for key groups
# y_p_ms <- yield_func %>%
#   filter(LongName %in% key_grps) %>%
#   ggplot(aes(x = depletion, y = yfrac))+
#   geom_point()+
#   geom_line()+
#   scale_x_reverse()+
#   geom_text(data = annotations %>% filter(LongName %in% key_grps),
#             aes(x = 0.5, y = 0.5, hjust=0.5, vjust=1,
#                 label=paste0('Depletion(%)=', round(depletion,2),
#                              '\n',
#                              'Yield fraction(%)=', round(yfrac, 2),
#                              '\n',
#                              'SSB(1000mt)=', round(ssb, 2),
#                              '\n',
#                              'Catch(1000mt)=', round(catch, 2))),
#             size = 3)+
#   theme_bw()+
#   labs(x = "Depletion", y = "Yield fraction")+
#   facet_wrap(~LongName)
# y_p_ms
