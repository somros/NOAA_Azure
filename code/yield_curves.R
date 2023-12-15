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
this_job <- "job20231207060427"
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
write.csv(atlantis_fmsy, "NOAA_Azure/data/f35_vector_PROXY.csv", row.names = F)

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
  labs(x = 'F as perceived by the model', y = '1000\'s of tons')+
  facet_grid2(LongName~type, scales = 'free', independent = 'all')+
  theme(strip.text.y = element_text(angle=0))
f_plot

# make a figure
t <- format(Sys.time(),'%Y-%m-%d %H-%M-%S')
ggsave(paste0('NOAA_Azure/results/figures/yield_curves',t,'.png'), f_plot, width = 8, height = 16)

# make figures for slides (break into two columns)
# # plot
# grp1 <- unique(f_df$LongName)[1:6]
# f_plot1 <- f_df %>%
#   filter(LongName %in% grp1) %>%
#   ggplot(aes(x = f, y = mt/1000))+
#   geom_line()+
#   geom_point(size = 2)+
#   geom_vline(data = fmsy %>% filter(LongName %in% grp1), aes(xintercept = FMSY, group = LongName), linetype = 'dashed', color = 'orange')+
#   geom_vline(data = atlantis_fmsy %>% filter(LongName %in% grp1), aes(xintercept = atlantis_fmsy, group = LongName), linetype = 'dashed', color = 'blue')+
#   geom_hline(data = b35 %>% filter(LongName %in% grp1), aes(yintercept = b35/1000, group = LongName), linetype = 'dashed', color = 'red')+
#   geom_hline(data = atlantis_fmsy %>% filter(LongName %in% grp1) %>% mutate(type = 'Biomass'), 
#              aes(yintercept = b_fmsy/1000, group = LongName), linetype = 'dashed', color = 'blue')+
#   geom_text(data = annotations %>% filter(LongName %in% grp1) %>% mutate(type = 'Biomass'), 
#             aes(x=Inf,y=Inf,hjust=1,vjust=1,label=paste0('Depletion=',depletion)), color = 'blue')+
#   theme_bw()+
#   scale_y_continuous(limits = c(0, NA))+
#   labs(x = 'F (as perceived by the model)', y = '1000\'s of tons')+
#   facet_grid2(LongName~type, scales = 'free', independent = 'all')+
#   theme(strip.text.y = element_text(angle=0))
# 
# grp2 <- unique(f_df$LongName)[7:12]
# f_plot2 <- f_df %>%
#   filter(LongName %in% grp2) %>%
#   ggplot(aes(x = f, y = mt/1000))+
#   geom_line()+
#   geom_point(size = 2)+
#   geom_vline(data = fmsy %>% filter(LongName %in% grp2), aes(xintercept = FMSY, group = LongName), linetype = 'dashed', color = 'orange')+
#   geom_vline(data = atlantis_fmsy %>% filter(LongName %in% grp2), aes(xintercept = atlantis_fmsy, group = LongName), linetype = 'dashed', color = 'blue')+
#   geom_hline(data = b35 %>% filter(LongName %in% grp2), aes(yintercept = b35/1000, group = LongName), linetype = 'dashed', color = 'red')+
#   geom_hline(data = atlantis_fmsy %>% filter(LongName %in% grp2) %>% mutate(type = 'Biomass'), 
#              aes(yintercept = b_fmsy/1000, group = LongName), linetype = 'dashed', color = 'blue')+
#   geom_text(data = annotations %>% filter(LongName %in% grp2) %>% mutate(type = 'Biomass'), 
#             aes(x=Inf,y=Inf,hjust=1,vjust=1,label=paste0('Depletion=',depletion)), color = 'blue')+
#   theme_bw()+
#   scale_y_continuous(limits = c(0, NA))+
#   labs(x = 'F (as perceived by the model)', y = '1000\'s of tons')+
#   facet_grid2(LongName~type, scales = 'free', independent = 'all')+
#   theme(strip.text.y = element_text(angle=0))
# 
# # make a figure
# t <- format(Sys.time(),'%Y-%m-%d %H-%M-%S')
# ggsave(paste0('NOAA_Azure/results/figures/yield_curves',t,'_1.png'), f_plot1, width = 6, height = 6)
# ggsave(paste0('NOAA_Azure/results/figures/yield_curves',t,'_2.png'), f_plot2, width = 6, height = 6)


