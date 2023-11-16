# Alberto Rovellini
# 11/15/2023
# Based on an idea from Grant Adams

# The single-species experiments have more information that we have not used.
# For each predator, we can plot the biomass of all other groups across the values of F
# F on the other groups is fixed and low (1/4 FMSY)

library(tidyverse)

# read data
# read in Groups.csv file
grp_path <- here('NOAA_Azure/data/GOA_Groups.csv') # functional groups
grps <- read.csv(grp_path)

# read maturity and selectivity information
mat <- read.csv("NOAA_Azure/data/age_at_mat.csv", header = T) # age at maturity
selex <- read.csv("NOAA_Azure/data/age_at_selex_new.csv", header = T) # age at selectivity

# read in lookup tables and f35 proxy values
f_lookup <- read.csv("NOAA_Azure/data/f_lookup.csv")

# list Tier 3 stocks 
t3_fg <- f_lookup %>% pull(species) %>% unique()

# list single-species F output files
filesdir <- "NOAA_Azure/results/pre-processing/job20231103012611/results/"
f_files <- list.files(filesdir)
f_files_df <- data.frame('file' = f_files,
                         idx = as.numeric(gsub("-","",substr(f_files,1,2)))) %>%
  arrange(idx)

# assuming anything that isn't the predator at hand is "prey"
pull_prey_biom <- function(predator){

  these_f_idx <- f_lookup %>% filter(species == predator) %>% pull(idx)
  these_files <- paste0(filesdir, f_files_df %>% filter(idx %in% these_f_idx) %>% pull(file))
  
  # loop over files and pull the biomage table from each
  f_biom_list <- list()
  for (i in 1:length(these_files)){
    
    this_file <- these_files[i]
    biomage <- readRDS(this_file)[[1]][[2]] # read biomage
    
    # get terminal biomass (all age classes combined)
    biomage_end <- biomage %>%
      slice_tail(n = 5) %>% # use last xxx years
      summarise(across(-"Time", ~ mean(.x, na.rm = TRUE))) %>%
      ungroup() %>%
      pivot_longer(everything(), names_to = 'Code.Age', values_to = 'biomass_mt') %>%
      separate_wider_delim(Code.Age, delim = '.', names = c('Code', 'Age')) %>%
      filter(Code %in% t3_fg) %>%
      group_by(Code) %>%
      summarise(biomass_mt = sum(biomass_mt)) %>%
      ungroup() %>%
      mutate(f = f_lookup %>% filter(species == predator, fidx == i) %>% pull(f))
    
    f_biom_list[[i]] <- biomage_end
  }
  
  f_biom <- bind_rows(f_biom_list) %>% 
    mutate(CodePred = predator) %>%
    rename(CodePrey = Code)
}

prey_biom <- lapply(t3_fg, pull_prey_biom) %>% bind_rows()

# join long names
prey_biom <- prey_biom %>%
  left_join(grps %>% select(Code, LongName), by = c("CodePrey"="Code")) %>%
  rename(PreyName = LongName) %>%
  left_join(grps %>% select(Code, LongName), by = c("CodePred"="Code")) %>%
  rename(PredName = LongName) 

prey_biom <- prey_biom %>%
  select(PredName, PreyName, f, biomass_mt) %>%
  left_join(tt %>% 
              filter(f == 0) %>% 
              select(PredName, PreyName, biomass_mt) %>%
              rename(biomass_f0 = biomass_mt), 
            by = c("PredName","PreyName")) %>%
  mutate(biomchange = biomass_mt / biomass_f0) %>%
  rowwise() %>%
  mutate(biomchange = ifelse(PreyName == PredName, NA, biomchange)) %>%
  ungroup()

# plot by predator being fished out
p <- prey_biom %>%
  ggplot(aes(x = f, y = biomchange, color = PreyName))+
  geom_line()+
  geom_point()+
  theme_bw()+
  facet_wrap(~PredName, scales = "free_x")

# overall confirmed fairly small effects. ATF release on POL is the biggest, still "only" x1.5
# fishing some species has no effect on the food web
# most importantly, there are no negative bottom-up effects of depleting prey
# this compounds the dynamics observed with the lack of propagation in the low LTL productivity scenarios
# they seem to either switch prey or not eat enough of each other, or both
# prioritize reparameterizing diets and consumption
ggsave("NOAA_Azure/results/figures/prey_biomass.png", p, width = 10, height = 6)

# plot by prey being affected
p1 <- prey_biom %>%
  ggplot(aes(x = f, y = biomchange, color = PredName))+
  geom_line()+
  geom_point()+
  theme_bw()+
  facet_wrap(~PreyName, scales = "free_x")
p1
