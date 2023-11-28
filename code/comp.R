# compare fixed ATF to base runs
library(tidyverse)

# depletion and catch curves
dep_base <- read.csv("NOAA_Azure/results/for_comp/base_experiment_ms_vs_ss.csv", header = T)
dep_atf <- read.csv("NOAA_Azure/results/for_comp/atf_fixed_ms_vs_ss.csv", header = T)

dep_all <- rbind(dep_base, dep_atf)

# plot
grp1 <- unique(dep_all$LongNamePlot)[1:6]
p1 <- dep_all %>%
  filter(experiment == "Multispecies") %>%
  filter(LongNamePlot %in% grp1) %>%
  ggplot(aes(x = f, y = mt/1000, color = batch))+
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

grp2 <- unique(dep_all$LongNamePlot)[7:12]
p2 <- dep_all %>%
  filter(experiment == "Multispecies") %>%
  filter(LongNamePlot %in% grp2) %>%
  ggplot(aes(x = f, y = mt/1000, color = batch))+
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
ggsave(paste0('NOAA_Azure/results/figures/COMP_yield_curves',t,'_MS_1.png'), p1, width = 6, height = 6)
ggsave(paste0('NOAA_Azure/results/figures/COMP_yield_curves',t,'_MS_2.png'), p2, width = 6, height = 6)

# unselected age classes
un_base <- read.csv("NOAA_Azure/results/for_comp/base_experiment_unselected.csv", header = T)
un_atf <- read.csv("NOAA_Azure/results/for_comp/atf_fixed_unselected.csv", header = T)

un_all <- rbind(un_base, un_atf)

# plot
p3 <- un_all %>%
  filter(experiment == "ms") %>%
  ggplot(aes(x = f, y = biomass_mt/1000, color = batch))+
  geom_line(linewidth = 1)+
  scale_color_viridis_d(begin = 0.2, end = 0.8)+
  geom_vline(data = fmsy, aes(xintercept = FMSY, group = LongName), linetype = 'dashed', color = 'orange')+
  geom_vline(data = atlantis_fmsy, aes(xintercept = atlantis_fmsy, group = LongName), linetype = 'dashed', color = 'blue')+
  theme_bw()+
  scale_y_continuous(limits = c(0, NA))+
  labs(x = 'F as perceived by the model', y = '1000\'s of tons')+
  facet_wrap(~LongName, scales = 'free', ncol = 3)+
  theme(strip.text.y = element_text(angle=0))
p3

ggsave(paste0("NOAA_Azure/results/figures/COMP_unselected.png"), p3, width = 9, height = 6)

# predator and forage fish biomass
top_base <- read.csv("NOAA_Azure/results/for_comp/base_experiment_top_pred.csv", header = T)
top_atf <- read.csv("NOAA_Azure/results/for_comp/atf_fixed_top_pred.csv", header = T)

top_all <- rbind(top_base, top_atf)

forage_base <- read.csv("NOAA_Azure/results/for_comp/base_experiment_forage.csv", header = T)
forage_atf <- read.csv("NOAA_Azure/results/for_comp/atf_fixed_forage.csv", header = T)

forage_all <- rbind(forage_base, forage_atf)

# plot (separate top preds and forage)
p4 <- top_all %>%
  filter(Code %in% c("KWT","KWR","WHT","WHH","WHB","WHG","DOL","SSL","PIN","BDF","BDI","BSF","BSI","SHD","SHP")) %>%
  ggplot(aes(x = mult, y = biomchange, color = batch))+
  geom_line(linewidth = 1)+
  scale_color_viridis_d(begin = 0.2, end = 0.8)+
  geom_hline(yintercept = 1, color = "red", linetype = "dashed")+
  theme_bw()+
  labs(x = "F35 permutation", y = "Change in biomass from unfished")+
  facet_wrap(~LongName)

p5 <- forage_all %>%
  filter(Code %in% c("CAP","SAN","EUL","HER","FOS")) %>%
  ggplot(aes(x = mult, y = biomchange, color = batch))+
  geom_line(linewidth = 1)+
  scale_color_viridis_d(begin = 0.2, end = 0.8)+
  geom_hline(yintercept = 1, color = "red", linetype = "dashed")+
  theme_bw()+
  labs(x = "F35 permutation", y = "Change in biomass from unfished")+
  facet_wrap(~LongName)

ggsave(paste0("NOAA_Azure/results/figures/COMP_other_top.png"), p4, width = 8, height = 6)
ggsave(paste0("NOAA_Azure/results/figures/COMP_other_forage.png"), p5, width = 8, height = 3)

# mortality

mort_base <- read.csv("NOAA_Azure/results/for_comp/base_experiment_rel_m.csv", header = T)
mort_atf <- read.csv("NOAA_Azure/results/for_comp/atf_fixed_rel_m.csv", header = T)

mort_all <- rbind(mort_base, mort_atf)

# plot
p6 <- mort_all %>%
  filter(mult > 0) %>%
  ggplot(aes(x = mult, y = rel_m, color = batch))+
  geom_line(linewidth = 1)+
  scale_color_viridis_d(begin = 0.2, end = 0.8)+
  geom_hline(yintercept = 1, color = "red", linetype = "dashed")+
  theme_bw()+
  labs(x = "F35 permutation", y = "Relative change in M (from mult = 0.2)")+
  facet_wrap(~LongName)
p6

ggsave(paste0("NOAA_Azure/results/figures/COMP_relative_m.png"), p6, width = 8, height = 4)
