# compare two sets of runs / setups
library(tidyverse)

run1 <- "fspb_beta_new"
run2 <- "beta10"

# make folder for plots
outdir <- paste0("NOAA_Azure/results/figures/comparisons/",run1,"_vs_",run2,"/")
dir.create(outdir, recursive = T)

# depletion and catch curves
dep1 <- read.csv(paste0("NOAA_Azure/results/for_comp/",run1,"_ms_vs_ss.csv"), header = T)
dep2 <- read.csv(paste0("NOAA_Azure/results/for_comp/",run2,"_ms_vs_ss.csv"), header = T)

dep_all <- rbind(dep1, dep2)

# plot
grp1 <- unique(dep_all$LongNamePlot)[1:6]
p1 <- dep_all %>%
  filter(experiment == "Multispecies") %>%
  filter(LongNamePlot %in% grp1) %>%
  ggplot(aes(x = f, y = mt/1000, color = batch))+
  geom_line(linewidth = 0.8)+
  scale_color_manual(values = c("blue","red"))+
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
  geom_line(linewidth = 0.8)+
  scale_color_manual(values = c("blue","red"))+
  geom_vline(data = fmsy %>% filter(LongNamePlot %in% grp2), aes(xintercept = FMSY, group = LongNamePlot), linetype = 'dashed', color = 'orange')+
  geom_vline(data = atlantis_fmsy %>% filter(LongNamePlot %in% grp2), aes(xintercept = atlantis_fmsy, group = LongNamePlot), linetype = 'dashed', color = 'blue')+
  geom_hline(data = b35 %>% filter(LongNamePlot %in% grp2), aes(yintercept = b35/1000, group = LongNamePlot), linetype = 'dashed', color = 'red')+
  theme_bw()+
  scale_y_continuous(limits = c(0, NA))+
  labs(x = 'F (as perceived by the model)', y = '1000\'s of tons')+
  facet_grid2(LongNamePlot~type, scales = 'free', independent = 'all')+
  theme(strip.text.y = element_text(angle=0))

# make a figure
ggsave(paste0(outdir,'COMP_yield_curves_MS_1.png'), p1, width = 6, height = 6)
ggsave(paste0(outdir,'COMP_yield_curves_MS_2.png'), p2, width = 6, height = 6)

# for joint meeting AFSC December 2023
grp3 <- c("Pacific\ncod", "Walleye\npollock", "Pacific\nOcean\nPerch", "Shallow-water\nflatfish")
f_plot3 <- dep_all %>%
  filter(experiment == "Multispecies") %>%
  filter(LongNamePlot %in% grp3) %>%
  ggplot(aes(x = f, y = mt/1000, color = batch))+
  geom_line(linewidth = 1)+
  scale_color_viridis_d(begin = 0.2, end = 0.8)+
  #geom_vline(data = fmsy %>% filter(LongNamePlot %in% grp2), aes(xintercept = FMSY, group = LongNamePlot), linetype = 'dashed', color = 'orange')+
  geom_vline(data = atlantis_fmsy %>% filter(LongNamePlot %in% grp3), aes(xintercept = atlantis_fmsy, group = LongNamePlot), linetype = 'dashed', color = 'blue')+
  geom_hline(data = b35 %>% filter(LongNamePlot %in% grp3), aes(yintercept = b35/1000, group = LongNamePlot), linetype = 'dashed', color = 'red')+
  theme_bw()+
  scale_y_continuous(limits = c(0, NA))+
  labs(x = 'F', y = '1000\'s of tons')+
  facet_grid2(LongNamePlot~type, scales = 'free', independent = 'all')+
  theme(strip.text.y = element_text(angle=0))
ggsave(paste0(outdir,'COMP_yield_curves_MS_3.png'), f_plot3, width = 7, height = 5)

# unselected age classes
un_base <- read.csv(paste0("NOAA_Azure/results/for_comp/",run1,"_unselected.csv"), header = T)
un_comp <- read.csv(paste0("NOAA_Azure/results/for_comp/",run2,"_unselected.csv"), header = T)

un_all <- rbind(un_base, un_comp)

# plot
p3 <- un_all %>%
  filter(experiment == "ms") %>%
  ggplot(aes(x = f, y = biomass_mt/1000, color = batch))+
  geom_line(linewidth = 0.8)+
  scale_color_manual(values = c("blue","red"))+
  geom_vline(data = fmsy, aes(xintercept = FMSY, group = LongName), linetype = 'dashed', color = 'orange')+
  geom_vline(data = atlantis_fmsy, aes(xintercept = atlantis_fmsy, group = LongName), linetype = 'dashed', color = 'blue')+
  theme_bw()+
  scale_y_continuous(limits = c(0, NA))+
  labs(x = 'F as perceived by the model', y = '1000\'s of tons')+
  facet_wrap(~LongName, scales = 'free', ncol = 3)+
  theme(strip.text.y = element_text(angle=0))
p3

ggsave(paste0(outdir,"COMP_unselected.png"), p3, width = 9, height = 6)

# predator and forage fish biomass
top_base <- read.csv(paste0("NOAA_Azure/results/for_comp/",run1,"_top_pred.csv"), header = T)
top_comp <- read.csv(paste0("NOAA_Azure/results/for_comp/",run2,"_top_pred.csv"), header = T)

top_all <- rbind(top_base, top_comp)

forage_base <- read.csv(paste0("NOAA_Azure/results/for_comp/",run1,"_forage.csv"), header = T)
forage_comp <- read.csv(paste0("NOAA_Azure/results/for_comp/",run2,"_forage.csv"), header = T)

forage_all <- rbind(forage_base, forage_comp)

# plot (separate top preds and forage)
p4 <- top_all %>%
  filter(Code %in% c("KWT","KWR","WHT","WHH","WHB","WHG","DOL","SSL","PIN","BDF","BDI","BSF","BSI","SHD","SHP")) %>%
  ggplot(aes(x = mult, y = biomchange, color = batch))+
  geom_line(linewidth = 0.8)+
  scale_color_manual(values = c("blue","red"))+
  geom_hline(yintercept = 1, linetype = "dashed")+
  theme_bw()+
  labs(x = "F35 permutation", y = "Change in biomass from unfished")+
  facet_wrap(~LongName)

p5 <- forage_all %>%
  filter(Code %in% c("CAP","SAN","EUL","HER","FOS")) %>%
  ggplot(aes(x = mult, y = biomchange, color = batch))+
  geom_line(linewidth = 0.8)+
  scale_color_manual(values = c("blue","red"))+
  geom_hline(yintercept = 1, linetype = "dashed")+
  theme_bw()+
  labs(x = "F35 permutation", y = "Change in biomass from unfished")+
  facet_wrap(~LongName)

ggsave(paste0(outdir,"COMP_other_top.png"), p4, width = 8, height = 6)
ggsave(paste0(outdir,"COMP_other_forage.png"), p5, width = 8, height = 3)

# mortality

mort_base <- read.csv(paste0("NOAA_Azure/results/for_comp/",run1,"_rel_m.csv"), header = T)
mort_comp <- read.csv(paste0("NOAA_Azure/results/for_comp/",run2,"_rel_m.csv"), header = T)

mort_all <- rbind(mort_base, mort_comp)

# plot
p6 <- mort_all %>%
  filter(mult > 0) %>%
  ggplot(aes(x = mult, y = rel_m, color = batch))+
  geom_line(linewidth = 0.8)+
  scale_color_manual(values = c("blue","red"))+
  geom_hline(yintercept = 1, linetype = "dashed")+
  theme_bw()+
  labs(x = "F35 permutation", y = "Relative change in M (from mult = 0.2)")+
  facet_wrap(~LongName)
p6

ggsave(paste0(outdir,"COMP_relative_m.png"), p6, width = 8, height = 4)

# SR

sr_base <- read.csv(paste0("NOAA_Azure/results/for_comp/",run1,"_sr.csv"), header = T)
sr_comp <- read.csv(paste0("NOAA_Azure/results/for_comp/",run2,"_sr.csv"), header = T)

sr_all <- rbind(sr_base, sr_comp)

# drop some cols
sr_all <- sr_all %>%
  dplyr::select(Name, depletion,rprop,batch)

# bring in pre-calibration SR curves
precal <- read.csv("NOAA_Azure/results/for_comp/SR_pre_calibration.csv", header = T)
precal <- precal %>%
  filter(Name %in% unique(sr_all$Name)) %>%
  mutate(depletion = SSB / S0,
         rprop = R / BHalpha) %>%
  dplyr::select(Name, depletion, rprop) %>%
  mutate(batch = "pre-calibration")

# bind
sr_all <- rbind(sr_all, precal)

p7 <- sr_all %>%
  #filter(batch != "pre-calibration") %>% # I want to see this one only for the 4F case
  ggplot(aes(x = depletion, y = rprop, color = batch))+
  geom_line(linewidth = 1)+
  scale_color_manual(values = c("blue","red","orange"))+
  geom_hline(yintercept = 0.5, linetype = "dashed")+
  #geom_vline(xintercept = 0.125, color = "red", linetype = "dashed")+
  theme_bw()+
  scale_x_continuous(limits = c(0,1))+
  scale_y_continuous(limits = c(0,1))+
  labs(x = "SSB/S0", y = "R/R0")+
  facet_wrap(~Name)
p7
  
ggsave(paste0(outdir,"COMP_SR.png"), p7, width = 8, height = 4)
