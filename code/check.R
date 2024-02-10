# check when the scalar on plankton is getting applied - plot full biomass series
# closest to base but with climate would be idx = 2 or 3, probably

res_30 <- readRDS("NOAA_Azure/results/f35/job20240129042511/results/30-result.rds")

biom <- res_30$`30`$biomage_30
catch <- res_30$`30`$catch_30

biom_long <- biom %>%
  pivot_longer(-Time, names_to = 'Code.Age', values_to = 'biomass_mt') %>%
  separate(Code.Age, into = c('Code', 'Age'), sep = '\\.') %>%
  filter(Code %in% c("EUP", "ZM", "PL", "POL","COD","ATF","CAP","SAN","EUL")) %>%
  group_by(Time, Code) %>%
  summarize(biomass_mt = sum(biomass_mt)) %>%
  ungroup() %>%
  mutate(Time = Time / 365)
  

biom_long %>%
  ggplot(aes(x = Time, y = biomass_mt))+
  geom_point()+
  geom_line()+
  facet_wrap(~Code, scales = "free")


# wrong forcings, congratulations