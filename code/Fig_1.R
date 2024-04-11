# Alberto Rovellini
# 4/11/2024
# This cript takes the harvest specification Table from AKFIN answers (1986-2024)
# We can probably drop Fig. S1.2 and 1.3 
library(readxl)
library(tidyverse)
library(viridis)

grps <- read.csv("data/GOA_Groups.csv")

specs <- read_excel("data/GOA_harvest specs_1986-2024.xlsx", 
                    sheet = 1,
                    na = "n/a",
                    n_max = 131)

# lots of cleaning to do
# drop asterisks and commas and turn to numeric
for(col in names(specs)){
  specs[[col]] <- gsub("\\*","", specs[[col]])
  specs[[col]] <- gsub(",","", specs[[col]])
}

# now handle column names
# pad years
colnames(specs) <- c("", "", rep(2024:1986, each = 3))
# collapse column names with the first row for pivot later
new_row <- rep(NA, ncol(specs))
for(i in 1:ncol(specs)){
  new_row[i] <- paste(specs[1,i], names(specs)[i], sep = "_")
  new_row[1:2] <- gsub("_","",new_row[1:2])
}

# set new colnames
colnames(specs) <- new_row
# drop old row 1
specs <- specs[-1,]

# now pad the species column
for(i in 1:nrow(specs)){
  if(is.na(specs[i,1])){
    specs[i,1] <- specs[i-1,1]
  }
}

# pivot longer, split spec and year, add Tier and Atlantis functional group
specs_long <- specs %>%
  pivot_longer(-c(Species, Area), names_to = "Spec_Year", values_to = "mt") %>%
  separate(Spec_Year, into = c("Spec", "Year"), sep = "_") %>%
  filter(Year > 1990, Area == "Total") %>%
  mutate(mt = as.numeric(mt))

species <- sort(unique(specs_long$Species))
key <- data.frame("Species" = species,
                  "Tier" = c(3,5,3,3,3,5,3,4,4,3,3,3,3,3,3,3,3,3,3),
                  "Code" = c("ATF","SKB","FFD","RFP","FHS","SKL","RFS","FFS","RFD","COD","POP","RFP","POL","REX","RFS","SBF","FFS","RFS","RFS"))

specs_long <- specs_long %>%
  left_join(key, by = "Species") %>%
  left_join(grps %>% select(Code, LongName), by = "Code")

# order factors
specs_long$Spec <- factor(specs_long$Spec, levels = c("OFL", "ABC", "TAC"))


# make a bar chart
new_fig1 <- specs_long %>%
  filter(Tier == 3) %>%
  group_by(Year, Spec, LongName) %>%
  summarise(mt = sum(mt, na.rm = T)) %>%
  ggplot(aes(x = Year, y = mt/1000, fill = LongName))+
  geom_bar(stat = "identity", position = "stack")+
  scale_fill_viridis_d()+
  geom_hline(yintercept = 800, linetype = "dashed", color = "red")+
  theme_bw()+
  scale_x_discrete(breaks = seq(1992,2024,2))+
  labs(x = "", y = "1000 mt", fill = "Stock")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  theme(legend.position="bottom",
        legend.spacing.x = unit(0.1, 'cm'))+
  guides(fill = guide_legend(nrow = 4))+
  facet_grid(~Spec)

ggsave("../Paper2OY/text/figures/fig1.png", new_fig1, width = 8, height = 4)

