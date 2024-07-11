library(dplyr)
library(readxl)
library(gridExtra)
library(grid)
library(ggplot2)
library(cowplot)

full.limit <- read.csv("/Users/michelle/Documents/DOSTraits/new_filtered_SA/Sulfate Assimilationprocessed_above75genomecompletion_2.csv")

tmp <- full.limit %>% 
  filter(MaxCompletion != 1 & (Step8==1 | Step9==1) & Step1==0 & Step2==0 & Step3==0 & Step5==0 & Step6==0)%>%
  group_by(Group) %>%
  summarise(n())

print(tmp[11:20,])
