setwd("/Users/michelle/Documents/DOSTraits/")
library(tidyr)
library(dplyr)
library(ggplot2)
# make figure one, using: a file with the max completion per genome, per pathway; the genome id to group name file; and the lifestyle

si_lifestyle <- read.csv("new_filtered_SA/Sulfate Assimilationprocessed_above75genomecompletion_wlifestyle.csv")
si_lifestyle <- si_lifestyle[,c(1,23)]
names(si_lifestyle)[1] <- "Genome"
names_ids <- read.csv("genome_summary_ELM_2.csv")
names_ids <- names_ids[,c(1,56)]
names(names_ids)[1] <- "Genome"
maxcompletion <- read.csv("kofam_results/summary_max_table.csv")
names(maxcompletion)[1] <- "Genome"
names_id_filtered <- names_ids[which(names_ids$Genome %in% maxcompletion$Genome),]

# combine our lifestyle and max completion tables
df1 <- merge(maxcompletion, si_lifestyle, by.x = "Genome", by.y = "Genome")
df2 <- gather(df1, key = "key", value = "value", Cysteine.Biosynthesis:S.adenosyl.L.homocysteine.degradation)
df2 <- df2 %>% filter(value ==1)
df3 <- df2 %>% group_by(key) %>% summarise(count = n()) %>% arrange(., count)
order <- rev(df3$key)
ggplot(df2, aes(x = factor(key, order))) + geom_bar() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("Pathway")
ggsave(plot = last_plot(), file = "/Users/michelle/Documents/GitHub/sulfurtraits/ocean_microbiomics_results/kofam_results/summary_barplot_1.pdf", unit = "in", height = 7, width = 11) 
ggplot(df2, aes(x = factor(key, order), fill = Lifestyle)) + geom_bar(position = "fill") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("Pathway")
ggsave(plot = last_plot(), file = "/Users/michelle/Documents/GitHub/sulfurtraits/ocean_microbiomics_results/kofam_results/summary_barplot_2.pdf", unit = "in", height = 7, width = 11) 
df4 <- df2 %>% filter(Lifestyle != "SemiAuxotroph")
df5 <- df4 %>% group_by(key) %>% summarise(count = n()) %>% arrange(., count)
order2 <- rev(df5$key)
ggplot(df4, aes(x = factor(key, order))) + geom_bar() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("Pathway")
ggsave(plot = last_plot(), file = "/Users/michelle/Documents/GitHub/sulfurtraits/ocean_microbiomics_results/kofam_results/summary_barplot_3.pdf", unit = "in", height = 7, width = 11) 
ggplot(df4, aes(x = factor(key, order2), fill = Lifestyle)) + geom_bar(position = "fill") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("Pathway")
ggsave(plot = last_plot(), file = "/Users/michelle/Documents/GitHub/sulfurtraits/ocean_microbiomics_results/kofam_results/summary_barplot_4.pdf", unit = "in", height = 7, width = 11) 

# let's take the max table and change the values a bit and plot as a heatmap just to figure out what we did before
maxcompletion_conservative <- maxcompletion %>% mutate(across(2:43,
                                \(x)if_else(x<1,0,1)))
maxcompletion_conservative1 <- maxcompletion_conservative[,2:43]
row.names(maxcompletion_conservative1) <- maxcompletion$Genome
maxcompletion_generous <- maxcompletion %>% mutate(across(2:43,
                                \(x)if_else(x>0,1,0)))
maxcompletion_generous1 <- maxcompletion_generous[,2:43]
row.names(maxcompletion_generous1) <- maxcompletion$Genome

library(pheatmap)
save_pheatmap_png <- function(x, filename, width=12, height=12) {
  pdf(filename, width = width, height = height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

plot <- pheatmap(maxcompletion_conservative1, cluster_rows = T, show_rownames =F, cluster_cols = T)
save_pheatmap_png(plot, "../../../DOSTraits/2025Results/summary_heatmap_1.pdf")

plot <- pheatmap(maxcompletion_generous1, cluster_rows = T, show_rownames =F, cluster_cols = T)
save_pheatmap_png(plot, "../../../DOSTraits/2025Results/summary_heatmap_2.pdf")

# now, using both generous and conservative, we can work on the rest of the analysis
# first, make the bar plots in a number of different versions
df1 <- merge(maxcompletion_generous, si_lifestyle, by.x = "Genome", by.y = "Genome")
df2 <- gather(df1, key = "key", value = "value", Cysteine.Biosynthesis:S.adenosyl.L.homocysteine.degradation)
df2 <- df2 %>% filter(value ==1)
df3 <- df2 %>% group_by(key) %>% summarise(count = n()) %>% arrange(., count)
order <- rev(df3$key)
ggplot(df2, aes(x = factor(key, order))) + geom_bar() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("Pathway")
ggsave(plot = last_plot(), file = "../../../DOSTraits/2025Results/summary_barplot_1_generous_present.pdf", unit = "in", height = 7, width = 11)
ggplot(df2, aes(x = factor(key, order), fill = Lifestyle)) + geom_bar(position = "fill") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("Pathway")
ggsave(plot = last_plot(), file = "../../../DOSTraits/2025Results/summary_barplot_2_generous_present.pdf", unit = "in", height = 7, width = 11) 
df4 <- df2 %>% filter(Lifestyle != "SemiAuxotroph")
df5 <- df4 %>% group_by(key) %>% summarise(count = n()) %>% arrange(., count)
order2 <- rev(df5$key)
ggplot(df4, aes(x = factor(key, order))) + geom_bar() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("Pathway")
ggsave(plot = last_plot(), file = "../../../DOSTraits/2025Results/summary_barplot_3_generous_present.pdf", unit = "in", height = 7, width = 11) 
ggplot(df4, aes(x = factor(key, order2), fill = Lifestyle)) + geom_bar(position = "fill") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("Pathway")
ggsave(plot = last_plot(), file = "../../../DOSTraits/2025Results/summary_barplot_4_generous_present.pdf", unit = "in", height = 7, width = 11) 

# we can compare the plot above to our original plots, which we'll make again
df1 <- merge(maxcompletion, si_lifestyle, by.x = "Genome", by.y = "Genome")
df2 <- gather(df1, key = "key", value = "value", Cysteine.Biosynthesis:S.adenosyl.L.homocysteine.degradation)
df2 <- df2 %>% filter(value ==1)
df3 <- df2 %>% group_by(key) %>% summarise(count = n()) %>% arrange(., count)
order <- rev(df3$key)
ggplot(df2, aes(x = factor(key, order))) + geom_bar() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("Pathway")
ggsave(plot = last_plot(), file = "../../../DOSTraits/2025Results/summary_barplot_1_present.pdf", unit = "in", height = 7, width = 11) 
ggplot(df2, aes(x = factor(key, order), fill = Lifestyle)) + geom_bar(position = "fill") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("Pathway")
ggsave(plot = last_plot(), file = "../../../DOSTraits/2025Results/summary_barplot_2_present.pdf", unit = "in", height = 7, width = 11) 
df4 <- df2 %>% filter(Lifestyle != "SemiAuxotroph")
df5 <- df4 %>% group_by(key) %>% summarise(count = n()) %>% arrange(., count)
order2 <- rev(df5$key)
ggplot(df4, aes(x = factor(key, order))) + geom_bar() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("Pathway")
ggsave(plot = last_plot(), file = "../../../DOSTraits/2025Results/summary_barplot_3_present.pdf", unit = "in", height = 7, width = 11) 
ggplot(df4, aes(x = factor(key, order2), fill = Lifestyle)) + geom_bar(position = "fill") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("Pathway")
ggsave(plot = last_plot(), file = "../../../DOSTraits/2025Results/summary_barplot_4_present.pdf", unit = "in", height = 7, width = 11) 

df2 <- gather(df1, key = "key", value = "value", Cysteine.Biosynthesis:S.adenosyl.L.homocysteine.degradation)
df2 <- df2 %>% filter(value ==0)
df3 <- df2 %>% group_by(key) %>% summarise(count = n()) %>% arrange(., count)
# order <- rev(df3$key)
ggplot(df2, aes(x = factor(key, order))) + geom_bar() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("Pathway")
ggsave(plot = last_plot(), file = "../../../DOSTraits/2025Results/summary_barplot_5_absent.pdf", unit = "in", height = 7, width = 11) 
ggplot(df2, aes(x = factor(key, order), fill = Lifestyle)) + geom_bar(position = "fill") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("Pathway")
ggsave(plot = last_plot(), file = "../../../DOSTraits/2025Results/summary_barplot_6_absent.pdf", unit = "in", height = 7, width = 11) 
df4 <- df2 %>% filter(Lifestyle != "SemiAuxotroph")
df5 <- df4 %>% group_by(key) %>% summarise(count = n()) %>% arrange(., count)
# order2 <- rev(df5$key)
ggplot(df4, aes(x = factor(key, order))) + geom_bar() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("Pathway")
ggsave(plot = last_plot(), file = "../../../DOSTraits/2025Results/summary_barplot_7_absent.pdf", unit = "in", height = 7, width = 11) 
ggplot(df4, aes(x = factor(key, order2), fill = Lifestyle)) + geom_bar(position = "fill") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("Pathway")
ggsave(plot = last_plot(), file = "../../../DOSTraits/2025Results/summary_barplot_8_absent.pdf", unit = "in", height = 7, width = 11) 

# lastly, plot the conservative one, which is the one I feel least comfortable using
df1 <- merge(maxcompletion_conservative, si_lifestyle, by.x = "Genome", by.y = "Genome")
df2 <- gather(df1, key = "key", value = "value", Cysteine.Biosynthesis:S.adenosyl.L.homocysteine.degradation)
df2 <- df2 %>% filter(value ==0)
df3 <- df2 %>% group_by(key) %>% summarise(count = n()) %>% arrange(., count)
# order <- rev(df3$key)
ggplot(df2, aes(x = factor(key, order))) + geom_bar() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("Pathway")
ggsave(plot = last_plot(), file = "../../../DOSTraits/2025Results/summary_barplot_1_conservative_absent.pdf", unit = "in", height = 7, width = 11) 
ggplot(df2, aes(x = factor(key, order), fill = Lifestyle)) + geom_bar(position = "fill") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("Pathway")
ggsave(plot = last_plot(), file = "../../../DOSTraits/2025Results/summary_barplot_2_conservative_absent.pdf", unit = "in", height = 7, width = 11) 
df4 <- df2 %>% filter(Lifestyle != "SemiAuxotroph")
df5 <- df4 %>% group_by(key) %>% summarise(count = n()) %>% arrange(., count)
# order2 <- rev(df5$key)
ggplot(df4, aes(x = factor(key, order))) + geom_bar() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("Pathway")
ggsave(plot = last_plot(), file = "../../../DOSTraits/2025Results/summary_barplot_3_conservative_absent.pdf", unit = "in", height = 7, width = 11) 
ggplot(df4, aes(x = factor(key, order2), fill = Lifestyle)) + geom_bar(position = "fill") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("Pathway")
ggsave(plot = last_plot(), file = "../../../DOSTraits/2025Results/summary_barplot_4_conservative_absent.pdf", unit = "in", height = 7, width = 11) 



# we need to visualize things relative to each kind of lifestyle, as a total
library(RColorBrewer)
df1 <- merge(maxcompletion_generous, si_lifestyle, by.x = "Genome", by.y = "Genome")
df2 <- gather(df1, key = "key", value = "value", Cysteine.Biosynthesis:S.adenosyl.L.homocysteine.degradation)
df2 <- df2 %>% mutate(across(4, \(x)if_else(x==0,"Absent","Present")))
df2$value <- factor(df2$value)
df2$Lifestyle <- factor(df2$Lifestyle)
df2$key <- factor(df2$key)
df3 <- df2 %>% filter(value == "Present") %>% group_by(key) %>% summarise(count = n()) %>% arrange(., count)
order <- rev(df3$key)
ggplot(df2, aes(x = factor(key, order), fill = value)) + geom_bar( position = "fill") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + xlab("Pathway") + 
  facet_grid(Lifestyle ~ .)+
  theme(plot.margin = unit(c(0,0,0,2), "cm"))
ggsave(plot = last_plot(), file = "../../../DOSTraits/2025Results/summary_barplot_1_relative_generous.pdf", unit = "in", height = 7, width = 11) 
ggplot(df2, aes(x = factor(key, order), fill = Lifestyle)) + geom_bar(position = "fill") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + xlab("Pathway") + 
  facet_grid(value ~ .)+
  theme(plot.margin = unit(c(0,0,0,2), "cm"))
ggsave(plot = last_plot(), file = "../../../DOSTraits/2025Results/summary_barplot_2_relative_generous.pdf", unit = "in", height = 7, width = 11) 
df2$value <- factor(df2$value)
df2$Lifestyle <- factor(df2$Lifestyle)
df2$key <- as.character(df2$key)
ggplot(df2, aes(x = interaction(Lifestyle,factor(key, order)), fill = value)) + geom_bar(aes(fill = interaction(value, Lifestyle)),position = "fill") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + xlab("Pathway") +
  scale_fill_manual(values = brewer.pal(6, "Paired"))+
  theme(plot.margin = unit(c(0,0,0,2), "cm"))
ggsave(plot = last_plot(), file = "../../../DOSTraits/2025Results/summary_barplot_3_relative_generous.pdf", unit = "in", height = 7, width = 11) 
ggplot(df2, aes(x = interaction(value,factor(key, order)), fill = Lifestyle)) + geom_bar(aes(fill = interaction(value, Lifestyle)),position = "fill") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + xlab("Pathway") +
  scale_fill_manual(values = brewer.pal(6, "Paired"))+
  theme(plot.margin = unit(c(0,0,0,2), "cm"))
ggsave(plot = last_plot(), file = "../../../DOSTraits/2025Results/summary_barplot_4_relative_generous.pdf", unit = "in", height = 7, width = 11) 

# use the conservative one
df1 <- merge(maxcompletion_conservative, si_lifestyle, by.x = "Genome", by.y = "Genome")
df2 <- gather(df1, key = "key", value = "value", Cysteine.Biosynthesis:S.adenosyl.L.homocysteine.degradation)
df2 <- df2 %>% mutate(across(4, \(x)if_else(x==0,"Absent","Present")))
df2$value <- factor(df2$value)
df2$Lifestyle <- factor(df2$Lifestyle)
df2$key <- factor(df2$key)
df3 <- df2 %>% filter(value == "Present") %>% group_by(key) %>% summarise(count = n()) %>% arrange(., count)
order <- rev(df3$key)
ggplot(df2, aes(x = factor(key, order), fill = value)) + geom_bar( position = "fill") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + xlab("Pathway") + 
  facet_grid(Lifestyle ~ .)+
  theme(plot.margin = unit(c(0,0,0,2), "cm"))
ggsave(plot = last_plot(), file = "../../../DOSTraits/2025Results/summary_barplot_1_relative_conservative.pdf", unit = "in", height = 7, width = 11) 
ggplot(df2, aes(x = factor(key, order), fill = Lifestyle)) + geom_bar(position = "fill") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + xlab("Pathway") + 
  facet_grid(value ~ .)+
  theme(plot.margin = unit(c(0,0,0,2), "cm"))
ggsave(plot = last_plot(), file = "../../../DOSTraits/2025Results/summary_barplot_2_relative_conservative.pdf", unit = "in", height = 7, width = 11) 
df2$value <- factor(df2$value)
df2$Lifestyle <- factor(df2$Lifestyle)
df2$key <- as.character(df2$key)
ggplot(df2, aes(x = interaction(Lifestyle,factor(key, order)), fill = value)) + geom_bar(aes(fill = interaction(value, Lifestyle)),position = "fill") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + xlab("Pathway") +
  scale_fill_manual(values = brewer.pal(6, "Paired"))+
  theme(plot.margin = unit(c(0,0,0,2), "cm"))
ggsave(plot = last_plot(), file = "../../../DOSTraits/2025Results/summary_barplot_3_relative_conservative.pdf", unit = "in", height = 7, width = 11) 
ggplot(df2, aes(x = interaction(value,factor(key, order)), fill = Lifestyle)) + geom_bar(aes(fill = interaction(value, Lifestyle)),position = "fill") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + xlab("Pathway") +
  scale_fill_manual(values = brewer.pal(6, "Paired"))+
  theme(plot.margin = unit(c(0,0,0,2), "cm"))
ggsave(plot = last_plot(), file = "../../../DOSTraits/2025Results/summary_barplot_4_relative_conservative.pdf", unit = "in", height = 7, width = 11) 


# use the regular one, which will just consider the zeros and the ones and nothing in between
df1 <- merge(maxcompletion, si_lifestyle, by.x = "Genome", by.y = "Genome")
df2 <- gather(df1, key = "key", value = "value", Cysteine.Biosynthesis:S.adenosyl.L.homocysteine.degradation)
df2 <- df2 %>% filter(value==0 | value==1)
df2 <- df2 %>% mutate(across(4, \(x)if_else(x==0,"Absent","Present")))
df2$value <- factor(df2$value)
df2$Lifestyle <- factor(df2$Lifestyle)
df2$key <- factor(df2$key)
df3 <- df2 %>% filter(value == "Present") %>% group_by(key) %>% summarise(count = n()) %>% arrange(., count)
order <- rev(df3$key)
ggplot(df2, aes(x = factor(key, order), fill = value)) + geom_bar( position = "fill") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + xlab("Pathway") + 
  facet_grid(Lifestyle ~ .)+
  theme(plot.margin = unit(c(0,0,0,2), "cm"))
ggsave(plot = last_plot(), file = "../../../DOSTraits/2025Results/summary_barplot_1_relative.pdf", unit = "in", height = 7, width = 11) 
ggplot(df2, aes(x = factor(key, order), fill = Lifestyle)) + geom_bar(position = "fill") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + xlab("Pathway") + 
  facet_grid(value ~ .) +
  theme(plot.margin = unit(c(0,0,0,2), "cm"))
ggsave(plot = last_plot(), file = "../../../DOSTraits/2025Results/summary_barplot_2_relative.pdf", unit = "in", height = 7, width = 11) 
df2$value <- factor(df2$value)
df2$Lifestyle <- factor(df2$Lifestyle)
df2$key <- as.character(df2$key)
ggplot(df2, aes(x = interaction(Lifestyle,factor(key, order)), fill = value)) + geom_bar(aes(fill = interaction(value, Lifestyle)),position = "fill") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + xlab("Pathway") +
  scale_fill_manual(values = brewer.pal(6, "Paired")) +
  theme(plot.margin = unit(c(0,0,0,2), "cm"))
ggsave(plot = last_plot(), file = "../../../DOSTraits/2025Results/summary_barplot_3_relative.pdf", unit = "in", height = 7, width = 11) 
ggplot(df2, aes(x = interaction(value,factor(key, order)), fill = Lifestyle)) + geom_bar(aes(fill = interaction(value, Lifestyle)),position = "fill") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + xlab("Pathway") +
  scale_fill_manual(values = brewer.pal(6, "Paired")) +
  theme(plot.margin = unit(c(0,0,0,2), "cm"))
ggsave(plot = last_plot(), file = "../../../DOSTraits/2025Results/summary_barplot_4_relative.pdf", unit = "in", height = 7, width = 11) 


# make next plot that incorporates lifestyle and taxa and mean sulfur protein content
s_content <- read.csv("avg_estimate_S_content.csv")
df6 <- merge(df1, names_ids, by.x = "Genome", by.y = "Genome")
df7 <- df6 %>% group_by(taxa_elm) %>% summarise(n = n())
df8 <- df6 %>% filter(Lifestyle == "Auxotroph") %>% group_by(taxa_elm)%>% summarise(n_aux = n())
df9 <- df6 %>% filter(Lifestyle == "Assimilator") %>% group_by(taxa_elm)%>% summarise(n_assim = n())
df7 <- merge(df7, df8, all.x = T)
df7 <- merge(df7, df9, all.x = T)
df7$n_aux_frac <- df7$n_aux/df7$n
df7$n_assim_frac <- df7$n_assim/df7$n
df7[is.na(df7)] <- 0
df7 <- merge(df7, s_content, by.x = "taxa_elm", by.y = "Genome")
# we can plot the assimilators vs s content
ggplot(df7, aes(y=n_aux_frac, x = s_content)) +geom_point(aes(color = taxa_elm)) + geom_smooth(method='lm', colour = "black")
ggsave(plot = last_plot(), file = "../../../DOSTraits/2025Results/summary_scontent_auxfrac_withlm.pdf", unit = "in", height = 7, width = 10) 
# can also plot the auxotrophs vs s content
ggplot(df7, aes(y=n_assim_frac, x = s_content)) +geom_point(aes(color = taxa_elm)) + geom_smooth(method='lm', colour = "black")
ggsave(plot = last_plot(), file = "../../../DOSTraits/2025Results/summary_scontent_assimfrac_withlm.pdf", unit = "in", height = 7, width = 10) 
ggplot(df7, aes(y=n_aux_frac, x = s_content)) +geom_point(aes(color = taxa_elm)) 
ggsave(plot = last_plot(), file = "../../../DOSTraits/2025Results/summary_scontent_auxfrac.pdf", unit = "in", height = 7, width = 10) 
# can also plot the auxotrophs vs s content
ggplot(df7, aes(y=n_assim_frac, x = s_content)) +geom_point(aes(color = taxa_elm)) 
ggsave(plot = last_plot(), file = "../../../DOSTraits/2025Results/summary_scontent_assimfrac.pdf", unit = "in", height = 7, width = 10) 
#gams
mod_gamlm <- gam(s_content ~ s(n_aux_frac), data = df7)
summary(mod_gamlm)
plot(mod_gamlm)
mod_gam1 <- gam(s_content ~ s(n_aux_frac, bs = "cr"), data = df7)
summary(mod_gam1)
plot(mod_gam1)
mod_gam2 <- gam(s_content ~ s(n_aux_frac), data = df7, method = "REML")
summary(mod_gam2)
plot(mod_gam2)
mod_gam3 <- gam(s_content ~ s(n_aux_frac, bs = "cr"), data = df7, method = "REML")
summary(mod_gam3)
plot(mod_gam3)
# use n aux frac as a function of s content, which I think is opposite what we want
mod_gamlm <- gam(n_aux_frac ~ s(s_content), data = df7)
summary(mod_gamlm)
plot(mod_gamlm)
mod_gam1 <- gam(n_aux_frac ~ s(s_content, bs = "cr"), data = df7)
summary(mod_gam1)
plot(mod_gam1)
mod_gam2 <- gam(n_aux_frac ~ s(s_content), data = df7, method = "REML")
summary(mod_gam2)
plot(mod_gam2)
mod_gam3 <- gam(n_aux_frac ~ s(s_content, bs = "cr"), data = df7, method = "REML")
summary(mod_gam3)
plot(mod_gam3)
# see if exponential
model <- lm(log(s_content) ~ n_aux_frac, data = df7)
summary(model)
plot(model)
xmodel <- lm(df7$n_aux_frac~df7$s_content)
xmodel
xmodel2 <- lm(n_aux_frac~poly(s_content, 2), data = df7)
summary(xmodel2)
# let's determine the spearman correlation
cor.test(df7$n_aux_frac, df7$s_content, method = "spearman")
# correlation: -0.404756
# p-value: 0.02941
prediction_interval <- exp(predict(
  model, 
  newdata = df7,
  interval="prediction",
  level = 0.95
))
plot(df7$n_aux_frac, df7$s_content, main="Exponential Regression", xlab="S Content", 
     ylab="Auxotroph Fraction", pch=19)
lines(df$n_aux_frac, prediction_interval[,1], col="red", lty=2)
lines(df$n_aux_frac, prediction_interval[,2], col="blue", lty=2)
lines(df$n_aux_frac, prediction_interval[,3], col="blue", lty=2)
legend("topright", legend="Exponential Regression", col="red", lwd=2)

#correlate max value with lifestyle and display
npathways <- ncol(maxcompletion_conservative1)
pathways <- names(maxcompletion_conservative1)
sapathways_coefficients <- as.data.frame(matrix(nrow = 3, ncol = npathways+1))
sapathways_pvalues <- as.data.frame(matrix(nrow = 3, ncol = npathways+1))
names(sapathways_coefficients) <- c("sa_lifestyle",pathways)
names(sapathways_pvalues) <- c("sa_lifestyle",pathways)
sapathways_coefficients[,1] <- c("Default","Conservative","Generous")
sapathways_pvalues[,1] <- c("Default","Conservative","Generous")
sa_lifestyles <- read.csv("../../../DOSTraits/new_filtered_SA/Sulfate Assimilationprocessed_above75genomecompletion_wlifestyle.csv", row.names = 1)[,c(18,22)]

for (name in pathways){
  sub1 <- maxcompletion_conservative[,which(names(maxcompletion_conservative)== name)]
  sub2 <- maxcompletion_generous[,which(names(maxcompletion_generous)== name)]
  sub3 <- maxcompletion[,which(names(maxcompletion)== name)]
  sub <- as.data.frame(cbind(sub3, sub1, sub2))
  names(sub) <- c("Default","Conservative","Generous")
  row.names(sub) <- maxcompletion$Genome
  lifestyles <- names(sub)
  for (lifestyle in lifestyles){
    sub4 <- as.data.frame(sub[,which(names(sub)== lifestyle)])
    row.names(sub4) <- row.names(sub)
    names(sub4) <- "sub4"
    temp <- merge(sub4,sa_lifestyles, by = "row.names")
    cor <- cor.test(temp$sub4, temp$MaxCompletion, method = 'spearman')
    sapathways_pvalues[which(lifestyles == lifestyle ),which(pathways ==  name)+1] <- cor[[3]]
    sapathways_coefficients[which(lifestyles == lifestyle ),which(pathways == name)+1] <- cor[[4]]
  }
}
sapathways_coefficients_long <- gather(sapathways_coefficients, value = "value", key = "pathway", `Cysteine.Biosynthesis`:`S.adenosyl.L.homocysteine.degradation`)
sapathways_coefficients_long$p_value <- NA
sapathways_coefficients_long$bold_value <- NA
for (row in 1:nrow(sapathways_coefficients_long)){
  a <- which(sapathways_pvalues[,1] == sapathways_coefficients_long[row,1])
  b <- which(names(sapathways_pvalues) == sapathways_coefficients_long[row,2])
  sapathways_coefficients_long[row,4] <- sapathways_pvalues[a,b]
  if (is.na(sapathways_coefficients_long$p_value[row] <= 0.05)){
    sapathways_coefficients_long$bold_value[row] <- NA
  } else if (sapathways_coefficients_long$p_value[row] <= 0.05){
    sapathways_coefficients_long$bold_value[row] <- 2
  } else {sapathways_coefficients_long$bold_value[row] <- 1}
}

sapathways_coefficients_long$value <- round(sapathways_coefficients_long$value, digits = 2)
p <- ggplot(data=sapathways_coefficients_long, aes(x=sa_lifestyle, y=pathway, fill=value))
p <- p + geom_tile() + geom_text(aes(label=value, fontface=bold_value), color="black", size = 3,hjust=0.5) + 
  scale_fill_gradientn(colors = c("blue", "white", "red"),limits=c(-1,1)) + theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(filename = "../../../DOSTraits/2025Results/heatmap_sa_lifestyle_correlations.pdf", plot = p, width=8, height=10)

#lastly, let's look for a cluster method
# for now, we'll use a tree that the group produced
library(ggtree)
tree1 <- read.tree("gtdb_newick_tree")
tree2 <- read.tree("gtdb_newick_tree_ar")
p1 <- ggplot(tree1, aes(x, y)) + geom_tree() + theme_tree()
# make original taxa list
taxa1 <- get_taxa_name(p1)
#make new taxa list
new_taxa1 <- taxa1
# get the list of genomes we will keep
genomes_in_p1 <- match(taxa1, names_id_filtered$Genome)
#replace all not in our genome list with a blank
new_taxa1[which(is.na(genomes_in_p1))] <- ""
#replace those in our list with the taxa group name
new_taxa1[which(!is.na(genomes_in_p1))] <- names_id_filtered$taxa_elm[genomes_in_p1[which(!is.na(genomes_in_p1))]]
new_taxa1
ggtree(tree1) + geom_tiplab(aes(label = new_taxa1))
d <- data.frame(label = taxa1, label2 = new_taxa1)
tree1_2 <- full_join(tree1, d, by = "label")
p1_2 <- ggplot(tree1_2, aes(x, y)) + geom_tree() + theme_tree()

p2 <- ggplot(tree2, aes(x, y)) + geom_tree() + theme_tree()
taxa2 <- get_taxa_name(p2)

genomes_in_p2 <- which(taxa2 %in% df1$Genome)
length(genomes_in_p1) + length(genomes_in_p2)

install.packages("seqinr")
library("seqinr")

# load the 16s sequence file
fasta_16s_file <- seqinr::read.fasta(file = "../../../DOSTraits/go_microbiomics-integrated-cpl50_ctn10-16S.fasta", seqtype = "DNA",as.string = TRUE, set.attributes = FALSE)
names(fasta_16s_file) <- gsub(names(fasta_16s_file), pattern = "gene_.*", replacement = "")
# load the scaffold membership file
table <- read.table("../../../DOSTraits/genomes-scaffolds-membership.tsv", sep = "\t")
table$V1 <- paste0(table$V1, "-")

# correlate the sequences with the correct genomes
head(names(fasta_16s_file))
head(match(names(fasta_16s_file),table$V1))
names(fasta_16s_file) <- table$V2[match(names(fasta_16s_file),table$V1)]

# eliminate duplicates?
length(which(!duplicated(names(fasta_16s_file))==TRUE))
fasta_16s_file <- fasta_16s_file[!duplicated(names(fasta_16s_file))]

#resave the file with just our genomes in it
our_genomes <- read.csv("../../../DOSTraits/new_filtered_SA/2024_June19_ERIN_SulfateAssimilationprocessed_above75genomecompletion.csv")
our_genomes <- our_genomes$Genome
all(names(fasta_16s_file) %in% our_genomes)
#FALSE
fasta_16s_file <- fasta_16s_file[which(names(fasta_16s_file) %in% our_genomes)]
write.fasta(sequences = fasta_16s_file, names = names(fasta_16s_file), nbchar = 80, file.out = "../../../DOSTraits/go_microbiomics-16S.fasta")


# let's remake figure 2
library("NatParksPalettes")
library(nationalparkcolors)
library(ggpubr)
library(patchwork)
values5 <- c(natparks.pals("Volcanoes", 5))
"#082544" "#1E547D" "#79668C" "#DE3C37" "#F2DC7E"
values3 <- park_palette("Badlands", 3)

c(natparks.pals("KingsCanyon", 3))
values3 <- c("#44637D", "#613921",  "#AAC9ED") 
c(natparks.pals("IguazuFalls", 6))
"#415521" "#97AD3D" "#4C3425" "#7F6552" "#5A8093" "#9FBAD3"

"#4C3425" "#415521" "#5A8093"

"#7F6552" "#97AD3D" "#9FBAD3"
values6 <- c(natparks.pals("Chamonix", 6))
"#008FF8" "#B6AA0D" "#E2C2A2" "#E23B0E" "#F2C621" "#196689"
values3 <- c("#196689","#E23B0E", "#F2C621" )
#dull:
values3 <- c("#51659C", "#B97069", "#857D83")
values3 <- c("#2F397A","#415521","#5A8093")
"#2F397A" "#415521" "#5A8093" "#9399bd" "#9fb6c2" "#c0d998"
values6 <- c(natparks.pals("Torres", 6))
"#2F397A" "#7391BD" "#894846" "#E9988C" "#535260" "#B7A7A6"
our_genomes <- read.csv("new_filtered_SA/2024_June19_ERIN_SulfateAssimilationprocessed_above75genomecompletion.csv")
table1 <- our_genomes %>% group_by(taxa_elm) %>% summarise(n = n())
table2 <- our_genomes %>% filter(Lifestyle == "auxotroph") %>% group_by(taxa_elm) %>% summarise(n_aux = n())
table1 <- merge(table1, table2,all = T)
table1$n_aux[is.na(table1$n_aux)] <- 0
table1$frac <- table1$n_aux/table1$n
table3 <- our_genomes %>% group_by(taxa_elm) %>% summarise(mean = mean(Mean.Completeness))
table1 <- merge(table1, table3,all = T)
order1 <- table1$taxa_elm[order(table1$frac)]
our_genomes$taxa_elm <- factor(our_genomes$taxa_elm, levels = order1)
our_genomes$Lifestyle <- factor(
  our_genomes$Lifestyle, levels =c("assimilator", "incomplete", "auxotroph"))
p1 <- ggplot(our_genomes, aes(x = taxa_elm )) + geom_bar(aes(fill = Lifestyle),position = "fill") + coord_flip() + theme(legend.position = "top") +
  scale_fill_manual(values = values3, breaks = c("assimilator", "auxotroph", "incomplete"), label = c("Assimilator","Auxotroph", "Incomplete")) + 
  theme_light() + ylab("Proportion of Genomes")+
  xlab("Taxonomic Group")
table1$taxa_elm <- factor(table1$taxa_elm, order1)
p2 <- ggplot(table1, aes(x = taxa_elm, y = n, fill = mean)) + geom_bar(stat = "identity") +
  coord_flip() + theme_light() + ylab("Number of Genomes per Taxonomic Group") +  xlab("") +
  scale_fill_gradient(low = "lightgrey", high = "black", name = "Mean Genome Completion" )
p3 <- ggarrange(p1, p2, align = "hv", legend = "top")
ggsave(filename = "2025Results/fig2.pdf", plot = p3, width=10, height=5)

p2 <- ggplot(table1, aes(x = taxa_elm, y = n, fill = mean)) + geom_bar(stat = "identity") +
  coord_flip() + theme_light() + ylab("Number of Genomes per Taxonomic Group") +  xlab("") +
  scale_fill_viridis_c(option = "mako", direction = -1 )
p3 <- ggarrange(p1, p2, align = "hv", legend = "top")
ggsave(filename = "2025Results/fig2_2.pdf", plot = p3, width=10, height=5)

p2 <- ggplot(table1, aes(x = taxa_elm, y = n, fill = mean)) + geom_bar(stat = "identity") +
  coord_flip() + theme_light() + ylab("Number of Genomes per Taxonomic Group") +  xlab("") +
  scale_fill_viridis_c( direction = -1 )
p3 <- ggarrange(p1, p2, align = "hv", legend = "top")
ggsave(filename = "2025Results/fig2_3.pdf", plot = p3, width=10, height=5)

p2 <- ggplot(table1, aes(x = taxa_elm, y = n, fill = mean)) + geom_bar(stat = "identity") +
  coord_flip() + theme_light() + ylab("Number of Genomes per Taxonomic Group") +  xlab("") +
  scale_fill_viridis_c(option = "plasma", direction = -1 )
p3 <- ggarrange(p1, p2, align = "hv", legend = "top")
ggsave(filename = "2025Results/fig2_4.pdf", plot = p3, width=10, height=5)

p2 <- ggplot(table1, aes(x = taxa_elm, y = n, fill = mean)) + geom_bar(stat = "identity") +
  coord_flip() + theme_light() + ylab("Number of Genomes per Taxonomic Group") +  xlab("") +
  scale_fill_viridis_c(option = "mako", direction = -1 )+ scale_y_log10(minor_breaks = NULL)
p3 <- ggarrange(p1, p2, align = "hv", legend = "top")
ggsave(filename = "2025Results/fig2_5.pdf", plot = p3, width=10, height=5)

p2 <- ggplot(table1, aes(x = taxa_elm, y = n, fill = mean)) + geom_bar(stat = "identity") +
  coord_flip() + theme_light() + ylab("Number of Genomes per Taxonomic Group") +  xlab("") +
  scale_fill_viridis_c( direction = -1 ) + scale_y_log10(minor_breaks = NULL)
p3 <- ggarrange(p1, p2, align = "hv", legend = "top")
ggsave(filename = "2025Results/fig2_6.pdf", plot = p3, width=10, height=5)

p2 <- ggplot(table1, aes(x = taxa_elm, y = n, fill = mean)) + geom_bar(stat = "identity") +
  coord_flip() + theme_light() + ylab("Number of Genomes per Taxonomic Group") +  xlab("") +
  scale_fill_viridis_c(option = "rocket", direction = -1 ) + scale_y_log10(minor_breaks = NULL)
p3 <- ggarrange(p1, p2, align = "hv", legend = "top")
ggsave(filename = "2025Results/fig2_7.pdf", plot = p3, width=10, height=5)

p2 <- ggplot(our_genomes, aes(x = taxa_elm, y = Mean.Completeness)) + geom_boxplot( outlier.size=0) +
  geom_jitter(position=position_jitter(0.2), size = 0.1)+
  coord_flip() + theme_light() + ylab("Genome Completion") +  xlab("")
p3 <- ggarrange(p1, p2, align = "hv", legend = "top")
ggsave(filename = "2025Results/fig2_8.pdf", plot = p3, width=10, height=5)



#let's play with figure 4, using the colors in values6
df1 <- merge(maxcompletion, si_lifestyle, by.x = "Genome", by.y = "Genome")
#remove those pathays we aren't interested in
names(df1)[which(names(df1) == "DOS.Misc.Degradation.2")] <- "Sulfite.Export"
names(df1)[which(names(df1) == "DOS.Misc.Degradation.3")] <- "Sulfite.Dehydrogenase"
df1[,which(names(df1) == "DOS.Misc.Degradation.4")] <- NULL
df1[,which(names(df1) == "Methionine.Salvage.Part.1")] <- NULL
df1[,which(names(df1) == "Methionine.Salvage.Part.2")] <- NULL
df1[,which(names(df1) == "Methionine.Salvage.Fusion.Enzyme")] <- NULL
df1[,which(names(df1) == "Methionine.Salvage.Anaerobic.Isoprenoid.Shunt")] <- NULL
df1[,which(names(df1) == "Methionine.Salvage.Anaerobic.Degradation.of.MTRu1P")] <- NULL

df2 <- gather(df1, key = "key", value = "value", Cysteine.Biosynthesis:S.adenosyl.L.homocysteine.degradation)
nrow(df2)/length(unique(df2$key))
#none removed here
df2 <- df2 %>% filter(value==0 | value==1)
((36*16059)-nrow(df2))
#removing 35,204
((36*16059)-nrow(df2))/(36*16059)
# which is 6% of our dataset
df2 <- df2 %>% mutate(across(4, \(x)if_else(x==0,"Absent","Present")))
df2$value <- factor(df2$value)
df2$Lifestyle <- factor(df2$Lifestyle)
df2$key <- factor(df2$key)
df3 <- df2 %>% filter(value == "Present") %>% group_by(key) %>% summarise(count = n()) %>% arrange(., count)
order <- rev(df3$key)
p4 <- ggplot(df2, aes(x = factor(key, order), fill = value)) + geom_bar( position = "fill") + 
  theme_linedraw() +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + xlab("Pathway") + 
  facet_grid(Lifestyle ~ .)+
  scale_fill_manual(values = values6[c(2,1)]) +
  theme(plot.margin = unit(c(0,0,0,2), "cm"))
p4
ggsave(plot = last_plot(), file = "../../../DOSTraits/2025Results/fig4_a.pdf", unit = "in", height = 20, width = 20) 
p4 <- ggplot(df2, aes(x = factor(key, order), fill = value)) + geom_bar( position = "fill") + 
  theme_linedraw() +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + xlab("Pathway") + 
  facet_grid(Lifestyle ~ .) +
  scale_fill_manual(values = values6[c(4,3)]) +
  theme(plot.margin = unit(c(0,0,0,2), "cm"))
p4
ggsave(plot = last_plot(), file = "../../../DOSTraits/2025Results/fig4_b.pdf", unit = "in", height = 20, width = 20) 
p4 <- ggplot(df2, aes(x = factor(key, order), fill = value)) + geom_bar( position = "fill") + 
  theme_linedraw()  + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + xlab("Pathway")+ 
  facet_grid(Lifestyle ~ .)+
  scale_fill_manual(values = values6[c(6,5)]) +
  theme(plot.margin = unit(c(0,0,0,2), "cm"))
p4
ggsave(plot = last_plot(), file = "../../../DOSTraits/2025Results/fig4_c.pdf", unit = "in", height = 20, width = 20) 


#ggsave(plot = last_plot(), file = "../../../DOSTraits/2025Results/summary_barplot_1_relative.pdf", unit = "in", height = 7, width = 11) 
#df1 <- merge(maxcompletion, si_lifestyle, by.x = "Genome", by.y = "Genome")
df2 <- gather(df1, key = "key", value = "value", Cysteine.Biosynthesis:S.adenosyl.L.homocysteine.degradation)
df2 <- df2 %>% filter(value ==1)
df3 <- df2 %>% group_by(key) %>% summarise(count = n()) %>% arrange(., count)
core_pathways <- df3$key[which(df3$count >= .75*16059)]
order <- rev(df3$key)
df2$core_value <- NA
df2$core_value[which(df2$key %in% core_pathways)] <- "Core"
df2$core_value[which(!(df2$key %in% core_pathways))] <- "Accessory"

p5 <- ggplot(df2, aes(x = factor(key, order))) + 
  geom_bar(fill = "black") + theme_linedraw() + theme(axis.text.x = element_text( angle = 45, vjust = 1, hjust=1)) + 
  xlab("Pathway")  
p5
ggsave(plot = last_plot(), file = "../../../DOSTraits/2025Results/fig4_d.pdf", unit = "in", height = 20, width = 20) 

df2 <- df1 %>% group_by(Lifestyle) %>% summarise(count = n())
p6 <- ggplot(df1, aes(x = Lifestyle)) + coord_flip() +
  geom_bar(fill = "black") + theme_linedraw() 
p6 
ggsave(plot = last_plot(), file = "../../../DOSTraits/2025Results/fig4_e.pdf", unit = "in", height = 20, width = 20) 

#### figure 5
# make tree with fastree
library(ggtree)
library(ggtreeExtra)
df2 <- merge(df1, names_id_filtered)
tree2 <- read.tree("2025Results/16s_newick_tree_3")
p9 <- ggtree(tree2,size =0.1) + geom_tiplab(size = 0)
gheatmap(p9, df5, offset=0.2, width=0.02, 
         colnames=TRUE,colnames_offset_x = 0, color = NA,
         hjust = 1, legend_title="Taxa", colnames_angle = 90, 
         colnames_offset_y = 0) + scale_fill_manual(values = tcols)
ggsave(plot = last_plot(), file = "2025Results/fig5_test.pdf", unit = "in", height = 40, width = 40, limitsize = FALSE) 
# okay
tree2 <- read.tree("2025Results/16s_newick_tree_4")
p9 <- ggtree(tree2,size =0.1) + geom_tiplab(size = 0)
gheatmap(p9, df5, offset=0.2, width=0.02, 
         colnames=TRUE,colnames_offset_x = 0, color = NA,
         hjust = 1, legend_title="Taxa", colnames_angle = 90, 
         colnames_offset_y = 0) + scale_fill_manual(values = tcols)
ggsave(plot = last_plot(), file = "2025Results/fig5_test1.pdf", unit = "in", height = 40, width = 40, limitsize = FALSE) 
# version 4 of the tree rerooted is best

tree1 <- read.tree("2025Results/16s_newick_tree_4")
p7 <- ggtree(tree1, layout = "fan", open.angle=15, size =0.1) + geom_tiplab(size = 0)
df3 <- as.data.frame(df2[,2:39])
row.names(df3) <- df2$Genome
df4 <- subset.data.frame(df3, select = "Lifestyle")
df5 <- subset.data.frame(df3, select = "taxa_elm")
df6 <- our_genomes
row.names(df6) <- df6$Genome
df6 <- subset.data.frame(df6, select = "Mean.Completeness")
df3 <- subset.data.frame(df3, select = -c(Lifestyle, taxa_elm))
df3 <- df3[,match(order, names(df3))]
length(get_taxa_name(p7))

write.csv(df3, "2025Results/fig5_annotation_pathways.csv", quote = F, row.names = T)
write.csv(df4, "2025Results/fig5_annotation_lifestyles.csv", quote = F, row.names = T)
write.csv(df5, "2025Results/fig5_annotation_taxa.csv", quote = F, row.names = T)
write.csv(df6, "2025Results/fig5_annotation_genomecompletion.csv", quote = F, row.names = T)

# colors to taxa names
cbind(tcols, unique(df5$taxa_elm)[order(unique(df5$taxa_elm))])

gheatmap(p7, df3, offset=0.2, width=0.6, 
         colnames=TRUE, hjust = 1, legend_title="Present or Absent", colnames_angle = 90, colnames_offset_y = 40) + scale_fill_gradient(low = "white", high = "black")
ggsave(plot = last_plot(), file = "2025Results/fig5_a.pdf", unit = "in", height = 30, width = 30) 
gheatmap(p7, df4, offset=0.2, width=0.6, 
         colnames=TRUE, colnames_offset_x = 0, legend_title="Lifestyle", colnames_angle = 90, colnames_offset_y = 40) + scale_fill_manual(values = c("#0016a6", "#629e02", "#0271a8"))
ggsave(plot = last_plot(), file = "2025Results/fig5_b.pdf", unit = "in", height = 30, width = 30) 
gheatmap(p7, df5, offset=0.2, width=0.6, 
         colnames=TRUE, legend_title="Taxa", colnames_angle = 90, colnames_offset_y = 40)
ggsave(plot = last_plot(), file = "2025Results/fig5_c.pdf", unit = "in", height = 30, width = 30) 

library(ggnewscale)
p1 <- gheatmap(p7, df3, offset=0.2, width=0.6, 
         colnames=TRUE, hjust = 1, legend_title="Present or Absent", colnames_angle = 90, colnames_offset_y = 40) + scale_fill_gradient(low = "white", high = "black")
p2 <- p1 + new_scale_fill()
p2 <- gheatmap(p2, df4, offset=3, width=0.6, 
         colnames=TRUE, colnames_offset_x = 0, legend_title="Lifestyle", colnames_angle = 90, colnames_offset_y = 40) + scale_fill_manual(values = c("#0016a6", "#629e02", "#0271a8"))
p3 <- p2 + new_scale_fill()
gheatmap(p3, df5, offset=5, width=0.6, 
         colnames=TRUE, legend_title="Taxa", colnames_angle = 90, colnames_offset_y = 40)
ggsave(plot = last_plot(), file = "2025Results/fig5_d.pdf", unit = "in", height = 40, width = 40) 

p1 <- gheatmap(p7, df5, offset=0.2, width=0.02, 
               colnames=TRUE,colnames_offset_x = 0, hjust = 1, legend_title="Taxa", colnames_angle = 90, colnames_offset_y = 0) + scale_fill_manual(values = tcols)
p2 <- p1 + new_scale_fill()
p2 <- gheatmap(p2, df4, offset=0.3, width=0.02, 
               colnames=TRUE, colnames_offset_x = 0, hjust = 1,legend_title="Lifestyle", colnames_angle = 90, colnames_offset_y = 0) + scale_fill_manual(values = c("#0016a6", "#629e02", "#0271a8"))
p3 <- p2 + new_scale_fill()
gheatmap(p3, df3, offset=.4, width=0.82, 
         colnames=TRUE,colnames_offset_x = 0,hjust = 1, legend_title="Present or Absent", colnames_angle = 90, colnames_offset_y = 0) + scale_fill_gradient(low = "white", high = "black")
ggsave(plot = last_plot(), file = "2025Results/fig5_e.pdf", unit = "in", height = 40, width = 40) 

#round
p1 <- gheatmap(p7, df5, offset=0.2, width=0.02, 
               colnames=TRUE,colnames_offset_x = 0, color = NA,
               hjust = 1, legend_title="Taxa", colnames_angle = 90, 
               colnames_offset_y = 0) + scale_fill_manual(values = tcols)
p2 <- p1 + new_scale_fill()
p2 <- gheatmap(p2, df4, offset=0.31, width=0.02, 
               colnames=TRUE, colnames_offset_x = 0, color = NA,
               hjust = 1,legend_title="Lifestyle", colnames_angle = 90, 
               colnames_offset_y = 0) + 
  scale_fill_manual(values = c("#0016a6", "#629e02", "#0271a8"))
p3 <- p2 + new_scale_fill()
p3 <- gheatmap(p3, df6, offset=0.42, width=0.02, 
               colnames=TRUE, colnames_offset_x = 0, color = NA,
               hjust = 1,legend_title="Genome Completion", colnames_angle = 90, 
               colnames_offset_y = 0) + scale_fill_viridis(direction = -1, option = "A")
p4 <- p3 + new_scale_fill()
gheatmap(p4, df3, offset=.53, width=0.82, 
         colnames=TRUE,colnames_offset_x = 0, color = NA,
         hjust = 1, legend_title="Present or Absent", colnames_angle = 90, 
         colnames_offset_y = 0) + 
  scale_fill_gradient(low = "white", high = "black")
ggsave(plot = last_plot(), file = "2025Results/fig5_f.pdf", unit = "in", height = 40, width = 40) 
gheatmap(p4, df3, offset=.53, width=0.82, 
         colnames=TRUE,colnames_offset_x = 0, color = NA,
         hjust = 1, legend_title="Present or Absent", colnames_angle = 90, 
         colnames_offset_y = 0) + 
  scale_fill_gradient(low = "white", high = "#A88205")
ggsave(plot = last_plot(), file = "2025Results/fig5_f_ochre.pdf", unit = "in", height = 40, width = 40) 
gheatmap(p4, df3, offset=.53, width=0.82, 
         colnames=TRUE,colnames_offset_x = 0, color = NA,
         hjust = 1, legend_title="Present or Absent", colnames_angle = 90, 
         colnames_offset_y = 0) + 
  scale_fill_gradient(low = "white", high = "#2F397A")
ggsave(plot = last_plot(), file = "2025Results/fig5_f_blue.pdf", unit = "in", height = 40, width = 40) 
gheatmap(p4, df3, offset=.53, width=0.82, 
         colnames=TRUE,colnames_offset_x = 0, color = NA,
         hjust = 1, legend_title="Present or Absent", colnames_angle = 90, 
         colnames_offset_y = 0) + 
  scale_fill_gradient(low = "white", high = "#415521")
ggsave(plot = last_plot(), file = "2025Results/fig5_f_green.pdf", unit = "in", height = 40, width = 40) 



#rectangular
p8 <- ggtree(tree1, size =0.1) + geom_tiplab(size = 0) + xlim(0,10) + ylim(-1,4460)
p1 <- gheatmap(p8, df5, offset=0.2, width=0.02, 
               colnames=TRUE,colnames_offset_x = 0, color = NA,
               hjust = 1, legend_title="Taxa", colnames_angle = 90, 
               colnames_offset_y = 0) + scale_fill_manual(values = tcols)
p2 <- p1 + new_scale_fill()
p2 <- gheatmap(p2, df4, offset=0.31, width=0.02, 
               colnames=TRUE, colnames_offset_x = 0, color = NA,
               hjust = 1,legend_title="Lifestyle", colnames_angle = 90, 
               colnames_offset_y = 0) + 
  scale_fill_manual(values = c("#0016a6", "#629e02", "#0271a8"))
p3 <- p2 + new_scale_fill()
p3 <- gheatmap(p3, df6, offset=0.42, width=0.02, 
               colnames=TRUE, colnames_offset_x = 0, color = NA,
               hjust = 1,legend_title="Genome Completion", colnames_angle = 90, 
               colnames_offset_y = 0) + scale_fill_viridis(direction = -1, option = "A")
p4 <- p3 + new_scale_fill()
gheatmap(p4, df3, offset=.53, width=0.82, 
         colnames=TRUE,colnames_offset_x = 0, color = NA,
         hjust = 1, legend_title="Present or Absent", colnames_angle = 90, 
         colnames_offset_y = 0) + 
  scale_fill_gradient(low = "white", high = "black")
ggsave(plot = last_plot(), file = "2025Results/fig5_g.pdf", unit = "in", height = 60, width = 30, limitsize = FALSE) 
gheatmap(p4, df3, offset=.53, width=0.82, 
         colnames=TRUE,colnames_offset_x = 0, color = NA,
         hjust = 1, legend_title="Present or Absent", colnames_angle = 90, 
         colnames_offset_y = 0) + 
  scale_fill_gradient(low = "white", high = "#A88205")
ggsave(plot = last_plot(), file = "2025Results/fig5_g_ochre.pdf", unit = "in", height = 40, width = 40) 
gheatmap(p4, df3, offset=.53, width=0.82, 
         colnames=TRUE,colnames_offset_x = 0, color = NA,
         hjust = 1, legend_title="Present or Absent", colnames_angle = 90, 
         colnames_offset_y = 0) + 
  scale_fill_gradient(low = "white", high = "#2F397A")
ggsave(plot = last_plot(), file = "2025Results/fig5_g_blue.pdf", unit = "in", height = 40, width = 40) 
gheatmap(p4, df3, offset=.53, width=0.82, 
         colnames=TRUE,colnames_offset_x = 0, color = NA,
         hjust = 1, legend_title="Present or Absent", colnames_angle = 90, 
         colnames_offset_y = 0) + 
  scale_fill_gradient(low = "white", high = "#415521")
ggsave(plot = last_plot(), file = "2025Results/fig5_g_green.pdf", unit = "in", height = 40, width = 40) 


# figure 3a
library(ggstatsplot)
data <- read.csv("20250816_protein.filt.sum.csv", row.names = 1)
s_content_df <- merge(x = our_genomes[,c(2,20,23)], y = data, by.x = "Genome", by.y = "genome")
ggplot(s_content_df, aes(x = factor(Lifestyle, levels = c("assimilator", "incomplete", "auxotroph")), y = perTot.mu.wt)) + 
  geom_violin(aes(color = factor(Lifestyle))) +theme_light()+ geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  scale_colour_manual(values = c(values3[1], values3[3], values3[2]))
ggsave(plot = last_plot(), file = "2025Results/fig3a_1.pdf", unit = "in", height = 7, width = 10) 

plt <- ggbetweenstats(
  data = s_content_df,
  x = Lifestyle,
  y = perTot.mu.wt, digits = 3,
  point.dodge.width = 0.8)
plt + labs(
  x = "Lifestyle", 
  title = "Distribution of S Proteome Content across S Lifestyles"
) + scale_colour_manual(values = c(values3[1], values3[3], values3[2]))
ggsave(plot = last_plot(), file = "2025Results/fig3a_2.pdf", unit = "in", height = 7, width = 10) 

s_content_df_temp <- s_content_df %>% filter(taxa_elm != "Alpha-Pelagibacter")
plt <- ggbetweenstats(
  data = s_content_df_temp,
  x = Lifestyle,
  y = perTot.mu.wt, digits = 3,
  point.dodge.width = 0.8)
plt + labs(
  x = "Lifestyle",
  title = "Distribution of S Proteome Content across S Lifestyles"
) + scale_colour_manual(values = c(values3[1], values3[3], values3[2]))
ggsave(plot = last_plot(), file = "2025Results/fig3a_3_sansSar11.pdf", unit = "in", height = 7, width = 10) 

s_content_df_temp <- s_content_df %>% filter(taxa_elm != "Alpha-Pelagibacter") %>% filter(taxa_elm != "Firmicutes")
plt <- ggbetweenstats(
  data = s_content_df_temp,
  x = Lifestyle,
  y = perTot.mu.wt, digits = 3,
  point.dodge.width = 0.8)
plt + labs(
  x = "Lifestyle",
  title = "Distribution of S Proteome Content across S Lifestyles"
) + scale_colour_manual(values = c(values3[1], values3[3], values3[2]))
ggsave(plot = last_plot(), file = "2025Results/fig3a_3_sansSar11Frim.pdf", unit = "in", height = 7, width = 10) 


s_content_df_temp <- s_content_df[grep(x = s_content_df$Genome, pattern = "_REFG_"),]
plt <- ggbetweenstats(
  data = s_content_df_temp,
  x = Lifestyle,
  y = perTot.mu.wt, digits = 3,
  point.dodge.width = 0.8)
plt + labs(
  x = "Lifestyle",
  title = "Distribution of S Proteome Content across S Lifestyles"
) + scale_colour_manual(values = c(values3[1], values3[3], values3[2]))
ggsave(plot = last_plot(), file = "2025Results/fig3a_4_onlyrefs.pdf", unit = "in", height = 7, width = 10) 

plt <- ggbetweenstats(
  data = s_content_df,
  x = Lifestyle,
  y = perTot.mu.wt,
  type = "r", digits = 3,
  point.dodge.width = 0.8)
plt + labs(
  x = "Lifestyle",
  title = "Distribution of S Proteome Content across S Lifestyles"
) + scale_colour_manual(values = c(values3[1], values3[3], values3[2]))
ggsave(plot = last_plot(), file = "2025Results/fig3a_2_robust.pdf", unit = "in", height = 7, width = 10) 

s_content_df_temp <- s_content_df %>% filter(taxa_elm != "Alpha-Pelagibacter")
plt <- ggbetweenstats(
  data = s_content_df_temp,
  x = Lifestyle,
  y = perTot.mu.wt,
  type = "r", digits = 3,
  point.dodge.width = 0.8)
plt + labs(
  x = "Lifestyle",
  title = "Distribution of S Proteome Content across S Lifestyles"
) + scale_colour_manual(values = c(values3[1], values3[3], values3[2]))
ggsave(plot = last_plot(), file = "2025Results/fig3a_3_sansSar11_robust.pdf", unit = "in", height = 7, width = 10) 

s_content_df_temp <- s_content_df %>% filter(taxa_elm != "Alpha-Pelagibacter") %>% filter(taxa_elm != "Firmicutes")
plt <- ggbetweenstats(
  data = s_content_df_temp,
  x = Lifestyle,
  y = perTot.mu.wt,
  type = "r", digits = 3,
  point.dodge.width = 0.8)
plt + labs(
  x = "Lifestyle",
  title = "Distribution of S Proteome Content across S Lifestyles"
) + scale_colour_manual(values = c(values3[1], values3[3], values3[2]))
ggsave(plot = last_plot(), file = "2025Results/fig3a_3_sansSar11Frim_robust.pdf", unit = "in", height = 7, width = 10) 


s_content_df_temp <- s_content_df[grep(x = s_content_df$Genome, pattern = "_REFG_"),]
plt <- ggbetweenstats(
  data = s_content_df_temp,
  x = Lifestyle,
  y = perTot.mu.wt,
  type = "r", digits = 3,
  point.dodge.width = 0.8)
plt + labs(
  x = "Lifestyle",
  title = "Distribution of S Proteome Content across S Lifestyles"
) + scale_colour_manual(values = c(values3[1], values3[3], values3[2]))
ggsave(plot = last_plot(), file = "2025Results/fig3a_4_onlyrefs_robust.pdf", unit = "in", height = 7, width = 10) 



# figure 3b
library(ggridges)
s_content_df$taxa_elm <- as.factor(s_content_df$taxa_elm)
s_content_df$Lifestyle <- as.factor(s_content_df$Lifestyle)
df1 <- s_content_df %>% 
  group_by(taxa_elm) %>%
  summarise(m = mean(perTot.mu.wt))
fig3_order <- df1$taxa_elm[order(df1$m)] 
ggplot(s_content_df, aes(x = perTot.mu.wt, y = factor(taxa_elm, fig3_order), fill = Lifestyle)) + 
  geom_density_ridges(alpha=0.6) +   theme(axis.title.x=element_blank()) +
  theme_ridges() + scale_fill_manual(values = c("#0016a6", "#0271a8", "#629e02"))
# ggplot(s_content_df, aes(x = perTot.mu.wt, y = MaxCompletion, color = taxa_elm)) + geom_point()
# ggplot(s_content_df, aes(x = perTot.mu.wt, y = taxa_elm, color = taxa_elm)) + geom_violin()
# ggplot(s_content_df, aes(x = perTot.mu.wt, y = taxa_elm, color = taxa_elm)) + geom_boxplot()
ggsave(plot = last_plot(), file = "2025Results/fig3b_1.pdf", unit = "in", height = 10, width = 10) 
 
s_content_df_manu <- s_content_df %>% filter(taxa_elm %in% c("Archaea-Nitrosopumilaceae", "Archaea-Poseidoniia MGIIa", "Archaea-Poseidoniia MGIIb", "Alpha-Puniceispirilla", "Archaea-Halobacteria", "Alpha-Pelagibacter", "Dadabacteria", "Actinobacteriota"))
df1 <- s_content_df_manu %>% 
  group_by(taxa_elm) %>%
  summarise(m = mean(perTot.mu.wt))
fig3_order <- df1$taxa_elm[order(df1$m)] 
ggplot(s_content_df_manu, aes(x = perTot.mu.wt, y = factor(taxa_elm, fig3_order), fill = Lifestyle)) + 
  geom_density_ridges(alpha=0.6) + 
  theme_ridges() + scale_fill_manual(values = c("#0016a6", "#0271a8", "#629e02"))
# ggplot(s_content_df, aes(x = perTot.mu.wt, y = MaxCompletion, color = taxa_elm)) + geom_point()
# ggplot(s_content_df, aes(x = perTot.mu.wt, y = taxa_elm, color = taxa_elm)) + geom_violin()
# ggplot(s_content_df, aes(x = perTot.mu.wt, y = taxa_elm, color = taxa_elm)) + geom_boxplot()
ggsave(plot = last_plot(), file = "2025Results/fig3b_2.pdf", unit = "in", height = 10, width = 10) 


for (taxa in unique(s_content_df$taxa_elm)){
  s_content_df_temp <- s_content_df %>% filter(taxa_elm == taxa)
  p1 <- ggbetweenstats(
    data = s_content_df_temp,
    x = Lifestyle,
    y = perTot.mu.wt,
    violin.args = list(width = 0),
    type = "p",
    conf.level = 0.99,
    title = "Parametric test",
    package = "ggsci",
    palette = "nrc_npg"
  )
  
  p2 <- ggbetweenstats(
    data = s_content_df_temp,
    x = Lifestyle,
    y = perTot.mu.wt,
    boxplot.args = list(width = 0),
    type = "np",
    conf.level = 0.99,
    title = "Non-parametric Test",
    package = "ggsci",
    palette = "uniform_startrek"
  )
  
  p3 <- ggbetweenstats(
    data = s_content_df_temp,
    x = Lifestyle,
    y = perTot.mu.wt,
    type = "r",
    conf.level = 0.99,
    title = "Robust Test",
    tr = 0.005,
    package = "wesanderson",
    palette = "Royal2",
    digits = 3
  )
  
  ## Bayes Factor for parametric t-test and boxviolin plot
  p4 <- ggbetweenstats(
    data = s_content_df_temp,
    x = Lifestyle,
    y = perTot.mu.wt,
    type = "bayes",
    violin.args = list(width = 0),
    boxplot.args = list(width = 0),
    point.args = list(alpha = 0),
    title = "Bayesian Test",
    package = "ggsci",
    palette = "nrc_npg"
  )
  
  ## combining the individual plots into a single plot
  combine_plots(
    list(p1, p2, p3, p4),
    plotgrid.args = list(nrow = 2L),
    annotation.args = list(
      title = paste0("Comparison of S content for lifesytles within ", taxa)
    )
  )
  ggsave(plot = last_plot(), file = paste0("2025Results/fig3b_", taxa, ".pdf"), unit = "in", height = 15, width = 15) 
}

for (taxa in unique(s_content_df$taxa_elm)[20:22]){
  s_content_df_temp <- s_content_df %>% filter(taxa_elm == taxa)
  p1 <- ggbetweenstats(
    data = s_content_df_temp,
    x = Lifestyle,
    y = perTot.mu.wt,
    violin.args = list(width = 0),
    type = "p",
    conf.level = 0.99,
    title = "Parametric test",
    package = "ggsci",
    palette = "nrc_npg"
  )
  
  p2 <- ggbetweenstats(
    data = s_content_df_temp,
    x = Lifestyle,
    y = perTot.mu.wt,
    boxplot.args = list(width = 0),
    type = "np",
    conf.level = 0.99,
    title = "Non-parametric Test",
    package = "ggsci",
    palette = "uniform_startrek"
  )
  
  p3 <- ggbetweenstats(
    data = s_content_df_temp,
    x = Lifestyle,
    y = perTot.mu.wt,
    type = "r",
    conf.level = 0.99,
    title = "Robust Test",
    tr = 0.005,
    package = "wesanderson",
    palette = "Royal2",
    digits = 3
  )
  
  ## Bayes Factor for parametric t-test and boxviolin plot
  p4 <- ggbetweenstats(
    data = s_content_df_temp,
    x = Lifestyle,
    y = perTot.mu.wt,
    type = "bayes",
    violin.args = list(width = 0),
    boxplot.args = list(width = 0),
    point.args = list(alpha = 0),
    title = "Bayesian Test",
    package = "ggsci",
    palette = "nrc_npg"
  )
  
  ## combining the individual plots into a single plot
  combine_plots(
    list(p1, p2, p3, p4),
    plotgrid.args = list(nrow = 2L),
    annotation.args = list(
      title = paste0("Comparison of S content for lifesytles within ", taxa)
    )
  )
  ggsave(plot = last_plot(), file = paste0("2025Results/fig3b_", taxa, ".pdf"), unit = "in", height = 15, width = 15) 
} 

for (taxa in unique(s_content_df$taxa_elm)[24:30]){
  s_content_df_temp <- s_content_df %>% filter(taxa_elm == taxa)
  p1 <- ggbetweenstats(
    data = s_content_df_temp,
    x = Lifestyle,
    y = perTot.mu.wt,
    violin.args = list(width = 0),
    type = "p",
    conf.level = 0.99,
    title = "Parametric test",
    package = "ggsci",
    palette = "nrc_npg"
  )
  
  p2 <- ggbetweenstats(
    data = s_content_df_temp,
    x = Lifestyle,
    y = perTot.mu.wt,
    boxplot.args = list(width = 0),
    type = "np",
    conf.level = 0.99,
    title = "Non-parametric Test",
    package = "ggsci",
    palette = "uniform_startrek"
  )
  
  p3 <- ggbetweenstats(
    data = s_content_df_temp,
    x = Lifestyle,
    y = perTot.mu.wt,
    type = "r",
    conf.level = 0.99,
    title = "Robust Test",
    tr = 0.005,
    package = "wesanderson",
    palette = "Royal2",
    digits = 3
  )
  
  ## Bayes Factor for parametric t-test and boxviolin plot
  p4 <- ggbetweenstats(
    data = s_content_df_temp,
    x = Lifestyle,
    y = perTot.mu.wt,
    type = "bayes",
    violin.args = list(width = 0),
    boxplot.args = list(width = 0),
    point.args = list(alpha = 0),
    title = "Bayesian Test",
    package = "ggsci",
    palette = "nrc_npg"
  )
  
  ## combining the individual plots into a single plot
  combine_plots(
    list(p1, p2, p3, p4),
    plotgrid.args = list(nrow = 2L),
    annotation.args = list(
      title = paste0("Comparison of S content for lifesytles within ", taxa)
    )
  )
  ggsave(plot = last_plot(), file = paste0("2025Results/fig3b_", taxa, ".pdf"), unit = "in", height = 15, width = 15) 
} 

s_content_df %>% filter(taxa_elm %in% unique(s_content_df$taxa_elm)[1:6]) %>% 
grouped_ggbetweenstats(
  ## arguments relevant for ggbetweenstats
  x = Lifestyle,
  y = perTot.mu.wt,
  grouping.var = taxa_elm,
  pairwise.display = "significant", ## display only significant pairwise comparisons
  p.adjust.method = "fdr", ## adjust p-values for multiple tests using this method
  # ggtheme = ggthemes::theme_tufte(),
  package = "ggsci",
  palette = "default_jco",
  ## arguments relevant for combine_plots
  annotation.args = list(title = "Comparison of S content for lifesytles between taxa 1-6"),
  plotgrid.args = list(ncol = 3)
)
ggsave(plot = last_plot(), file = paste0("2025Results/fig3b_taxa1.pdf"), unit = "in", height = 20, width = 30)

s_content_df %>% filter(taxa_elm %in% unique(s_content_df$taxa_elm)[7:12]) %>% 
  grouped_ggbetweenstats(
    ## arguments relevant for ggbetweenstats
    x = Lifestyle,
    y = perTot.mu.wt,
    grouping.var = taxa_elm,
    pairwise.display = "significant", ## display only significant pairwise comparisons
    p.adjust.method = "fdr", ## adjust p-values for multiple tests using this method
    # ggtheme = ggthemes::theme_tufte(),
    package = "ggsci",
    palette = "default_jco",
    ## arguments relevant for combine_plots
    annotation.args = list(title = "Comparison of S content for lifesytles between taxa 7-12"),
    plotgrid.args = list(ncol = 3)
  )
ggsave(plot = last_plot(), file = paste0("2025Results/fig3b_taxa2.pdf"), unit = "in", height = 20, width = 30)

s_content_df %>% filter(taxa_elm %in% unique(s_content_df$taxa_elm)[13:18]) %>% 
  grouped_ggbetweenstats(
    ## arguments relevant for ggbetweenstats
    x = Lifestyle,
    y = perTot.mu.wt,
    grouping.var = taxa_elm,
    pairwise.display = "significant", ## display only significant pairwise comparisons
    p.adjust.method = "fdr", ## adjust p-values for multiple tests using this method
    # ggtheme = ggthemes::theme_tufte(),
    package = "ggsci",
    palette = "default_jco",
    ## arguments relevant for combine_plots
    annotation.args = list(title = "Comparison of S content for lifesytles between taxa 7-12"),
    plotgrid.args = list(ncol = 3)
  )
ggsave(plot = last_plot(), file = paste0("2025Results/fig3b_taxa3.pdf"), unit = "in", height = 20, width = 30)

s_content_df %>% filter(taxa_elm %in% unique(s_content_df$taxa_elm)[19:24]) %>% 
  grouped_ggbetweenstats(
    ## arguments relevant for ggbetweenstats
    x = Lifestyle,
    y = perTot.mu.wt,
    grouping.var = taxa_elm,
    pairwise.display = "significant", ## display only significant pairwise comparisons
    p.adjust.method = "fdr", ## adjust p-values for multiple tests using this method
    # ggtheme = ggthemes::theme_tufte(),
    package = "ggsci",
    palette = "default_jco",
    ## arguments relevant for combine_plots
    annotation.args = list(title = "Comparison of S content for lifesytles between taxa 7-12"),
    plotgrid.args = list(ncol = 3)
  )
ggsave(plot = last_plot(), file = paste0("2025Results/fig3b_taxa4.pdf"), unit = "in", height = 20, width = 30)

s_content_df %>% filter(taxa_elm %in% unique(s_content_df$taxa_elm)[25:30]) %>% 
  grouped_ggbetweenstats(
    ## arguments relevant for ggbetweenstats
    x = Lifestyle,
    y = perTot.mu.wt,
    grouping.var = taxa_elm,
    pairwise.display = "significant", ## display only significant pairwise comparisons
    p.adjust.method = "fdr", ## adjust p-values for multiple tests using this method
    # ggtheme = ggthemes::theme_tufte(),
    package = "ggsci",
    palette = "default_jco",
    ## arguments relevant for combine_plots
    annotation.args = list(title = "Comparison of S content for lifesytles between taxa 7-12"),
    plotgrid.args = list(ncol = 3)
  )
ggsave(plot = last_plot(), file = paste0("2025Results/fig3b_taxa5.pdf"), unit = "in", height = 20, width = 30)

s_content_df_manu <- s_content_df %>% filter(taxa_elm %in% c("Gemmatimonadota", "Archaea-Poseidoniia MGIIa", "Archaea-Poseidoniia MGIIb", "Planctomycetota", "Archaea-Halobacteria", "SAR324"))
df1 <- s_content_df_manu %>% 
  group_by(taxa_elm) %>%
  summarise(m = mean(perTot.mu.wt))
fig3_order <- df1$taxa_elm[order(df1$m)] 
ggplot(s_content_df_manu, aes(x = perTot.mu.wt, y = factor(taxa_elm, fig3_order), fill = Lifestyle)) + 
  geom_density_ridges(alpha=0.6, scale = 1) + 
  theme_ridges() + scale_fill_manual(values = c("#0016a6", "#0271a8", "#629e02"))
# ggplot(s_content_df, aes(x = perTot.mu.wt, y = MaxCompletion, color = taxa_elm)) + geom_point()
# ggplot(s_content_df, aes(x = perTot.mu.wt, y = taxa_elm, color = taxa_elm)) + geom_violin()
# ggplot(s_content_df, aes(x = perTot.mu.wt, y = taxa_elm, color = taxa_elm)) + geom_boxplot()

ggsave(plot = last_plot(), file = "2025Results/fig3b_3.pdf", unit = "in", height = 10, width = 10) 

s_content_df_manu <- s_content_df %>% filter(taxa_elm %in% c( "Alpha-Pelagibacter", "Gamma-SAR92", "Alpha-Rhodobacter"))
df1 <- s_content_df_manu %>% 
  group_by(taxa_elm) %>%
  summarise(m = mean(perTot.mu.wt))
fig3_order <- df1$taxa_elm[order(df1$m)] 
ggplot(s_content_df_manu, aes(x = perTot.mu.wt, y = factor(taxa_elm, fig3_order), fill = Lifestyle)) + 
  geom_density_ridges(alpha=0.6, scale = 1) + 
  theme_ridges() + scale_fill_manual(values = c("#0016a6", "#0271a8", "#629e02"))
# ggplot(s_content_df, aes(x = perTot.mu.wt, y = MaxCompletion, color = taxa_elm)) + geom_point()
# ggplot(s_content_df, aes(x = perTot.mu.wt, y = taxa_elm, color = taxa_elm)) + geom_violin()
# ggplot(s_content_df, aes(x = perTot.mu.wt, y = taxa_elm, color = taxa_elm)) + geom_boxplot()

ggsave(plot = last_plot(), file = "2025Results/fig3b_4_model.pdf", unit = "in", height = 10, width = 10) 


df1 <- s_content_df %>% 
  group_by(taxa_elm) %>%
  summarise(m = mean(perTot.mu.wt))
fig3_order <- df1$taxa_elm[order(df1$m)] 
ggplot(s_content_df, aes(x = perTot.mu.wt, y = factor(taxa_elm, fig3_order), fill = Lifestyle)) + 
  geom_density_ridges(alpha=0.6, scale = 1) + 
  theme_ridges() + scale_fill_manual(values = c("#0016a6", "#0271a8", "#629e02"))
# ggplot(s_content_df, aes(x = perTot.mu.wt, y = MaxCompletion, color = taxa_elm)) + geom_point()
# ggplot(s_content_df, aes(x = perTot.mu.wt, y = taxa_elm, color = taxa_elm)) + geom_violin()
# ggplot(s_content_df, aes(x = perTot.mu.wt, y = taxa_elm, color = taxa_elm)) + geom_boxplot()

ggsave(plot = last_plot(), file = "2025Results/fig3b_5.pdf", unit = "in", height = 10, width = 10) 

s_content_df_manu <- s_content_df %>% filter(taxa_elm %in% c( "Chloroflexota", "Alpha-Other", "Alpha-Rhodobacter"))
df1 <- s_content_df_manu %>% 
  group_by(taxa_elm) %>%
  summarise(m = mean(perTot.mu.wt))
fig3_order <- df1$taxa_elm[order(df1$m)] 
ggplot(s_content_df_manu, aes(x = perTot.mu.wt, y = factor(taxa_elm, fig3_order), fill = Lifestyle)) + 
  geom_density_ridges(alpha=0.6, scale = 1) + 
  theme_ridges() + scale_fill_manual(values = c("#0016a6", "#0271a8", "#629e02"))
# ggplot(s_content_df, aes(x = perTot.mu.wt, y = MaxCompletion, color = taxa_elm)) + geom_point()
# ggplot(s_content_df, aes(x = perTot.mu.wt, y = taxa_elm, color = taxa_elm)) + geom_violin()
# ggplot(s_content_df, aes(x = perTot.mu.wt, y = taxa_elm, color = taxa_elm)) + geom_boxplot()

ggsave(plot = last_plot(), file = "2025Results/fig3b_6_statdifferent.pdf", unit = "in", height = 10, width = 10) 

# figure 3a new
library(ggstatsplot)
data <- read.csv("20250905_protein.filt.sum.csv", row.names = 1)
s_content_df <- merge(x = our_genomes[,c(2,20,23)], y = data, by.x = "Genome", by.y = "genome")
ggplot(s_content_df, aes(x = factor(Lifestyle, levels = c("assimilator", "incomplete", "auxotroph")), y = avgS)) + 
  geom_violin(aes(color = factor(Lifestyle))) +theme_light()+ geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  scale_colour_manual(values = c(values3[1], values3[3], values3[2]))
ggsave(plot = last_plot(), file = "2025Results/newProteinFilter_09052025/fig3a_1.pdf", unit = "in", height = 7, width = 10) 

plt <- ggbetweenstats(
  data = s_content_df,
  x = Lifestyle,
  y = avgS, digits = 3,
  point.dodge.width = 0.8)
plt + labs(
  x = "Lifestyle", 
  title = "Distribution of S Proteome Content across S Lifestyles"
) + scale_colour_manual(values = c(values3[1], values3[3], values3[2]))
ggsave(plot = last_plot(), file = "2025Results/newProteinFilter_09052025/fig3a_2.pdf", unit = "in", height = 7, width = 10) 

s_content_df_temp <- s_content_df %>% filter(taxa_elm != "Alpha-Pelagibacter")
plt <- ggbetweenstats(
  data = s_content_df_temp,
  x = Lifestyle,
  y = avgS, digits = 3,
  point.dodge.width = 0.8)
plt + labs(
  x = "Lifestyle",
  title = "Distribution of S Proteome Content across S Lifestyles"
) + scale_colour_manual(values = c(values3[1], values3[3], values3[2]))
ggsave(plot = last_plot(), file = "2025Results/newProteinFilter_09052025/fig3a_3_sansSar11.pdf", unit = "in", height = 7, width = 10) 

s_content_df_temp <- s_content_df %>% filter(taxa_elm != "Alpha-Pelagibacter") %>% filter(taxa_elm != "Firmicutes")
plt <- ggbetweenstats(
  data = s_content_df_temp,
  x = Lifestyle,
  y = avgS, digits = 3,
  point.dodge.width = 0.8)
plt + labs(
  x = "Lifestyle",
  title = "Distribution of S Proteome Content across S Lifestyles"
) + scale_colour_manual(values = c(values3[1], values3[3], values3[2]))
ggsave(plot = last_plot(), file = "2025Results/newProteinFilter_09052025/fig3a_3_sansSar11Frim.pdf", unit = "in", height = 7, width = 10) 


s_content_df_temp <- s_content_df[grep(x = s_content_df$Genome, pattern = "_REFG_"),]
plt <- ggbetweenstats(
  data = s_content_df_temp,
  x = Lifestyle,
  y = avgS, digits = 3,
  point.dodge.width = 0.8)
plt + labs(
  x = "Lifestyle",
  title = "Distribution of S Proteome Content across S Lifestyles"
) + scale_colour_manual(values = c(values3[1], values3[3], values3[2]))
ggsave(plot = last_plot(), file = "2025Results/newProteinFilter_09052025/fig3a_4_onlyrefs.pdf", unit = "in", height = 7, width = 10) 

plt <- ggbetweenstats(
  data = s_content_df,
  x = Lifestyle,
  y = avgS,
  type = "r", digits = 3,
  point.dodge.width = 0.8)
plt + labs(
  x = "Lifestyle",
  title = "Distribution of S Proteome Content across S Lifestyles"
) + scale_colour_manual(values = c(values3[1], values3[3], values3[2]))
ggsave(plot = last_plot(), file = "2025Results/newProteinFilter_09052025/fig3a_2_robust.pdf", unit = "in", height = 7, width = 10) 

s_content_df_temp <- s_content_df %>% filter(taxa_elm != "Alpha-Pelagibacter")
plt <- ggbetweenstats(
  data = s_content_df_temp,
  x = Lifestyle,
  y = avgS,
  type = "r", digits = 3,
  point.dodge.width = 0.8)
plt + labs(
  x = "Lifestyle",
  title = "Distribution of S Proteome Content across S Lifestyles"
) + scale_colour_manual(values = c(values3[1], values3[3], values3[2]))
ggsave(plot = last_plot(), file = "2025Results/newProteinFilter_09052025/fig3a_3_sansSar11_robust.pdf", unit = "in", height = 7, width = 10) 

s_content_df_temp <- s_content_df %>% filter(taxa_elm != "Alpha-Pelagibacter") %>% filter(taxa_elm != "Firmicutes")
plt <- ggbetweenstats(
  data = s_content_df_temp,
  x = Lifestyle,
  y = avgS,
  type = "r", digits = 3,
  point.dodge.width = 0.8)
plt + labs(
  x = "Lifestyle",
  title = "Distribution of S Proteome Content across S Lifestyles"
) + scale_colour_manual(values = c(values3[1], values3[3], values3[2]))
ggsave(plot = last_plot(), file = "2025Results/newProteinFilter_09052025/fig3a_3_sansSar11Frim_robust.pdf", unit = "in", height = 7, width = 10) 


s_content_df_temp <- s_content_df[grep(x = s_content_df$Genome, pattern = "_REFG_"),]
plt <- ggbetweenstats(
  data = s_content_df_temp,
  x = Lifestyle,
  y = avgS,
  type = "r", digits = 3,
  point.dodge.width = 0.8)
plt + labs(
  x = "Lifestyle",
  title = "Distribution of S Proteome Content across S Lifestyles"
) + scale_colour_manual(values = c(values3[1], values3[3], values3[2]))
ggsave(plot = last_plot(), file = "2025Results/newProteinFilter_09052025/fig3a_4_onlyrefs_robust.pdf", unit = "in", height = 7, width = 10) 

#new figure 3b
s_content_df$taxa_elm <- as.factor(s_content_df$taxa_elm)
s_content_df$Lifestyle <- as.factor(s_content_df$Lifestyle)
df1 <- s_content_df %>% 
  group_by(taxa_elm) %>%
  summarise(m = mean(avgS))
fig3_order <- df1$taxa_elm[order(df1$m)] 
ggplot(s_content_df, aes(x = avgS, y = factor(taxa_elm, fig3_order), fill = Lifestyle)) + 
  geom_density_ridges(alpha=0.6, scale =1) +   theme(axis.title.x=element_blank()) +
  theme_ridges() + scale_fill_manual(values = c("#0016a6", "#0271a8", "#629e02"))
# ggplot(s_content_df, aes(x = perTot.mu.wt, y = MaxCompletion, color = taxa_elm)) + geom_point()
# ggplot(s_content_df, aes(x = perTot.mu.wt, y = taxa_elm, color = taxa_elm)) + geom_violin()
# ggplot(s_content_df, aes(x = perTot.mu.wt, y = taxa_elm, color = taxa_elm)) + geom_boxplot()
ggsave(plot = last_plot(), file = "2025Results/newProteinFilter_09052025/fig3b_1.pdf", unit = "in", height = 20, width = 10) 

ggplot(s_content_df, aes(x = avgS, y = factor(taxa_elm, fig3_order), fill = Lifestyle)) + 
  geom_density_ridges(alpha=0.6, stat = "binline", scale = 1) +   theme(axis.title.x=element_blank()) +
  theme_ridges() + scale_fill_manual(values = c("#0016a6", "#0271a8", "#629e02"))
ggsave(plot = last_plot(), file = "2025Results/newProteinFilter_09052025/fig3b_2.pdf", unit = "in", height = 20, width = 10) 

ggplot(s_content_df) + 
  geom_density_ridges(aes(x = avgS, y = factor(taxa_elm, fig3_order), fill = Lifestyle, height = after_stat(count)),
                                            stat = "density", scale = 1, alpha = 0.6) +   
  theme(axis.title.x=element_blank()) +
  theme_ridges() + scale_fill_manual(values = c("#0016a6", "#0271a8", "#629e02"))
ggsave(plot = last_plot(), file = "2025Results/newProteinFilter_09052025/fig3b_3.pdf", unit = "in", height = 20, width = 10) 

ggplot(s_content_df) + 
  geom_density_ridges(aes(x = avgS, y = factor(taxa_elm, fig3_order), fill = Lifestyle, height = after_stat(count)),
                      stat = "binline", scale = 1, alpha = 0.6) +   
  theme(axis.title.x=element_blank()) +
  theme_ridges() + scale_fill_manual(values = c("#0016a6", "#0271a8", "#629e02"))
ggsave(plot = last_plot(), file = "2025Results/newProteinFilter_09052025/fig3b_4.pdf", unit = "in", height = 20, width = 10) 

s_content_df %>% 
  ggplot(aes(x = avgS, fill = Lifestyle)) + 
  geom_histogram(alpha=0.6, position = 'identity') + facet_wrap(taxa_elm~., ncol = 1, scales = "free_y") + scale_fill_manual(values = c("#0016a6", "#0271a8", "#629e02"))
ggsave(plot = last_plot(), file = paste0("2025Results/newProteinFilter_09052025/fig3b_v2.pdf"), unit = "in", height = 40, width = 10) 


for (taxa in unique(s_content_df$taxa_elm)){
  s_content_df %>% filter(taxa_elm == taxa) %>% ggplot(aes(x = avgS, fill = Lifestyle)) + geom_histogram(alpha=0.6,position = 'identity') + facet_wrap(~Lifestyle)
  ggsave(plot = last_plot(), file = paste0("2025Results/newProteinFilter_09052025/fig3b_",taxa,".pdf"), unit = "in", height = 5, width = 10) 
}

# S auxoptrophs : oligotrophs and copiotrophs

fig3_order2 <- c("Alpha-Pelagibacter","Gamma-SAR86","SAR324","Patescibacteria", "Firmicutes","Alpha-Puniceispirilla", "Marinisomatota")
s_content_df %>% filter(taxa_elm %in% c("Alpha-Pelagibacter","Gamma-SAR86","Patescibacteria", "Firmicutes","SAR324","Alpha-Puniceispirilla", "Marinisomatota")) %>% 
  ggplot(aes(x = avgS, y = factor(taxa_elm, fig3_order), fill = Lifestyle)) + 
  geom_density_ridges(alpha=0.6, scale =1) +   theme(axis.title.x=element_blank()) +
  theme_ridges() + scale_fill_manual(values = c("#0016a6", "#0271a8", "#629e02"))
ggsave(plot = last_plot(), file = "2025Results/newProteinFilter_09052025/fig3b_1_Sauxotrophs.pdf", unit = "in", height = 10, width = 10) 

s_content_df %>% filter(taxa_elm %in% c("Alpha-Pelagibacter","Gamma-SAR86","Patescibacteria", "Firmicutes","SAR324","Alpha-Puniceispirilla", "Marinisomatota")) %>% 
  ggplot( aes(x = avgS, y = factor(taxa_elm, fig3_order2), fill = Lifestyle)) + 
  geom_density_ridges(alpha=0.6, stat = "binline", scale = 1) +   theme(axis.title.x=element_blank()) +
  theme_ridges() + scale_fill_manual(values = c("#0016a6", "#0271a8", "#629e02"))
ggsave(plot = last_plot(), file = "2025Results/newProteinFilter_09052025/fig3b_2_Sauxotrophs.pdf", unit = "in", height = 10, width = 10) 

s_content_df %>% filter(taxa_elm %in% c("Alpha-Pelagibacter","Gamma-SAR86","Patescibacteria", "Firmicutes","SAR324","Alpha-Puniceispirilla", "Marinisomatota")) %>% ggplot() + 
  geom_density_ridges(aes(x = avgS, y = factor(taxa_elm, fig3_order), fill = Lifestyle, height = after_stat(count)),
                      stat = "density", scale = 1, alpha = 0.6) +   
  theme(axis.title.x=element_blank()) +
  theme_ridges() + scale_fill_manual(values = c("#0016a6", "#0271a8", "#629e02"))
ggsave(plot = last_plot(), file = "2025Results/newProteinFilter_09052025/fig3b_3_Sauxotrophs.pdf", unit = "in", height = 10, width = 10) 

 s_content_df %>% filter(taxa_elm %in% c("Alpha-Pelagibacter","Gamma-SAR86","Patescibacteria", "Firmicutes","SAR324","Alpha-Puniceispirilla", "Marinisomatota")) %>% ggplot() + 
  geom_density_ridges(aes(x = avgS, y = factor(taxa_elm, fig3_order2), fill = Lifestyle, height = after_stat(count)),
                      stat = "binline", scale = 1, alpha = 0.6) +   
  theme(axis.title.x=element_blank()) +
  theme_ridges() + scale_fill_manual(values = c("#0016a6", "#0271a8", "#629e02"))
ggsave(plot = last_plot(), file = "2025Results/newProteinFilter_09052025/fig3b_4_Sauxotrophs.pdf", unit = "in", height = 10, width = 10) 
 

#playing with figure 4 again
#first, make a values6
values6 <- c(values3, "#9399bd", "#9fb6c2", "#c0d998")
df1 <- merge(maxcompletion, si_lifestyle, by.x = "Genome", by.y = "Genome")
#remove those pathays we aren't interested in
names(df1)[which(names(df1) == "DOS.Misc.Degradation.2")] <- "Sulfite.Export"
names(df1)[which(names(df1) == "DOS.Misc.Degradation.3")] <- "Sulfite.Dehydrogenase"
df1[,which(names(df1) == "DOS.Misc.Degradation.4")] <- NULL
df1[,which(names(df1) == "Methionine.Salvage.Part.1")] <- NULL
df1[,which(names(df1) == "Methionine.Salvage.Part.2")] <- NULL
df1[,which(names(df1) == "Methionine.Salvage.Fusion.Enzyme")] <- NULL
df1[,which(names(df1) == "Methionine.Salvage.Anaerobic.Isoprenoid.Shunt")] <- NULL
df1[,which(names(df1) == "Methionine.Salvage.Anaerobic.Degradation.of.MTRu1P")] <- NULL

df2 <- gather(df1, key = "key", value = "value", Cysteine.Biosynthesis:S.adenosyl.L.homocysteine.degradation)
nrow(df2)/length(unique(df2$key))
#none removed here
df2 <- df2 %>% filter(value==0 | value==1)
((36*16059)-nrow(df2))
#removing 35,204
((36*16059)-nrow(df2))/(36*16059)
# which is 6% of our dataset
df2 <- df2 %>% mutate(across(4, \(x)if_else(x==0,"Absent","Present")))
df2$value <- factor(df2$value)
df2$Lifestyle <- factor(df2$Lifestyle)
df2$key <- factor(df2$key)
df3 <- df2 %>% filter(value == "Present") %>% group_by(key) %>% summarise(count = n()) %>% arrange(., count)
order <- rev(df3$key)
p4 <- ggplot(df2, aes(x = factor(key, order), fill = value)) + geom_bar() + 
  theme_linedraw() +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + xlab("Pathway") + 
  facet_grid(Lifestyle ~ .)+
  scale_fill_manual(values = values6[c(6,2)]) +
  theme(plot.margin = unit(c(0,0,0,2), "cm")) + theme(axis.text.x = element_text(size = 14))
p4
ggsave(plot = last_plot(), file = "2025Results/fig4Varations/fig4_a.pdf", unit = "in", height = 20, width = 20) 
p4 <- ggplot(df2, aes(x = factor(key, order), fill = value)) + geom_bar() + 
  theme_linedraw() +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + xlab("Pathway") + 
  facet_grid(Lifestyle ~ .)+
  scale_fill_manual(values = values6[c(5,3)]) +
  theme(plot.margin = unit(c(0,0,0,2), "cm")) + theme(axis.text.x = element_text(size = 14))
p4
ggsave(plot = last_plot(), file = "2025Results/fig4Varations/fig4_b.pdf", unit = "in", height = 20, width = 20) 
p4 <- ggplot(df2, aes(x = factor(key, order), fill = value)) + geom_bar() + 
  theme_linedraw() +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + xlab("Pathway") + 
  facet_grid(Lifestyle ~ .)+
  scale_fill_manual(values = values6[c(4,1)]) +
  theme(plot.margin = unit(c(0,0,0,2), "cm")) + theme(axis.text.x = element_text(size = 14))
p4
ggsave(plot = last_plot(), file = "2025Results/fig4Varations/fig4_c.pdf", unit = "in", height = 20, width = 20) 
df1 %>% group_by(Lifestyle) %>% summarise(count = n())

fd_IV <- read.csv("2025Results/IndicatorValues.csv", row.names = 1)
df_Indgroup <- read.csv("2025Results/IndicatorGroups.csv", row.names = 1)[,1:2][duplicated(read.csv("2025Results/IndicatorGroups.csv", row.names = 1)[,1:2]),]
fd_IV$ind.vari[which(!(fd_IV$ind.vari %in% order))]
order[which(!(order %in% fd_IV$ind.vari))]
#Cysteine.Biosynthesis Sulfite.Export        Sulfite.Dehydrogenase
tmppaths <- c("Cysteine.Biosynthesis","Sulfite.Export", "Sulfite.Dehydrogenase")
tmpdf <- as.data.frame(cbind(tmppaths, c(NA, NA, NA), c(NA, NA, NA)))
names(tmpdf) <- names(fd_IV)
fd_IV <- fd_IV %>% filter(ind.vari %in% order)
fd_IV <- rbind(fd_IV, tmpdf)
tmpdf <- as.data.frame(cbind(tmppaths, c(NA, NA, NA)))
names(tmpdf) <- names(df_Indgroup)
df_Indgroup <- rbind(df_Indgroup, tmpdf)
fd_IV <- merge(fd_IV, df_Indgroup, by = "ind.vari")
fd_IV$value <- as.numeric(fd_IV$value)
ggplot(fd_IV, aes(x = factor(ind.vari, order), y = value, shape = key, color = indicator.group)) + 
  geom_point(size = 4) +
  scale_color_manual(values = c('#2F397A','#415521','#A88205')) +
  ylim(0,1) + theme_linedraw()+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.margin = unit(c(0,0,0,2), "cm")) 
ggsave(plot = last_plot(), file = "2025Results/fig4Varations/fig4_d1.pdf", unit = "in", height = 5, width = 20)   
ggplot(fd_IV, aes(x = factor(ind.vari, order), y = value, shape = key)) + 
  geom_point(color = "black", size = 4) +
  ylim(0,1) + theme_linedraw()+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.margin = unit(c(0,0,0,2), "cm")) 
ggsave(plot = last_plot(), file = "2025Results/fig4Varations/fig4_d2.pdf", unit = "in", height = 5, width = 20)   

# work on methionine salvage figure
# need to import the processed data for part 1 (MSP1), part 2 (MSP2), fusion enzyme (MSFE), anaerobic isoprenoid shunt (MSIS), anaerobic degradation (MSAD) 
MSP1 <- read.csv("kofam_results/Methionine Salvage Part 1processed_above75genomecompletion.csv", row.names = 1)
MSP2 <- read.csv("kofam_results/Methionine Salvage Part 2processed_above75genomecompletion.csv", row.names = 1)
MSFE <- read.csv("kofam_results/Methionine Salvage Fusion Enzymeprocessed_above75genomecompletion.csv", row.names = 1)
MSIS <- read.csv("kofam_results/Methionine Salvage Anaerobic Isoprenoid Shuntprocessed_above75genomecompletion.csv", row.names = 1)
MSAD <- read.csv("kofam_results/Methionine Salvage Anaerobic Degradation of MTRu1Pprocessed_above75genomecompletion.csv", row.names = 1)
names(MSP1)[1:7] <- paste0("MSP1:", names(MSP1)[1:7])
names(MSP2)[1:6] <- paste0("MSP2:", names(MSP2)[1:6])
names(MSFE)[1:2] <- paste0("MSFE:", names(MSFE)[1:2])
names(MSIS)[1:3] <- paste0("MSIS:", names(MSIS)[1:3])
names(MSAD)[1] <- paste0("MSAD:", names(MSAD)[1])
# bind the steps together
MS <- merge(MSP1[1:7],MSP2[1:6],by=0)
MS <- merge(MS,MSFE[1:2], by.x = "Row.names", by.y=0)
MS <- merge(MS,MSIS[1:3], by.x = "Row.names", by.y=0)
MS <- merge(MS,MSAD[c(1,4,5,6)], by.x = "Row.names", by.y=0)
# pull the genomes of interest
# alteromonas mac: MARD_SAMN10134660_REFG_MMP10134660
MS_alt <- MS[which(MS$Row.names == "MARD_SAMN10134660_REFG_MMP10134660"),]
MS_alt <- MS_alt[,2:20]
# SAR11: GORG_SAMEA6069913_SAGS_AG893F11
MS_sar11 <- MS[which(MS$Row.names == "GORG_SAMEA6069913_SAGS_AG893F11"),]
MS_sar11 <- MS_sar11[,2:20]
# ruegaria: MARD_SAMN05444358_REFG_MMP05444358
MS_ruegaria <- MS[which(MS$Row.names == "MARD_SAMN05444358_REFG_MMP05444358"),]
MS_ruegaria <- MS_ruegaria[,2:20]
# get averages
# ref genomes
MS_refg_means <- colMeans(MS[grep("REFG", MS$Row.names), 2:20])
# all
MS_allmeans <- colMeans(MS[,2:20])
# bind the set together to make the heatmap, keeping the order as is
order_MSsteps <- names(MS[2:20])
MS_summary <- rbind(MS_allmeans,MS_refg_means,MS_alt,MS_sar11,MS_ruegaria)
MS_names <- c("allmeans","refg_means","alt","sar11","ruegaria")
MS_summary <- cbind(MS_names, MS_summary)
MS_summary <- gather(MS_summary, key = "step", value = "value", 'MSP1:Step1':'MSAD:Step1')
MS_summary$MS_names <- as.factor(MS_summary$MS_names)
MS_summary$step <- as.factor(MS_summary$step)
ggplot(MS_summary, aes(x = factor(MS_names, levels = c("allmeans","refg_means","alt","sar11","ruegaria")), y = factor(step, levels = order_MSsteps), fill = value)) + 
  geom_tile() + scale_y_discrete(limits=rev) +
  scale_fill_viridis(option = "A")
ggsave(plot = last_plot(), file = "2025Results/metsal/heatmap_v1.pdf", unit = "in", height = 5, width = 5)   

