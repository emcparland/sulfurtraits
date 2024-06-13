setwd("/pool001/demers/")
library("tidyr")
library("dplyr")
library("ggplot2")
library("ggridges")
library("gridExtra")

results <- read.table("DOSTraits/new_filtered_SA/2024_0531_sulf_out_filt.txt")
# get the KOs
sulfateKOs <- read.csv("DOSTraits/Sulfate.csv")
sulfateKOs$StepKOs <- paste0("Step", sulfateKOs$Step, ":", sulfateKOs$KO)
genomessums <- read.csv("DOSTraits/genome_summary_ELM.csv")
genomesscaffolds <- read.table("DOSTraits/genomes-scaffolds-membership.tsv")
genomesscaffolds$V1 <- paste0(genomesscaffolds$V1, "-")
logics <- read.csv("DOSTraits/SulfateLogic.csv",header = 1)

#enter for the pathway
#ps <- length(unique(sulfateKOs$Pathway.Number))
# for (p in c(1)){
#   sulfateKO <- sulfateKOs %>%
#     filter(Pathway.Number == p)
#   
#   # make these KOs OR steps the columns in a df and filter the results to those KOs
#   df <- as.data.frame(matrix(ncol=length(sulfateKO$StepKOs), nrow=0))
#   names(df) <- sulfateKO$StepKOs
#   result <- results %>% filter(KO %in% sulfateKO$KO)
#   
#   #enter for the genome
#   for (genome in 1:nrow(genomessums)){
#     print(genome)
#     # we can go through the scaffolds membership file by the second column, for which we will use the genomes summary file
#     genomescaffolds <- genomesscaffolds %>%
#       filter(V2 == genomessums$Genome[genome])
#     
#     # filter the hits in the all filtered hits file to just the one genome
#     result1 <- result %>% 
#       dplyr::filter(genome %in% genomescaffolds$V2)
#     
#     # enter for that KO
#     for (KO in 1:length(sulfateKO$KO)){
#       # filter to one column at a time in the already filtered df
#       # add value to that row/column combination
#       if (sulfateKO$KO[KO] %in% result1$KO){
#         df[genome,KO]<-1
#       } else {
#         df[genome,KO]<-0
#       }
#     }
#     
#     # populate row with that genome name
#     row.names(df)[genome] <- genomessums$Genome[genome]
#   }
#   write.csv(x = df, file = paste0("DOSTraits/new_filtered_SA/", unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), "Raw.csv"))
# }

df <- read.csv("/pool001/demers/DOSTraits/new_filtered_SA/new_filtered_SA_Sulfate_AssimilationRaw.csv", row.names = 1)
df <- df[match(genomessums$Genome, row.names(df)),]

df$Group <- genomessums$taxa_elm
df$Completion <- genomessums$Mean.Completeness
p <- 1
sulfateKO <- sulfateKOs %>%
  filter(Pathway.Number == p)
# want to open the above file and see the max of each step
steps <- na.omit(unique(sulfateKO$Step))
dfsteps <- as.data.frame(matrix(ncol=length(steps), nrow =nrow(df)))
#loop through each step
for (step in steps){
  dfstep <- as.data.frame(df[,grep(pattern = paste0("Step", step,"."), x = names(df))])
  if (ncol(dfstep) >1){
    dfstep[, "max"] <- apply(dfstep[,], 1, max)
  } else {
    dfstep$max <- dfstep[,1]
  }
  dfsteps[,which(steps == step)]<- dfstep$max
  names(dfsteps)[which(steps == step)]<- paste0("Step",step)
}

row.names(dfsteps)<- row.names(df)
# then we will check the logic using the logic file and get the max of that
logic <- logics %>%
  filter(Pathway.Number == p)
for (r in unique(logic$Subpathway)){
  subpathway <- paste0("Step", strsplit(x = logic[r,4], split = ",")[[1]])
  dfsteptotal <- dfsteps %>%
    select(subpathway)
  if (ncol(dfsteptotal) >1){
    dfsteptotal$total <- apply(dfsteptotal[,], 1, sum)
  } else {
    dfsteptotal$total <- dfsteptotal[,1]
  }
  dfsteptotal$fraction <- dfsteptotal$total/length(subpathway)
  dfsteps <- cbind(dfsteps, dfsteptotal$fraction)
  names(dfsteps)[which(names(dfsteps)=="dfsteptotal$fraction")] <- paste0("Subpathway", r)
}
if (length(grep(x = names(dfsteps),pattern="Subpathway")) >1 ){
  dfsteps$MaxCompletion <- apply(dfsteps[,grep(x = names(dfsteps),pattern="Subpathway")],1,max)
} else {
  dfsteps$MaxCompletion <- dfsteps[,grep(x = names(dfsteps),pattern="Subpathway")]
}

# split into the different groups and then plot as box plots per group
dfsteps$Group <- genomessums$taxa_elm
dfsteps$Completion <- genomessums$Mean.Completeness
dfsteps$Group <- as.factor(dfsteps$Group)

#get the types of genomes here
dfsteps$Type <- "MAG"
dfsteps$Type[grep("_REFG_", row.names(dfsteps))] <- "RefGenome"
dfsteps$Type[grep("_SAGS_", row.names(dfsteps))] <- "SAG"

#save this file
write.csv(dfsteps, file = paste0("DOSTraits/new_filtered_SA/", unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), "processed.csv"), quote = F)

dfsteps50 <- dfsteps %>%
  filter(Completion >= 50)
dfsteps75 <- dfsteps %>%
  filter(Completion >= 75)
write.csv(dfsteps50, file = paste0("DOSTraits/new_filtered_SA/", unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), "processed_above50genomecompletion.csv"), quote = F)
write.csv(dfsteps75, file = paste0("DOSTraits/new_filtered_SA/", unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), "processed_above75genomecompletion.csv"), quote = F)

# plot and save
plot <- ggplot(dfsteps, aes(x=Group, y=MaxCompletion)) + geom_point(position = "jitter", aes(colour = Completion, shape=Type), alpha=.2, size = 2) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_color_viridis_c() +
  geom_boxplot(fill="white", alpha=0) +
  ggtitle(paste0(unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]),": Genome Completion Above 25%"))
ggsave(filename = paste0("DOSTraits/new_filtered_SA/", unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), ".pdf"), plot = plot, width = 50, height = 25, limitsize = FALSE)


# do the above box plots with higher completion
plot <- ggplot(dfsteps50, aes(x=Group, y=MaxCompletion)) + geom_point(position = "jitter", aes(colour = Completion, shape=Type), alpha=.2, size = 2) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_color_viridis_c() +
  geom_boxplot(fill="white", alpha=0) +
  ggtitle(paste0(unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]),": Genome Completion Above 50%"))
ggsave(filename = paste0("DOSTraits/new_filtered_SA/", unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), "_above50genomecompletion.pdf"), plot = plot, width = 50, height = 25, limitsize = FALSE)

plot <- ggplot(dfsteps75, aes(x=Group, y=MaxCompletion)) + geom_point(position = "jitter", aes(colour = Completion, shape=Type), alpha=.2, size = 2) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_color_viridis_c() +
  geom_boxplot(fill="white", alpha=0) +
  ggtitle(paste0(unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]),": Genome Completion Above 75%"))
ggsave(filename = paste0("DOSTraits/new_filtered_SA/", unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), "_above75genomecompletion.pdf"), plot = plot, width = 50, height = 25, limitsize = FALSE)

grouplevels <- unique(s_dfsteps$Group)
p_tab <- tableGrob(logic[,2:4], rows = NULL)
# also make a heatmap with the groups, with the heatmap color being the average value per gene, per group
s_dfsteps <- dfsteps %>%
  group_by(Group) %>%
  dplyr::summarise_all("mean")
s_dfsteps <- gather( s_dfsteps[,1:(ncol(s_dfsteps)-2)], key = "key", value = "value", Step1:MaxCompletion)
plot <- ggplot(s_dfsteps, aes(x=factor(key, level=unique(s_dfsteps$key)), y = factor(Group, level=grouplevels), fill = value)) +
  geom_tile() + scale_fill_viridis_c(option = "rocket", direction = -1) +
  labs(title = paste0(unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), ": Step Summaries"),
       x = "Step or Calculation",
       y = "Group")
plot <-  ggdraw() + draw_plot(p_tab, x = 0.85, width = .15) + draw_plot(plot, width = 0.85)
ggsave(filename = paste0("DOSTraits/new_filtered_SA/", unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), "_StepsHeatmap.pdf"), plot = plot, width = 30, height = 20, limitsize = FALSE)

s_dfsteps50 <- dfsteps50 %>%
  group_by(Group) %>%
  dplyr::summarise_all("mean")
s_dfsteps50 <- gather( s_dfsteps50[,1:(ncol(s_dfsteps50)-2)], key = "key", value = "value", Step1:MaxCompletion)
plot <- ggplot(s_dfsteps50, aes(x=factor(key, level=unique(s_dfsteps50$key)), y = factor(Group, level=grouplevels), fill = value)) +
  geom_tile() + scale_fill_viridis_c(option = "rocket", direction = -1) +
  labs(title = paste0(unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), ": Step Summaries for Genomes above 50% Complete"),
       x = "Step or Calculation",
       y = "Group")
plot <-  ggdraw() + draw_plot(p_tab, x = 0.85, width = .15) + draw_plot(plot, width = 0.85)
ggsave(filename = paste0("DOSTraits/new_filtered_SA/", unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), "_StepsHeatmap_genomecompletionabove50.pdf"), plot = plot, width = 30, height = 20, limitsize = FALSE)

s_dfsteps75 <- dfsteps75 %>%
  group_by(Group) %>%
  dplyr::summarise_all("mean")
s_dfsteps75 <- gather( s_dfsteps75[,1:(ncol(s_dfsteps75)-2)], key = "key", value = "value", Step1:MaxCompletion)
plot <- ggplot(s_dfsteps75, aes(x=factor(key, level=unique(s_dfsteps75$key)), y = factor(Group, level=grouplevels), fill = value)) +
  geom_tile() + scale_fill_viridis_c(option = "rocket", direction = -1) +
  labs(title = paste0(unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), ": Step Summaries for Genomes above 75% Complete"),
       x = "Step or Calculation",
       y = "Group")
plot <-  ggdraw() + draw_plot(p_tab, x = 0.85, width = .15) + draw_plot(plot, width = 0.85)
ggsave(filename = paste0("DOSTraits/new_filtered_SA/", unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), "_StepsHeatmap_genomecompletionabove75.pdf"), plot = plot, width = 30, height = 20, limitsize = FALSE)

for (p in c(1)){
  sulfateKO <- sulfateKOs %>%
    filter(Pathway.Number == p)
  #save this file
  dfsteps <- read.csv(file = paste0("DOSTraits/new_filtered_SA/", unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), "processed.csv"), row.names = 1)
  dfsteps50 <- read.csv( file = paste0("DOSTraits/new_filtered_SA/", unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), "processed_above50genomecompletion.csv"), row.names = 1)
  dfsteps75 <- read.csv( file = paste0("DOSTraits/new_filtered_SA/", unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), "processed_above75genomecompletion.csv"), row.names = 1)
  
  length(which(dfsteps$MaxCompletion == 0))
  #19110 down to 16727
  length(which(dfsteps50$MaxCompletion == 0))
  #16437 down to 14247
  length(which(dfsteps75$MaxCompletion == 0))
  #6439 down to 5466
  
  # pull and plot a scatter plot of genome completion vs pathway completion for known auxotrophs (Sar86, SAR11, halobacteria)
  auxotrophs <- dfsteps %>%
    filter(Group == "Gamma-SAR86" | Group == "Alpha-Pelagibacter" | Group == "Archaea-Halobacteria")
  plot <- ggplot(auxotrophs, aes(x=Completion, y=MaxCompletion)) + geom_point(position = "jitter", aes(colour = Group, shape=Type), alpha=0.3, size =2) +
    ggtitle(paste0(unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]),": Pathway Completion as a function of Genome Completion for Known Auxotroph"))
  ggsave(filename = paste0("DOSTraits/new_filtered_SA/", unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), "_knownauxotrophs.pdf"), plot = plot, width = 10, height = 10, limitsize = FALSE)
  
  # pull and plot in the form of a bar plot the number of genomes in each group, with the in each bar being the type of genome
  dfsteps$Lifestyle <- "SemiAuxotroph"
  dfsteps$Lifestyle[which(dfsteps$MaxCompletion == 0)] <- "Auxotroph"
  dfsteps$Lifestyle[which(dfsteps$MaxCompletion == 1)] <- "Assimilator"
  refgenomes <- dfsteps %>%
    filter(Type == "RefGenome")
  notrefgenomes <- dfsteps %>%
    filter(Type == "SAG" | Type == "MAG")
  
  plot <- ggplot(dfsteps, aes(x=Group, fill=factor(Lifestyle, level=c("Assimilator", "SemiAuxotroph", "Auxotroph")))) + 
    geom_histogram(position = "stack", stat="count") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggtitle(paste0(unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]),": All Genomes Lifestyles"))
  ggsave(filename = paste0("DOSTraits/new_filtered_SA/", unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), "_allgenomesLifestyle.pdf"), plot = plot, width = 20, height = 20, limitsize = FALSE)
  
  plot <- ggplot(refgenomes, aes(x=Group, fill=factor(Lifestyle, level=c("Assimilator", "SemiAuxotroph", "Auxotroph")))) + 
    geom_histogram(position = "stack", stat="count") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggtitle(paste0(unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]),": Reference Genomes Lifestyles"))
  ggsave(filename = paste0("DOSTraits/new_filtered_SA/", unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), "_refgenomesLifestyle.pdf"), plot = plot, width = 20, height = 20, limitsize = FALSE)
  
  plot <- ggplot(notrefgenomes, aes(x=Group, fill=factor(Lifestyle, level=c("Assimilator", "SemiAuxotroph", "Auxotroph")))) + 
    geom_histogram(position = "stack", stat="count") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggtitle(paste0(unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]),": Not Reference Genomes Lifestyles"))
  ggsave(filename = paste0("DOSTraits/new_filtered_SA/", unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), "_notrefgenomesLifestyle.pdf"), plot = plot, width = 20, height = 20, limitsize = FALSE)
  
  ## do the heatmaps of jsut reference genomes
  refgenomes50 <- dfsteps50 %>%
    filter(Type == "RefGenome")
  refgenomes75 <- dfsteps75 %>%
    filter(Type == "RefGenome")
  
  # also make a ref heatmap with the groups, with the heatmap color being the average value per gene, per group
  s_refgenomes <- refgenomes %>%
    group_by(Group) %>%
    dplyr::summarise_all("mean")
  s_refgenomes <- gather( s_refgenomes[,1:(ncol(s_refgenomes)-2)], key = "key", value = "value", Step1:MaxCompletion)
  plot <- ggplot(s_refgenomes, aes(x=factor(key, level=unique(s_refgenomes$key)), y = factor(Group, level=grouplevels), fill = value)) +
    geom_tile() + scale_fill_viridis_c(option = "rocket", direction = -1) +
    labs(title = paste0(unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), ": Step Summaries"),
         x = "Step or Calculation",
         y = "Group")
  plot <-  ggdraw() + draw_plot(p_tab, x = 0.85, width = .15) + draw_plot(plot, width = 0.85) 
  ggsave(filename = paste0("DOSTraits/new_filtered_SA/", unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), "_RefGenomes_StepsHeatmap.pdf"), plot = plot, width = 30, height = 20, limitsize = FALSE)
  
  s_refgenomes50 <- refgenomes50 %>%
    group_by(Group) %>%
    dplyr::summarise_all("mean")
  s_refgenomes50 <- gather( s_refgenomes50[,1:(ncol(s_refgenomes50)-2)], key = "key", value = "value", Step1:MaxCompletion)
  plot <- ggplot(s_refgenomes50, aes(x=factor(key, level=unique(s_refgenomes50$key)), y = factor(Group, level=grouplevels), fill = value)) +
    geom_tile() + scale_fill_viridis_c(option = "rocket", direction = -1) +
    labs(title = paste0(unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), ": Step Summaries for Genomes above 50% Complete"),
         x = "Step or Calculation",
         y = "Group")
  plot <-  ggdraw() + draw_plot(p_tab, x = 0.85, width = .15) + draw_plot(plot, width = 0.85) 
  ggsave(filename = paste0("DOSTraits/new_filtered_SA/", unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), "_RefGenomes_StepsHeatmap_genomecompletionabove50.pdf"), plot = plot, width = 30, height = 20, limitsize = FALSE)
  
  s_refgenomes75 <- refgenomes75 %>%
    group_by(Group) %>%
    dplyr::summarise_all("mean")
  s_refgenomes75 <- gather( s_refgenomes75[,1:(ncol(s_refgenomes75)-2)], key = "key", value = "value", Step1:MaxCompletion)
  plot <- ggplot(s_refgenomes75, aes(x=factor(key, level=unique(s_refgenomes75$key)), y = factor(Group, level=grouplevels), fill = value)) +
    geom_tile() + scale_fill_viridis_c(option = "rocket", direction = -1) +
    labs(title = paste0(unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), ": Step Summaries for Genomes above 75% Complete"),
         x = "Step or Calculation",
         y = "Group")
  plot <-  ggdraw() + draw_plot(p_tab, x = 0.85, width = .15) + draw_plot(plot, width = 0.85) 
  ggsave(filename = paste0("DOSTraits/new_filtered_SA/", unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), "_RefGenomes_StepsHeatmap_genomecompletionabove75.pdf"), plot = plot, width = 30, height = 20, limitsize = FALSE)
  
  # here we want to plot the combinations
  # import the logic for the combinations here, which specify the steps present, then not, then the missing steps in each combination
  # STEPS_NOT_MISSINGSTEPS
  combinations <- read.csv("DOSTraits/SulfateCombinations.csv")
  combinationsprocessed <- as.data.frame(matrix(ncol=length(combinations$Combinations), nrow = nrow(dfsteps)))
  row.names(combinationsprocessed) <- row.names(dfsteps)
  names(combinationsprocessed) <- combinations$Combinations
  combinationsprocessed[,] <- 0
  # Going to use dfsteps, loop through the columns so we can sum by the steps
  for (col in 1:ncol(combinationsprocessed)){
    combo <- names(combinationsprocessed)[col]
    combo <- strsplit(combo, split = "NOT")
    need <- combo[[1]][1]
    miss <- combo[[1]][2]
    if (miss == " "){
      miss <- ""
      misscols <- as.data.frame(matrix(nrow=nrow(dfsteps), ncol=1))
      row.names(misscols) <- row.names(dfsteps)
      names(misscols) <- "total"
      misscols$total <- 0
    } else {
      misscols <- dfsteps %>%
        select(paste0("Step",strsplit(miss, split=",")[[1]]))
      if (ncol(misscols) >1){
        misscols$total <- apply(misscols[,], 1, sum)
      } else if (ncol(misscols) == 1){
        misscols$total <- misscols[,1]
      }
    }
    
    if (need == " "){
      need <- ""
      needcols <- as.data.frame(matrix(nrow=nrow(dfsteps), ncol=1))
      row.names(needcols) <- row.names(dfsteps)
      names(needcols) <- "total"
      needcols$total <- 1
    } else {
      needcols <- dfsteps %>%
        select(paste0("Step",strsplit(need, split=",")[[1]]))
      if (ncol(needcols) >1){
        needcols$total <- apply(needcols[,], 1, sum)
        needcols$total <- (needcols$total)/length(strsplit(need, split=",")[[1]])
      } else if (ncol(needcols) == 1){
        needcols$total <- needcols[,1]
      } 
    }
    combinationsprocessed[which(misscols$total == 0 & needcols$total == 1),col] <- 1
  }
  combinationsprocessed$Group <- genomessums$taxa_elm
  s_combinationsprocessed <- combinationsprocessed %>%
    group_by(Group) %>%
    dplyr::summarise_all("mean")
  s_combinationsprocessed <- gather( s_combinationsprocessed, key = "key", value = "value", `2,5,6,8NOT `:`9NOT1,4,6`)
  plot <- ggplot(s_combinationsprocessed, aes(x=factor(key, level=unique(s_combinationsprocessed$key)), y = Group, fill = value)) +
    geom_tile() + scale_fill_viridis_c(option = "rocket", direction = -1) +
    labs(title = paste0(unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), ": Step Summaries for Genomes above 75% Complete"),
         x = "Step or Calculation",
         y = "Group")
  plot <-  ggdraw() + draw_plot(p_tab, x = 0.85, width = .15) + draw_plot(plot, width = 0.85) 
  ggsave(filename = paste0("DOSTraits/new_filtered_SA/", unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), "_combinations.pdf"), plot = plot, width = 80, height = 20, limitsize = FALSE)
  write.csv(file = paste0("DOSTraits/new_filtered_SA/", unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), "_combinations.csv"), x= combinationsprocessed)
  
  ## do relative bar plot here, for all completion, all ref genomes, 75% completion and 50% completion
  plot <- ggplot(dfsteps, aes(x=Group, fill=factor(Lifestyle, level=c("Assimilator", "SemiAuxotroph", "Auxotroph")))) + 
    geom_histogram(position = "fill", stat="count") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggtitle(paste0(unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]),": All Genomes Lifestyles"))
  ggsave(filename = paste0("DOSTraits/new_filtered_SA/", unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), "_Lifestyle.pdf"), plot = plot, width = 15, height = 10, limitsize = FALSE)
  
  
  plot <- ggplot(refgenomes, aes(x=Group, fill=factor(Lifestyle, level=c("Assimilator", "SemiAuxotroph", "Auxotroph")))) + 
    geom_histogram(position = "fill", stat="count") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggtitle(paste0(unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]),": All Genomes Lifestyles"))
  ggsave(filename = paste0("DOSTraits/new_filtered_SA/", unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), "_refgenomesLifestyle.pdf"), plot = plot, width = 15, height = 10, limitsize = FALSE)
  
  dfsteps50$Lifestyle <- "SemiAuxotroph"
  dfsteps50$Lifestyle[which(dfsteps50$MaxCompletion == 0)] <- "Auxotroph"
  dfsteps50$Lifestyle[which(dfsteps50$MaxCompletion == 1)] <- "Assimilator"
  plot <- ggplot(dfsteps50, aes(x=Group, fill=factor(Lifestyle, level=c("Assimilator", "SemiAuxotroph", "Auxotroph")))) + 
    geom_histogram(position = "fill", stat="count") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggtitle(paste0(unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]),": All Genomes (above 50% complete) Lifestyles"))
  ggsave(filename = paste0("DOSTraits/new_filtered_SA/", unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), "_Lifestyle_above50.pdf"), plot = plot, width = 15, height = 10, limitsize = FALSE)
  
  dfsteps75$Lifestyle <- "SemiAuxotroph"
  dfsteps75$Lifestyle[which(dfsteps75$MaxCompletion == 0)] <- "Auxotroph"
  dfsteps75$Lifestyle[which(dfsteps75$MaxCompletion ==1)] <- "Assimilator"
  plot <- ggplot(dfsteps75, aes(x=Group, fill=factor(Lifestyle, level=c("Assimilator", "SemiAuxotroph", "Auxotroph")))) + 
    geom_histogram(position = "fill", stat="count") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggtitle(paste0(unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]),": All Genomes (above 75% complete) Lifestyles"))
  ggsave(filename = paste0("DOSTraits/new_filtered_SA/", unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), "_Lifestyle_above75.pdf"), plot = plot, width = 15, height = 10, limitsize = FALSE)
  
  # we removed the genus work that was here as instead of using the max from each genus as the representative for that genus and then averaging over the groups, we are just using a different hmm score filter so we'll still use each genome
}
length(results$scaffold_gene[which(results$KO == "K13811")])

dupes <- results %>%
  filter(scaffold_gene %in% (results$scaffold_gene[which(results$KO == "K13811")])) %>%
  group_by(scaffold_gene) %>%
  dplyr::summarise(count = n()) %>%
  group_by(count)
which(dupes$scaffold_gene == 4)
unique(dupes$KO)

plot <- ggplot(dfsteps, aes(x=Group, fill=factor(Lifestyle, level=c("Assimilator", "SemiAuxotroph", "Auxotroph")))) + 
  geom_histogram(position = "stack", stat="count") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle(paste0(unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]),": All Genomes Lifestyles")) +
  ylim(0,5500)
ggsave(filename = paste0("DOSTraits/new_filtered_SA/", unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), "_stack_Lifestyle_allgenomes.pdf"), plot = plot, width = 15, height = 10, limitsize = FALSE)

plot <- ggplot(dfsteps50, aes(x=Group, fill=factor(Lifestyle, level=c("Assimilator", "SemiAuxotroph", "Auxotroph")))) + 
  geom_histogram(position = "stack", stat="count") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle(paste0(unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]),": All Genomes (above 50% complete) Lifestyles")) +
  ylim(0,5500)
ggsave(filename = paste0("DOSTraits/new_filtered_SA/", unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), "_stack_Lifestyle_above50.pdf"), plot = plot, width = 15, height = 10, limitsize = FALSE)

plot <- ggplot(dfsteps75, aes(x=Group, fill=factor(Lifestyle, level=c("Assimilator", "SemiAuxotroph", "Auxotroph")))) + 
  geom_histogram(position = "stack", stat="count") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle(paste0(unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]),": All Genomes (above 75% complete) Lifestyles")) +
  ylim(0,5500)
ggsave(filename = paste0("DOSTraits/new_filtered_SA/", unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), "_stack_Lifestyle_above75.pdf"), plot = plot, width = 15, height = 10, limitsize = FALSE)


