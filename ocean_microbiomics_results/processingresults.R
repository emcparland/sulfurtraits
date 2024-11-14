setwd("/pool001/demers/")
library("tidyr")
library("dplyr")
library("ggplot2")
library("ggridges")

# results <- read.table("DOSTraits/all_filtered_results.txt")
# results <- results %>%
#   mutate(V2 = gsub(V2, pattern = "gene_.*", replacement = ""))
#write.table(x=results, "DOSTraits/formatted_all_filtered_results.txt", row.names = F, quote = F, sep = "\t", col.names = F)


# results <- read.table(file = "DOSTraits/formatted_all_filtered_results.txt", sep ="\t", header = F, quote = "")
# get the KOs
sulfateKOs <- read.csv("DOSTraits/Sulfate.csv")
sulfateKOs$StepKOs <- paste0("Step", sulfateKOs$Step, ":", sulfateKOs$KO)
genomessums <- read.csv("DOSTraits/genome_summary_ELM.csv")
genomesscaffolds <- read.table("DOSTraits/genomes-scaffolds-membership.tsv")
genomesscaffolds$V1 <- paste0(genomesscaffolds$V1, "-")

# assign the genome ids to the results file
# results$V8 <- genomesscaffolds$V2[match(results$V2, genomesscaffolds$V1)]
# write.table(x = results, file = "DOSTraits/formatted_all_filtered_results_withgenomes.txt", sep = "\t", quote = F, col.names = T)

# results <- read.table("DOSTraits/kofam_results/formatted_all_filtered_results_withgenomes.txt", sep = "\t", header = T, quote = "")
# 
# # add new KOs:
# 
# results <- read.table("DOSTraits/all_filtered_results_addKO.txt")
# results <- results %>%
#   mutate(V2 = gsub(V2, pattern = "gene_.*", replacement = ""))
# write.table(x=results, "DOSTraits/formatted_all_filtered_results_addKO.txt", row.names = F, quote = F, sep = "\t", col.names = F)
# 
# 
# results <- read.table(file = "DOSTraits/formatted_all_filtered_results_addKO.txt", sep ="\t", header = F, quote = "")
# # get the KOs
# genomesscaffolds <- read.table("DOSTraits/genomes-scaffolds-membership.tsv")
# genomesscaffolds$V1 <- paste0(genomesscaffolds$V1, "-")
# 
# # assign the genome ids to the results file
# results$V8 <- genomesscaffolds$V2[match(results$V2, genomesscaffolds$V1)]
# write.table(x = results, file = "DOSTraits/kofam_results/formatted_all_filtered_results_addKO_withgenomes.txt", sep = "\t", quote = F, col.names = T)
# 
# # add hmms:
# 
# results <- read.table("DOSTraits/allcustom_filtered.txt")
# results <- results %>%
#   mutate(V2 = gsub(V1, pattern = "gene_.*", replacement = ""))
# write.table(x=results, "DOSTraits/formatted_allcustom_filtered.txt", row.names = F, quote = F, sep = "\t", col.names = F)
# 
# 
# results <- read.table(file = "DOSTraits/formatted_allcustom_filtered.txt", sep ="\t", header = F, quote = "")
# # get the KOs
# genomesscaffolds <- read.table("DOSTraits/genomes-scaffolds-membership.tsv")
# genomesscaffolds$V1 <- paste0(genomesscaffolds$V1, "-")
# 
# # assign the genome ids to the results file
# results$V20 <- genomesscaffolds$V2[match(results$V2, genomesscaffolds$V1)]
# write.table(x = results, file = "DOSTraits/kofam_results/formatted_allcustom_filtered_withgenomes.txt", sep = "\t", quote = F, col.names = T)
# 
# #concatenate results files
# results <- read.table("DOSTraits/kofam_results/formatted_all_filtered_results_withgenomes.txt", sep = "\t", header = T, quote = "")
# results1 <- read.table("DOSTraits/kofam_results/formatted_all_filtered_results_addKO_withgenomes.txt", sep = "\t", header = T, quote = "")
# results2 <- read.table("DOSTraits/kofam_results/formatted_allcustom_filtered_withgenomes.txt", sep = "\t", header = T, quote = "")
# # from these three files, we really just want the genome and the KO/hmm name
# names(results)[c(8,3)] <- c("Genome", "GeneName")
# names(results1)[c(8,3)] <- c("Genome", "GeneName")
# names(results2)[c(20,3)] <- c("Genome", "GeneName")
# results3 <- rbind(results[,c(8,3)], results1[,c(8,3)], results2[,c(20,3)])
# #save
# write.table(x = results3, file = "DOSTraits/kofam_results/formatted_allcustom_and_KOs_filtered_withgenomes.txt", sep = "\t", quote = F, col.names = T)

# # add new KOs from Oct 2, 2024:
# 
# results <- read.table("DOSTraits/all_filtered_results_addKO_Oct02.txt")
# results <- results %>%
#   mutate(V2 = gsub(V2, pattern = "gene_.*", replacement = ""))
# write.table(x=results, "DOSTraits/formatted_all_filtered_results_addKO_Oct02.txt", row.names = F, quote = F, sep = "\t", col.names = F)
# 
# 
# results <- read.table(file = "DOSTraits/formatted_all_filtered_results_addKO_Oct02.txt", sep ="\t", header = F, quote = "")
# # get the KOs
# genomesscaffolds <- read.table("DOSTraits/genomes-scaffolds-membership.tsv")
# genomesscaffolds$V1 <- paste0(genomesscaffolds$V1, "-")
# 
# # assign the genome ids to the results file
# results$V8 <- genomesscaffolds$V2[match(results$V2, genomesscaffolds$V1)]
# write.table(x = results, file = "DOSTraits/kofam_results/formatted_all_filtered_results_addKO_Oct02_withgenomes.txt", sep = "\t", quote = F, col.names = T)
# 
# #concatenate results files
# results <- read.table("DOSTraits/kofam_results/formatted_allcustom_and_KOs_filtered_withgenomes.txt", sep = "\t", header = T, quote = "")
# results1 <- read.table("DOSTraits/kofam_results/formatted_all_filtered_results_addKO_Oct02_withgenomes.txt", sep = "\t", header = T, quote = "")
# # from these three files, we really just want the genome and the KO/hmm name
# 
# names(results1)[c(8,3)] <- c("Genome", "GeneName")
# results3 <- rbind(results[,c(1,2)], results1[,c(8,3)])
# #save
# write.table(x = results3, file = "DOSTraits/kofam_results/formatted_allcustom_and_KOs_filtered_withgenomes_Oct02.txt", sep = "\t", quote = F, col.names = T)

results <- read.table("DOSTraits/kofam_results/formatted_allcustom_and_KOs_filtered_withgenomes_Oct02.txt", sep = "\t", header = T, quote = "")

#enter for the pathways
#ps <- c(1:length(na.omit(unique(sulfateKOs$Pathway.Number))))
for (p in c(8,9,10,11,12)){
  sulfateKO <- sulfateKOs %>%
    filter(Pathway.Number == p)

  # make these KOs OR steps the columns in a df and filter the results to those KOs
  df <- as.data.frame(matrix(ncol=length(sulfateKO$StepKOs), nrow=0))
  names(df) <- sulfateKO$StepKOs
  result <- results %>% filter(GeneName %in% sulfateKO$KO)

  #enter for the genome
  for (g in 1:nrow(genomessums)){
    # filter the hits in the all filtered hits file to just the one genome
    result1 <- result %>% filter(Genome == genomessums$Genome[g])

    # enter for that KO
    for (k in 1:length(sulfateKO$KO)){
      # filter to one column at a time in the already filtered df
      # add value to that row/column combination
      if (sulfateKO$KO[k] %in% result1$GeneName){
        df[g,k]<-1
      } else {
        df[g,k]<-0
      }
    }

    # populate row with that genome name
    row.names(df)[g] <- genomessums$Genome[g]
  }
  write.csv(x = df, file = paste0("DOSTraits/kofam_results/", unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), "Raw.csv"))

  df$Group <- genomessums$taxa_elm
  df$Completion <- genomessums$Mean.Completeness
}


#visualize
#where necessary, uncomment lines to make additional figures

logics <- read.csv("DOSTraits/SulfateLogic.csv",header = 1)
#1:length(na.omit(unique(sulfateKOs$Pathway.Number)))
for (p in c(1:length(na.omit(unique(sulfateKOs$Pathway.Number))))){
  sulfateKO <- sulfateKOs %>%
    filter(Pathway.Number == p)
  # want to open the above file and see the max of each step
  df <- read.csv(paste0("DOSTraits/kofam_results/",unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]),"Raw.csv"), row.names = 1)
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
  # write.csv(dfsteps, file = paste0("DOSTraits/kofam_results/", unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), "processed.csv"), quote = F)
  # 
  # dfsteps50 <- dfsteps %>%
  #   filter(Completion >= 50)
  dfsteps75 <- dfsteps %>%
    filter(Completion >= 75)
  # write.csv(dfsteps50, file = paste0("DOSTraits/kofam_results/", unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), "processed_above50genomecompletion.csv"), quote = F)
  write.csv(dfsteps75, file = paste0("DOSTraits/kofam_results/", unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), "processed_above75genomecompletion.csv"), quote = F)

  # plot and save
  # plot <- ggplot(dfsteps, aes(x=Group, y=MaxCompletion)) + geom_point(position = "jitter", aes(colour = Completion, shape=Type), alpha=.2, size = 2) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  #   scale_color_viridis_c() +
  #   geom_boxplot(fill="white", alpha=0) +
  #   ggtitle(paste0(unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]),": Genome Completion Above 25%"))
  # ggsave(filename = paste0("DOSTraits/kofam_results/", unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), ".pdf"), plot = plot, width = 50, height = 25, limitsize = FALSE)
  
  #some QC below
  # dfstepstemp <- dfsteps %>%
  #   filter(Group == "Pelagibacterales")
  # plot <- ggplot(dfstepstemp, aes(x=MaxCompletion)) +
  #   geom_histogram(bins = 100) +
  #   ggtitle(paste0(unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]),": Pelagibacterales"))
  # ggsave(filename = paste0("DOSTraits/kofam_results/", unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), "Pelagibacterales.pdf"), plot = plot, width = 15, height = 15, limitsize = FALSE)
  # 
  # dfstepstemp <- dfsteps %>%
  #   filter(Group == "Alphaproteobacteria")
  # plot <- ggplot(dfstepstemp, aes(x=MaxCompletion)) +
  #   geom_histogram(bins = 100) +
  #   ggtitle(paste0(unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]),": Alphaproteobacteria"))
  # ggsave(filename = paste0("DOSTraits/kofam_results/", unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), "Alphaproteobacteria.pdf"), plot = plot, width = 15, height = 15, limitsize = FALSE)
  # 
  # dfstepstemp <- dfsteps %>%
  #   filter(Group == "Bacilli")
  # plot <- ggplot(dfstepstemp, aes(x=MaxCompletion)) +
  #   geom_histogram(bins = 100) +
  #   ggtitle(paste0(unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]),": Bacilli"))
  # ggsave(filename = paste0("DOSTraits/kofam_results/", unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), "Bacilli.pdf"), plot = plot, width = 15, height = 15, limitsize = FALSE)
  # 
  # dfstepstemp <- dfsteps %>%
  #   filter(Group == "Cyanobacteriia")
  # plot <- ggplot(dfstepstemp, aes(x=MaxCompletion)) +
  #   geom_histogram(bins = 100) +
  #   ggtitle(paste0(unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]),": Cyanobacteriia"))
  # ggsave(filename = paste0("DOSTraits/kofam_results/", unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), "Cyanobacteriia.pdf"), plot = plot, width = 15, height = 15, limitsize = FALSE)
  # 
  # dfstepstemp <- dfsteps75 %>%
  #   filter(Group == "Pelagibacterales")
  # plot <- ggplot(dfstepstemp, aes(x=MaxCompletion)) +
  #   geom_histogram(bins = 100) +
  #   ggtitle(paste0(unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]),": Pelagibacterales Above 75% Complete"))
  # ggsave(filename = paste0("DOSTraits/kofam_results/", unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), "Pelagibacteralesabove75.pdf"), plot = plot, width = 15, height = 15, limitsize = FALSE)
  # 
  # dfstepstemp <- dfsteps75 %>%
  #   filter(Group == "Alphaproteobacteria")
  # plot <- ggplot(dfstepstemp, aes(x=MaxCompletion)) +
  #   geom_histogram(bins = 100) +
  #   ggtitle(paste0(unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]),": Alphaproteobacteria Above 75% Complete"))
  # ggsave(filename = paste0("DOSTraits/kofam_results/", unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), "Alphaproteobacteriaabove75.pdf"), plot = plot, width = 15, height = 15, limitsize = FALSE)
  # 
  # dfstepstemp <- dfsteps75 %>%
  #   filter(Group == "Bacilli")
  # plot <- ggplot(dfstepstemp, aes(x=MaxCompletion)) +
  #   geom_histogram(bins = 100) +
  #   ggtitle(paste0(unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]),": Bacilli Above 75% Complete"))
  # ggsave(filename = paste0("DOSTraits/kofam_results/", unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), "Bacilliabove75.pdf"), plot = plot, width = 15, height = 15, limitsize = FALSE)
  # 
  # dfstepstemp <- dfsteps75 %>%
  #   filter(Group == "Cyanobacteriia")
  # plot <- ggplot(dfstepstemp, aes(x=MaxCompletion)) +
  #   geom_histogram(bins = 100) +
  #   ggtitle(paste0(unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]),": Cyanobacteriia Above 75% Complete"))
  # ggsave(filename = paste0("DOSTraits/kofam_results/", unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), "Cyanobacteriiaabove75.pdf"), plot = plot, width = 15, height = 15, limitsize = FALSE)
  # 
  # #plot pathway completion as a function of genome completion
  # plot <- ggplot(dfsteps, aes(x=Completion, y=MaxCompletion)) + geom_point(position = "jitter", aes(shape=Type, colour=Type), alpha=.2, size = 1) +
  #   ggtitle(paste0(unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]),": Max Completion as a function of Genome Completion"))
  # ggsave(filename = paste0("DOSTraits/kofam_results/", unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), "_genomeCompletion.pdf"), plot = plot, width = 20, height = 20, limitsize = FALSE)
  # 
  # #further plot pathway completion as a function of genome completion
  # plot <- ggplot(dfsteps, aes(x=Completion, y=MaxCompletion, group = MaxCompletion)) +
  #   geom_density_ridges(aes(point_color = Type, point_fill = Type, point_shape = Type), point_alpha = 1, jittered_points = TRUE) +
  #   ggtitle(paste0(unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]),": Max Completion as a function of Genome Completion")) +
  #   theme_ridges()
  # ggsave(filename = paste0("DOSTraits/kofam_results/", unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), "_genomeCompletionHistograms.pdf"), plot = plot, width = 30, height = 30, limitsize = FALSE)
  # QC Completed above
  
  # do the above box plots with higher completion
  # plot <- ggplot(dfsteps50, aes(x=Group, y=MaxCompletion)) + geom_point(position = "jitter", aes(colour = Completion, shape=Type), alpha=.2, size = 2) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  #   scale_color_viridis_c() +
  #   geom_boxplot(fill="white", alpha=0) +
  #   ggtitle(paste0(unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]),": Genome Completion Above 50%"))
  # ggsave(filename = paste0("DOSTraits/kofam_results/", unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), "_above50genomecompletion.pdf"), plot = plot, width = 50, height = 25, limitsize = FALSE)

  plot <- ggplot(dfsteps75, aes(x=Group, y=MaxCompletion)) + geom_point(position = "jitter", aes(colour = Completion, shape=Type), alpha=.2, size = 2) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_color_viridis_c() +
    geom_boxplot(fill="white", alpha=0) +
    ggtitle(paste0(unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]),": Genome Completion Above 75%"))
  ggsave(filename = paste0("DOSTraits/kofam_results/", unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), "_above75genomecompletion.pdf"), plot = plot, width = 50, height = 25, limitsize = FALSE)

  p_tab <- tableGrob(logic[,2:4], rows = NULL)
  # also make a heatmap with the groups, with the heatmap color being the average value per gene, per group
  # s_dfsteps <- dfsteps %>%
  #   group_by(Group) %>%
  #   dplyr::summarise_all("mean")
  # s_dfsteps <- gather( s_dfsteps[,1:(ncol(s_dfsteps)-2)], key = "key", value = "value", Step1:MaxCompletion)
  # plot <- ggplot(s_dfsteps, aes(x=factor(key, level=unique(s_dfsteps$key)), y = Group, fill = value)) +
  #   geom_tile() + scale_fill_viridis_c(option = "rocket", direction = -1) +
  #   labs(title = paste0(unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), ": Step Summaries"),
  #        x = "Step or Calculation",
  #        y = "Group")
  # plot <-  ggdraw() + draw_plot(p_tab, x = 0.85, width = .15) + draw_plot(plot, width = 0.85) 
  # ggsave(filename = paste0("DOSTraits/kofam_results/", unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), "_StepsHeatmap.pdf"), plot = plot, width = 30, height = 20, limitsize = FALSE)
  # 
  # s_dfsteps50 <- dfsteps50 %>%
  #   group_by(Group) %>%
  #   dplyr::summarise_all("mean")
  # s_dfsteps50 <- gather( s_dfsteps50[,1:(ncol(s_dfsteps50)-2)], key = "key", value = "value", Step1:MaxCompletion)
  # plot <- ggplot(s_dfsteps50, aes(x=factor(key, level=unique(s_dfsteps50$key)), y = Group, fill = value)) +
  #   geom_tile() + scale_fill_viridis_c(option = "rocket", direction = -1) +
  #   labs(title = paste0(unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), ": Step Summaries for Genomes above 50% Complete"),
  #        x = "Step or Calculation",
  #        y = "Group")
  # plot <-  ggdraw() + draw_plot(p_tab, x = 0.85, width = .15) + draw_plot(plot, width = 0.85) 
  # ggsave(filename = paste0("DOSTraits/kofam_results/", unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), "_StepsHeatmap_genomecompletionabove50.pdf"), plot = plot, width = 30, height = 20, limitsize = FALSE)

  s_dfsteps75 <- dfsteps75 %>%
    group_by(Group) %>%
    dplyr::summarise_all("mean")
  s_dfsteps75 <- gather( s_dfsteps75[,1:(ncol(s_dfsteps75)-2)], key = "key", value = "value", Step1:MaxCompletion)
  plot <- ggplot(s_dfsteps75, aes(x=factor(key, level=unique(s_dfsteps75$key)), y = Group, fill = value)) +
    geom_tile() + scale_fill_viridis_c(option = "rocket", direction = -1, limits = c(0,1)) +
    labs(title = paste0(unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), ": Step Summaries for Genomes above 75% Complete"),
         x = "Step or Calculation",
         y = "Group")
  plot <-  ggdraw() + draw_plot(p_tab, x = 0.85, width = .15) + draw_plot(plot, width = 0.85) 
  ggsave(filename = paste0("DOSTraits/kofam_results/", unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), "_StepsHeatmap_genomecompletionabove75.pdf"), plot = plot, width = 30, height = 20, limitsize = FALSE)

}

#write the empty max table
maxtable <- data.frame(matrix(nrow=nrow(read.csv(paste0("DOSTraits/kofam_results/",unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == 2)]),"processed_above75genomecompletion.csv"), row.names = 1)),ncol=(length(na.omit(unique(sulfateKOs$Pathway.Number)))-1)))

for (p in 2:length(na.omit(unique(sulfateKOs$Pathway.Number)))){
  sulfateKO <- sulfateKOs %>%
    filter(Pathway.Number == p)
  # want to open the above file and see the max of each step
  df <- read.csv(paste0("DOSTraits/kofam_results/",unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]),"processed_above75genomecompletion.csv"), row.names = 1)
  maxtable[,p-1]<- df$MaxCompletion
  names(maxtable)[p-1]<- unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)])
}
 
row.names(maxtable)<- row.names(df)

write.csv(x = maxtable, file = "DOSTraits/kofam_results/summary_max_table.csv", row.names = T)
groups <- data.frame(df$Group)
row.names(groups)<- row.names(df)
groups$df.Group <- as.factor(groups$df.Group)
names(groups) <- c("Group")

library(pheatmap)

save_pheatmap_png <- function(x, filename, width=12, height=12) {
  pdf(filename, width = width, height = height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
plot <- pheatmap(maxtable, cluster_rows = T, show_rownames =F)
save_pheatmap_png(plot, "DOSTraits/kofam_results/summary_heatmap_1.pdf")

my_colour  = c(`Acidobacteriota` = "#c19aff", `Actinobacteriota` = "#ff9ac2", `Alpha−Other` = "#db6ad9", `Alpha−Pelagibacter` = "#7e2799", `Alpha−Puniceispirilla` = "#ff70d1", `Alpha−Rhodobacter` ="#df8eff",
                              `Alpha−Rhodobacter−LowGC` = "380c60", `Archaea−Halobacteria` = "#54db76", `Archaea−Nitrosopumilaceae` = "#01c887", `Archaea−Other` = "#004f02", `Archaea−Poseidoniia MGIIa` = "#008b4d", 
                              `Archaea−Poseidoniia MGIIb` = "#3c5c00", `Bacteria−Other` = "#757500", `Bacteroidota`="#e89126", `Chloroflexota`="#c8de7a", `Cyano−Other`="#01eacc", `Cyano−Pro`="#0170d5", `Cyano−Syn`="#93b6ff", `Dadabacteria`="#d25523", 
                              `Firmicutes` = "#74508f", `Gamma−Alteromonas`="#ee3f74", `Gamma−Other`="#870013", `Gamma−SAR86`="#ff609c", `Gamma−SAR92`="#e54c44", `Gemmatimonadota`="#ffab4f", `Marinisomatota`="#8c0050", `Other`="#ffa35e", 
                              `Patescibacteria`="#610020", `Planctomycetota`="#ffa87c", `SAR324`="#754014", `Verrucomicrobiota`="#a56500")
my_colour  = list(Group = c("#c19aff",  "#ff9ac2", "#db6ad9", "#7e2799","#ff70d1","#df8eff", "380c60", "#54db76", "#01c887",  "#004f02", "#008b4d", 
                "#3c5c00", "#757500", "#e89126", "#c8de7a", "#01eacc","#0170d5", "#93b6ff", "#d25523", 
                "#74508f", "#ee3f74", "#870013", "#ff609c", "#e54c44", "#ffab4f", "#8c0050", "#ffa35e", 
               "#610020", "#ffa87c", "#754014", "#a56500"))

plot <- pheatmap(maxtable, cluster_rows = T, annotation_row = groups, show_rownames =F)
save_pheatmap_png(plot, "DOSTraits/kofam_results/summary_heatmap_2.pdf")
my_colour  = list(Group = c(Acidobacteriota = "#c19aff", Actinobacteriota = "#ff9ac2", AlphaOther = "#db6ad9", AlphaPelagibacter = "#7e2799", AlphaPuniceispirilla = "#ff70d1", AlphaRhodobacter ="#df8eff",
               AlphaRhodobacterLowGC = "#380c60", ArchaeaHalobacteria = "#54db76", ArchaeaNitrosopumilaceae = "#01c887", ArchaeaOther = "#004f02", ArchaeaPoseidoniiaMGIIa = "#008b4d", 
               ArchaeaPoseidoniiaMGIIb = "#3c5c00", BacteriaOther = "#757500", Bacteroidota="#e89126", Chloroflexota="#c8de7a", CyanoOther="#01eacc", CyanoPro="#0170d5", CyanoSyn="#93b6ff", Dadabacteria="#d25523", 
               Firmicutes = "#74508f", GammaAlteromonas="#ee3f74", GammaOther="#870013", GammaSAR86="#ff609c", GammaSAR92="#e54c44", Gemmatimonadota="#ffab4f", Marinisomatota="#8c0050", Other="#ffa35e", 
               Patescibacteria="#610020", Planctomycetota="#ffa87c", SAR324="#754014", Verrucomicrobiota="#a56500"))
groups$Group <- gsub(pattern = "-", replacement = "", x = groups$Group)
groups$Group <- gsub(pattern = " ", replacement = "", x = groups$Group)
plot <- pheatmap(maxtable, cluster_rows = T, annotation_row = groups, show_rownames =F, annotation_colors = my_colour)
save_pheatmap_png(plot, "DOSTraits/kofam_results/summary_heatmap_3.pdf")

SA_status <- read.csv("DOSTraits/new_filtered_SA/2024_June19_ERIN_SulfateAssimilationprocessed_above75genomecompletion.csv", row.names = 2)
SA_status_1 <- data.frame(SA_status$Lifestyle)
row.names(SA_status_1) <- row.names(SA_status)
SA_status_1$SA_status.Lifestyle <- as.factor(SA_status_1$SA_status.Lifestyle)
plot <- pheatmap(maxtable, cluster_rows = T, annotation_row = SA_status_1, show_rownames =F)
save_pheatmap_png(plot, "DOSTraits/kofam_results/summary_heatmap_4.pdf")

my_colour2  = list(Group = c(Acidobacteriota = "#c19aff", Actinobacteriota = "#ff9ac2", AlphaOther = "#db6ad9", AlphaPelagibacter = "#7e2799", AlphaPuniceispirilla = "#ff70d1", AlphaRhodobacter ="#df8eff",
                            AlphaRhodobacterLowGC = "#380c60", ArchaeaHalobacteria = "#54db76", ArchaeaNitrosopumilaceae = "#01c887", ArchaeaOther = "#004f02", ArchaeaPoseidoniiaMGIIa = "#008b4d", 
                            ArchaeaPoseidoniiaMGIIb = "#3c5c00", BacteriaOther = "#757500", Bacteroidota="#e89126", Chloroflexota="#c8de7a", CyanoOther="#01eacc", CyanoPro="#0170d5", CyanoSyn="#93b6ff", Dadabacteria="#d25523", 
                            Firmicutes = "#74508f", GammaAlteromonas="#ee3f74", GammaOther="#870013", GammaSAR86="#ff609c", GammaSAR92="#e54c44", Gemmatimonadota="#ffab4f", Marinisomatota="#8c0050", Other="#ffa35e", 
                            Patescibacteria="#610020", Planctomycetota="#ffa87c", SAR324="#754014", Verrucomicrobiota="#a56500"),
                  Lifestyle = c(assimilator = "black", auxotroph = "white", incomplete = "grey"))
SA_status_2 <- data.frame(SA_status_1[match(row.names(groups), row.names(SA_status_1)),])
row.names(SA_status_2) <- row.names(SA_status_1)[match(row.names(groups), row.names(SA_status_1))]
names(SA_status_2) <- c("Lifestyle")
annotations <- cbind(groups, SA_status_2)
plot <- pheatmap(maxtable, cluster_rows = T, annotation_row = annotations, show_rownames =F, annotation_colors = my_colour2)
save_pheatmap_png(plot, "DOSTraits/kofam_results/summary_heatmap_5.pdf")

# add oxidation state per pathway
oxstates <- as.data.frame(numeric(42))
row.names(oxstates) <- names(maxtable)
names(oxstates) <- c("Oxidation State")
ox <- c(-2, 6,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2)
oc <- c(ox, -2,-2,4,4,-2,4,4,4,4,4,4,4,-2,6,6,4,-2,-2,-2,-2,-2,-2,-1,0,0,4,4,4,4,-2)
oxstates$`Oxidation State` <- as.factor(oc)
my_colour3  = list(Group = c(Acidobacteriota = "#c19aff", Actinobacteriota = "#ff9ac2", AlphaOther = "#db6ad9", AlphaPelagibacter = "#7e2799", AlphaPuniceispirilla = "#ff70d1", AlphaRhodobacter ="#df8eff",
                             AlphaRhodobacterLowGC = "#380c60", ArchaeaHalobacteria = "#54db76", ArchaeaNitrosopumilaceae = "#01c887", ArchaeaOther = "#004f02", ArchaeaPoseidoniiaMGIIa = "#008b4d", 
                             ArchaeaPoseidoniiaMGIIb = "#3c5c00", BacteriaOther = "#757500", Bacteroidota="#e89126", Chloroflexota="#c8de7a", CyanoOther="#01eacc", CyanoPro="#0170d5", CyanoSyn="#93b6ff", Dadabacteria="#d25523", 
                             Firmicutes = "#74508f", GammaAlteromonas="#ee3f74", GammaOther="#870013", GammaSAR86="#ff609c", GammaSAR92="#e54c44", Gemmatimonadota="#ffab4f", Marinisomatota="#8c0050", Other="#ffa35e", 
                             Patescibacteria="#610020", Planctomycetota="#ffa87c", SAR324="#754014", Verrucomicrobiota="#a56500"),
                   Lifestyle = c(assimilator = "black", auxotroph = "white", incomplete = "grey"),
                   `Oxidation State` = c(`-2` = "darkcyan", `-1`= "cyan", `0`="white", `4`="brown1", `6`="brown4"))

plot <- pheatmap(maxtable, cluster_rows = T, annotation_row = annotations, show_rownames =F, annotation_col = oxstates, annotation_colors = my_colour3)
save_pheatmap_png(plot, "DOSTraits/kofam_results/summary_heatmap_6.pdf")

# going to eliminate some of the pathways (all transport except isethionate and N-acetyltaaurine, plus some empty ones, plus fusion enzyme)
# here we have to make a new maxtable:
maxtable2 <- maxtable[,-c(4,9,12,14,15,18,24,25,27,29,39)]
maxtable3 <- maxtable[,-c(4,7,8,9,10,11,12,14,15,18,24,25,27,29,39)]
metsaltable <- maxtable[,c(7,8,9,10,11)]
plot <- pheatmap(maxtable2, cluster_rows = T, annotation_row = annotations, show_rownames =F, annotation_col = oxstates, annotation_colors = my_colour3)
save_pheatmap_png(plot, "DOSTraits/kofam_results/summary_heatmap_7.pdf")
plot <- pheatmap(maxtable3, cluster_rows = T, annotation_row = annotations, show_rownames =F, annotation_col = oxstates, annotation_colors = my_colour3)
save_pheatmap_png(plot, "DOSTraits/kofam_results/summary_heatmap_8.pdf")
plot <- pheatmap(metsaltable, cluster_rows = T, annotation_row = annotations, show_rownames =F, annotation_col = oxstates, annotation_colors = my_colour3)
save_pheatmap_png(plot, "DOSTraits/kofam_results/summary_heatmap_9.pdf",  width=6, height=12)

# do the above figures but simiplifying any in between value to 0.5 to see if we get a better signal out of the individual genomes
mypalettelength <- 50
mycolor <- viridis::rocket(mypalettelength,direction = -1)
mybreaks <- c(seq(min(maxtable), 0.001, length.out=ceiling(mypalettelength/2)),0.5, seq(0.999, max(maxtable), length.out=floor(mypalettelength/2)))
plot <- pheatmap(maxtable, cluster_rows = T, annotation_row = annotations, show_rownames =F, annotation_col = oxstates, annotation_colors = my_colour3, breaks = mybreaks, color = mycolor)
save_pheatmap_png(plot, "DOSTraits/kofam_results/summary_heatmap_10.pdf")

plot <- pheatmap(maxtable2, cluster_rows = T, annotation_row = annotations, show_rownames =F, annotation_col = oxstates, annotation_colors = my_colour3, breaks = mybreaks, color = mycolor)
save_pheatmap_png(plot, "DOSTraits/kofam_results/summary_heatmap_11.pdf")

plot <- pheatmap(maxtable3, cluster_rows = T, annotation_row = annotations, show_rownames =F, annotation_col = oxstates, annotation_colors = my_colour3, breaks = mybreaks, color = mycolor)
save_pheatmap_png(plot, "DOSTraits/kofam_results/summary_heatmap_12.pdf",  width=12, height=12)

plot <- pheatmap(metsaltable, cluster_rows = T, annotation_row = annotations, show_rownames =F, annotation_col = oxstates, annotation_colors = my_colour3, breaks = mybreaks, color = mycolor)
save_pheatmap_png(plot, "DOSTraits/kofam_results/summary_heatmap_13.pdf",  width=6, height=12)
# going to get the averages per taxanomic group with all pathways, from maxtable
maxtablesimp <- maxtable
maxtablesimp$Groups <- groups[match( row.names(maxtable), row.names(groups)),]
mean_maxtablesimp <- maxtablesimp %>%
  group_by(Groups) %>%
  dplyr::summarise_all("mean")
maxtable2simp <- maxtable2
maxtable2simp$Groups <- groups[match( row.names(maxtable2), row.names(groups)),]
mean_maxtable2simp <- maxtable2simp %>%
  group_by(Groups) %>%
  dplyr::summarise_all("mean")
maxtable3simp <- maxtable3
maxtable3simp$Groups <- groups[match( row.names(maxtable3), row.names(groups)),]
mean_maxtable3simp <- maxtable3simp %>%
  group_by(Groups) %>%
  dplyr::summarise_all("mean")
metsaltablesimp <- metsaltable
metsaltablesimp$Groups <- groups[match( row.names(metsaltable), row.names(groups)),]
mean_metsaltablesimp <- metsaltablesimp %>%
  group_by(Groups) %>%
  dplyr::summarise_all("mean")
mean_maxtablesimp <- as.data.frame(mean_maxtablesimp)
mean_maxtable2simp <- as.data.frame(mean_maxtable2simp)
mean_maxtable3simp <- as.data.frame(mean_maxtable3simp)
mean_metsaltablesimp <- as.data.frame(mean_metsaltablesimp)
row.names(mean_maxtablesimp)<- mean_maxtablesimp$Groups
row.names(mean_maxtable2simp)<- mean_maxtable2simp$Groups
row.names(mean_maxtable3simp)<- mean_maxtable3simp$Groups
row.names(mean_metsaltablesimp)<- mean_metsaltablesimp$Groups
mean_maxtablesimp <- mean_maxtablesimp[,2:43]
mean_maxtable2simp <- mean_maxtable2simp[,2:33]
mean_maxtable3simp <- mean_maxtable3simp[,2:28]
mean_metsaltablesimp <- mean_metsaltablesimp[,2:6]

SA_status3 <- as.data.frame(cbind( SA_status$MaxCompletion, SA_status$taxa_elm))
names(SA_status3) <- c("MaxCompletion", "Group")
SA_status3$MaxCompletion <- as.numeric(SA_status3$MaxCompletion)
row.names(SA_status3) <- row.names(SA_status)
mean_SA_status <- SA_status3 %>%
  group_by(Group) %>%
  dplyr::summarise_all("mean")
groupnames <- mean_SA_status$Group
mean_SA_status <- as.data.frame(mean_SA_status$MaxCompletion)
groupnames <- gsub(x=groupnames, "-", "")
groupnames <- gsub(x=groupnames, " ", "")
row.names(mean_SA_status) <- groupnames
names(mean_SA_status) <- c("MaxCompletionAverage")
my_colour4  = list( MaxCompletionAverage = c( `0` = "white", `1` = "black"),
                   `Oxidation State` = c(`-2` = "darkcyan", `-1`= "cyan", `0`="white", `4`="brown1", `6`="brown4"))
plot <- pheatmap(mean_maxtablesimp, cluster_rows = T, annotation_row = mean_SA_status, show_rownames =T, annotation_col = oxstates, annotation_colors = my_colour4)
save_pheatmap_png(plot, "DOSTraits/kofam_results/summary_heatmap_14.pdf")

# averages per taxanomic group with only the important pathways for summary, from maxtable2, maxtable3, and metsaltable
plot <- pheatmap(mean_maxtable2simp, cluster_rows = T, annotation_row = mean_SA_status, show_rownames =T, annotation_col = oxstates, annotation_colors = my_colour4)
save_pheatmap_png(plot, "DOSTraits/kofam_results/summary_heatmap_15.pdf")
plot <- pheatmap(mean_maxtable3simp, cluster_rows = T, annotation_row = mean_SA_status, show_rownames =T, annotation_col = oxstates, annotation_colors = my_colour4)
save_pheatmap_png(plot, "DOSTraits/kofam_results/summary_heatmap_16.pdf")
plot <- pheatmap(mean_metsaltablesimp, cluster_rows = T, annotation_row = mean_SA_status, show_rownames =T, annotation_col = oxstates, annotation_colors = my_colour4)
save_pheatmap_png(plot, "DOSTraits/kofam_results/summary_heatmap_17.pdf")

plot <- pheatmap(mean_maxtablesimp, cluster_rows = T, annotation_row = mean_SA_status, show_rownames =T, annotation_col = oxstates, annotation_colors = my_colour4, color = mycolor)
save_pheatmap_png(plot, "DOSTraits/kofam_results/summary_heatmap_18.pdf")
plot <- pheatmap(mean_maxtable2simp, cluster_rows = T, annotation_row = mean_SA_status, show_rownames =T, annotation_col = oxstates, annotation_colors = my_colour4, color = mycolor)
save_pheatmap_png(plot, "DOSTraits/kofam_results/summary_heatmap_19.pdf")
plot <- pheatmap(mean_maxtable3simp, cluster_rows = T, annotation_row = mean_SA_status, show_rownames =T, annotation_col = oxstates, annotation_colors = my_colour4, color = mycolor)
save_pheatmap_png(plot, "DOSTraits/kofam_results/summary_heatmap_20.pdf")
plot <- pheatmap(mean_metsaltablesimp, cluster_rows = T, annotation_row = mean_SA_status, show_rownames =T, annotation_col = oxstates, annotation_colors = my_colour4, color = mycolor)
save_pheatmap_png(plot, "DOSTraits/kofam_results/summary_heatmap_21.pdf")

# work on the figure from the paper:


#extra figures for sulfate assimilation
# for (p in c(1)){
#   sulfateKO <- sulfateKOs %>%
#     filter(Pathway.Number == p)
#   #save this file
#   dfsteps <- read.csv(file = paste0("DOSTraits/kofam_results/", unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), "processed.csv"), row.names = 1)
#   dfsteps50 <- read.csv( file = paste0("DOSTraits/kofam_results/", unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), "processed_above50genomecompletion.csv"), row.names = 1)
#   dfsteps75 <- read.csv( file = paste0("DOSTraits/kofam_results/", unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), "processed_above75genomecompletion.csv"), row.names = 1)
#   
#   length(which(dfsteps$MaxCompletion == 0))
#   #19110
#   length(which(dfsteps50$MaxCompletion == 0))
#   #16437
#   length(which(dfsteps75$MaxCompletion == 0))
#   #6439
#   
#   # pull and plot a scatter plot of genome completion vs pathway completion for known auxotrophs (Sar86, SAR11, halobacteria)
#   auxotrophs <- dfsteps %>%
#     filter(Group == "SAR86" | Group == "Pelagibacterales" | Group == "Halobacteria")
#   plot <- ggplot(auxotrophs, aes(x=Completion, y=MaxCompletion)) + geom_point(position = "jitter", aes(colour = Group, shape=Type), alpha=0.3, size =2) +
#     ggtitle(paste0(unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]),": Pathway Completion as a function of Genome Completion for Known Auxotroph"))
#   ggsave(filename = paste0("DOSTraits/kofam_results/", unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), "_knownauxotrophs.pdf"), plot = plot, width = 10, height = 10, limitsize = FALSE)
#   
#   # pull and plot in the form of a bar plot the number of genomes in each group, with the in each bar being the type of genome
#   dfsteps$Lifestyle <- "SemiAuxotroph"
#   dfsteps$Lifestyle[which(dfsteps$MaxCompletion == 0)] <- "Auxotroph"
#   dfsteps$Lifestyle[which(dfsteps$MaxCompletion == 1)] <- "Assimilator"
#   refgenomes <- dfsteps %>%
#     filter(Type == "RefGenome")
#   notrefgenomes <- dfsteps %>%
#     filter(Type == "SAG" | Type == "MAG")
#   plot <- ggplot(dfsteps, aes(x=Group, fill=Lifestyle)) + 
#     geom_histogram(position = "stack", stat="count") +
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#     ggtitle(paste0(unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]),": All Genomes Lifestyles"))
#   ggsave(filename = paste0("DOSTraits/kofam_results/", unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), "_allgenomesLifestyle.pdf"), plot = plot, width = 20, height = 20, limitsize = FALSE)
#   
#   plot <- ggplot(refgenomes, aes(x=Group, fill=Lifestyle)) + 
#     geom_histogram(position = "stack", stat="count") +
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#     ggtitle(paste0(unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]),": Reference Genomes Lifestyles"))
#   ggsave(filename = paste0("DOSTraits/kofam_results/", unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), "_refgenomesLifestyle.pdf"), plot = plot, width = 20, height = 20, limitsize = FALSE)
#   
#   plot <- ggplot(notrefgenomes, aes(x=Group, fill=Lifestyle)) + 
#     geom_histogram(position = "stack", stat="count") +
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#     ggtitle(paste0(unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]),": Not Reference Genomes Lifestyles"))
#   ggsave(filename = paste0("DOSTraits/kofam_results/", unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), "_notrefgenomesLifestyle.pdf"), plot = plot, width = 20, height = 20, limitsize = FALSE)
#   
#   ## do the heatmaps of jsut reference genomes
#   refgenomes50 <- dfsteps50 %>%
#     filter(Type == "RefGenome")
#   refgenomes75 <- dfsteps75 %>%
#     filter(Type == "RefGenome")
#   
#   # also make a heatmap with the groups, with the heatmap color being the average value per gene, per group
#   s_refgenomes <- refgenomes %>%
#     group_by(Group) %>%
#     dplyr::summarise_all("mean")
#   s_refgenomes <- gather( s_refgenomes[,1:(ncol(s_refgenomes)-2)], key = "key", value = "value", Step1:MaxCompletion)
#   plot <- ggplot(s_refgenomes, aes(x=factor(key, level=unique(s_refgenomes$key)), y = Group, fill = value)) +
#     geom_tile() + scale_fill_viridis_c(option = "rocket", direction = -1) +
#     labs(title = paste0(unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), ": Step Summaries"),
#          x = "Step or Calculation",
#          y = "Group")
#   plot <-  ggdraw() + draw_plot(p_tab, x = 0.85, width = .15) + draw_plot(plot, width = 0.85) 
#   ggsave(filename = paste0("DOSTraits/kofam_results/", unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), "_RefGenomes_StepsHeatmap.pdf"), plot = plot, width = 30, height = 20, limitsize = FALSE)
#   
#   s_refgenomes50 <- refgenomes50 %>%
#     group_by(Group) %>%
#     dplyr::summarise_all("mean")
#   s_refgenomes50 <- gather( s_refgenomes50[,1:(ncol(s_refgenomes50)-2)], key = "key", value = "value", Step1:MaxCompletion)
#   plot <- ggplot(s_refgenomes50, aes(x=factor(key, level=unique(s_refgenomes50$key)), y = Group, fill = value)) +
#     geom_tile() + scale_fill_viridis_c(option = "rocket", direction = -1) +
#     labs(title = paste0(unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), ": Step Summaries for Genomes above 50% Complete"),
#          x = "Step or Calculation",
#          y = "Group")
#   plot <-  ggdraw() + draw_plot(p_tab, x = 0.85, width = .15) + draw_plot(plot, width = 0.85) 
#   ggsave(filename = paste0("DOSTraits/kofam_results/", unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), "_RefGenomes_StepsHeatmap_genomecompletionabove50.pdf"), plot = plot, width = 30, height = 20, limitsize = FALSE)
#   
#   s_refgenomes75 <- refgenomes75 %>%
#     group_by(Group) %>%
#     dplyr::summarise_all("mean")
#   s_refgenomes75 <- gather( s_refgenomes75[,1:(ncol(s_refgenomes75)-2)], key = "key", value = "value", Step1:MaxCompletion)
#   plot <- ggplot(s_refgenomes75, aes(x=factor(key, level=unique(s_refgenomes75$key)), y = Group, fill = value)) +
#     geom_tile() + scale_fill_viridis_c(option = "rocket", direction = -1) +
#     labs(title = paste0(unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), ": Step Summaries for Genomes above 75% Complete"),
#          x = "Step or Calculation",
#          y = "Group")
#   plot <-  ggdraw() + draw_plot(p_tab, x = 0.85, width = .15) + draw_plot(plot, width = 0.85) 
#   ggsave(filename = paste0("DOSTraits/kofam_results/", unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), "_RefGenomes_StepsHeatmap_genomecompletionabove75.pdf"), plot = plot, width = 30, height = 20, limitsize = FALSE)
#   
#   # here we want to plot the combinations
#   # import the logic for the combinations here, which specify the steps present, then not, then the missing steps in each combination
#   # STEPS_NOT_MISSINGSTEPS
#   combinations <- read.csv("DOSTraits/SulfateCombinations.csv")
#   combinationsprocessed <- as.data.frame(matrix(ncol=length(combinations$Combinations), nrow = nrow(dfsteps)))
#   row.names(combinationsprocessed) <- row.names(dfsteps)
#   names(combinationsprocessed) <- combinations$Combinations
#   combinationsprocessed[,] <- 0
#   # Going to use dfsteps, loop through the columns so we can sum by the steps
#   for (col in 1:ncol(combinationsprocessed)){
#     combo <- names(combinationsprocessed)[col]
#     combo <- strsplit(combo, split = "NOT")
#     need <- combo[[1]][1]
#     miss <- combo[[1]][2]
#     if (miss == " "){
#       miss <- ""
#       misscols <- as.data.frame(matrix(nrow=nrow(dfsteps), ncol=1))
#       row.names(misscols) <- row.names(dfsteps)
#       names(misscols) <- "total"
#       misscols$total <- 0
#     } else {
#       misscols <- dfsteps %>%
#         select(paste0("Step",strsplit(miss, split=",")[[1]]))
#       if (ncol(misscols) >1){
#         misscols$total <- apply(misscols[,], 1, sum)
#       } else if (ncol(misscols) == 1){
#         misscols$total <- misscols[,1]
#       }
#     }
#     
#     if (need == " "){
#       need <- ""
#       needcols <- as.data.frame(matrix(nrow=nrow(dfsteps), ncol=1))
#       row.names(needcols) <- row.names(dfsteps)
#       names(needcols) <- "total"
#       needcols$total <- 1
#     } else {
#       needcols <- dfsteps %>%
#         select(paste0("Step",strsplit(need, split=",")[[1]]))
#       if (ncol(needcols) >1){
#         needcols$total <- apply(needcols[,], 1, sum)
#         needcols$total <- (needcols$total)/length(strsplit(need, split=",")[[1]])
#       } else if (ncol(needcols) == 1){
#         needcols$total <- needcols[,1]
#       } 
#     }
#     combinationsprocessed[which(misscols$total == 0 & needcols$total == 1),col] <- 1
#   }
#   combinationsprocessed$Group <- genomessums$taxa_elm
#   s_combinationsprocessed <- combinationsprocessed %>%
#     group_by(Group) %>%
#     dplyr::summarise_all("mean")
#   s_combinationsprocessed <- gather( s_combinationsprocessed, key = "key", value = "value", `2,5,6,8NOT `:`9NOT3,5,6`)
#   plot <- ggplot(s_combinationsprocessed, aes(x=factor(key, level=unique(s_combinationsprocessed$key)), y = Group, fill = value)) +
#     geom_tile() + scale_fill_viridis_c(option = "rocket", direction = -1) +
#     labs(title = paste0(unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), ": Step Summaries for Genomes above 75% Complete"),
#          x = "Step or Calculation",
#          y = "Group")
#   plot <-  ggdraw() + draw_plot(p_tab, x = 0.85, width = .15) + draw_plot(plot, width = 0.85) 
#   ggsave(filename = paste0("DOSTraits/kofam_results/", unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), "_combinations.pdf"), plot = plot, width = 60, height = 20, limitsize = FALSE)
#   write.csv(file = paste0("DOSTraits/kofam_results/", unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), "_combinations.csv"), x= combinationsprocessed)
#   
#   ## do relative bar plot here, for all completion and 75%
#   plot <- ggplot(dfsteps, aes(x=Group, fill=factor(Lifestyle, level=c("Assimilator", "SemiAuxotroph", "Auxotroph")))) + 
#     geom_histogram(position = "fill", stat="count") +
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#     ggtitle(paste0(unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]),": All Genomes Lifestyles"))
#   ggsave(filename = paste0("DOSTraits/kofam_results/", unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), "_Lifestyle.pdf"), plot = plot, width = 15, height = 10, limitsize = FALSE)
#   
#   
#   plot <- ggplot(refgenomes, aes(x=Group, fill=factor(Lifestyle, level=c("Assimilator", "SemiAuxotroph", "Auxotroph")))) + 
#     geom_histogram(position = "fill", stat="count") +
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#     ggtitle(paste0(unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]),": All Genomes Lifestyles"))
#   ggsave(filename = paste0("DOSTraits/kofam_results/", unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), "_refgenomesLifestyle.pdf"), plot = plot, width = 15, height = 10, limitsize = FALSE)
#   
#   dfsteps50$Lifestyle <- "SemiAuxotroph"
#   dfsteps50$Lifestyle[which(dfsteps50$MaxCompletion == 0)] <- "Auxotroph"
#   dfsteps50$Lifestyle[which(dfsteps50$MaxCompletion == 1)] <- "Assimilator"
#   plot <- ggplot(dfsteps50, aes(x=Group, fill=factor(Lifestyle, level=c("Assimilator", "SemiAuxotroph", "Auxotroph")))) + 
#     geom_histogram(position = "fill", stat="count") +
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#     ggtitle(paste0(unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]),": All Genomes (above 50% complete) Lifestyles"))
#   ggsave(filename = paste0("DOSTraits/kofam_results/", unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), "_Lifestyle_above50.pdf"), plot = plot, width = 15, height = 10, limitsize = FALSE)
#   
#   dfsteps75$Lifestyle <- "SemiAuxotroph"
#   dfsteps75$Lifestyle[which(dfsteps75$MaxCompletion == 0)] <- "Auxotroph"
#   dfsteps75$Lifestyle[which(dfsteps75$MaxCompletion >=.75)] <- "Assimilator"
#   plot <- ggplot(dfsteps75, aes(x=Group, fill=factor(Lifestyle, level=c("Assimilator", "SemiAuxotroph", "Auxotroph")))) + 
#     geom_histogram(position = "fill", stat="count") +
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#     ggtitle(paste0(unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]),": All Genomes (above 75% complete) Lifestyles"))
#   ggsave(filename = paste0("DOSTraits/kofam_results/", unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), "_Lifestyle_above75.pdf"), plot = plot, width = 15, height = 10, limitsize = FALSE)
#   
#   
#   dfsteps$Genus <- genomessums$Genus
#   dfsteps_genus <- dfsteps[,c(1:9,21)]
#   dfsteps_genus <- dfsteps_genus %>%
#     group_by(Genus) %>%
#     dplyr::summarise_all("max")
#   # redo the logic, now just using the max per genus as our step values
#   logic <- logics %>%
#     filter(Pathway.Number == p)
#   for (r in unique(logic$Subpathway)){
#     subpathway <- paste0("Step", strsplit(x = logic[r,4], split = ",")[[1]])
#     dfsteptotal <- dfsteps_genus %>%
#       select(subpathway)
#     if (ncol(dfsteptotal) >1){
#       dfsteptotal$total <- apply(dfsteptotal[,], 1, sum)
#     } else {
#       dfsteptotal$total <- dfsteptotal[,1]
#     }
#     dfsteptotal$fraction <- dfsteptotal$total/length(subpathway)
#     dfsteps_genus <- cbind(dfsteps_genus, dfsteptotal$fraction)
#     names(dfsteps_genus)[which(names(dfsteps_genus)=="dfsteptotal$fraction")] <- paste0("Subpathway", r)
#   }
#   if (length(grep(x = names(dfsteps_genus),pattern="Subpathway")) >1 ){
#     dfsteps_genus$MaxCompletion <- apply(dfsteps_genus[,grep(x = names(dfsteps_genus),pattern="Subpathway")],1,max)
#   } else {
#     dfsteps_genus$MaxCompletion <- dfsteps_genus[,grep(x = names(dfsteps_genus),pattern="Subpathway")]
#   }
#   
#   # split into the different groups and then plot as box plots per group
#   for (genus in unique(dfsteps_genus$Genus)){
#     group <- filter(genomessums, Genus == genus)$taxa_elm[1]
#     dfsteps_genus$Group[which(dfsteps_genus$Genus == genus)] <- group
#   }
#   dfsteps_genus$Group <- as.factor(dfsteps_genus$Group)
#   
#   #get the types of genomes here
#   dfsteps_genus$Type <- "MAG"
#   dfsteps_genus$Type[grep("_REFG_", row.names(dfsteps_genus))] <- "RefGenome"
#   dfsteps_genus$Type[grep("_SAGS_", row.names(dfsteps_genus))] <- "SAG"
#   
#   write.csv(dfsteps_genus, file = paste0("DOSTraits/kofam_results/", unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), "processed_genus.csv"), quote = F)
#   
#   # also make a heatmap with the groups, with the heatmap color being the average value per gene, per group
#   s_dfsteps_genus <- dfsteps_genus %>%
#     group_by(Group) %>%
#     dplyr::summarise_all("mean")
#   s_dfsteps_genus <- gather( s_dfsteps_genus[,1:(ncol(s_dfsteps_genus)-1)], key = "key", value = "value", Step1:MaxCompletion)
#   plot <- ggplot(s_dfsteps_genus, aes(x=factor(key, level=unique(s_dfsteps_genus$key)), y = Group, fill = value)) +
#     geom_tile() + scale_fill_viridis_c(option = "rocket", direction = -1) +
#     labs(title = paste0(unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), ": Step Summaries"),
#          x = "Step or Calculation",
#          y = "Group")
#   plot <-  ggdraw() + draw_plot(p_tab, x = 0.85, width = .15) + draw_plot(plot, width = 0.85) 
#   ggsave(filename = paste0("DOSTraits/kofam_results/", unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), "_StepsHeatmap_genus.pdf"), plot = plot, width = 30, height = 20, limitsize = FALSE)
#   
#   dfsteps_genus$Lifestyle <- "SemiAuxotroph"
#   dfsteps_genus$Lifestyle[which(dfsteps_genus$MaxCompletion == 0)] <- "Auxotroph"
#   dfsteps_genus$Lifestyle[which(dfsteps_genus$MaxCompletion == 1)] <- "Assimilator"
#   plot <- ggplot(dfsteps_genus, aes(x=Group, fill=factor(Lifestyle, level=c("Assimilator", "SemiAuxotroph", "Auxotroph")))) + 
#     geom_histogram(position = "fill", stat="count") +
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#     ggtitle(paste0(unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]),": All Genomes Lifestyles"))
#   ggsave(filename = paste0("DOSTraits/kofam_results/", unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), "_allgenomesLifestyle_bygenus.pdf"), plot = plot, width = 15, height = 10, limitsize = FALSE)
#   
#   # lastly make the combinations summaries
#   dfsteps_genus <- dfsteps_genus[1:(nrow(dfsteps_genus)-1),]
#   combinations <- read.csv("DOSTraits/SulfateCombinations.csv")
#   combinationsprocessed <- as.data.frame(matrix(ncol=length(combinations$Combinations), nrow = nrow(dfsteps_genus)))
#   row.names(combinationsprocessed) <- dfsteps_genus$Genus
#   names(combinationsprocessed) <- combinations$Combinations
#   combinationsprocessed[,] <- 0
#   # Going to use dfsteps_genus, loop through the columns so we can sum by the steps
#   for (col in 1:ncol(combinationsprocessed)){
#     combo <- names(combinationsprocessed)[col]
#     combo <- strsplit(combo, split = "NOT")
#     need <- combo[[1]][1]
#     miss <- combo[[1]][2]
#     if (miss == " "){
#       miss <- ""
#       misscols <- as.data.frame(matrix(nrow=nrow(dfsteps_genus), ncol=1))
#       row.names(misscols) <- dfsteps_genus$Genus
#       names(misscols) <- "total"
#       misscols$total <- 0
#     } else {
#       misscols <- dfsteps_genus %>%
#         select(paste0("Step",strsplit(miss, split=",")[[1]]))
#       if (ncol(misscols) >1){
#         misscols$total <- apply(misscols[,], 1, sum)
#       } else if (ncol(misscols) == 1){
#         misscols$total <- misscols[,1]
#       }
#     }
#     
#     if (need == " "){
#       need <- ""
#       needcols <- as.data.frame(matrix(nrow=nrow(dfsteps_genus), ncol=1))
#       row.names(needcols) <- row.names(dfsteps_genus)
#       names(needcols) <- "total"
#       needcols$total <- 1
#     } else {
#       needcols <- dfsteps_genus %>%
#         select(paste0("Step",strsplit(need, split=",")[[1]]))
#       if (ncol(needcols) >1){
#         needcols$total <- apply(needcols[,], 1, sum)
#         needcols$total <- (needcols$total)/length(strsplit(need, split=",")[[1]])
#       } else if (ncol(needcols) == 1){
#         needcols$total <- needcols[,1]
#       } 
#     }
#     combinationsprocessed[which(misscols$total == 0 & needcols$total == 1),col] <- 1
#   }
#   
#   combinationsprocessed$Group <- dfsteps_genus$Group
#   s_combinationsprocessed <- combinationsprocessed %>%
#     group_by(Group) %>%
#     dplyr::summarise_all("mean")
#   s_combinationsprocessed <- gather( s_combinationsprocessed, key = "key", value = "value", `2,5,6,8NOT `:`9NOT3,5,6`)
#   plot <- ggplot(s_combinationsprocessed, aes(x=factor(key, level=unique(s_combinationsprocessed$key)), y = Group, fill = value)) +
#     geom_tile() + scale_fill_viridis_c(option = "rocket", direction = -1) +
#     labs(title = paste0(unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), ": Step Summaries for Genomes above 75% Complete"),
#          x = "Step or Calculation",
#          y = "Group")
#   plot <-  ggdraw() + draw_plot(p_tab, x = 0.85, width = .15) + draw_plot(plot, width = 0.85) 
#   ggsave(filename = paste0("DOSTraits/kofam_results/", unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), "_combinations_bygenus.pdf"), plot = plot, width = 60, height = 20, limitsize = FALSE)
#   write.csv(file = paste0("DOSTraits/kofam_results/", unique(sulfateKOs$Broad.pathway.description[which(sulfateKOs$Pathway.Number == p)]), "_combinations_bygenus.csv"), x= combinationsprocessed)
#   
#   
# }
