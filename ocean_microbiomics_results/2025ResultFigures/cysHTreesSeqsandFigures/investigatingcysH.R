library(tidyr)
library(dplyr)
library(ggplot2)
library("seqinr")

# load the large sequence file
cysH_genes <- seqinr::read.fasta(file = "genes-redundant.faa", seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
# subset the large seq file to just those genes designated as cysH
cysh_gene_list <- read.table("new_filtered_SA/cysH_hits_filt.txt")
all(names(cysH_genes) %in% cysh_gene_list$V5)
#FALSE
cysH_genes <- cysH_genes[which(names(cysH_genes) %in% cysh_gene_list$V5)]
all(names(cysH_genes) %in% cysh_gene_list$V5)
write.fasta(sequences = cysH_genes, names = names(cysH_genes), nbchar = 80, file.out = "cysH_genes_notgenomefiltered.fasta")

# correlate the sequences with the correct genomes
head(names(cysH_genes))
head(match(names(cysH_genes),cysh_gene_list$V5))
names(cysH_genes) <- cysh_gene_list$V3[match(names(cysH_genes),cysh_gene_list$V5)]

# eliminate duplicates?
length(which(!duplicated(names(cysH_genes))==TRUE))
cysH_genes <- cysH_genes[!duplicated(names(cysH_genes))]

#resave the file with just our genomes in it
our_genomes <- read.csv("new_filtered_SA/2025_1022_MAD_SulfateAssimilationprocessed_above75genomecompletion.csv")
all(names(cysH_genes) %in% our_genomes$Genome)
#FALSE
cysH_genes <- cysH_genes[which(names(cysH_genes) %in% our_genomes$Genome)]
write.fasta(sequences = cysH_genes, names = names(cysH_genes), nbchar = 80, file.out = "cysH_genes.fasta")
