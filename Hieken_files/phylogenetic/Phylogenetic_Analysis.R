### Phylogenetic Analysis ###

#Load libraries
library(dada2)

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("DECIPHER")

library(DECIPHER)
library(phangorn)
library(ggplot2)
library(phyloseq)

#Phylogenetic Tree
## Extract sequences from DADA2 output
sequences<-getSequences(seqtab.nochim)
names(sequences)<-sequences

## Run Sequence Alignment (MSA) using DECIPHER
alignment <- AlignSeqs(DNAStringSet(sequences), anchor=NA)

## Change sequence alignment output into a phyDat structure
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")

## Create distance matrix
dm <- dist.ml(phang.align)

## Perform Neighbor joining
treeNJ <- NJ(dm) # Note, tip order != sequence order

## Internal maximum likelihood
fit = pml(treeNJ, data=phang.align)

## negative edges length changed to 0!
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))

#Import into phyloseq
map <- import_qiime_sample_data("~/Desktop/HMetadata_Tissues.txt")
##Assigning species
taxa.plus <- addSpecies(taxa, "/media/mbb/Sidras_Projects/Hieken_paper/fastq_files/primer_folder-nogz/silva.nr_v138.align", verbose=TRUE)
colnames(taxa.plus) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus","Species")
unname(taxa.plus)
##0 out of 6111 were assigned to the species level.
##Of which 0 had genera consistent with the input table.
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               tax_table(taxa.plus),phy_tree(fitGTR$tree))

## Merge PhyloSeq object with map
ps <- merge_phyloseq(ps,map)
ps

## Will need to root tree to calculate phylogeny based 
## diversity metrics (like Unifrac)

set.seed(711)
phy_tree(ps) <- root(phy_tree(ps), sample(taxa_names(ps), 1), resolve.root = TRUE)
is.rooted(phy_tree(ps))
#plot(phy_tree(ps))


#For UniFrac Plots
##Unweighted UniFrac
UniFrac(ps) #Weighted default is False
unweighted <- UniFrac(ps)
plot(unweighted)

##Weighted UniFrac
UniFrac(ps, weighted = TRUE)
weighted <- UniFrac(ps, weighted = TRUE)
plot(weighted)


