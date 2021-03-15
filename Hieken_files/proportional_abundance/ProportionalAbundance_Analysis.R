### Proportional Abundance Plots ###

library(phyloseq)

#First go through this set of commands, so to use the correct ps file as input for making
#the phylum, genus, and family plots.
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),
               sample_data(samdf), # originally sample_data
               tax_table(taxa))
ps <- prune_samples(sample_names(ps) != "Mock", ps)
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray", color="env_biome")
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)

#Using the ps file from above to generate the Phylum, Family, and Genus plots

ps <- tax_glom(ps, "Phylum")
ps0 <- transform_sample_counts(ps, function(x) x / sum(x))
ps1 <- merge_samples(ps0, "final_dx")
ps2 <- transform_sample_counts(ps1, function(x) x / sum(x))
plot_bar(ps2, fill="Phylum")

ps <- tax_glom(ps, "Family")
ps0 <- transform_sample_counts(ps, function(x) x / sum(x))
ps1 <- merge_samples(ps0, "final_dx")
ps2 <- transform_sample_counts(ps1, function(x) x / sum(x))
plot_bar(ps2, fill="Family")

                               
ps <- tax_glom(ps, "Genus")
ps0 <- transform_sample_counts(ps, function(x) x / sum(x))
ps1 <- merge_samples(ps0, "final_dx")
ps2 <- transform_sample_counts(ps1, function(x) x / sum(x))
plot_bar(ps2, fill="Genus")
