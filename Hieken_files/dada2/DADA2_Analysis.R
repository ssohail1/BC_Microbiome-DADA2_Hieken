### Note: Reverse reads were of very low quality, thus, not included in the analysis

# Set the working directory
library(dada2)
setwd("/media/mbb/Sidras_Projects/Hieken_paper/fastq_files/primer_folder-nogz/") #change to directory where files are located
dir()
path <- "/media/mbb/Sidras_Projects/Hieken_paper/fastq_files/primer_folder-nogz/"
list.files(path)
fnFs <- sort(list.files(path, pattern="_1.fastq_primer.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# To get the quality profile of the forward reads:
plotQualityProfile(fnFs[1:2])

# Assigning filenames for the filtered fastq.gz files
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_primer_filt.fastq"))
names(filtFs) <- sample.names

##minLen = 50 so that there are at least 50 nucleotides in all samples for the assignTaxonomy command
#truncLen is omitted as too many reads were being removed
out <- filterAndTrim(fnFs, filtFs,
                     maxN=0, maxEE= 3, truncQ=2, minLen = 50, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)


head(out)

# Learn the Error rates
errF <- learnErrors(filtFs, multithread=TRUE)
##195538305 total bases in 726311 reads from 2 samples will be used for learning the error rates.
plotErrors(errF, nominalQ=TRUE)

derepFs <- derepFastq(filtFs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names

# Applying the core sample inference algorithm to the filtered and trimmed sequence data.
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
#Sample 1 - 326354 reads in 87589 unique sequences.
#Sample 2 - 399957 reads in 93695 unique sequences.
#Sample 3 - 402883 reads in 90223 unique sequences.
#Sample 4 - 316201 reads in 79845 unique sequences.
#Sample 5 - 335077 reads in 66301 unique sequences.
#Sample 6 - 118973 reads in 48812 unique sequences.
#Sample 7 - 310267 reads in 67065 unique sequences.
#Sample 8 - 249862 reads in 64940 unique sequences.
#Sample 9 - 406087 reads in 79326 unique sequences.
#Sample 10 - 92598 reads in 42555 unique sequences.
#Sample 11 - 897202 reads in 183649 unique sequences.
#Sample 12 - 331370 reads in 91834 unique sequences.
#Sample 13 - 260231 reads in 51813 unique sequences.
#Sample 14 - 470699 reads in 102064 unique sequences.
#Sample 15 - 97625 reads in 40733 unique sequences.
#Sample 16 - 229642 reads in 57765 unique sequences.
#Sample 17 - 349734 reads in 89184 unique sequences.
#Sample 18 - 112871 reads in 59353 unique sequences.
#Sample 19 - 36062 reads in 16303 unique sequences.
#Sample 20 - 265094 reads in 65282 unique sequences.
#Sample 21 - 363439 reads in 87404 unique sequences.
#Sample 22 - 88018 reads in 37561 unique sequences.
#Sample 23 - 164147 reads in 32769 unique sequences.
#Sample 24 - 649398 reads in 207874 unique sequences.
#Sample 25 - 285799 reads in 73100 unique sequences.
#Sample 26 - 327734 reads in 58796 unique sequences.
#Sample 27 - 107387 reads in 56610 unique sequences.
#Sample 28 - 110617 reads in 52483 unique sequences.
#Sample 29 - 390287 reads in 85668 unique sequences.
#Sample 30 - 48309 reads in 23595 unique sequences.
#Sample 31 - 157149 reads in 101645 unique sequences.
#Sample 32 - 399103 reads in 105753 unique sequences.
#Sample 33 - 91741 reads in 46741 unique sequences.
#Sample 34 - 97103 reads in 63258 unique sequences.
#Sample 35 - 424581 reads in 79479 unique sequences.
#Sample 36 - 142243 reads in 88422 unique sequences.
#Sample 37 - 73210 reads in 34469 unique sequences.
#Sample 38 - 103189 reads in 50305 unique sequences.
#Sample 39 - 135728 reads in 85793 unique sequences.
#Sample 40 - 172999 reads in 107066 unique sequences.
#Sample 41 - 44066 reads in 26624 unique sequences.
#Sample 42 - 30995 reads in 15843 unique sequences.
#Sample 43 - 103335 reads in 64328 unique sequences.
#Sample 44 - 167338 reads in 108497 unique sequences.
#Sample 45 - 74690 reads in 45944 unique sequences.
#Sample 46 - 86832 reads in 38666 unique sequences.
#Sample 47 - 175963 reads in 107363 unique sequences.
#Sample 48 - 111814 reads in 64940 unique sequences.
#Sample 49 - 60126 reads in 26548 unique sequences.
#Sample 50 - 107084 reads in 67299 unique sequences.
#Sample 51 - 90515 reads in 41600 unique sequences.
#Sample 52 - 27487 reads in 13014 unique sequences.
#Sample 53 - 88929 reads in 61421 unique sequences.
#Sample 54 - 148413 reads in 95407 unique sequences.
#Sample 55 - 80409 reads in 45490 unique sequences.
#Sample 56 - 219131 reads in 136290 unique sequences.
#Sample 57 - 234684 reads in 74137 unique sequences.
#Sample 58 - 158067 reads in 101533 unique sequences.
#Sample 59 - 63126 reads in 36226 unique sequences.
#Sample 60 - 35213 reads in 19350 unique sequences.
#Sample 61 - 114872 reads in 74778 unique sequences.
#Sample 62 - 123884 reads in 79037 unique sequences.
#Sample 63 - 104758 reads in 62316 unique sequences.
#Sample 64 - 136162 reads in 88063 unique sequences.
#Sample 65 - 113804 reads in 50806 unique sequences.
#Sample 66 - 36211 reads in 19009 unique sequences.
#Sample 67 - 120208 reads in 77527 unique sequences.
#Sample 68 - 362191 reads in 96822 unique sequences.
#Sample 69 - 214070 reads in 134905 unique sequences.
#Sample 70 - 34191 reads in 20725 unique sequences.
#Sample 71 - 21739 reads in 10595 unique sequences.
#Sample 72 - 132427 reads in 83697 unique sequences.
#Sample 73 - 470232 reads in 261322 unique sequences.
#Sample 74 - 23740 reads in 14323 unique sequences.
#Sample 75 - 35267 reads in 15680 unique sequences.
#Sample 76 - 156898 reads in 98103 unique sequences.
#Sample 77 - 29819 reads in 16762 unique sequences.
#Sample 78 - 34877 reads in 17957 unique sequences.
#Sample 79 - 333156 reads in 80207 unique sequences.
#Sample 80 - 128346 reads in 84187 unique sequences.
#Sample 81 - 156644 reads in 102040 unique sequences.
#Sample 82 - 30892 reads in 17240 unique sequences.
#Sample 83 - 131752 reads in 84058 unique sequences.
#Sample 84 - 129898 reads in 79488 unique sequences.
#Sample 85 - 30353 reads in 18569 unique sequences.
#Sample 86 - 30878 reads in 15876 unique sequences.
#Sample 87 - 180152 reads in 112071 unique sequences.
#Sample 88 - 91561 reads in 53463 unique sequences.
#Sample 89 - 32523 reads in 17582 unique sequences.
#Sample 90 - 83521 reads in 34962 unique sequences.
#Sample 91 - 77539 reads in 51303 unique sequences.
#Sample 92 - 162603 reads in 103393 unique sequences.
#Sample 93 - 54556 reads in 31672 unique sequences.
#Sample 94 - 31613 reads in 15457 unique sequences.
#Sample 95 - 375393 reads in 228714 unique sequences.
#Sample 96 - 122459 reads in 64522 unique sequences.
#Sample 97 - 55610 reads in 24888 unique sequences.
#Sample 98 - 111116 reads in 69112 unique sequences.

dadaFs[[1]]
#dada-class: object describing DADA2 denoising results
#77 sequence variants were inferred from 87589 input unique sequences.
#Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16

# Constructing sequence table
seqtab <- makeSequenceTable(dadaFs)
dim(seqtab)
##[1]   98 8204

table(nchar(getSequences(seqtab)))

#50   51   52   54   55   57   59   61   62   63   66   67   68   69   70   71   72   73 
#1    2    1    1    1    1    2    5    2    1    1    3   12    1    7    1    4    3 
#74   75   76   77   78   79   80   81   82   83   85   86   87   88   89   92   93   94 
#129   47   12    6   16    4   17    5    9   33    4    2    1    1    2    2    1    1 
#96   97   98  100  101  103  106  107  113  116  119  129  130  139  142  143  146  148 
#9    1   12    2    1    5    1    2    1    2    1    2    1    1    1    1    1    1 
#149  152  153  155  156  164  172  176  182  186  187  188  189  195  199  208  209  214 
#1    1    1    1    1    1    1    1    1    1    3    1    1    1    1    1    4    2 
#223  227  232  235  236  239  241  244  246  251  256  258  259  261  262  263  264  265 
#1    1    2    1    1    2    1    1    1    3    2    4    3    1    1    2    3    1 
#266  267  268  269  270  271  272  273  274  275  276  278  279  280  281  282  283  284 
#3    5    2   10    3    4    2    5    6    3    6    2    4    7   15   12   37   24 
#285  286  287  288  289  290  291  292  293  294  295  296  297  298  299  300  301 
#58   29   38   61   90   75  153  107   60   92  109  394  582  114  178 1751 3722


# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
#Identified 2093 bimeras out of 8204 input sequences.
#minFoldParentOverAbundance not required here

sum(seqtab.nochim)/sum(seqtab)
##[1] 0.9734543  ~ should be close to 1

rowSums(seqtab.nochim)
rowSums(seqtab)

# Track reads thru pipeline
getN <- function(x) sum(getUniques(x))

# According to Dada2 tutorial, make the following modification to the track command: 
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) 
# with getN(dadaFs)
track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "nonchim")
rownames(track) <- sample.names
head(track)


# Assign Taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "/media/mbb/Sidras_Projects/Hieken_paper/fastq_files/primer_folder-nogz/silva.nr_v138.align", multithread = TRUE)
## Went through!

taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)


# Did not evaluate accuracy with mock data as no mock samples in our data

# Phyloseq
## Need to download phyloseq first
#if (!requireNamespace("BiocManager", quietly = TRUE))
  # install.packages("BiocManager")
#  BiocManager::install("phyloseq")
## Then load library
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")

# Steps to create a metadata (samdf file) - have metadata so can read in to samdf below
theme_set(theme_bw())
samples.out <- rownames(seqtab.nochim)
Sample_Name <- sapply(strsplit(samples.out, "D"), `[`, 1)
gender <- substr(Sample_Name,1,1)
subject <- substr(Sample_Name,2,999)
day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
samdf <- data.frame(sample_name=Sample_Name, Name=gender, Day=day)
samdf$When <- "Early"
samdf$When[samdf$Day>100] <- "Late"
rownames(samdf) <- samples.out

# Made different samdf files for different graphs specific for the treatment types
samdf <- read.table("~/Desktop/H_Metadata.txt", header=TRUE)
samdf <- read.table("~/Desktop/HMetadata_Tissues.txt", header=TRUE)

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),
               sample_data(samdf), # originally sample_data
               tax_table(taxa))
ps <- prune_samples(sample_names(ps) != "Mock", ps)
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps
plot_richness(ps, x = "env_biome", measures=c("Shannon", "Simpson"), color = "final_dx")
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray", color="env_biome")
plot_ordination(ps.prop, ord.nmds.bray, color="env_biome", title="Bray NMDS")
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x = "env_biome", fill = "Family") + facet_wrap(~final_dx, scales = "free_x")
plot_bar(ps.top20, x = "env_biome", fill = "Phylum") + facet_wrap(~final_dx, scales = "free_x")
plot_bar(ps.top20, x = "env_biome", fill = "Genus") + facet_wrap(~final_dx, scales = "free_x") #978 873

# End of dada2 Tutorial
