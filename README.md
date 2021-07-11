# BC_Microbiome-DADA2_Hieken
## Sidra Sohail 
## Analyzing the Breast Cancer Microbiome
### DADA2 Analysis of the study by [Hieken et al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4971513/)
#### Background
The microbiome is a collection of bacteria that reside in tissues where in cancer it has been implicated in a variety of tissues, and specifically distinct microbial communities have been linked to different types of cancers [[Hieken et al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4971513/)]. The association of microbial communities to human health has prompted researchers to investigate the microbial ecology of various tissues leading to a surplus of microbial data. This in turn has led to the development of various software and tools that can analyze microbial data and disentangle the biological variation from amplicon sequencing errors [[Callahan et al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4927377/)]. Here, I analyze data from _The Microbiome of Aseptically Collected Human Breast Tissue in Benign and Malignant Disease_ study by [[Hieken et al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4971513/)] using DADA2 [[Callahan et al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4927377/)], where the dataset can be downloaded from the [SRA](https://www.ncbi.nlm.nih.gov/sra) using this accession number PRJNA335375.
#### Packages and software to install before running
R: Bioconductor, dada2, phyloseq, Biostrings, ggplot2

#### Files in repo
- PCoA folder:
  - Starting_PCoA_Analysis.R: Code for starting the PCoA analysis
- dada2 folder:
  - DADA2_Analysis.R: Rscript for running [DADA2 pipeline](https://benjjneb.github.io/dada2/tutorial_1_8.html)
- phylogenetic folder:
  - Phylogenetic_Analysis.R: Code for phylogenetic analyses where the final product is a tree that is visualized through weighted and unweighted UniFrac plots.
- proportional_abundance folder:
  - ProportionalAbundance_Analysis.R: Code to visualize the phyloseq object in the form of abundance plots at the phylum, family, and genus levels.

#### How to use
Clone repository into personal directory using this command,  
`git clone https://github.com/ssohail1/BC_Microbiome-DADA2_Hieken.git`

To move into BC_Microbiome-DADA2_Hieken directory use `cd`,  
`cd BC_Microbiome-DADA2_Hieken`

