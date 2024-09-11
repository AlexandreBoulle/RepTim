<h1 align="center">RepTim</h1>


&nbsp;


## Presentation
RepTim is an R package for fast (< 2 minutes with test datasets) and efficient DNA Replication Timing analysis.

<div align="justify"> The package aims to perform the differential analysis between replication timing profiles (biological conditions) in order to detect significantly modified regions. Then, it localizes genes in regions of interest (e.g. modified regions) and tests if biological pathways are enriched. It does also a Fisher's exact test to find transcription factors for wich more targets (genes) are found in the modified regions in comparison with non-modified regions. </div>


&nbsp;


## Installation

### Installation with R commands

Open an R console or RStudio and use this code:

```
if (!require("devtools", character.only = TRUE)){
  install.packages("devtools", dependencies = TRUE)
  library("devtools", character.only = TRUE)
}

devtools::install_github("AlexandreBoulle/RepTim")
```

### Installation from source code

#### Linux / macOS installation
Open a terminal and use these commands:

```
git clone https://github.com/AlexandreBoulle/RepTim.git
tar -czvf RepTim.tar.gz ./RepTim/
```

Open an R console or RStudio and use this command:

```
install.packages("/your_path/Reptim.tar.gz", repo = NULL, type = "source")
```

#### Windows installation

* Click on the green "Code" button (top right) and choose "Download ZIP"
* Uncompress the folder and then compress it to ".tar.gz" using 7-Zip or another tool (compress to ".tar" and then compress ".tar" folder to ".gz")
* Open an R console or RStudio and use this command:

```
install.packages("/your_path/RepTim-main.tar.gz", repo = NULL, type = "source")
```


&nbsp;


## Package using

### 1/8: Load libraries

```
if (!require("enrichR", character.only = TRUE)){
  install.packages("enrichR")
  library("enrichR", character.only = TRUE)
}

if (!require("org.Hs.eg.db", character.only = TRUE)){
  BiocManager::install("org.Hs.eg.db")
  library("org.Hs.eg.db", character.only = TRUE)
}

if (!require("tftargets", character.only = TRUE)){
  devtools::install_github("slowkow/tftargets")
  library("tftargets", character.only = TRUE)
}

library(stringr)
library(tftargets)
library(org.Hs.eg.db)
library(enrichR)
library(RepTim)
```

### 2/8: Set a path to load bedgraph files and to write results

```
setwd("/your_path/")
```

### 3/8: Load input data (databases)

For human (GRCh38.p13), you can load databases directly from the package. \
If you need another genome version or genomes of other species, it is necessary to use your own databases:
* A list of known transcription factors (TF) for a species (don't write a header in the file)
* A dataframe containing at least 4 columns: "Chromosome", "Gene_name", "Gene_start_bp" and "Gene_end_bp"
* A dataframe containing 2 columns: "Chromosome" and "Seq_length"

```
# List of human TF
fpath.tf <- system.file("extdata", "Human_TF.txt", package = "RepTim")
tf <- read.table(fpath.tf, header = FALSE, sep = "\t")[, 1]

# Ensembl database: gene names with positions
fpath.ensembldb <- system.file("extdata", "BioMart_Ensembl_Grch38p13_02092022.txt", package = "RepTim")
ensembldb <- read.table(fpath.ensembldb, header = TRUE, sep = "\t")
ensembldb.genetype <- ensembldb[ensembldb$Gene_type == "protein_coding", ]

# Chromosome Size
fpath.chr <- system.file("extdata", "Human_chr_size.txt", package = "RepTim")
size.chr.table <- read.table(fpath.chr, header = TRUE, sep = "\t")
```

### 4/8: Load input data (biological conditions = Replication Timing profiles)

To load this data, you must have a "Bedgraph_files" folder in your directory. \
This folder needs to contain bedgraph files with the same type of name ("condition.bedGraph" and "condition_Loess.bedGraph"). \
Example: if condition = "54" so the file name are "54.bedGraph" and "54_Loess.bedGraph"

```
# WT
cond1.WT <- load.data(name.experiment = "5967")[[1]]
cond1.loess.WT <- load.data(name.experiment = "5967")[[2]]
cond2.WT <- load.data(name.experiment = "MP21")[[1]]
cond2.loess.WT <- load.data(name.experiment = "MP21")[[2]]

# NA early
cond1.NA.early <- load.data(name.experiment = "54")[[1]]
cond1.loess.NA.early <- load.data(name.experiment = "54")[[2]]
cond2.NA.early <- load.data(name.experiment = "70")[[1]]
cond2.loess.NA.early <- load.data(name.experiment = "70")[[2]]

# NA late
cond1.NA.late <- load.data(name.experiment = "672")[[1]]
cond1.loess.NA.late <- load.data(name.experiment = "672")[[2]]
cond2.NA.late <- load.data(name.experiment = "673")[[1]]
cond2.loess.NA.late <- load.data(name.experiment = "673")[[2]]
```

### 5/8: Merge replicates

```
# WT
cond.WT <- replicate.merging(list(cond1.WT, cond2.WT))
cond.loess.WT <- replicate.merging(list(cond1.loess.WT, cond2.loess.WT))
# NA early
cond.NA.early <- replicate.merging(list(cond1.NA.early, cond2.NA.early))
cond.loess.NA.early <- replicate.merging(list(cond1.loess.NA.early, cond2.loess.NA.early))
# NA late
cond.NA.late <- replicate.merging(list(cond1.NA.late, cond2.NA.late))
cond.loess.NA.late <- replicate.merging(list(cond1.loess.NA.late, cond2.loess.NA.late))
```

### 6/8: Choose Thresholds

**NOTE 1**: The variables "per.dist.detect" and "per.dist.elong" are values between 0 and 100 (= a percentage). \

**NOTE 2**: If "per.dist.detect" is equal to 90% so it means we keep chromosome regions where the distance between profiles (curves) are in the 9th decile of the distance distribution (distances calculated between all replicates). \

**NOTE 3**: If "per.dist.elong" is equal to 10% so it means we extend discovered regions as long as the distance between profiles is higher than the 1st decile.

```
pval <- 1e-3
per.dist.detect <- 90
per.dist.elong <- 10
organism <- "org.Hs.eg"
```

### 7/8: Choose conditions to compare

```
cond1 <- cond.WT
cond2 <- cond.NA.late
cond1.loess <- cond.loess.WT
cond2.loess <- cond.loess.NA.late
```

### 8/8: RepTim functions

**NOTE 1**: Result folders are written directly in your directory. \

**NOTE 2**: "Files_START-R_Viewer" folder contains ".SRV" files that can be loaded into the START-R Viewer shiny application in order to visualize result graphs (curve of Replication Timing profiles, annotations to see modified regions and their status "advanced" or "delayed"). \

**NOTE 3**: To obtain more informations about parameters used in R functions, you can run this R command:

```
?modified.regions.detection
```

```
list.cond <- list(cond1.WT, cond2.WT, cond1.NA.early, cond2.NA.early, cond1.NA.late, cond2.NA.late)
comparison <- "WT_vs_NA-late_pval-1e-3_dist-detect-90_dist-elong-10"

res <- modified.regions.detection(list.cond, cond1, cond2, cond1.loess, cond2.loess, pval, per.dist.detect, per.dist.elong, comparison)
df.nmregions <- non.modified.regions(res[[1]], size.chr.table, comparison)
df.genes <- gene.annotation(ensembldb.genetype, comparison)
df.tf.targets <- targets.tf.enrichment(comparison, tf, organism)
df.biopath <- enrich.biopathways(comparison)
```
