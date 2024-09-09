<h1 align="center">RepTim</h1>


&nbsp;


## Presentation
RepTim is an R package for fast and efficient DNA replication timing analysis.

The package aims to perform the differential analysis between replication timing profiles (biological conditions) in order to detect modified regions. 
Then, it localizes genes in regions of interest (e.g. modified regions) and tests if biological pathways are enriched.
It does also a Fisher's exact test to find transcription factors for wich more genes (targets) are found in the modified regions in comparison with non-modified regions.


&nbsp;


## Installation

### Installation with R commands

Open an R console or RStudio and use this code :

```
load.install.package <- lapply(
  c("devtools"),
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

devtools::install_github("AlexandreBoulle/RepTim")
```

### Installation from source code

#### Linux / macOS installation
Open a terminal and use these commands :

```
git clone https://github.com/AlexandreBoulle/RepTim.git
tar -czvf RepTim.tar.gz ./RepTim/
```

Open an R console or RStudio and use this command :

```
install.packages("/your_path/Reptim.tar.gz", repo = NULL, type = "source")
```

#### Windows installation

* Click on the green "Code" button (top right) and choose "Download ZIP"
* Uncompress the folder and then compress it to ".tar.gz" using 7-Zip or another tool (compress to ".tar" and then compress ".tar" folder to ".gz")
* Open an R console or RStudio and use this command :

```
install.packages("/your_path/RepTim-main.tar.gz", repo = NULL, type = "source")
```


&nbsp;


## Package using

### 1/8 : Load libraries

```
library(stringr)
library(tftargets)
library(org.Hs.eg.db)
library(enrichR)
library(RepTim)
```

### 2/8 : Set a path to load bedgraph files and to write results

```
setwd("/your_path/")
```

### 3/8 : Load input data (databases)

```
# List of human TF
fpath.tf <- system.file("extdata", "Human_TF.txt", package = "RepTim")
tf <- read.table(fpath.tf, header = FALSE, sep = "\t")[, 1]

# Ensembl database : gene names with positions
fpath.ensembldb <- system.file("extdata", "BioMart_Ensembl_Grch38p13_02092022.txt", package = "RepTim")
ensembldb <- read.table(fpath.ensembldb, header = TRUE, sep = "\t")
ensembldb.genetype <- ensembldb[ensembldb$Gene_type == "protein_coding", ]

# Chromosome Size
fpath.chr <- system.file("extdata", "Human_chr_size.txt", package = "RepTim")
size.chr.table <- read.table(fpath.chr, header = TRUE, sep = "\t")
```

### 4/8 : Load input data (biological conditions = Replication Timing profiles)

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

### 5/8 : Merge replicates

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

### 6/8 : Choose Thresholds

```
pval <- 1e-3
per.dist.detect <- 90
per.dist.elong <- 10
organism <- "org.Hs.eg"
```

### 7/8 : Choose conditions to compare

```
cond1 <- cond.WT
cond2 <- cond.NA.late
cond1.loess <- cond.loess.WT
cond2.loess <- cond.loess.NA.late
```

### 8/8 : RepTim functions

```
list.cond <- list(cond1.WT, cond2.WT, cond1.NA.early, cond2.NA.early, cond1.NA.late, cond2.NA.late)
comparison <- "WT_vs_NA-late_pval-1e-3_dist-detect-90_dist-elong-10"

res <- modified.regions.detection(list.cond, cond1, cond2, cond1.loess, cond2.loess, pval, per.dist.detect, per.dist.elong, comparison)
df.nmregions <- non.modified.regions(res[[1]], size.chr.table, comparison)
df.genes <- gene.annotation(ensembldb.genetype, comparison)
df.tf.targets <- targets.tf.enrichment(comparison, tf, organism)
df.biopath <- enrich.biopathways(comparison)
```
