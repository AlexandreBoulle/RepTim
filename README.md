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
* Open a R console or RStudio and use this command :

```
install.packages("/your_path/RepTim-main.tar.gz", repo = NULL, type = "source")
```


&nbsp;


## Package using
