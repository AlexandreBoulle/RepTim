####################
#      REPTIM      #
# Alexandre BOULLE #
####################


#########################
# FIND GENES IN REGIONS #
#########################


#' @import stringr

#' @param table.regions as the dataframe containing modified regions
#' @param gene.database as the database containing gene positions
find.genes.in.regions <- function(table.regions, gene.database){
  # Initialize columns
  table.regions$GENES.BODY <- ""
  table.regions$GENES.OVERLAPPING <- ""
  for (i in 1:nrow(table.regions)){
    # Genomic coordinates of different regions
    chr <- unlist(str_split(table.regions$CHR[i], "r"))[2]
    start <- table.regions$START[i]
    end <- table.regions$END[i]
    # Find genes within these regions of interest
    ind.body <- which(gene.database$Chromosome == chr & gene.database$Gene_start_bp >= start & gene.database$Gene_end_bp <= end)
    genes.body <- sort(unique(gene.database[ind.body, "Gene_name"]))
    genes.body <- genes.body[genes.body != ""]
    gene.body.list <- paste(genes.body, collapse = " ; ")
    # Find overlapping genes
    if (length(ind.body) > 0){
      subset.db <- gene.database[-ind.body, ]
    }else{
      subset.db <- gene.database
    }
    ind.overlap1 <- which(subset.db$Chromosome == chr & subset.db$Gene_start_bp <= start & subset.db$Gene_end_bp >= end)
    ind.overlap2 <- which(subset.db$Chromosome == chr & subset.db$Gene_start_bp >= start & subset.db$Gene_start_bp <= end & subset.db$Gene_end_bp >= end)
    ind.overlap3 <- which(subset.db$Chromosome == chr & subset.db$Gene_start_bp <= start & subset.db$Gene_end_bp >= start & subset.db$Gene_end_bp <= end)
    genes.overlap <- sort(unique(subset.db[c(ind.overlap1, ind.overlap2, ind.overlap3), "Gene_name"]))
    genes.overlap <- genes.overlap[genes.overlap != ""]
    gene.overlap.list <- paste(genes.overlap, collapse = " ; ")
    # Add gene names in the appropriate column
    table.regions$GENES.BODY[i] <- gene.body.list
    table.regions$GENES.OVERLAPPING[i] <- gene.overlap.list
  }
  return(table.regions)
}

#' @title Gene detection
#' @description Find genes in regions of interest
#' @param ensembldb.genetype the database containing gene name and gene positions
#' @param comparison the comparison name used in the result folder
#' @export
gene.annotation <- function(ensembldb.genetype, comparison){
  for (status in c("advanced", "delayed")){
    differential.result <- read.table(str_glue("Differential-Analysis_{comparison}/Modified_regions/{status}_regions.txt"), header = TRUE, sep = "\t")
    df.non.modified.reduce <- read.table(str_glue("Differential-Analysis_{comparison}/Non-modified_regions/Reduced_non-modified_regions_for-{status}.txt"), header = TRUE, sep = "\t")
    if (dim(differential.result)[1] > 0){
      dir.create(str_glue("Differential-Analysis_{comparison}/Genes_in_modified_regions/"))
      differential.result.genes <- find.genes.in.regions(table.regions = differential.result, gene.database = ensembldb.genetype)
      write.table(differential.result.genes, str_glue("Differential-Analysis_{comparison}/Genes_in_modified_regions/Genes_in_Regions_{status}.txt"), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
    }
    if (dim(df.non.modified.reduce)[1] > 0 ){
      df.non.modified.reduce.genes <- find.genes.in.regions(table.regions = df.non.modified.reduce, gene.database = ensembldb.genetype)
      write.table(df.non.modified.reduce.genes, str_glue("Differential-Analysis_{comparison}/Non-modified_regions/Genes_in_reduced-non-modified-regions_{status}.txt"), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
    }
  }
  return(list(differential.result.genes, df.non.modified.reduce.genes))
}

