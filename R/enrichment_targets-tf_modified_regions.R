####################
#      REPTIM      #
# Alexandre BOULLE #
####################


###########################
# ENRICHMENT OF TF TARGET #
###########################


#' @import annotate
#' @importFrom stats fisher.test
#' @import stringr


#' @param tf as all transcription factors known for an organism
#' @param genes.modified as genes localized in modified regions
#' @param genes.non.modified as genes localized in non-modified regions
#' @param organism as the name of the organism
tf.enrichment <- function(tf, genes.modified, genes.non.modified, organism){
  # List of transcription factor targets
  list.target <- list(ENCODE[tf],
                      ITFP[tf],
                      Marbach2016[tf],
                      TRED[tf],
                      TRRUST[tf])
  # Initialize a dataframe
  res <- data.frame(TF = rep(NA, length(tf)),
                    n.1.1 = rep(NA, length(tf)),
                    n.1.2 = rep(NA, length(tf)),
                    n.2.1 = rep(NA, length(tf)),
                    n.2.2 = rep(NA, length(tf)),
                    pv = rep(NA, length(tf)))
  res.targets <- data.frame(TF = rep(NA, length(tf)), Targets = rep(NA, length(tf)))
  res.targets.ideogram <- data.frame(TF = rep(NA, length(tf)), Targets = rep(NA, length(tf)))
  for (i in 1:length(tf)){
    res[i, 1] <- tf[i]
    res.targets[i, 1] <- tf[i]
    res.targets.ideogram[i, 1] <- tf[i]
    if (!is.null(list.target[[1]][[i]])){
      targets.ENCODE <- unname(getSYMBOL(list.target[[1]][[i]], data = organism))
    }else{
      targets.ENCODE <- ""
    }
    if (!is.null(list.target[[2]][[i]])){
      targets.ITFP <- list.target[[2]][[i]]
    }else{
      targets.ITFP <- ""
    }
    if (!is.null(list.target[[3]][[i]])){
      targets.Marbach2016 <- list.target[[3]][[i]]
    }else{
      targets.Marbach2016 <- ""
    }
    if (!is.null(list.target[[4]][[i]])){
      targets.TRED <- unname(getSYMBOL(list.target[[4]][[i]], data = organism))
    }else{
      targets.TRED <- ""
    }
    if (!is.null(list.target[[5]][[i]])){
      targets.TRRUST <- list.target[[5]][[i]]
    }else{
      targets.TRRUST <- ""
    }
    targets <- sort(unique(c(targets.ENCODE, targets.ITFP, targets.Marbach2016, targets.TRED, targets.TRRUST)))
    ind.remove <- which(targets %in% "")
    targets <- targets[-ind.remove]
    genes.in.targets <- genes.modified[genes.modified %in% targets]
    if (length(genes.in.targets) == 0){
      res.targets[i, 2] <- ""
      res.targets.ideogram[i, 2] <- ""
    }else{
      res.targets[i, 2] <- paste(sort(genes.in.targets), collapse = " ; ")
      res.targets.ideogram[i, 2] <- paste(sort(genes.in.targets), collapse = ",")
    }
    # Number of genes located in modified regions and which are TF targets
    n.1.1 <- length(which(genes.modified %in% targets))
    # Number of genes located in non-modified regions and which are TF targets
    n.1.2 <- length(which(genes.non.modified %in% targets))
    # Number of genes located in modified regions and which aren't TF targets
    n.2.1 <- length(which(!(genes.modified %in% targets)))
    # Number of genes located in non-modified regions and which aren't TF targets
    n.2.2 <- length(which(!(genes.non.modified %in% targets)))
    # Enrichment analysis : Fisher test
    stat <- fisher.test(matrix(c(n.1.1, n.1.2, n.2.1, n.2.2), nrow = 2, byrow = TRUE), alternative = "greater")
    res[i, 2] <- n.1.1
    res[i, 3] <- n.1.2
    res[i, 4] <- n.2.1
    res[i, 5] <- n.2.2
    res[i, 6] <- stat$p.value
  }
  res <- res[order(res[, "pv"]),]
  res$odds.ratio <- (res$n.1.1 / res$n.1.2) / (res$n.2.1 / res$n.2.2)
  res$padj <- p.adjust(res$pv, method = "fdr")
  res <- res[, c("TF", "n.1.1", "n.1.2", "n.2.1", "n.2.2", "odds.ratio", "pv", "padj")]
  return(list(res, res.targets, res.targets.ideogram))
}


#' @title Transcription factor enrichment analysis
#' @description Perform a fisher test to find the transcription factors with the most gene targets in modified regions
#' @param comparison the name of the differential analysis performed between profiles
#' @param tf all transcription factors known for an organism
#' @param organism the name of the organism
#' @export
targets.tf.enrichment <- function(comparison, tf, organism){
  list.tf.adv <- list()
  list.tf.del <- list()
  dir.create(str_glue("Differential-Analysis_{comparison}/TF_enrichment_analysis/"))
  # for (status in c("advanced", "delayed")){
  for (status in c("advanced", "delayed")){
    # 1) MODIFIED REGIONS
    differential.result <- read.table(str_glue("Differential-Analysis_{comparison}/Genes_in_modified_regions/Genes_in_Regions_{status}.txt"), header = TRUE, sep = "\t")
    genes.body.modified <- unique(unlist(str_split(differential.result$GENES.BODY, " ; ")))
    genes.overlap.modified <- unique(unlist(str_split(differential.result$GENES.OVERLAPPING, " ; ")))
    genes.modified <- unique(c(genes.body.modified, genes.overlap.modified))
    # Remove "" string recovered from empty lines
    genes.modified <- genes.modified[genes.modified != ""]
    # 2) NON-MODIFIED REGIONS
    df.non.modified.remove <- read.table(str_glue("Differential-Analysis_{comparison}/Non-modified_regions/Genes_in_reduced-non-modified-regions_{status}.txt"), header = TRUE, sep = "\t")
    genes.body.non.modified <- unique(unlist(str_split(df.non.modified.remove$GENES.BODY, " ; ")))
    genes.overlap.non.modified <- unique(unlist(str_split(df.non.modified.remove$GENES.OVERLAPPING, " ; ")))
    # Remove "" string recovered from empty lines
    genes.non.modified <- unique(c(genes.body.non.modified, genes.overlap.non.modified))
    genes.non.modified <- genes.non.modified[genes.non.modified != ""]
    # 3) ENRICHMENT ANALYSIS
    res.tf <- tf.enrichment(tf, genes.modified, genes.non.modified, organism)[[1]]
    res.targets <- tf.enrichment(tf, genes.modified, genes.non.modified, organism)[[2]]
    res.targets.ideogram <- tf.enrichment(tf, genes.modified, genes.non.modified, organism)[[3]]
    if (status == "advanced"){
      list.tf.adv[[length(list.tf.adv) + 1]] <- res.tf
      list.tf.adv[[length(list.tf.adv) + 1]] <- res.targets
      list.tf.adv[[length(list.tf.adv) + 1]] <- res.targets.ideogram
    }else{
      list.tf.del[[length(list.tf.del) + 1]] <- res.tf
      list.tf.del[[length(list.tf.del) + 1]] <- res.targets
      list.tf.del[[length(list.tf.del) + 1]] <- res.targets.ideogram
    }
    write.table(res.tf, str_glue("Differential-Analysis_{comparison}/TF_enrichment_analysis/Targets-TF_Enrichment_{status}-regions.txt"), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
    write.table(res.targets, str_glue("Differential-Analysis_{comparison}/TF_enrichment_analysis/Targets-TF_{status}-regions.txt"), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
    write.table(res.targets.ideogram, str_glue("Differential-Analysis_{comparison}/TF_enrichment_analysis/Targets-TF_{status}-regions_ideogram-viewer.txt"), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
  }
  return(list(list.tf.adv, list.tf.del))
}
