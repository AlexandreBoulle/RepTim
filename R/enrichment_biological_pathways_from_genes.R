####################
#      REPTIM      #
# Alexandre BOULLE #
####################


#################################
# BIOLOGICAL PATHWAY ENRICHMENT #
#################################


#' @import enrichR
#' @import stringr

#' @title Biological pathway enrichment
#' @description Find biological pathway enrichment from genes localized in modified regions
#' @param comparison the comparison name used in the result folder
#' @export
enrich.biopathways <- function(comparison){
  dir.create(str_glue("Differential-Analysis_{comparison}/Gene_enrichment_analysis_in_modified_regions/"))
  list.df.adv <- list()
  list.df.del <- list()
  for (status in c("advanced", "delayed")){
    differential.result <- read.table(str_glue("Differential-Analysis_{comparison}/Genes_in_modified_regions/Genes_in_Regions_{status}.txt"), header = TRUE, sep = "\t")
    genes.body.modified <- unique(unlist(str_split(differential.result$GENES.BODY, " ; ")))
    genes.overlap.modified <- unique(unlist(str_split(differential.result$GENES.OVERLAPPING, " ; ")))
    genes.modified <- unique(c(genes.body.modified, genes.overlap.modified))
    # Remove "" string recovered from empty lines
    genes.modified <- genes.modified[genes.modified != ""]
    dbs <- c("GO_Molecular_Function_2018", "GO_Cellular_Component_2018", "GO_Biological_Process_2018")
    enriched <- enrichr(genes.modified, dbs)
    if (status == "advanced"){
      list.df.adv[[length(list.df.adv) + 1]] <- enriched$GO_Biological_Process_2018
      list.df.adv[[length(list.df.adv) + 1]] <- enriched$GO_Cellular_Component_2018
      list.df.adv[[length(list.df.adv) + 1]] <- enriched$GO_Molecular_Function_2018
    }else{
      list.df.del[[length(list.df.del) + 1]] <- enriched$GO_Biological_Process_2018
      list.df.del[[length(list.df.del) + 1]] <- enriched$GO_Cellular_Component_2018
      list.df.del[[length(list.df.del) + 1]] <- enriched$GO_Molecular_Function_2018
    }
    write.table(enriched$GO_Biological_Process_2018, str_glue("Differential-Analysis_{comparison}/Gene_enrichment_analysis_in_modified_regions/Enrichment_{status}-regions_GO_Biological_Process_2018.txt"),
                col.names = TRUE,
                row.names = FALSE,
                sep = "\t",
                quote = FALSE)
    write.table(enriched$GO_Cellular_Component_2018, str_glue("Differential-Analysis_{comparison}/Gene_enrichment_analysis_in_modified_regions/Enrichment_{status}-regions_GO_Cellular_Component_2018.txt"),
                col.names = TRUE,
                row.names = FALSE,
                sep = "\t",
                quote = FALSE)
    write.table(enriched$GO_Molecular_Function_2018, str_glue("Differential-Analysis_{comparison}/Gene_enrichment_analysis_in_modified_regions/Enrichment_{status}-regions_GO_Molecular_Function_2018.txt"),
                col.names = TRUE,
                row.names = FALSE,
                sep = "\t",
                quote = FALSE)
  }
  return(list(list.df.adv, list.df.del))
}


