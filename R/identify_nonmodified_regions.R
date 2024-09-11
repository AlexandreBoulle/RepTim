####################
#      REPTIM      #
# Alexandre BOULLE #
####################


#############################
# KEEP NON-MODIFIED REGIONS #
#############################


#' @import stringr

#' @param differential.result as the dataframe containing modified regions detection after the comparison of two conditions
#' @param size.chr.table as the dataframe containing the size of chromosomes
#' @param name as the comparison name used in the result folder
find.non.modified.regions <- function(differential.result, size.chr.table, name){
  all.chr <- c()
  start <- c()
  end <- c()
  for (i in 1:dim(size.chr.table)[1]){
    chr = size.chr.table$Chromosome[i]
    # region.chr <- differential.result %>% filter(CHR == chr)
    region.chr <- differential.result[differential.result$CHR == chr, ]
    region.chr <- region.chr[order(region.chr$START), ]
    # print(region.chr)
    # Check if it exists different regions for the chromosome of interest
    if (dim(region.chr)[1] > 0){
      for (j in 1:dim(region.chr)[1]){
        # Find the first non-modified region
        if (j == 1 & region.chr$START[j] > 0){
          all.chr = c(all.chr, region.chr$CHR[j])
          start <- c(start, 0)
          end <- c(end, (region.chr$START[j] - 1))
        }
        # Find the non-modified regions between different regions
        if ((j > 1) && (region.chr$END[j-1] - region.chr$START[j] != 0)){
          all.chr = c(all.chr, region.chr$CHR[j])
          start <- c(start, region.chr$END[j-1] + 1)
          end <- c(end, region.chr$START[j] - 1)
        }
        # Find the last non-modified region
        if (j == dim(region.chr)[1] & region.chr$END[j] < size.chr.table$Seq_length[i]){
          all.chr = c(all.chr, region.chr$CHR[j])
          start = c(start, region.chr$END[j] + 1)
          end = c(end, size.chr.table$Seq_length[i])
        }
      }
    }
  }
  df.non.modified <- data.frame(CHR = all.chr, START = start, END = end, STATUS = rep("non_modified", length(all.chr)))
  write.table(df.non.modified, str_glue("Differential-Analysis_{name}/Non-modified_regions/Non-modified_regions.txt"), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
  return(df.non.modified)
}

#' @param res.adv as the dataframe containing only the advanced regions
#' @param res.del as the dataframe containing only the delayed regions
#' @param df.non.modified as the dataframe containing all non-modified regions
#' @param name as the comparison name used in the result folder
reduce.size.non.modified.regions <- function(res.adv, res.del, df.non.modified, comparison){
  list.df <- list()
  for (status in c("advanced", "delayed")){
    if (status == "advanced"){
      differential.result <- res.adv
    }else{
      differential.result <- res.del
    }
    columns <- c("CHR", "START", "END", "STATUS")
    df.non.modified.reduce <- data.frame(matrix(nrow = 0, ncol = length(columns)))
    colnames(df.non.modified.reduce) <- columns
    all.chr.df <- unique(differential.result$CHR)
    for (chr in all.chr.df){
      # Consider a chromosome
      res.chr <- differential.result[differential.result$CHR == chr, ]
      nbr.regions <- dim(res.chr)[1]
      df.nm.chr <- df.non.modified[df.non.modified$CHR == chr, ]
      df.nm.chr$SIZE <- abs(df.nm.chr$END - df.nm.chr$START)
      df.nm.chr.sort <- df.nm.chr[order(df.nm.chr$SIZE, decreasing = TRUE), ]
      # Keep the same non-modified region number as modified regions
      df.nm.select.regions <- df.nm.chr.sort[1:nbr.regions, ]
      df.nm.sort <- df.nm.select.regions[order(df.nm.select.regions$START, decreasing = FALSE), columns]
      # Reduce the size of non-modified regions to have an equivalent size between non-modified and modified regions
      percentage <- sum(res.chr$END - res.chr$START) / sum(df.nm.sort$END - df.nm.sort$START)
      remove.bp <- round( ( (df.nm.sort$END - df.nm.sort$START) - ((df.nm.sort$END - df.nm.sort$START) * percentage) ) / 2 )
      df.nm.subset <- data.frame(CHR = df.nm.sort$CHR, START = (df.nm.sort$START + remove.bp), END = (df.nm.sort$END - remove.bp), STATUS = rep("non_modified", length(df.nm.sort$CHR)))
      df.non.modified.reduce <- rbind(df.non.modified.reduce, df.nm.subset)
    }
    print(status)
    print("Modified regions (bp number)")
    print(sum(differential.result$END - differential.result$START))
    print("Non-modified regions (bp number)")
    print(sum(df.non.modified.reduce$END - df.non.modified.reduce$START))
    list.df[[length(list.df) + 1]] <- df.non.modified.reduce
    write.table(df.non.modified.reduce, str_glue("Differential-Analysis_{comparison}/Non-modified_regions/Reduced_non-modified_regions_for-{status}.txt"), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
  }
  return(list.df)
}

#' @title Non-modified regions generation
#' @description Method to recover non-modified regions of the same size as modified regions
#' @param differential.result the dataframe containing modified regions between conditions of interest
#' @param size.chr.table the dataframe containing the size of chromosomes
#' @param comparison the comparison name used in the result folder
#' @export
non.modified.regions <- function(differential.result, size.chr.table, comparison){
  dir.create(str_glue("Differential-Analysis_{comparison}/Non-modified_regions/"))
  res.adv <- differential.result[differential.result$STATUS == "ADVANCED", ]
  res.del <- differential.result[differential.result$STATUS == "DELAYED", ]
  write.table(res.adv, str_glue("Differential-Analysis_{comparison}/Modified_regions/advanced_regions.txt"), col.names = TRUE, row.names= FALSE, quote = FALSE, sep = "\t")
  write.table(res.del, str_glue("Differential-Analysis_{comparison}/Modified_regions/delayed_regions.txt"), col.names = TRUE, row.names= FALSE, quote = FALSE, sep = "\t")
  df.non.modified <- find.non.modified.regions(differential.result, size.chr.table, comparison)
  df.non.modified.reduce <- reduce.size.non.modified.regions(res.adv, res.del, df.non.modified, comparison)
  return(df.non.modified.reduce)
}




