####################
#      REPTIM      #
# Alexandre BOULLE #
####################

# Methods of analysis for replication timing profiles

#' @import stringr
#' @importFrom stats na.omit p.adjust quantile wilcox.test
#' @importFrom utils read.table write.table

#' @title Load data
#' @description recover values of replication timing profile for each condition of interest
#' @param name.experiment the name of a condition
#' @export
load.data <- function(name.experiment){
  cond <- read.table(str_glue("Bedgraph_files/{name.experiment}.bedGraph", header = FALSE))
  cond.loess <- read.table(str_glue("Bedgraph_files/{name.experiment}_Loess.bedGraph", header = FALSE))
  colnames(cond) <- c("CHR", "START", "END", "Intensity")
  colnames(cond.loess) <- c("CHR", "START", "END", "Intensity")
  return(list(cond, cond.loess))
}

#' @title Merge samples
#' @description Merge replicates for a condition of interest
#' @param list.replicates a list containing dataframes of replicate samples
#' @export
replicate.merging <- function(list.replicates){
  for (i in 1:(length(list.replicates))){
    if (i == 1){
      intensities <- list.replicates[[i]]$Intensity
    }else{
      intensities <- cbind(intensities, list.replicates[[i]]$Intensity)
    }
  }
  condition <- data.frame(cbind(list.replicates[[1]]$CHR, list.replicates[[1]]$START, list.replicates[[1]]$END, apply(intensities, 1, mean)))
  condition[, 2:4] <- lapply(condition[, 2:4], as.numeric)
  colnames(condition) <- c("CHR", "START", "END", "Intensity")
  return(condition)
}

#' @param condition1 as a condition of interest
#' @param condition2 as the other condition of interest
euclidean.distance <- function(condition1, condition2){
  # Calculate distance between curves
  dist.euclid <- sqrt((condition1[, 2] - condition2[, 2])**2 + (condition1[, 4] - condition2[, 4])**2)
  return(dist.euclid)
}

#' @param list.conditions as the list of conditions in order to calculate the distance between different conditions
euclidean.distance.multiple.conditions <- function(list.conditions){
  dist.euclid <- c()
  for (i in 1:length(list.conditions)){
    for (j in i:length(list.conditions)){
      if (j < length(list.conditions)){
        dist.euclid <- c(dist.euclid, euclidean.distance(condition1 = list.conditions[[i]], condition2 = list.conditions[[j+1]]))
      }
    }
  }
  return(dist.euclid)
}

#' @param condition1 as the first condition
#' @param condition2 as the second condition
#' @param threshold.pvalue as the pvalue threshold
differential.analysis.1 <- function(condition1, condition2, threshold.pvalue){
  condition1$pvalue <- NA
  min.pos <- c()
  max.pos <- c()
  chr <- c()
  k <- 1
  for (i in 1:nrow(condition1)){
    if(i >= 2 && abs(condition1[i, 2] - condition1[(i - 1), 2]) > 50*10**3){
      # Keep all positions of a region where the distance is greater than the threshold
      subset <- condition1[k:(i-1), ]
      # Find all positions for which the p-value is lower than the threshold
      ind.pval <- which(subset$pvalue < threshold.pvalue)
      # Keep the chromosome
      chr <- c(chr, subset[1, 1])
      # Keep the first position of the region of interest
      min.pos <- c(min.pos, subset[1, 2])
      # Keep the last position for which p-value is lower than the threshold
      # If p-value is greater than the threshold then NA appears
      max.pos <- c(max.pos, subset[max(ind.pval), 2])
      # The position of k corresponds to the beginning of the new region
      k <- i
    }
    # Test if the difference is significant
    test <- wilcox.test(condition1[k:i, 4], condition2[k:i, 4], alternative = "two.sided")
    condition1$pvalue[i] <- test$p.value
  }
  # Check positions and p-values in the table
  # write.table(condition1, "condition1_pvalue.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
  # Create a dataframe of the results
  df.differential <- data.frame(CHR = chr, START = min.pos, END = max.pos)
  df.differential <- na.omit(df.differential)
  ind.non.zero <- which((df.differential$START - df.differential$END) != 0)
  df.differential <- df.differential[ind.non.zero, ]
  return(df.differential)
}

#' @param condition1
#' @param condition2
#' @param df.differential.1 as the result of the first differential analysis step
#' @param threshold.pvalue as the p-value threshold
#' @param threshold.dist as the distance threshold to perform an elongation of different regions
differential.analysis.2 <- function(condition1, condition2, df.differential.1, threshold.pvalue, threshold.dist){
  min.pos <- c()
  max.pos <- c()
  chr.list <- c()
  dist.euclid <- euclidean.distance(condition1, condition2)
  for (i in 1:nrow(df.differential.1)){
    start <- df.differential.1[i, "START"]
    end <- df.differential.1[i, "END"]
    chr <- df.differential.1[i, "CHR"]
    ind.chr <- which(condition1[, 1] == chr)
    dist.chr <- dist.euclid[ind.chr]
    subset.cond1 <- condition1[condition1[, 1] == chr, ]
    subset.cond2 <- condition2[condition2[, 1] == chr, ]
    ind.start <- which(subset.cond1[, 2] == start)
    ind.end <- which(subset.cond1[, 2] == end)
    sign.diff.start <- sign(subset.cond1[ind.start, "Intensity"] - subset.cond2[ind.start, "Intensity"])
    sign.diff.end <- sign(subset.cond1[ind.end, "Intensity"] - subset.cond2[ind.end, "Intensity"])
    test <- wilcox.test(subset.cond1[ind.start:ind.end, 4], subset.cond2[ind.start:ind.end, 4], alternative = "two.sided")
    pval <- test$p.value
    stop.start <- 0
    stop.end <- 0
    while((pval < threshold.pvalue) & (stop.start == 0) | (pval < threshold.pvalue) & (stop.end == 0)){
      # If the distance between points is greater than the threshold (median, mean, deciles, ...) then we continue
      # Move to the left
      if (ind.start-1 >= 1){
        sign.diff.decrease <- sign(subset.cond1[(ind.start - 1), "Intensity"] - subset.cond2[(ind.start - 1), "Intensity"])
        if (dist.chr[ind.start - 1] >= threshold.dist & sign.diff.decrease == sign.diff.start){
          ind.start <- ind.start - 1
        }else{
          stop.start <- stop.start + 1
        }
      }else{
        stop.start <- stop.start + 1
      }
      # Move to the right
      if (ind.end+1 <= nrow(subset.cond1)){
        sign.diff.increase <- sign(subset.cond1[(ind.end + 1), "Intensity"] - subset.cond2[(ind.end + 1), "Intensity"])
        if (dist.chr[ind.end + 1] >= threshold.dist & sign.diff.increase == sign.diff.end){
          ind.end <- ind.end + 1
        }else{
          stop.end <- stop.end + 1
        }
      }else{
        stop.end <- stop.end + 1
      }
      # If the distance between points is greater than the threshold (median, mean, deciles, ...)
      # then we can perform a Wilcoxon test to know the significance of the difference
      if (stop.start == 0 | stop.end == 0){
        test <- wilcox.test(subset.cond1[ind.start:ind.end, 4], subset.cond2[ind.start:ind.end, 4], alternative = "two.sided")
        pval <- test$p.value
      }
    }
    # Recover the genomic coordinates
    chr.list <- c(chr.list, chr)
    if (pval >= threshold.pvalue & stop.start == 0 & stop.end != 0){
      min.pos <- c(min.pos, subset.cond1[(ind.start + 1), 2])
      max.pos <- c(max.pos, subset.cond1[ind.end, 2])
    }else if (pval >= threshold.pvalue & stop.start != 0 & stop.end == 0){
      min.pos <- c(min.pos, subset.cond1[ind.start, 2])
      max.pos <- c(max.pos, subset.cond1[(ind.end - 1), 2])
    }else if (pval >= threshold.pvalue & stop.start == 0 & stop.end == 0){
      min.pos <- c(min.pos, subset.cond1[(ind.start + 1), 2])
      max.pos <- c(max.pos, subset.cond1[(ind.end - 1), 2])
    }else {
      min.pos <- c(min.pos, subset.cond1[ind.start, 2])
      max.pos <- c(max.pos, subset.cond1[ind.end, 2])
    }
  }
  # Create a table with the new results
  df.differential.2 <- data.frame(CHR = chr.list, START = min.pos, END = max.pos)
  df.differential.2 <- na.omit(df.differential.2)
  ind.non.zero <- which((df.differential.2$START - df.differential.2$END) != 0)
  df.differential.2 <- df.differential.2[ind.non.zero, ]
  return(df.differential.2)
}

#' @param condition1 as the first condition
#' @param condition2 as the second condition
#' @param df.differential as the dataframe containing the positions of modified regions
#' @param threshold.pvalue as the p-value threshold
add.status.differential <- function(condition1, condition2, df.differential, threshold.pvalue){
  df.differential$STATUS <- ""
  df.differential$PVALUE <- ""
  line.remove <- c()
  for (i in 1:nrow(df.differential)){
    chr <- df.differential$CHR[i]
    start <- df.differential$START[i]
    end <- df.differential$END[i]
    subset.chr <- condition1[condition1$CHR == chr, ]
    ind.start <- which(subset.chr$START == start)
    ind.end <- which(subset.chr$END == end)
    test.less <- wilcox.test(condition1[condition1$CHR == chr, "Intensity"][ind.start:ind.end],
                             condition2[condition2$CHR == chr, "Intensity"][ind.start:ind.end],
                             alternative = "less")
    test.greater <- wilcox.test(condition1[condition1$CHR == chr, "Intensity"][ind.start:ind.end],
                                condition2[condition2$CHR == chr, "Intensity"][ind.start:ind.end],
                                alternative = "greater")
    if (test.less$p.value <= threshold.pvalue){
      df.differential$STATUS[i] <- "ADVANCED"
      df.differential$PVALUE[i] <- test.less$p.value
    }else if (test.greater$p.value <= threshold.pvalue){
      df.differential$STATUS[i] <- "DELAYED"
      df.differential$PVALUE[i] <- test.greater$p.value
    }else if (test.greater$p.value > threshold.pvalue){
      line.remove <- c(line.remove, i)
    }else if (test.greater$p.value > threshold.pvalue){
      line.remove <- c(line.remove, i)
    }
  }
  if (length(line.remove) > 0){
    df.differential <- df.differential[-line.remove, ]
  }
  return(df.differential)
}

#' @param df.differential as the dataframe containing the different regions
#' @param i as the indice corresponding to a line in the dataframe
#' @param chr as the chromosome
#' @param start as the start position of modified region
#' @param end as the end position of modified region
#' @param status as the status ("advanced" or delayed") of the modified region
replace.positions <- function(df.differential, i, chr, start, end, status){
  ind1 <- which(df.differential$CHR == chr & df.differential$START >= start & df.differential$END >= end & df.differential$START <= end & df.differential$STATUS == status)
  ind2 <- which(df.differential$CHR == chr & df.differential$START <= start & df.differential$END <= end & df.differential$END >= start & df.differential$STATUS == status)
  ind3 <- which(df.differential$CHR == chr & df.differential$START <= start & df.differential$END >= end & df.differential$STATUS == status)
  # Recover "start" and "end" positions from regions with an intersection
  start.overlap <- c(df.differential[ind1, "START"], df.differential[ind2, "START"], df.differential[ind3, "START"])
  end.overlap <- c(df.differential[ind1, "END"], df.differential[ind2, "END"], df.differential[ind3, "END"])
  # Keep the smallest "start" and the bigger "end"
  min.pos <- min(c(start, start.overlap))
  max.pos <- max(c(end, end.overlap))
  # Remove the smallest regions  = regions included in the larger region (min.pos to max.pos)
  ind <- unique(ind1, ind2, ind3)
  # Avoid to remove the current indice
  df.differential[ind[ind != i], "START"] <- 0
  df.differential[ind[ind != i], "END"] <- 0
  # Keep the largest region : replace "start" and "end" positions with the new values
  df.differential[i, "START"] <- min.pos
  df.differential[i, "END"] <- max.pos
  return(list(df.differential, min.pos, max.pos))
}

#' @param df.differential as the dataframe containing the modified regions linked with a status
fusion.common.regions <- function(df.differential){
  for (i in 1:nrow(df.differential)){
    chr  <- df.differential[i, "CHR"]
    start <- df.differential[i, "START"]
    end <- df.differential[i, "END"]
    status <- df.differential[i, "STATUS"]
    df.differential <- replace.positions(df.differential, i, chr, start, end, status)[[1]]
    min.pos <- replace.positions(df.differential, i, chr, start, end, status)[[2]]
    max.pos <- replace.positions(df.differential, i, chr, start, end, status)[[3]]
    # If "min.pos" and "max.pos" values are different of "start" and "end" positions then we continue because it could exist other intersections
    while (min.pos != start | max.pos != end){
      start <- min.pos
      end <- max.pos
      df.differential <- replace.positions(df.differential, i, chr, start, end, status)[[1]]
      min.pos <- replace.positions(df.differential, i, chr, start, end, status)[[2]]
      max.pos <- replace.positions(df.differential, i, chr, start, end, status)[[3]]
    }
  }
  # Keep only lines of interest (= positions other than 0)
  df.differential <- df.differential[df.differential$END != 0, ]
  return(df.differential)
}

#' @param cond1 as the first condition
#' @param cond2 as the second condition
#' @param pvalue as the p-value threshold
#' @param indice.decile as the positions for which we have a sufficient distance between conditions
#' @param min.dist.euclid as the distance to stop the elongation of regions
differential.analysis <- function(cond1, cond2, pvalue, indice.decile, min.dist.euclid, name){
  # STEP 1 : detect regions with differences between curves
  df.diff.1 <- differential.analysis.1(condition1 = cond1[indice.decile, ], condition2 = cond2[indice.decile, ], threshold.pvalue = pvalue)
  # percentage.1 <- percentage.difference(condition = cond1, df.differential = df.diff.1)
  # Check if different regions were detected
  if (dim(df.diff.1)[1] > 0){
    # STEP 2 : increase the size of the detected regions
    df.diff.2 <- differential.analysis.2(condition1 = cond1, condition2 = cond2, df.differential.1 = df.diff.1, threshold.pvalue = pvalue, threshold.dist = min.dist.euclid)
    # percentage.2 <- percentage.difference(condition = cond1, df.differential = df.diff.2)
    # Check if it remains different regions
    if (dim(df.diff.2)[1] > 0){
      # STEP 3 : add "ADVANCED" or "DELAYED" annotation
      df.diff.3 <- add.status.differential(condition1 = cond1, condition2 = cond2, df.differential = df.diff.2, threshold.pvalue = pvalue)
      # percentage.3 <- percentage.difference(condition = cond1, df.differential = df.diff.3)
      # Check if it remains different regions
      if (dim(df.diff.3)[1] > 0){
        # Fusion of regions sharing an intersection
        df.diff.unique <- fusion.common.regions(df.differential = df.diff.3)
        # Update p-values
        df.diff.unique.pval <- add.status.differential(condition1 = cond1, condition2 = cond2, df.differential = df.diff.unique[, 1:3], threshold.pvalue = pvalue)
        # Adjust p-values with FDR method
        df.diff.unique.pval$ADJUSTED_PVALUE <- p.adjust(df.diff.unique.pval$PVALUE, method = "fdr")
        df.diff.final <- df.diff.unique.pval[df.diff.unique.pval$ADJUSTED_PVALUE < pvalue, ]
        percentage.final <- percentage.difference(condition = cond1, df.differential = df.diff.unique)
        dir.create(str_glue("Differential-Analysis_{comparison}/Modified_regions/"))
        write.table(df.diff.final, str_glue("Differential-Analysis_{name}/Modified_regions/Modified_regions_fdr-padj.txt"), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t", append = T)
        write.table(percentage.final, str_glue("Differential-Analysis_{name}/Modified_regions/Differential_percentage.txt"), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t", append = T)
        return(list(df.diff.final, percentage.final))
      }else{
        return("No differences between conditions")
      }
    }else{
      return("No differences between conditions")
    }
  }else{
    return("No differences between conditions")
  }
}

#' @param condition as a condition of the study allowing to recover chromosome positions
#' @param df.differential as a dataframe containing positions of modified regions
percentage.difference <- function(condition, df.differential){
  size.genome <- 0
  chrs <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",
            "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
            "chr20", "chr21", "chr22", "chrX")
  df.percent <- data.frame(Chromosome = character(length(chrs) + 1), Percentage = numeric(length(chrs) + 1))
  ind <- 1
  for (chr in chrs){
    cond.chr <- condition[condition$CHR == chr, ]
    # Chromosome size corresponds to the last position
    size.chr <- cond.chr$END[nrow(cond.chr)]
    # Add each chromosome size to calculate the genome size
    size.genome <- size.genome + size.chr
    percent.diff.chr <- (sum(df.differential[df.differential$CHR == chr, "END"] - df.differential[df.differential$CHR == chr, "START"]) / size.chr) * 100
    df.percent$Chromosome[ind] <- chr
    df.percent$Percentage[ind] <- percent.diff.chr
    ind <- ind + 1
  }
  percent.diff.all <- (sum(df.differential$END - df.differential$START) / size.genome) * 100
  df.percent$Chromosome[ind] <- "Total"
  df.percent$Percentage[ind] <- percent.diff.all
  return(df.percent)
}

#' @param differential.1 the dataframe containing results of the first differential analysis performed between conditions
#' @param differential.2 the dataframe containing results of the second differential analysis performed between conditions
size.common.region <- function(differential.1, differential.2){
  sum.total <- 0
  nbr.regions <- 0
  ind.not.find <- c()
  for (i in 1:nrow(differential.1)){
    sum.case.0 <- 0
    sum.case.1 <- 0
    sum.case.2 <- 0
    sum.case.3 <- 0
    sum.case.4 <- 0
    chr <- differential.1$CHR[i]
    start <- differential.1$START[i]
    end <- differential.1$END[i]
    status <- differential.1$STATUS[i]
    # Case 0 : regions have the sames start and end positions
    ind <- which(differential.2$CHR == chr & differential.2[, "START"] == start & differential.2[, "END"] == end & differential.2$STATUS == status)
    if (length(ind) > 0){
      sum.case.0 <- sum(abs(end - start))
      # Other case : regions are different but we test if it exists an intersection
    }else{
      ind1 <- which(differential.2$CHR == chr & differential.2$START >= start & differential.2$END >= end & differential.2$START <= end & differential.2$STATUS == status)
      if (length(ind1) > 0){
        sum.case.1 <- sum(abs(end - differential.2[ind1, "START"]))
      }
      ind2 <- which(differential.2$CHR == chr & differential.2$START <= start & differential.2$END <= end & differential.2$END >= start & differential.2$STATUS == status)
      if (length(ind2) > 0){
        sum.case.2 <- sum(abs(differential.2[ind2, "END"] - start))
      }
      ind3 <- which(differential.2$CHR == chr & differential.2$START <= start & differential.2$END >= end & differential.2$STATUS == status)
      if (length(ind3) > 0){
        sum.case.3 <- sum(abs(end - start))
      }
      ind4 <- which(differential.2$CHR == chr & differential.2$START >= start & differential.2$END <= end & differential.2$STATUS == status)
      if (length(ind4) > 0){
        sum.case.4 <- sum(abs(differential.2[ind4, "END"] - differential.2[ind4, "START"]))
      }
    }
    if ((sum.case.0 + sum.case.1 + sum.case.2 + sum.case.3 + sum.case.4) > 0){
      nbr.regions <- nbr.regions + 1
    }else{
      ind.not.find <- c(ind.not.find, i)
    }
    sum.total <- sum.total + sum.case.0 + sum.case.1 + sum.case.2 + sum.case.3 + sum.case.4
  }
  return(list(sum.total, nbr.regions, ind.not.find))
}

#' @param condition1 the first condition (values from quantile normalization)
#' @param condition2 the second condition (values from quantile normalization)
#' @param condition1.loess the first condition (values from loess smoothing)
#' @param condition2.loess the second condition (values from loess smoothing)
#' @param df.differential as the result (dataframe) of differential analysis performed between two conditions
#' @param name as the name of the comparison
files.startr.viewer <- function(condition1, condition2, condition1.loess, condition2.loess, df.differential, name){
  dir.create(str_glue("Differential-Analysis_{name}/Files_START-R_Viewer/"))
  chrs <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",
            "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
            "chr20", "chr21", "chr22", "chrX")
  for (chr in chrs){
    cond1.chr <- condition1[condition1$CHR == chr, ]
    cond2.chr <- condition2[condition2$CHR == chr, ]
    positions <- cond1.chr$START
    intensities.cond1 <- cond1.chr$Intensity
    intensities.cond2 <- cond2.chr$Intensity
    cond1.loess.chr <- condition1.loess[condition1.loess$CHR == chr, ]
    cond2.loess.chr <- condition2.loess[condition2.loess$CHR == chr, ]
    loess.value1 <- cond1.loess.chr$Intensity
    loess.value2 <- cond2.loess.chr$Intensity
    # Cond1
    df.chr.cond1 <- data.frame(Position = positions, Intensity = intensities.cond1, Name = rep("Exp 1", nrow(cond1.chr)))
    # Cond1 Loess
    df.chr.cond1.loess <- data.frame(Position = positions, Intensity = loess.value1, Name = rep("Smooth Exp 1", nrow(cond1.chr)))
    # Cond2
    df.chr.cond2 <- data.frame(Position = positions, Intensity = intensities.cond2, Name = rep("Exp 2", nrow(cond2.chr)))
    # Cond2 Loess
    df.chr.cond2.loess <- data.frame(Position = positions, Intensity = loess.value2, Name = rep("Smooth Exp 2", nrow(cond2.chr)))
    df.chr <- rbind(df.chr.cond1, df.chr.cond2, df.chr.cond1.loess, df.chr.cond2.loess)
    # min.int is the ordinate value where the annotation ("Advanced" or "Delayed" bar) is positioned
    min.int <- round(min(intensities.cond1, intensities.cond2) - 0.5)
    df.diff.chr <- df.differential[df.differential$CHR == chr, ]
    # Add positions where we have detected differences between profiles
    if (dim(df.diff.chr)[1] > 0){
      for (j in 1:nrow(df.diff.chr)){
        if (df.diff.chr$STATUS[j] == "ADVANCED"){
          status <- "Advanced"
        }else {
          status <- "Delayed"
        }
        df.chr <- rbind(df.chr, c(df.diff.chr$START[j], min.int, status))
        df.chr <- rbind(df.chr, c(df.diff.chr$END[j], min.int, status))
        df.chr <- rbind(df.chr, c(NA, NA, status))
      }
    }
    # Add a line before colnames to allow reading in START-R Viewer
    cat(str_glue("loess	NA	0	NA	hg38	Loess	{chr}"), file = str_glue("Differential-Analysis_{name}/Files_START-R_Viewer/{chr}.SRV"), sep = "\n")
    write.table(df.chr, str_glue("Differential-Analysis_{name}/Files_START-R_Viewer/{chr}.SRV"), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t", append = T)
  }
}

#' @title Differential analysis of replication timing profiles
#' @description Determine modified regions between conditions of interest
#' @param list.cond a list of all conditions of the study
#' @param cond1 the first condition (values from quantile normalization)
#' @param cond2 the second condition (values from quantile normalization)
#' @param cond1.loess the first condition (values from loess smoothing)
#' @param cond2.loess the second condition (values from loess smoothing)
#' @param pval the p-value threshold
#' @param max.q.dist the distance at which potential modified regions are detected
#' @param min.q.dist the distance at which the elongation of the regions stops
#' @param comparison the name of comparison to build result folder
#' @export
modified.regions.detection <- function(list.cond, cond1, cond2, cond1.loess, cond2.loess, pval, max.q.dist, min.q.dist, comparison){
  dir.create(str_glue("Differential-Analysis_{comparison}/"))
  all.dist <- euclidean.distance.multiple.conditions(list.cond)
  dist.cond1.vs.cond2 <- euclidean.distance(condition1 = cond1, condition2 = cond2)
  # ind.decile <- which(dist.cond1.vs.cond2 >= quantile(dist.cond1.vs.cond2, probs = seq(0.1, 1, by = 0.05))[str_glue("{max.q.dist}%")])
  ind.decile <- which(dist.cond1.vs.cond2 >= quantile(all.dist, probs = seq(0, 1, by = 0.05))[str_glue("{max.q.dist}%")])
  min.dist <- as.numeric(quantile(dist.cond1.vs.cond2, probs = seq(0.1, 1, by = 0.05))[str_glue("{min.q.dist}%")])
  # Results of differential analysis performed between profiles
  results.diff <- differential.analysis(cond1 = cond1, cond2 = cond2,
                                        pvalue = pval,
                                        indice.decile = ind.decile,
                                        min.dist.euclid = min.dist,
                                        name = comparison)
  # Files to obtain profiles
  if (is.list(results.diff)){
    files.startr.viewer(condition1 = cond1, condition2 = cond2,
                        condition1.loess = cond1.loess, condition2.loess = cond2.loess,
                        df.differential = results.diff[[1]], name = comparison)
  }
  return(results.diff)
}




