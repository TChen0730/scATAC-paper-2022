


#### INSTRUCTIONS

## This function takes two granges list objects (grlist1, grlist2): 
# (1) To test the 15-state enrichment for a set of peaks: one list for 15-state model (a list with 15 elements, one element for one state), the other list for cluster-specific peaks (each element in the list correspond to the DA peaks for each cluster).
# (2) To test the 15-state enrichment for each peak: one list for 15-state model (the same as (1)), the other list for list of peaks (each peak is a list element).

## chrom.size.dir: Should include file name! File was downloaded from https://hgdownload-test.gi.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes


#### ENRICHMENT FUNCTION
peak_enrich <- function(grlist1, grlist2, chrom.size.dir = chrom.size.dir) {
  
  # calculate total bases in whole genome 
  chrom.size <- read.table(chrom.size.dir, stringsAsFactors = F)
  chrom.size <- chrom.size[1:24, ]
  whole.genome.bases <- sum(chrom.size[, 2])
  
  # output
  output <- matrix(0, nrow = length(grlist1)*length(grlist2), ncol = 5)
  colnames(output) <- c("grList1", "grList2", "FE", "pvalue", "overlap_bases")
  
  # iterator
  output.i <- 1
  
  for (i in seq_along(grlist1)) {
    print(i)
    for (j in seq_along(grlist2)) {
      # test over-representation
      x <- GenomicRanges::intersect(grlist1[[i]], grlist2[[j]]) %>% width() %>% sum() %>% as.numeric() # hitInSample 
      n <- as.numeric(sum(width(grlist2[[j]])))    # hitInPop
      m <- as.numeric(whole.genome.bases - sum(width(grlist2[[j]])))    # failInPop
      k <- as.numeric(sum(width(grlist1[[i]])))      # sampleSize
      
      # fold enrichment
      FE <- as.numeric(x*(m+n)/(k*n))
      
      # p value
      enr.p <- phyper(x-1,n,m,k,lower.tail = F)
      
      # save output
      output[output.i, ] <- c(names(grlist1)[i], names(grlist2)[j], FE, enr.p, x)
      output.i <- output.i + 1
    }
  }
  
  # add adjusted p
  output <- data.frame(output, stringsAsFactors = F)
  output$adj.p <- p.adjust(as.numeric(output$pvalue), method = "BH")
  
  return(output)
}