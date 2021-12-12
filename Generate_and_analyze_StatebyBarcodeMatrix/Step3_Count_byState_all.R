
# packages
library(utils)
library(stats)
library(Seurat)
library(Signac)
library(tidyr)
library(dplyr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# list files
files <- list.files(path = "../FragmentsbyCell/", pattern = "\\.txt$")

# prepare 15-state granges
## specify the epigenome
epi.genome <- "E124"
state15 <- read.table(paste0("~/4_scATACseq_AnnotationFile/Roadmap_15StateModel/", epi.genome, "_15_coreMarks_hg38lift_mnemonics.bed"), stringsAsFactors = F)
state15 <- state15[state15$V1 != "chrM", ]

# convert 15 state model (in peaks) into list of granges
state <- unique(sort(state15[, 4]))
state15.gr.list <- list()
for (i in seq_along(state)) {
  tmp <- state15[state15[, 4] == state[i], ]
  state15.gr.list[[i]] <- GRanges(tmp[, 1], IRanges(tmp[, 2], tmp[, 3]))
}
names(state15.gr.list) <- state

# output matrix
res <- matrix(0, nrow = 15, ncol = length(files))
rownames(res) <- names(state15.gr.list)
colnames(res) <- gsub("\\.txt", "", files)

# loop for files: this step will take long (hours to 1day) for >10k cells
for (n in seq_along(files)) {
  
  # cut site ranges
  fragments <- read.table(paste0("../FragmentsbyCell/", files[n]), stringsAsFactors = F)
  colnames(fragments) <- c("chr", "start", "end", "barcode", "duplicates")
  
  # prepare cut site granges
  fragments.long <- pivot_longer(fragments, cols = c("start", "end"), names_to = "cutsite", values_to = "position")
  cutsite.gr <- GRanges(fragments.long$chr, IRanges(fragments.long$position, fragments.long$position))
  
  # count cutsite for each state
  res[, n] <- sapply(state15.gr.list, function(x) GenomicRanges::intersect(cutsite.gr, x) %>% length())
  
  print(n) # check how many cells have been processed
  
}

# save
saveRDS(res, file = "CutSites_count_byState_all.rds")



