

### Extract the cell barcodes of interest, and generate 40 (customized n) barcodes in a text file to accelerate the speed of step 2


library(Seurat)
library(Signac)

# read
pbmc <- readRDS("signac_object.rds")

# output barcode
barcode <- data.frame(rownames(pbmc@meta.data), stringsAsFactors = F)
colnames(barcode) <- NULL

# barcode dir to save the outputs
folder.path <- "Barcode"
dir.create(path2<-file.path(folder.path), showWarnings = FALSE)

# bin 
cell.per.bin <- 40
bin_n <- nrow(barcode)/cell.per.bin

for (i in 1:ceiling(bin_n)) {
  if (i == ceiling(bin_n)) {
    tmp <- barcode[((i-1)*cell.per.bin+1):nrow(barcode), ]
    write.table(tmp, file = paste0("Barcode/BarcodeList_", i, ".txt"), row.names = F, sep = "", quote = F)
  } else {
    tmp <- barcode[((i-1)*cell.per.bin+1):(cell.per.bin*i), ]
    write.table(tmp, file = paste0("Barcode/BarcodeList_", i, ".txt"), row.names = F, sep = "", quote = F)
  }
}

