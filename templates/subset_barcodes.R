# $1 barcodes file
# $2 key
# $5 number of cells for this hashtag

args <- commandArgs(trailingOnly = TRUE)
barcodes<-args[1]
key<-args[2]
n<-args[3]
hashtag<-args[4]

library(tidyverse)
library(magrittr)

bcs<-read_table(barcodes, col_names = "barcode")
print(length(bcs$barcode))
set.seed(1)
bcs_sub<-bcs$barcode[sample(seq_along(bcs$barcode),n)]
write_tsv(tibble(bcs_sub),file=paste("Hashtag",hashtag,"_", key,"_",n,"_sub.tsv",sep=""),col_names=FALSE)