# source("https://bioconductor.org/biocLite.R")
# biocLite("RTCGA")
# biocLite("RTCGAToolbox")

library("Homo.sapiens")
locations <- genes(Homo.sapiens, columns = "SYMBOL")
locations <- as.data.frame(locations)
locations = subset(locations, is.na(locations$SYMBOL) == F)
locations$chr = gsub("chr", "", locations$seqnames)
locations$skat <- paste0(locations$SYMBOL, " ", locations$chr, ":", locations$start, "-", locations$end)

write.table(locations$skat, "hg19.ucsc.gene.set", row.names = F, col.names = F, quote = F, sep = "")

