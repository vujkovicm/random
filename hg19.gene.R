# http://bioconductor.org/packages/release/bioc/vignettes/RTCGAToolbox/inst/doc/RTCGAToolbox-vignette.html
# source("https://bioconductor.org/biocLite.R")
# biocLite("RTCGA")
# biocLite("RTCGAToolbox")

library("Homo.sapiens")
locations <- genes(Homo.sapiens, columns = "SYMBOL")
locations <- as.data.frame(locations)
locations <- subset(locations, is.na(locations$SYMBOL) == F)
locations$chr <- gsub("chr", "", locations$seqnames)
locations$skat <- paste0(locations$SYMBOL, " ", locations$chr, ":", locations$start, "-", locations$end)

write.table(locations$skat, "hg19.ucsc.gene.set", row.names = F, col.names = F, quote = F, sep = "")

# rvtest --inVcf ../data/hg19/ceu_hg19.vcf.gz \
# --kernel skat \
# --setFile ../hg19.ucsc.gene.set \
# --pheno ../data/phe.sample \
# --covar ../data/phe.sample \
# --pheno-name phe \
# --covar-name sex,age,bmi,P1 \
# --out ../output/phe
# ### MAF < 0.1
