# Manhattan plot for copy number data

setwd("T:\\Aplenc_New\\Projects\\CopyNumbers_Germline\\Data\\SNP_info")
source("T:\\Aplenc\\SampleSets\\AAML_0531\\GWAS\\sgpt\\qqman.r")

###################################
# IMPORT AML & Control CNV files  #
###################################
#(remove waves & telo/centromeres)

cnv  <- read.table("T:\\Aplenc_New\\Projects\\CopyNumbers_Germline\\Data\\cnvPartition\\array_all_clean.txt", header = T, sep = "\t", stringsAsFactors = F)
# cnv <- cnv[order(cnv$chr, cnv$start),]
chop <- read.table("T:\\Aplenc_New\\Projects\\CopyNumbers_Germline\\Data\\cnvPartition\\chop_clean.txt", header = T, sep = "\t", stringsAsFactors = F)

###########################
# IMPORT CHROMOSOME MAPS  #
###########################

for (i in 1:22)
{
	df <- read.table(paste("chr", i, "_template.txt", sep=""), header = T, sep = "\t", stringsAsFactors = F) ;
	colnames(df)[1] 	<- "CHR"
	colnames(df)[2] 	<- "BP"
	colnames(df)[3] 	<- "SNP"
	df$P			<- 0 ;
	df$Pgain		<- 0 ;
	df$LossCase 	<- 0 ;
	df$LossControl	<- 0 ;
	df$GainCase 	<- 0 ;
	df$GainControl 	<- 0 ;
	name <- paste("chr", i, sep = "") ;
	assign(name, df) ;
}

#########################
# AML CASE COUNT        #
#########################

/* you can do this more efficiently, not loop through entire chromosome, that's not necessary

for (i in 1:dim(cnv)[1])
{
	# first select the correct chromosome for internal count revision
	chr <- get(paste("chr", cnv$chr[i], sep=""))
	
	# loop through the entire chromosome
	for (j in cnv$start[i]:cnv$end[i])
	{
		if (chr$BP[j] >= cnv$start[i])
		{
			if (chr$BP[j] <= cnv$end[i])
			{
				if(cnv$value[i] < 2 )
				{
					chr$LossCase[j] <- chr$LossCase[j] + 1 ;
				}
				if(cnv$value[i] > 2 )
				{
					chr$GainCase[j] <- chr$GainCase[j] + 1 ;
				}	
			}
		}
	# save new internal count to external data frame
	assign(paste("chr", cnv$chr[i], sep=""), chr)
	}
}

#########################
# CHOP CONTROL COUNT    #
#########################

for (i in 1:dim(chop)[1])
{
	# first select the correct chromosome for internal count revision
	chr <- get(paste("chr", chop$chr[i], sep=""))
	print(i) ;

	# loop through the entire chromosome
	for (j in 1:dim(chr)[1])
	{
		if (chr$BP[j] >= chop$start_bp[i])
		{
			if (chr$BP[j] <= chop$end_bp[i])
			{
				if(chop$value[i] < 2 )
				{
					chr$LossControl[j] <- chr$LossControl[j] + 1 ;
				}
				if(chop$value[i] > 2 )
				{
					chr$GainControl[j] <- chr$GainControl[j] + 1 ;
				}	
			}
		}
	# save new internal count to external data frame
	assign(paste("chr", chop$chr[i], sep=""), chr)
	}
}

##############################
# ChiX + Manhattan + Export  #
##############################

for (i in 1:3)
{
	df <- get(paste("chr", i, sep=""))
	# Chi Square
	#for (j in 1:dim(df)[1])
	#{
	#	df$P[j] <- chisq.test(rbind(c(df$LossCase[j], length(unique(cnv$regno)) - df$LossCase[j]), c(df$LossControl[j], 2026 - df$LossControl[j])))$p.value ;
	#	df$Pgain[j] <- chisq.test(rbind(c(df$GainCase[j], length(unique(cnv$regno)) - df$GainCase[j]), c(df$GainControl[j], 2026 - df$GainControl[j])))$p.value ;
	#}
	# manhattan plot
	png(file = paste("chr", i, "_cnv_manhattan.png", sep = ""), width = 900, height = 350)
	manhattan(df, colors = c("black", "#666666", "#CC6600"), pch = 20, main = paste ("Germline CNV study: Chr ", i, sep = ""),  genomewideline = F, suggestiveline = F)
	dev.off()
	# qq-plot
	png(file = paste("chr", i, "_cnv_qq.png", sep = ""), width = 300, height = 300)
	qq(df$P)
	dev.off()
	# export table
 	# write.table(df, paste("chr", i, "_cnv_counts.txt", sep = ""), sep = "\t", row.names = F, col.names = T, quote = F) ;
	assign(paste("chr", i, sep = ""), df)
}

for (i in 1:22)
{
	df <- read.table(paste("chr", i, "_cnv_counts.txt", sep=""), header = T, sep = "\t", stringsAsFactors = F) ;
	png(file = paste("chr", i, "_cnv_manhattan.png", sep = ""), width = 900, height = 350)
	manhattan(df, colors = c("black", "#666666", "#CC6600"), pch = 20, main = paste ("Germline CNV study: Chr ", i, sep = ""),  genomewideline = F, suggestiveline = F)
	dev.off()
	png(file = paste("chr", i, "_cnv_qq.png", sep = ""), width = 300, height = 300)
	qq(df$P)
	dev.off()
}
