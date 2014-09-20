setwd("T:\\my_directory")

# Chromosome seperated SNP reference file format
#
# Chr	Position
# 1	82154
# 1	534247
# 1	565374
# 1	569624
# 1	723918
# 1	729632

# cnvPartition output [clean of LOH and telomeric CNVs]
#
# regno	chr	start		end	      	size	value
# 1	5	49447784	49548900	101116	1
# 1	12	37999925	38410357	410432	1
# 2	17	34443830	34815551	371721	3
# 3	8	24972808	24990418	17610	0
# 3	16	34508448	34756258	247810	3
# 4	8	39237843	39385979	148136	1

# CNV output file (remove waves & telo/centromeres)
cnv  <- read.table("my_cnvPartition_output.txt", header = T, sep = "\t", stringsAsFactors = F)

######################################################## 
# IMPORT chromosome-seperated SNP reference files      #
######################################################## 

for (i in 1:22)
{
	df <- read.table(paste("chr", i, "_template.txt", sep=""), header = T, sep = "\t", stringsAsFactors = F) ;
	df$Loss <- 0 ; 				# init Loss
	df$Gain <- 0 ; 				# init Gain
	name <- paste("chr", i, sep="") ; 	# varying df name [1/2]
	assign(name, df) ;			# varying df name [2/2] 
}

######################################################## 
# TRANFORM  chromosome-seperated SNP reference files  #
######################################################## 

for (i in 1:dim(cnv)[1])
{
	# first select the correct chromosome for internal gain/loss addition
	df <- get(paste("chr", cnv$chr[i], sep="")) # varying df name [1/2]
	
	# loop through the entire chromosome
	for (j in 1:dim(df)[1])
	{
		if (df[j] >= cnv$start[i])
		{
			if (df[j] <= cnv$end[i])
			{
				# add 1 to Loss
				if(cnv$value[i] < 2 )
				{
					df[j] <- df[j] + 1 ;
				}
				# add 1 to Gain
				if(cnv$value[i] > 2 )
				{
					df[j] <- df[j] + 1 ;
				}	
			}
		}
	# save new internal count to external data frame
	assign(paste("chr", cnv$chr[i], sep=""), df) # varying df name [1/2]
	}
}

######################################################## 
# EXPORT  chromosome-seperated SNP reference files     #
######################################################## 

for (i in 1:22)
{
	df <- get(paste("chr", i, sep=""))			# varying df name [1/1]
 	file <- paste("chr", i, "_count.txt", sep="") ;		# varying filename [1/1]
	write.table(df, file, sep = "\t", row.names = F, col.names = T, quote = F) ;
}
