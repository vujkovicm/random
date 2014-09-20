setwd("T:\\ymy_directory")

# Chromosome seperated SNP reference file format
#
# Chr	Position
# 1	82154
# 1	534247
# 1	565374
# 1	569624
# 1	723918
# 1	729632

# cnvPartition output
# regno chr	start	    end	      size	  value
# 1	    5	  49447784	49548900	101116	1
# 1	    12	37999925	38410357	410432	1
# 2	    17	34443830	34815551	371721	3
# 3	    8	  24972808	24990418	17610	  0
# 3	    16	34508448	34756258	247810	3
# 4	    8	  39237843	39385979	148136	1

# CNV output file (remove waves & telo/centromeres)
cnv  <- read.table("my_cnvPartition_output.txt", header = T, sep = "\t", stringsAsFactors = F)

######################################################## 
# IMPORT chromosome-seperated SNP reference files      #
######################################################## 

# chromosome-seperated SNP reference files
# init count

for (i in 1:22)
{
	df <- read.table(paste("chr", i, "_template.txt", sep=""), header = T, sep = "\t", stringsAsFactors = F) ;
	df$Loss <- 0 ; # init Loss
	df$Gain <- 0 ; # init Gain
	name <- paste("chr", i, sep="") ;
	assign(name, df) ;
}

###############
# TRANFORM    #
###############

# count by looping through cnv-file 
# for (i in 1:dim(cnv)[1])

for (i in 1:1)
{
	# first select the correct chromosome for internal count revision
	chr <- get(paste("chr", cnv$chr[i], sep=""))
	print(head(chr)) ;
	
	# loop through the entire chromosome
	for (j in 1:dim(chr)[1])
	{
		if (chr$Position[j] >= cnv$start[i])
		{
			if (chr$Position[j] <= cnv$end[i])
			{
				if(cnv$value[i] < 2 )
				{
					print (paste(j, ": ", chr$Position[j], sep = "")) ;
					chr$Loss[j] <- chr$Loss[j] + 1 ;
				}
				if(cnv$value[i] > 2 )
				{
					chr$Gain[j] <- chr$Gain[j] + 1 ;
				}	
			}
		}
	# save new internal count to external data frame
	assign(paste("chr", cnv$chr[i], sep=""), chr)
	}
}

##############
# EXPORT     #
##############

for (i in 1:22)
{
	df <- get(paste("chr", i, sep=""))
 	file <- paste("chr", i, "_count.txt", sep="") ;
	write.table(df, file, sep = "\t", row.names = F, col.names = T, quote = F) ;
}
