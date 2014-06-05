 ###################
#		      # 
 #  makeCircos     #
#		      #
 ###################

# invoke the R script file with makeCircos function by 
source("c:\\stat\\myScripts\\makeCircosFormat.txt", echo=T)

# directory name, with double slashes in windows environment
setwd("C:\\stat\\myData\\") 

# data file
# data file must include variable names  
# data file must include patient identifiers in the first row
# the output will be stored in the same directory as the input file
infile <- "myInputFile.txt"
outfile <- "myCircosFile.txt"

# call the function
makeCircos (infile, outfile)


#
#
#
# This is the function
# Don't change!
#
#
#
#

makeCircos <- function(infile, outfile)
{
	df <- read.table(infile, T);
	df.in <- df[ , 2:dim(df)[2]];
	df.circos <- data.frame(VAR1 = numeric(dim(df.in)[2]), stringsAsFactors = FALSE);
	for(i in 1:dim(df.in)[2]) 
	{
		row.names(df.circos)[i] <- names(df.in)[i];
		df.circos[ , 1] <- 0;
		for (j in 1:i)
		{ 
			if(!is.na(table(df.in[[i]], df.in[[j]])[4]))
			{	
				df.circos[j, i] <- table(df.in[[i]], df.in[[j]])[4];
			}
			if(is.na(table(df.in[[i]], df.in[[j]])[4]))
			{	
				df.circos[j, i] <- 0;
			}
			if (df.circos[j, i] <= 1 ) 
			{
				df.circos[j, i] <- 0;
			}
		}
		df.circos[i, i] <- 0;
		colnames(df.circos)[i]  <- names(df.in)[i];
		for (j in 1:dim(df.in)[2])
		{ 
			if (is.na(df.circos[j, i])) 
			{
				df.circos[j, i] <- 0;
			}
		}
	}
	write.table(df.circos, outfile, row.names = T, col.names = T, quote = F, sep ="\t")
}



