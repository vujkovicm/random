 ###################
#		      # 
 #  makeCircos     #
#		      #
 ###################

# invoke R function
source("c:\\stat\\myScripts\\makeCircosFormat.R", echo=T)

# set working directory
setwd("C:\\stat\\myData\\") 

# data file
# include variable names and patient identifiers in the first row
infile <- "myInputFile.txt"

# output is stored in wd
outfile <- "myCircosFile.txt"

# call the function
makeCircos (infile, outfile)


# 
# table rearrangement
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



