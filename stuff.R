# don't restrict R rounding of p-value to 2.2E-16 (default: 53 binary digits accuracy)
.Machine$double.eps # 2.220446e-16
.Machine$double.xmin # 2.225074e-308
# soluation
.Machine$double.eps = .Machine$double.xmin

# systematically rename column names
# remove everything after _ in column name
names(df) = gsub(pattern = "_*", replacement = "", x = names(df))

# remove everyting after dash (-) 
# e.g. 1:1520-1520 > 1:1520
df$CHRCBP = gsub(pattern = "(-).*", replacement = "", x = df$Location)

# keep position (e.g. 1520)
df$BP = gsub(pattern = "*.(-)", replacement = "", x = df$CHRCBP)

# keep chr (e.g. 1) 
df$CHR = gsub(pattern = "(:).*", replacement = "", x = df$CHRCBP)
