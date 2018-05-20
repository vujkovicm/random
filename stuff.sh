# print duplicate lines
awk 'seen[$0]++' filename 

#split chrX:123 into two columns and convert alleles to uppercase
awk '{split($1, a, ":"); print substr(a[1], 4), a[2], toupper($2), toupper($3)}' filename 

# grep exact match from list in file
grep -w -f lookup.snps filename

# inverse grep on multiple patterns
grep -v -E 'AGE|SEX|PC' filename
