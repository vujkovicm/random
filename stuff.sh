# print duplicate lines
awk 'seen[$0]++' filename 

#split chrX:123 into two columns
awk '{split($1, a, ":"); print substr(a[1], 4), a[2], toupper($2), toupper($3)}' filename 
