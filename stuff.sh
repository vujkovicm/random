# stating the obvious
bgzip -c file.vcf > file.vcf.gz
tabix -p vcf file.vcf.gz   

# extract patients from vcf
bcftools view -Ov --samples-file ids.keeps file.vcf.gz -o new.file.vcf

# single variant association [binary coded as 1/2]
rvtest --inVcf file.vcf.gz \
  --single score,wald \
  --pheno file.sample \
  --covar file.sample \
  --pheno-name phenotype \
  --covar-name COV1,COV2,COVN \
  --minMAF 0.0001 \
  --out file

# print duplicate lines
awk 'seen[$0]++' filename 

#split chrX:123 into two columns and convert alleles to uppercase
awk '{split($1, a, ":"); print substr(a[1], 4), a[2], toupper($2), toupper($3)}' filename 

# grep exact match from list in file
grep -w -f lookup.snps filename

# inverse grep on multiple patterns
grep -v -E 'AGE|SEX|PC' filename

# replace header
sed -i '1s/.*/NEW HEADER COLS/g' file.txt

# extract from gzipped text-file, disregarding header info (48 lines)
awk 'NR>48 {if($16 < 0.0000002) print $0}' <(gzip -dc file.gz)

# iterations to variables
chr=$1
i=$2
from=$((3000000*$i+1))
to=$((3000000*($i+1)))

# get left unique (ref check)
grep -Fxv -f mvp.x.snps g1k.snps > mvp.x.uniq




