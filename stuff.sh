# stating the obvious
bgzip -c file.vcf > file.vcf.gz
tabix -p vcf file.vcf.gz   

# extract patients from vcf
bcftools view -Ov --samples-file ids.keeps file.vcf.gz -o new.file.vcf

# single variant association
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
