# stating the obvious
bgzip -c file.vcf > file.vcf.gz
tabix -p vcf file.vcf.gz   

# extract SNPs from imputed file (impute2.gz)
plink --gen filename.gz  \
 --extract select.snp \
  --make-bed \
  --out filename \
  --sample filename.sample \
  --allow-extra-chr

# extract patients from vcf
bcftools view -Ov --samples-file ids.keeps file.vcf.gz -o new.file.vcf

# extract variant from vcf
bcftools query -R gene.bed \
 -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/nhomref\t%INFO/nhet\t%INFO/nhomvar\n' \
 -S samples.keep \
 -o output.txt \
  file.vcf.gz
  
 bcftools query -f'%CHROM\t%POS\t%REF\t%ALT\n' vcffile.vcf.gz > output.bed

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

# replace new line with comma
paste -sd, file.in > file.out

# push environemental variales into awk
awk -v "gene=$1" '{if($8 == gene) print $0}' file.in > file.out

# extract from gzipped text-file, disregarding header info (48 lines)
awk 'NR>48 {if($16 < 0.0000002) print $0}' <(gzip -dc file.gz)

# iterations to variables
chr=$1
i=$2
from=$((3000000*$i+1))
to=$((3000000*($i+1)))

# get left unique (ref check)
grep -Fxv -f mvp.x.snps g1k.snps > mvp.x.uniq

# impute2 to plink
plink2 --gen filename.gz --make-pgen --sort-vars --out filename --allow-extra-chr --sample pheno.sample
sed -i "s/---/19/g" filename.pvar
echo -e "#CHROM\tPOS\tID\tREF\tALT" > filename.new
awk '{if(NR>1) print $1"\t"$2"\tchr"$1":"$2"\t"$4"\t"$5}' filename.pvar >> filename.new
mv filename.pvar filename.tmp
mv filename.new filename.pvar
plink2 --pfile filename --make-bed --out filename
plink2 --bfile filename --extract ../extract.snp --recode A --out subset

# plink2 --gen chr6/chr_6_30e6_35e6.gz   --make-pgen --out chr6/filename  --allow-extra-chr --sample ../../../out/GWAS2/MI/GWAS2.MI.sample --sort-vars
# plink2 --gen chr19/chr_19_15e6_20e6.gz --make-pgen --out chr19/filename --allow-extra-chr --sample ../../../out/GWAS2/MI/GWAS2.MI.sample --sort-vars
# plink2 --gen chr7/chr_7_135e6_140e6.gz --make-pgen --out chr7/filename  --allow-extra-chr --sample ../../../out/GWAS2/MI/GWAS2.MI.sample --sort-vars
# plink2 --gen chr1/chr_1_195e6_200e6.gz --make-pgen --out chr1/filename  --allow-extra-chr --sample ../../../out/GWAS2/MI/GWAS2.MI.sample --sort-vars
# plink2 --gen chr3/chr_3_85e6_90e6.gz   --make-pgen --out chr3/filename  --allow-extra-chr --sample ../../../out/GWAS2/MI/GWAS2.MI.sample --sort-vars
# plink2 --gen chr5/chr_5_0e6_5e6.gz     --make-pgen --out chr5/filename  --allow-extra-chr --sample ../../../out/GWAS2/MI/GWAS2.MI.sample --sort-vars
# plink2 --gen chr4/chr_4_55e6_60e6.gz   --make-pgen --out chr4/filename  --allow-extra-chr --sample ../../../out/GWAS2/MI/GWAS2.MI.sample --sort-vars 

# second round
# plink2 --gen chr19/chr_19_50e6_55e6.gz --make-pgen --out chr19/filename --allow-extra-chr --sample ../../../out/GWAS2/MI/GWAS2.MI.sample --sort-vars 

# EXPORT VCF AS HARDCALLS in PLINK2
plink2 --vcf filename.vcf.gz dosage=DS --hard-call-threshold 0.499999 --extract common.snps --export vcf vcf-dosage=DS --out filename.common

