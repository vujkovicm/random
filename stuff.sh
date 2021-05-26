# Genisis, where is my SNP?
count=1
echo -n "" > pfiles.txt
while IFS= read -r line;
do
 chr="$(cut -d' ' -f2 <<< $line)"
 pos="$(cut -d' ' -f3 <<< $line)"
 grep -l :$pos: /data/data1/mvp003/mvp_imputed/Release4_PGEN/chr$chr/*.pvar >> pfiles.txt
 (( count++ ))
done < snps.txt

# download data from BOX folder
lftp -p 990 -e "open username@upenn.edu@ftp.box.com; mirror -R /local/path/to/folder box_name; exit"

# manipulate gzippe file
gzip -dc infile | awk '{print $0}' > outfile 

# selectively kill jobs
bkill `bjobs -p -o jobid | grep -v ^JOBID | tr '\n' ' '`
bjobs -w | grep 'PEND' | awk '{print $1}' | xargs bkill

# remove all comments, # from file
grep -o '^[^#]*' file > new.file

# quick grep based on list
grep -w -E 'first|second|last' textfile.txt > subset.txt

# run shell script without killing it due to inactive time (genisis stuff)
# keep running on the background
# save output every 60 seconds to nohup.out
nohup watch -n60 'bash script.sh' &
bash script &

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

# split filename into base and extension
base="${file%.*}"
ext="${filename##*.}"

# print duplicate lines
awk 'seen[$0]++' filename 

# replace first occurence in file
sed -i '1,/pattern/{s/pattern/newpattern/;}' filename.txt
 
# concatenate fields from 2 files
# paste 6 and 7 column from file 2 seperated by comma with entire of file1
awk 'FNR==NR {a[FNR""]= $6 FS $7; next}{print $0, a[FNR""]}' file2 file1 > file.out

# merge two files based on key identifierd (untion)
awk 'FNR==NR{a[$1]=$2 FS $3; next}{ print $0, a[$1]}' file2 file1 > file.out

# replace missing with value in specific column usign awk
awk '{ sub("NA", -9, $5); print }' file.in > file.out

# copy files f
rom aws s3 bucket
aws s3 cp s3://bucket/dir/ . --recursive

#split chrX:123 into two columns and convert alleles to uppercase
awk '{split($1, a, ":"); print substr(a[1], 4), a[2], toupper($2), toupper($3)}' filename 

# grep exact match from list in file
grep -w -f lookup.snps filename

# find duplicate value in row 4, print line
awk '{++a[$4]; b[$4]=b[$4]? b[$4] ORS $0:$0} END {for(i in a) {if(a[i]>1) {print b[i]}}}' filename

# inverse grep on multiple patterns
grep -v -E 'AGE|SEX|PC' filename

# remove files from a specific date (sugen - direct log)
rm `ls -all * | grep "Jan 11" | awk '{print $9}'`

# replace header
sed -i '1s/.*/NEW HEADER COLS/g' file.txt

# insert new header
sed -i '1 i\HEADER COLS TXT' file.txt

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

# if-else statement
if [ $POP == "TRANS" ]
then
        REF="EUR"
else
        REF=$POP
fi

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

