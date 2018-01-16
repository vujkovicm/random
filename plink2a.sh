plink2a --vcf filename.dose.vcf.gz --export bgen-1.1 --out filename

# sample file 
# phenotype: 0/1
# categories: 1/2 (cannot contain zero's)

for i in {1..22}
do
        bsub -q normal -n 24 -M 140000 -e ../output/chr$i/msg$i.e -o ../output/chr$i/msg$i.o plink2 \
        --bgen ../imp/ukb_imp_chr${i}_v2.bgen \
        --sample ../sample/BBUK.T2D.v4.sample \
        --allow-extra-chr \
        --glm \
        --ci 0.95 --threads 24 --missing-code NA \
        --out ../output/chr${i}/BBUK.T2D.chr${i}.out
done
