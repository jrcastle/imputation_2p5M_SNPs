#!/bin/sh
FILE=$1
echo "Converting lgen to ped ..."
plink --lfile $FILE --out ${FILE}_recode --recode
cp $FILE.map ${FILE}_recode.map

echo "Converting ped to bed ..."
plink --file ${FILE}_recode --out ${FILE}_recode --make-bed

echo "Splitting by chromosome ..."
perl split_by_chromosome.pl ${FILE}_recode ${FILE}_recode_Chr

echo "Converting to vcfgz ..."
for i in {1..22}
do
    echo "Processing Chromosome $i ..."
    plink --file ${FILE}_recode_Chr$i --out ${FILE}_recode_Chr$i --map ${FILE}_recode.map --recode vcf
    vcf-sort ${FILE}_recode_Chr$i.vcf | bgzip -c > ${FILE}_recode_Chr$i.vcf.gz
done