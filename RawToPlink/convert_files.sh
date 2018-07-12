#!/bin/sh
SECONDS=0 
FILE=$1
WORKING_DIR=$PWD
HRC_DIR='/home/jrca253/HRC/'
OUT_DIR='VCFGZ'
cd $WORKING_DIR

# Convert and first-pass QC
echo -e "\n########################## Converting lgen to ped, bed, and freq ##########################\n"
echo -e "\n########### plink_1_9 --lfile $FILE --out ${FILE}_semiclean  --geno 0.1 --maf 0.005 --mind 0.1 --recode --make-bed --freq ###########\n"
plink_1_9 --lfile $FILE --out ${FILE}_semiclean  --geno 0.1 --maf 0.005 --mind 0.1 --recode --make-bed --freq

# SEX CHECK
echo -e "\n########################## SEX CHECK ##########################\n"
echo -e "\n########### plink_1_9 --file ${FILE}_semiclean --out ${FILE}_semiclean --split-x \"hg19\" --make-bed ###########\n"
plink_1_9 --file ${FILE}_semiclean --out ${FILE}_semiclean --split-x 'hg19' --make-bed 

echo -e "\n########### plink_1_9 --bfile ${FILE}_semiclean --out ${FILE}_semiclean --check-sex ###########\n"
plink_1_9 --bfile ${FILE}_semiclean --out ${FILE}_semiclean --check-sex

echo -e "\n########### grep \"PROBLEM\" ${FILE}_semiclean.sexcheck |  awk '{print $1,$2}' > ambiguous_sex_samples.txt ###########\n"
grep "PROBLEM" ${FILE}_semiclean.sexcheck |  awk '{print $1,$2}' > ambiguous_sex_samples.txt

echo -e "\n########### plink_1_9 --bfile ${FILE}_semiclean --remove ambiguous_sex_samples.txt -out ${FILE}_semiclean2 --make-bed ###########\n"
plink_1_9 --bfile ${FILE}_semiclean --remove ambiguous_sex_samples.txt -out ${FILE}_semiclean2 --make-bed

echo -e "\n########### plink_1_9 --bfile ${FILE}_semiclean2 --recode --out ${FILE}_semiclean2 ###########\n"
plink_1_9 --bfile ${FILE}_semiclean2 --recode --out ${FILE}_semiclean2

echo -e "\n########### plink_1_9 --bfile ${FILE}_semiclean2 --out ${FILE}_semiclean2 --freq ###########\n"
plink_1_9 --bfile ${FILE}_semiclean2 --out ${FILE}_semiclean2 --freq


# HRC CHECK
echo -e "\n########################## Performing HRC Imputation preparation checks ##########################\n"
echo -e "\n########### perl ${HRC_DIR}HRC-1000G-check-bim-NoReadKey.pl -r ${HRC_DIR}HRC.r1-1.GRCh37.wgs.mac5.sites.tab -b ${FILE}_semiclean2.bim -f ${FILE}_semiclean2.frq -h ###########\n"
perl ${HRC_DIR}HRC-1000G-check-bim-NoReadKey.pl -r ${HRC_DIR}HRC.r1-1.GRCh37.wgs.mac5.sites.tab -b ${FILE}_semiclean2.bim -f ${FILE}_semiclean2.frq -h
echo -e "\n########### bash Run-plink.sh ###########\n"
bash Run-plink.sh

echo -e "\n########################## Converting to vcfgz ##########################\n"
for i in {1..22}
do
    echo -e "\n########################## Processing Chromosome $i ##########################\n"
    echo -e "\n########### plink_1_9 --bfile ${FILE}_semiclean2-updated-chr$i --recode tab --out ${FILE}_clean-chr$i --make-bed --geno 0.1 --maf 0.005 --hwe 0.00005 ###########\n"
    plink_1_9 --bfile ${FILE}_semiclean2-updated-chr$i --recode tab --out ${FILE}_clean-chr$i --make-bed --geno 0.1 --maf 0.005 --hwe 0.00005
    echo -e "\n########### plink_1_9 --file ${FILE}_clean-chr$i --out ${FILE}_clean-chr$i --recode vcf --geno 0.1 --maf 0.005 --hwe 0.00005 ###########\n"
    plink_1_9 --file ${FILE}_clean-chr$i --out ${FILE}_clean-chr$i --recode vcf --geno 0.1 --maf 0.005 --hwe 0.00005
    echo -e "\n########### vcf-sort ${FILE}_clean-chr$i.vcf | bgzip -c > ${FILE}_clean-chr$i.vcf.gz ###########\n"
    vcf-sort ${FILE}_clean-chr$i.vcf | bgzip -c > ${FILE}_clean-chr$i.vcf.gz
done

if [ -d $OUT_DIR ]; then
    rm -r $OUT_DIR
fi

echo -e "\n########### mkdir $OUT_DIR ###########\n"
mkdir $OUT_DIR
echo -e"\n########### mv *.vcf.gz ${OUT_DIR}/. ###########\n"
mv *.vcf.gz ${OUT_DIR}/.

#echo -e "\n########################## CLEANING UP ##########################\n"
#rm *semiclean*
#rm Run-plink.sh
#rm ambiguous_sex_samples.txt

echo "Process completed in $((SECONDS/60)) minutes"