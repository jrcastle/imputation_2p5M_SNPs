#!/bin/sh
SECONDS=0 
FILE=$1
WORKING_DIR=$PWD
HRC_DIR='/home/jrca253/HRC/'
OUT_DIR='VCFGZ'
cd $WORKING_DIR

# Convert and first-pass QC
echo -e "\n########################## Converting lgen to ped, bed, and freq ##########################\n"
plink_1_9 --lfile $FILE --out ${FILE}_semiclean  --geno 0.1 --maf 0.005 --mind 0.1 --recode --make-bed --freq

# SEX CHECK
echo -e "\n########################## Sex check ##########################\n"
#plink_1_9 --file ${FILE}_semiclean --out ${FILE}_semiclean --split-x 'hg19' --make-bed
plink_1_9 --bfile ${FILE}_semiclean --out ${FILE}_semiclean --check-sex 
grep "PROBLEM" ${FILE}_semiclean.sexcheck |  awk '{print $1,$2}' > ambiguous_sex_samples.txt
plink_1_9 --bfile ${FILE}_semiclean --remove ambiguous_sex_samples.txt -out ${FILE}_semiclean2 --make-bed
plink --bfile ${FILE}_semiclean2 --recode --out ${FILE}_semiclean2
plink --bfile ${FILE}_semiclean2 --out ${FILE}_semiclean2 --freq

# HRC CHECK
echo -e "\n########################## Performing HRC Imputation preparation checks ##########################\n"
perl ${HRC_DIR}HRC-1000G-check-bim-NoReadKey.pl -r ${HRC_DIR}HRC.r1-1.GRCh37.wgs.mac5.sites.tab -b ${FILE}_semiclean2.bim -f ${FILE}_semiclean2.frq -h
bash Run-plink.sh

echo -e "\n########################## Converting to vcfgz ##########################\n"
for i in {1..22}
do
    echo -e "\n########################## Processing Chromosome $i ##########################\n"
    cp ${FILE}_recode.map ${FILE}_semiclean2-updated-chr$i.map
    plink_1_9 --bfile ${FILE}_semiclean2-updated-chr$i --recode --tab --out ${FILE}_semiclean2-updated-chr$i --geno 0.1 --maf 0.005 --hwe 0.00005  
    plink_1_9 --file ${FILE}_semiclean2-updated-chr$i --out ${FILE}_semiclean2-updated-chr$i --recode vcf --geno 0.1 --maf 0.005 --hwe 0.00005  
    vcf-sort ${FILE}_semiclean2-updated-chr$i.vcf | bgzip -c > ${FILE}_semiclean2-updated-chr$i.vcf.gz
done

#echo "Cleaning up ..."
#cd $WORKING_DIR
#rm ${FILE}_recode_chr*.bed
#rm ${FILE}_recode_chr*.bim
#rm ${FILE}_recode_chr*.fam
#rm ${FILE}_recode_chr*.log
#rm ${FILE}_recode_chr*.map
#rm ${FILE}_recode_chr*.ped
#rm ${FILE}_recode_chr*.vcf
#rm ${FILE}_recode.*

if[ -d $OUT_DIR ]; then
    rm -r $OUT_DIR
fi

mkdir $OUT_DIR
mv *.vcf.gz ${OUT_DIR}/.

echo "Process completed in $((SECONDS/60)) minutes"