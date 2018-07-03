#!/bin/sh
SECONDS=0 
FILE=$1
WORKING_DIR=$PWD
HRC_DIR='/home/jrca253/HRC/'
cd $WORKING_DIR

echo -e "########################## Converting lgen to ped ##########################\n"
plink_1_9 --lfile $FILE --out ${FILE}_recode --snps-only just-acgt --recode

echo -e "########################## Converting ped to bed ##########################\n"
plink_1_9 --file ${FILE}_recode --out ${FILE}_recode --make-bed

echo -e "########################## Making frq file ##########################\n"
plink_1_9 --bfile ${FILE}_recode --out ${FILE}_recode --freq

echo -e "########################## Performing HRC Imputation preparation checks ##########################\n"
perl ${HRC_DIR}HRC-1000G-check-bim-NoReadKey.pl -r ${HRC_DIR}HRC.r1-1.GRCh37.wgs.mac5.sites.tab -b ${FILE}_recode.bim -f ${FILE}_recode.frq -h
bash Run-plink.sh

echo -e "########################## Converting to vcfgz ##########################\n"
for i in {1..22}
do
    echo -e "########################## Processing Chromosome $i ##########################\n"
    cp ${FILE}_recode.map ${FILE}_recode-updated-chr$i.map
    plink_1_9 --bfile ${FILE}_recode-updated-chr$i --recode --tab --out ${FILE}_recode-updated-chr$i
    plink_1_9 --file ${FILE}_recode-updated-chr$i --out ${FILE}_recode-updated-chr$i --recode vcf
    vcf-sort ${FILE}_recode-updated-chr$i.vcf | bgzip -c > ${FILE}_recode-updated-chr$i.vcf.gz
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

echo "Process completed in $((SECONDS/60)) minutes"