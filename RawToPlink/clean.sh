#!/bin/sh

# Cuts
MAF=0.01 
GENO=0.05
MIND=0.05
LDCUT=0.2
HWE=0.000001

# Directories/timing
SECONDS=0 
FILE=$1
WORKING_DIR=$PWD
HRC_DIR='/home/jrca253/HRC/'
OUT_DIR='VCFGZ'

cd $WORKING_DIR

# 0. FILE CONVERSION AND SNP FILTERING PART 1
echo -e "\n########################## 0. FILE CONVERSION AND SNP FILTERING PART 1 ##########################\n"
echo -e "\n########### plink --lfile $FILE --out ${FILE}_snpqc1  --geno $GENO --maf $MAF --recode --make-bed --freq ###########\n"
plink --lfile $FILE --out ${FILE}_snpqc1  --geno $GENO --maf $MAF --recode --make-bed --freq

# 1. FILTER SAMPLES WITH EXCESSIVE MISSINGNESS
echo -e "\n########################## 1. FILTER SAMPLES WITH EXCESSIVE MISSINGNESS  ##########################\n"
echo -e "\n########### plink --file ${FILE}_snpqc1 --out ${FILE}_samplefilter1  --mind $MIND --recode --make-bed --freq ###########\n"
plink --file ${FILE}_snpqc1 --out ${FILE}_samplefilter1  --mind $MIND --recode --make-bed --freq

# 2. LD PRUNE
echo -e "\n########################## 2. LD PRUNE ##########################\n"
echo -e "\n########### plink --file ${FILE}_samplefilter1 --out snps_to_ld_prune --indep-pairwise 50 5 $LDCUT ###########\n"
plink --file ${FILE}_samplefilter1 --out snps_to_ld_prune --indep-pairwise 50 5 $LDCUT
echo -e "\n########### plink --file ${FILE}_samplefilter1 --extract snps_to_ld_prune.in --make-bed --freq --out ${FILE}_samplefilter1_ldprune ###########\n"
plink --file ${FILE}_samplefilter1 --extract snps_to_ld_prune.prune.in --recode --make-bed --freq --out ${FILE}_samplefilter1_ldprune 

# 3. CUT BASED ON INBREEDING COEFFICIENT
echo -e "\n########################## 3. CUT BASED ON INBREEDING COEFFICIENT ##########################\n"
echo -e "\n########### plink --file ${FILE}_samplefilter1_ldprune --het --out inbreed_coefficient ###########\n"
plink --file ${FILE}_samplefilter1_ldprune --het --out inbreed_coefficient 
echo -e "\n########### awk 'sqrt($6*$6)>0.1 {print ;}' inbreed_coefficient.het | awk {'print $1,$2'} > inbreed_samples_to_remove.txt  ###########\n"
awk 'sqrt($6*$6)>0.1 {print ;}' inbreed_coefficient.het | awk {'print $1,$2'} > inbreed_samples_to_remove.txt
echo -e "\n########### plink --file ${FILE}_samplefilter1 --remove inbreed_samples_to_remove.txt -out ${FILE}_samplefilter2 --make-bed --recode --freq ###########\n"
plink --file ${FILE}_samplefilter1 --remove inbreed_samples_to_remove.txt --out ${FILE}_samplefilter2 --make-bed --recode --freq
echo -e "\n########### plink --file ${FILE}_samplefilter1_ldprune --remove inbreed_samples_to_remove.txt --out ${FILE}_samplefilter2_ldprune --make-bed --recode --freq ###########\n"
plink --file ${FILE}_samplefilter1_ldprune --remove inbreed_samples_to_remove.txt --out ${FILE}_samplefilter2_ldprune --make-bed --recode --freq 

# 4. CUT BASED ON IBD COEFFICIENT
echo -e "\n########################## 4. CUT BASED ON IBD COEFFICIENT ##########################\n"
echo -e "\n########### plink --file ${FILE}_samplefilter2_ldprune --genome --min 0.1 --out ibd_coeff ###########\n"
plink --file ${FILE}_samplefilter2_ldprune --genome --min 0.1 --out ibd_coeff
echo -e "\n########### awk {'print $1,$3'} ibd_coeff.genome > ibd.txt ###########\n"
awk {'print $1,$3'} ibd_coeff.genome > ibd.txt
echo -e "\n########### plink --file ${FILE}_samplefilter2 --missing --out missingess_report ###########\n"
plink --file ${FILE}_samplefilter2 --missing --out missingess_report
echo -e "\n########### awk {'print $1,$6'} missingess_report.imiss > miss.txt ###########\n"
awk {'print $1,$6'} missingess_report.imiss > miss.txt
echo -e "\n########### python detSamplesToIBDCut.py ###########\n"
python detSamplesToIBDCut.py
echo -e "\n########### plink --file ${FILE}_samplefilter2 --remove ibd_samples_to_cut.txt --out ${FILE}_samplefilter3 --make-bed --recode --freq ###########\n"
plink --file ${FILE}_samplefilter2 --remove ibd_samples_to_cut.txt --out ${FILE}_samplefilter3 --make-bed --recode --freq 

# 5. DETERMINE SAMPLES TO CUT BASED ON SEX MISMATCH
echo -e "\n########################## 5. DETERMINE SAMPLES TO CUT BASED ON SEX MISMATCH ##########################\n"
echo -e "\n########### plink --file ${FILE}_samplefilter3 --out ${FILE}_samplefilter3 --split-x 'hg19' --make-bed ###########\n"
plink --file ${FILE}_samplefilter3 --out ${FILE}_samplefilter3 --split-x 'hg19' --make-bed 
echo -e "\n########### plink --bfile ${FILE}_samplefilter3 --out ${FILE}_samplefilter3 --check-sex  ###########\n"
plink --bfile ${FILE}_samplefilter3 --out ${FILE}_samplefilter3 --check-sex
echo -e "\n########### grep \"PROBLEM\" ${FILE}_samplefilter3.sexcheck |  awk '{print $1,$2}' > ambiguous_sex_samples.txt ###########\n"
grep "PROBLEM" ${FILE}_samplefilter3.sexcheck |  awk '{print $1,$2}' > ambiguous_sex_samples.txt
echo -e "\n########### plink --bfile ${FILE}_semiclean --remove ambiguous_sex_samples.txt -out ${FILE}_semiclean2 --make-bed ###########\n"
plink --bfile ${FILE}_samplefilter3 --remove ambiguous_sex_samples.txt -out ${FILE}_snpqc1_sampleqc --make-bed --recode --freq

# 6. CLEAN UP
echo -e "\n########################## 6. CLEAN UP ##########################\n"
rm *samplefilter*
rm ibd.txt 
rm miss.txt
rm *.log


# 7. SNP FILTERING PART 2
echo -e "\n########################## 7. SNP FILTERING PART 2  ##########################\n" 
echo -e "\n########### plink --file ${FILE}_snpqc1_sampleqc --out ${FILE}_snpqc1_sampleqc_snpqc2 --make-bed --recode --hwe $HWE ##########\n"
plink --file ${FILE}_snpqc1_sampleqc --hwe $HWE --out ${FILE}_snpqc1_sampleqc_snpqc2 --make-bed --recode
echo -e "\n########### plink --file ${FILE}_snpqc1_sampleqc_snpqc2  --out ${FILE}_snpqc1_sampleqc_snpqc2 --freq ###########\n"
plink --file ${FILE}_snpqc1_sampleqc_snpqc2  --out ${FILE}_snpqc1_sampleqc_snpqc2 --freq


# 8. HRC CHECK
echo -e "\n########################## 8. HRC CHECK ##########################\n"
echo -e "\n########### perl ${HRC_DIR}HRC-1000G-check-bim-NoReadKey.pl -r ${HRC_DIR}HRC.r1-1.GRCh37.wgs.mac5.sites.tab -b ${FILE}_snpqc1_sampleqc_snpqc2.bim -f ${FILE}_snpqc1_sampleqc_snpqc2.frq -h ###########\n"
perl ${HRC_DIR}HRC-1000G-check-bim-NoReadKey.pl -r ${HRC_DIR}HRC.r1-1.GRCh37.wgs.mac5.sites.tab -b ${FILE}_snpqc1_sampleqc_snpqc2.bim -f ${FILE}_snpqc1_sampleqc_snpqc2.frq -h
echo -e "\n########### bash Run-plink.sh ###########\n"
bash Run-plink.sh

# 9. CONVERT VCFGZ
echo -e "\n########################## 9. CONVERT TO VCFGZ ##########################\n"
for i in {1..22}
do
    echo -e "\n########################## Processing Chromosome $i ##########################\n"
    echo -e "\n########### plink --bfile ${FILE}_semiclean2-updated-chr$i --recode tab --out ${FILE}_clean-chr$i --make-bed ###########\n"
    plink --bfile ${FILE}_snpqc1_sampleqc_snpqc2-updated-chr$i --recode tab --out ${FILE}_clean-chr$i --make-bed
    echo -e "\n########### plink --file ${FILE}_clean-chr$i --out ${FILE}_clean-chr$i --recode vcf ###########\n"
    plink --file ${FILE}_clean-chr$i --out ${FILE}_clean-chr$i --recode vcf
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