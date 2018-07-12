#!/bin/sh
SECONDS=0 
WORK_DIR="/home/jrca253/imputation_2p5M_SNPs/RawToPlink/"
DATA_DIR="/home/jrca253/DATA/EA11101_2011-09-28_FinalReport_1_to_16/PLINK_FILES/"

echo -e "#################### makeFam.py ####################\n"
cd $WORK_DIR
python -u makeFAM.py

echo -e "#################### Using sed to replace \"_\" with \"-\" in each fam file ####################\n"
cd $DATA_DIR

echo "sed -i \"s/_/-/g\" EUR_EA11101_2011-09-28.fam"
sed -i 's/_/-/g' EUR_EA11101_2011-09-28.fam
echo "sed -i \"s/_/-/g\" AFR_EA11101_2011-09-28.fam"
sed -i 's/_/-/g' AFR_EA11101_2011-09-28.fam
echo "sed -i \"s/_/-/g\" ASN_EA11101_2011-09-28.fam"
sed -i 's/_/-/g' ASN_EA11101_2011-09-28.fam
echo "sed -i \"s/_/-/g\" OTR_EA11101_2011-09-28.fam"
sed -i 's/_/-/g' OTR_EA11101_2011-09-28.fam

echo -e "#################### makeMAP.py ####################\n"
cd $WORK_DIR
python -u makeMAP.py

echo -e "#################### makeLGEN.py ####################\n"
cd $WORK_DIR
python -u makeLGEN.py

echo -e "#################### Using sed to replace \"_\" with \"-\" in each lgen file ####################\n"
cd $DATA_DIR

echo "sed -i \"s/_/-/g\" EUR_EA11101_2011-09-28.lgen"
sed -i 's/_/-/g' EUR_EA11101_2011-09-28.lgen
echo "sed -i \"s/_/-/g\" AFR_EA11101_2011-09-28.lgen"
sed -i 's/_/-/g' AFR_EA11101_2011-09-28.lgen
echo "sed -i \"s/_/-/g\" ASN_EA11101_2011-09-28.lgen"
sed -i 's/_/-/g' ASN_EA11101_2011-09-28.lgen
echo "sed -i \"s/_/-/g\" OTR_EA11101_2011-09-28.lgen"
sed -i 's/_/-/g' OTR_EA11101_2011-09-28.lgen

echo "makePLINK.sh process completed in $((SECONDS/3600)) hours" 