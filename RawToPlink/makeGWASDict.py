
#!/opt/anaconda/anaconda2/bin/python
import os
import pandas as pd
import numpy as np
import pickle
import time

DATA_DIR  = '/mnt/DATA/EA11101_2011-09-28/EA11101_2011-09-28/EA11101_2011-09-28_FinalReport_1_to_16/PLINK_FILES/'
GWAS_file = '/mnt/DATA/EA11101_2011-09-28/EA11101_2011-09-28/snp144.txt'
out_file  = DATA_DIR + "GWAS_dict.pkl"


start_time = time.time()

##### LOAD GWAS CSV #####
if False:
    print "Loading GWAS file ..."
    df = pd.read_table(
        GWAS_file,
        sep = '\t',
        header = None,
        names = ['Chr', 'BasePairPos', 'SNP'],
        usecols = [1,3,4]
        #nrows = 10 #5000000
    )



    ##### STRIP "chr" from CROMOSOME COLUMN #####
    df['Chr'] = df['Chr'].map( lambda x: x.lstrip('chr') )

                            

    ##### REORDER #####
    cols = ['SNP', 'Chr', 'BasePairPos']
    df = df[cols]
#ENDIF


##### FIND NON-UNIQUE SNPS #####
print "Loading GWAS file ..."
df = pd.read_table(
    GWAS_file,
    sep = '\t',
    header = None,
    names = ['Chr', 'BP End', 'SNP'],
    usecols = [1,3,4]
    #nrows = 10
)

print str(df.shape[0]) + " SNPs are in this file"

print "Remmoving Chromosome names other than \"chr#\" ..."
nbefore = df.shape[0]
df.drop(df[ df.Chr.str.contains("_") == True ].index, inplace=True)
nafter = df.shape[0]
nrm = nbefore-nafter
print "Removed %.0f SNPs" %(nrm)


print "Finding non-unique SNPs ..."
repeat_snps = df.loc[df.duplicated(["SNP"], keep='first')]['SNP'].unique().tolist()
repeat_chrs = df.loc[df.duplicated(["SNP"], keep='first')]['Chr'].unique().tolist()
print "%.0f non-unique SNPs found!" % (len(repeat_snps))
print "SNPs found on the following chromosomes: "
print repeat_chrs

iterations = 10
if len(repeat_snps) < iterations:
    iterations = len(repeat_snps)
#ENDIF


for i in range(0, iterations):
    print "Occurances of SNP " + str(repeat_snps[i])
    print df.loc[ df['SNP'] == str(repeat_snps[i]) ]
    print "\n\n"
#ENDFOR

end_time = time.time()
run_time = (end_time - start_time) / 60
print "Process completed in %.0f minutes" % (run_time)

exit()


##### SAVE RELEVANT COLUMNS AS A DICTIONARY #####
if os.path.isfile( out_file ):
    command = 'rm ' + out_file
    print 'Removing old GWAS_dict.pkl ...'
    os.system(command)
#ENDIF

print "Saving GWAS as a dictionary ..."
GWAS_dict = df.set_index('SNP').T.to_dict('list')
output = open(out_file, 'wb')
pickle.dump(GWAS_dict, output)

end_time = time.time()
run_time = end_time - start_time
print "Process completed in %.0f seconds" % (run_time) 


##### END SCRIPT #####

##### UNIQUE CHROMOSOMES IN FILE #####
#['chr1' 'chr10' 'chr11' 'chr11_gl000202_random' 'chr12' 'chr13' 'chr14'
# 'chr15' 'chr16' 'chr17' 'chr17_ctg5_hap1' 'chr17_gl000203_random'
# 'chr17_gl000204_random' 'chr17_gl000205_random' 'chr17_gl000206_random'
# 'chr18' 'chr18_gl000207_random' 'chr19' 'chr19_gl000208_random'
# 'chr19_gl000209_random' 'chr1_gl000191_random' 'chr1_gl000192_random'
# 'chr2' 'chr20' 'chr21' 'chr21_gl000210_random' 'chr22' 'chr3' 'chr4'
# 'chr4_ctg9_hap1' 'chr4_gl000193_random' 'chr4_gl000194_random' 'chr5'
# 'chr6' 'chr6_apd_hap1' 'chr6_cox_hap2' 'chr6_dbb_hap3' 'chr6_mann_hap4'
# 'chr6_mcf_hap5' 'chr6_qbl_hap6' 'chr6_ssto_hap7' 'chr7'
# 'chr7_gl000195_random' 'chr8' 'chr8_gl000196_random'
# 'chr8_gl000197_random' 'chr9' 'chr9_gl000198_random'
# 'chr9_gl000199_random' 'chr9_gl000200_random' 'chr9_gl000201_random'
# 'chrM' 'chrUn_gl000211' 'chrUn_gl000212' 'chrUn_gl000213'
# 'chrUn_gl000214' 'chrUn_gl000215' 'chrUn_gl000216' 'chrUn_gl000217'
# 'chrUn_gl000218' 'chrUn_gl000219' 'chrUn_gl000220' 'chrUn_gl000221'
# 'chrUn_gl000222' 'chrUn_gl000223' 'chrUn_gl000224' 'chrUn_gl000225'
# 'chrUn_gl000226' 'chrUn_gl000227' 'chrUn_gl000228' 'chrUn_gl000229'
# 'chrUn_gl000230' 'chrUn_gl000231' 'chrUn_gl000232' 'chrUn_gl000233'
# 'chrUn_gl000234' 'chrUn_gl000235' 'chrUn_gl000236' 'chrUn_gl000237'
# 'chrUn_gl000238' 'chrUn_gl000239' 'chrUn_gl000240' 'chrUn_gl000241'
# 'chrUn_gl000242' 'chrUn_gl000243' 'chrUn_gl000244' 'chrUn_gl000245'
# 'chrUn_gl000246' 'chrUn_gl000247' 'chrUn_gl000248' 'chrUn_gl000249'
# 'chrX' 'chrY']
