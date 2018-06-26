#!/opt/anaconda/anaconda2/bin/python
import os
import pandas as pd
import numpy as np
import pickle
import time

start_time = time.time()

TEST = 0
DATA_DIR = '/mnt/DATA/EA11101_2011-09-28/EA11101_2011-09-28/EA11101_2011-09-28_FinalReport_1_to_16/'
files = [
    DATA_DIR + 'EA11101_2011-09-28_FinalReport1.txt',
    DATA_DIR + 'EA11101_2011-09-28_FinalReport2.txt',
    DATA_DIR + 'EA11101_2011-09-28_FinalReport3.txt',
    DATA_DIR + 'EA11101_2011-09-28_FinalReport4.txt',
    DATA_DIR + 'EA11101_2011-09-28_FinalReport5.txt',
    DATA_DIR + 'EA11101_2011-09-28_FinalReport6.txt',
    DATA_DIR + 'EA11101_2011-09-28_FinalReport7.txt',
    DATA_DIR + 'EA11101_2011-09-28_FinalReport8.txt',
    DATA_DIR + 'EA11101_2011-09-28_FinalReport9.txt',
    DATA_DIR + 'EA11101_2011-09-28_FinalReport10.txt',
    DATA_DIR + 'EA11101_2011-09-28_FinalReport11.txt',
    DATA_DIR + 'EA11101_2011-09-28_FinalReport12.txt',
    DATA_DIR + 'EA11101_2011-09-28_FinalReport13.txt',
    DATA_DIR + 'EA11101_2011-09-28_FinalReport14.txt',
    DATA_DIR + 'EA11101_2011-09-28_FinalReport15.txt',
    DATA_DIR + 'EA11101_2011-09-28_FinalReport16.txt'
]

out_file       = DATA_DIR + 'PLINK_FILES/EA11101_2011-09-28.map'
manifest_file  = DATA_DIR + 'PLINK_FILES/manifest_dict_MERGED.pkl'
rsConverter    = DATA_DIR + 'PLINK_FILES/rsConverterDict.pkl'
GWAS_file      = DATA_DIR + 'PLINK_FILES/GWAS_dict.txt'
unmatched_file = DATA_DIR + "PLINK_FILES/unmatchedSNPs.txt"

if TEST:
    files          = ['test.txt', 'test2.txt']
    out_file       = 'EA11101_2011-09-28.map'
    unmatched_file = 'unmatchedSNPs.txt'


if not os.path.isfile( GWAS_file ):
    print 'GWAS dictionary not found. Make the GWAS dictionary file first and then run this script.'
    exit()
#ENDIF 



##### LOAD GWAS DICTIONARY #####
print "Loading GWAS dictionary at " + GWAS_file
chunksize = 1000000
GWAS_dict = {}

i = 0
for chunk in pd.read_table(GWAS_file, sep = '\t', header = 0, dtype = {'SNP': object, 'Chr': np.dtype('S2'), 'BasePairPos': np.uint32}, chunksize=chunksize):
    pct_cpt = 100.*float(i)/149.
    print "Loading chunk %i \t %.1f%% complete ..." % (i, pct_cpt)
    print chunk.info(memory_usage='deep')
    exit()
    GWAS_dict.update(chunk.set_index('SNP').T.to_dict('list'))
    i = i + 1
#ENDFOR



##### LOAD RS CONVERTER DICTIONARY #####
print "Loading rs converter dictionary at " + rsConverter
pkl_file2 = open(rsConverter, 'rb')
rs_converter_dict = pickle.load(pkl_file2)
pkl_file2.close()



##### INITIALIZE THE SNP ARRAY
all_snps_array = np.array( [] )



for f in files:

    ##### LOAD DATA FILES #####
    print "Loading " + f + " ..."
    df = pd.read_table(
        f,
        sep = '\t',
        skiprows = range(0,10),
        header = 0,
        usecols = [0]
    )



    ##### CONVERT NON-RS IDs to RS IDs #####
    print "Converting to rsIDs ..."
    df['SNP Name'] = df['SNP Name'].map(rs_converter_dict).fillna(df['SNP Name'])



    ##### GET UNIQUE SNPS #####
    print "Getting list of unique SNPs ..."
    snp_list = df['SNP Name'].unique()
    all_snps_array = np.union1d(all_snps_array, snp_list)


    ##### FREE MEMORY #####
    del df
    del snp_list
        
#ENDFOR



##### CREATE NEW DATAFRAME WITH UNIQUE SNPS
df_snps = pd.DataFrame(all_snps_array, columns = ["SNP"])
del all_snps_array


##### MATCH SNP TO CHROMOSOME AND BP POSITION #####
print "Matching SNP to chromosome and base-pair position ..."
df_snps['tmp'] = df_snps['SNP'].map(GWAS_dict)
del GWAS_dict
rows_before = df.shape[0]



##### FIND UNMATCHED SNPS #####
print "Finding unmatched SNPs"
df2 = df_snps[df_snps['tmp'].isnull()]
df2 = df2[['SNP']]

df_snps.dropna(inplace=True)
rows_after = df_snps.shape[0]
print "%.0f SNPS ARE UNMATCHED!" % (rows_before - rows_after)

print "List of Unmatched SNPs saved in " + unmatched_file
df2.to_csv(
    unmatched_file,
    sep = '\t',
    header = False,
    index = False
)
del df2


##### SPLIT THE DICTIONARY MAPPING INTO TWO COLUMNS #####
df_snps[['Chromosome','BasePairPos']] = pd.DataFrame(df_snps.tmp.values.tolist(), index= df_snps.index)

df_snps.drop("tmp",
        axis = 1,
        inplace = True
)



##### ADD GENETIC DISTANCE (MORGANS) AND SET TO 0 #####
df_snps['GeneticDist'] = 0



##### REORGANIZE COLUMNS AND WRITE #####
reordered_cols = [
    'Chromosome',
    'SNP',
    'GeneticDist',
    'BasePairPos'
]

df_snps = df_snps[ reordered_cols ]

print "Writing MAP file to " + out_file
df_snps.to_csv(
    out_file,
    sep = '\t',
    header = False,
    index = False
)


end_time = time.time()
run_time = (end_time - start_time) / 60 / 60
print "Process completed in %.1f hours" % (run_time)
