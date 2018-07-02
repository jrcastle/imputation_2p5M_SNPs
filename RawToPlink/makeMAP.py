#!/opt/anaconda/anaconda2/bin/python
import os
import pandas as pd
import numpy as np
import pickle
import time

start_time = time.time()

TEST = 0
DATA_DIR = '/home/jrca253/DATA/EA11101_2011-09-28_FinalReport_1_to_16/'
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
manifest_file  = DATA_DIR + 'PLINK_FILES/manifest_dict.txt'
unmatched_file = DATA_DIR + "PLINK_FILES/unmatchedSNPs.txt"

if TEST:
    files          = ['test.txt']
    out_file       = 'EA11101_2011-09-28.map'
    unmatched_file = 'unmatchedSNPs.txt'



##### LOAD MANIFEST DICTIONARY #####
print "Loading Manifest Dictionary at " + manifest_file
manifest_df = pd.read_table(
    manifest_file,
    sep = '\t',
    header = 0,
    dtype = {'Chr': object, 'MapInfo': int}
)

manifest_dict = manifest_df.set_index('Name').T.to_dict('list')
del manifest_df



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
df_snps['tmp'] = df_snps['SNP'].map(manifest_dict)
rows_before = df_snps.shape[0]



##### FIND UNMATCHED SNPS #####
print "Finding unmatched SNPs"
df2 = df_snps.loc[ df_snps['tmp'].isnull(), 'SNP' ]

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
print "Splitting the dictionary mapping into two columns ..." 
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
