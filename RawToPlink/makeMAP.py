
#!/opt/anaconda/anaconda2/bin/python
import os
import pandas as pd
import numpy as np
import pickle
import time

start_time = time.time()

TEST          = False
DATA_DIR='/mnt/DATA/EA11101_2011-09-28/EA11101_2011-09-28/EA11101_2011-09-28_FinalReport_1_to_16/'
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
unmatched_file = DATA_DIR + "PLINK_FILES/unmatchedSNPs.txt"

if TEST:
    files          = ['test.txt', 'test2.txt']
    out_file       = 'EA11101_2011-09-28.map'
    unmatched_file = 'unmatchedSNPs.txt'


if not os.path.isfile( manifest_file ):
    print 'Manifest dictionary not found. Make the manifest dictionary file first and then run this script.'
    exit()
#ENDIF 

##### LOAD MANIFEST DICTIONARY #####
print "Loading manifest dictionary at " + manifest_file
pkl_file = open(manifest_file, 'rb')
manifest_dict = pickle.load(pkl_file)


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
df = pd.DataFrame(all_snps_array, columns = ["SNP"])



##### MATCH SNP TO CHROMOSOME AND BP POSITION #####
print "Matching SNP to chromosome and base-pair position ..."
df['tmp'] = df['SNP'].map(manifest_dict)
rows_before = df.shape[0]



##### FIND UNMATCHED SNPS #####
print "Finding unmatched SNPs"
df2 = df[df['tmp'].isnull()]
df2 = df2[['SNP']]

df.dropna(inplace=True)
rows_after = df.shape[0]
print "%.0f SNPS ARE UNMATCHED!" % (rows_before - rows_after)

print "List of Unmatched SNPs saved in " + unmatched_file
df2.to_csv(
    unmatched_file,
    sep = '\t',
    header = False,
    index = False
)



##### SPLIT THE DICTIONARY MAPPING INTO TWO COLUMNS #####
df[['Chromosome','BasePairPos']] = pd.DataFrame(df.tmp.values.tolist(), index= df.index)

df.drop("tmp",
        axis = 1,
        inplace = True
)



##### ADD GENETIC DISTANCE (MORGANS) AND SET TO 0 #####
df['GeneticDist'] = 0



##### REORGANIZE COLUMNS AND WRITE #####
reordered_cols = [
    'Chromosome',
    'SNP',
    'GeneticDist',
    'BasePairPos'
]

df = df[ reordered_cols ]

print "Writing MAP file to " + out_file
df.to_csv(
    out_file,
    sep = '\t',
    header = False,
    index = False
)


end_time = time.time()
run_time = (end_time - start_time) / 60 / 60
print "Process completed in %.1f hours" % (run_time)
