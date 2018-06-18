import pandas as pd
import numpy as np
import pickle


manifest_file = '/mnt/DATA/EA11101_2011-09-28/EA11101_2011-09-28/InfiniumOmni2-5-8v1-4_A1.csv'
#manifest_file = 'TESTMANIFEST.csv'

##### LOAD MANIFEST CSV #####
cols = ['SNP', 'Chr', 'MapInfo']

print "Loading manifest file..."
df = pd.read_table(
    manifest_file,
    sep = ',',
    skiprows = range(0,7),
    header = 0,
    usecols = cols,
    dtype = {'SNP': str, 'Chr': str, 'MapInfo': int}
)

##### KEEP ONLY SNPS WITH RS NUMBERS #####
df = df.drop(df[ df['SNP'].str.contains("rs") == False ].index)


##### SAVE RELEVANT COLUMNS AS A DICTIONARY #####
manifest_dict = df.set_index('SNP').T.to_dict('list')
output = open('manifest_dict.pkl', 'wb')

print "Saving manifest as a dictionary..."
pickle.dump(manifest_dict, output)

###### READ PICKLE DICTIONARY
#pkl_file = open('manifest_dict.pkl', 'rb')
#dict = pickle.load(pkl_file)
