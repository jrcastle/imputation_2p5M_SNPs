#!/opt/anaconda/anaconda2/bin/python
import os
import pandas as pd
import numpy as np
import pickle

TEST          = False
DATA_DIR      = '/mnt/DATA/EA11101_2011-09-28/EA11101_2011-09-28/EA11101_2011-09-28_FinalReport_1_to_16/PLINK_FILES/'
if TEST:
    DATA_DIR = "./"
#ENDIF

#manifest_file = '/mnt/DATA/EA11101_2011-09-28/EA11101_2011-09-28/HumanOmni25-8v1-2_A1.csv'
manifest_file = '/mnt/DATA/EA11101_2011-09-28/EA11101_2011-09-28/InfiniumOmni2-5-8v1-3_A1.csv'
#manifest_file = '/mnt/DATA/EA11101_2011-09-28/EA11101_2011-09-28/InfiniumOmni2-5-8v1-4_A1.csv' 
out_file      = DATA_DIR + "manifest_dict.pkl"




##### LOAD MANIFEST CSV #####
cols = ['Name', 'Chr', 'MapInfo']

print "Loading manifest file ..."
df = pd.read_table(
    manifest_file,
    sep = ',',
    skiprows = range(0,7),
    header = 0,
    usecols = cols,
    dtype = {'Chr': object, 'MapInfo': object}
)

##### KEEP ONLY SNPS WITH RS NUMBERS #####
#df.drop(df[ df['Name'].str.contains("rs") == False ].index, inplace = True)



##### SAVE RELEVANT COLUMNS AS A DICTIONARY #####
if os.path.isfile( out_file ):
    command = 'rm ' + out_file
    print 'Removing old manifest_dict.pkl ...'
    os.system(command)
#ENDIF

print "Saving manifest as a dictionary ..."
manifest_dict = df.set_index('Name').T.to_dict('list')
output = open(out_file, 'wb')
pickle.dump(manifest_dict, output)
