import os
import pandas as pd
import numpy as np
import pickle


manifest_file = '/mnt/DATA/EA11101_2011-09-28/EA11101_2011-09-28/InfiniumOmni2-5-8v1-4_A1.csv'
#manifest_file = 'TESTMANIFEST.csv'

##### LOAD MANIFEST CSV #####
cols = ['Name', 'Chr', 'MapInfo']

print "Loading manifest file ..."
df = pd.read_table(
    manifest_file,
    sep = ',',
    skiprows = range(0,7),
    header = 0,
    usecols = cols,
    dtype = {'Chr': object}
)

##### KEEP ONLY SNPS WITH RS NUMBERS #####
df = df.drop(df[ df['Name'].str.contains("rs") == False ].index)

print df.head()

##### SAVE RELEVANT COLUMNS AS A DICTIONARY #####
if os.path.isfile( "manifest_dict.pkl" ):
    command = 'rm manifest_dict.pkl'
    print 'Removing old manifest_dict.pkl ...'
    os.system(command)
#ENDIF

print "Saving manifest as a dictionary ..."
manifest_dict = df.set_index('Name').T.to_dict('list')
output = open('manifest_dict.pkl', 'wb')
pickle.dump(manifest_dict, output)
