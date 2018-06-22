#!/opt/anaconda/anaconda2/bin/python
import os
import pandas as pd
import numpy as np
import pickle

out_file    = '/mnt/DATA/EA11101_2011-09-28/EA11101_2011-09-28/EA11101_2011-09-28_FinalReport_1_to_16/PLINK_FILES/rsConverterDict.pkl'
rsIDFile1p2 = '/mnt/DATA/EA11101_2011-09-28/EA11101_2011-09-28/IlluminaSupportFiles/HumanOmni2-5-8-v1-2-A-b138-rsIDs.txt'
rsIDFile1p3 = '/mnt/DATA/EA11101_2011-09-28/EA11101_2011-09-28/IlluminaSupportFiles/InfiniumOmni2-5-8v1-3_A1_b144_rsids.txt'
rsIDFile1p4 = '/mnt/DATA/EA11101_2011-09-28/EA11101_2011-09-28/IlluminaSupportFiles/InfiniumOmni2-5-8v1-4_A1_b150_rsids.txt'



##### LOAD FILES #####
print "Loading " + rsIDFile1p2 + " ..."
df1p2 = pd.read_table(rsIDFile1p2,
                      sep = '\t',
                      header = 0
)

print "Loading " + rsIDFile1p3 + " ..."
df1p3 = pd.read_table(rsIDFile1p3,
                      sep = '\t',
                      header = 0
)

print "Loading " + rsIDFile1p4 + " ..."
df1p4 = pd.read_table(rsIDFile1p4,
                      sep = '\t',
                      header = 0
)



##### DROP UNMATCHED SNPs #####
df1p2.drop(df1p2[ df1p2['RsID'] == "." ].index, inplace = True)
df1p3.drop(df1p3[ df1p3['RsID'] == "." ].index, inplace = True)
df1p4.drop(df1p4[ df1p4['RsID'] == "." ].index, inplace = True)



##### APPEND #####
print "Merging ..."
df = pd.concat( [df1p2, df1p3, df1p4], ignore_index = True )
del df1p2
del df1p3
del df1p4



##### DROP DUPLICATES (kgp maps that HAVE NOT changed between versions) #####
print "Dropping duplicates ..."
df.drop_duplicates(keep = 'last', inplace = True)



##### IF KGP MAPPING CHANGED BETWEEN VERSIONS, KEEP THE MOST RECENT #####
df.drop_duplicates(['Name'], keep='last',inplace = True) 



##### SAVE TO DICTIONARY #####
if os.path.isfile( out_file ):
    command = 'rm ' + out_file
    print 'Removing old dictionary ...'
    os.system(command)
#ENDIF

out_dict = df.set_index('Name').T.to_dict('list')

##### CONVERT DICTIONARY ELEMENTS FROM LIST TO STRING #####
for key in out_dict.keys():
    out_dict[key] = out_dict[key][0]
#ENDFOR

print "Saving as a dictionary ..."
output = open(out_file, 'wb')
pickle.dump(out_dict, output)
