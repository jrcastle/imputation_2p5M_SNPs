#!/opt/anaconda/anaconda2/bin/python
import os
import pandas as pd
import numpy as np
import pickle

DATA_DIR      = '/home/jrca253/DATA/IlluminaSupportFiles/'
OUT_DIR       = '/home/jrca253/DATA/EA11101_2011-09-28_FinalReport_1_to_16/PLINK_FILES/' 
manifest_file = DATA_DIR + 'humanomni2.5-8v1_c.csv'
out_file      = OUT_DIR + "manifest_dict.txt"



##### LOAD MANIFEST CSV #####
cols = ['Name', 'Chr', 'MapInfo']

print "Loading manifest file v1.0 ..."
df = pd.read_table(
    manifest_file,
    sep = ',',
    skiprows = range(0,7),
    header = 0,
    usecols = cols,
    dtype = {'Chr': object, 'MapInfo': object}
)



##### DROP NON-SNPs FROM NAME COLUMN #####
print "Dropping NON-SNPs from name column ..."
df.drop( df.loc[ df['Name'] == 'Extension'].index, inplace=True )
df.drop( df.loc[ df['Name'] == 'Hybridization'].index, inplace=True )
df.drop( df.loc[ df['Name'] == 'Non-Polymorphic'].index, inplace=True )
df.drop( df.loc[ df['Name'] == 'Non-Specific Binding'].index, inplace=True )
df.drop( df.loc[ df['Name'] == 'Staining'].index, inplace=True )
df.drop( df.loc[ df['Name'] == 'Stringency'].index, inplace=True )

##### DROP DUPLICATES AND NANS #####
print "Dropping NaNs ..."
df.dropna( inplace=True )


##### DROP ROWS WHERE CHR=0, or BP=0 #####
print "Dropping SNPs where Chr = 0 and/or BP = 0 ..."
df.drop( df.loc[ df['Chr'] == '0'].index, inplace=True )
df.drop( df.loc[ df['MapInfo'] == '0'].index, inplace=True )



##### FIND STILL DUPLICATED SNPs #####
#print "Checking to see if there are still duplicate SNP names ... "
#df_dup = df.loc[ df.duplicated(subset='Name', keep=False) ].copy()
#df_dup.sort_values('Name', inplace=True)
#
#print "Saving list ..."
#df_dup.to_csv("manifest_duplicates.txt",
#              sep = '\t',
#              index = False
#)
#exit()



##### SAVE RELEVANT COLUMNS AS A DICTIONARY #####
if os.path.isfile( out_file ):
    command = 'rm ' + out_file
    print 'Removing old manifest_dict.txt ...'
    os.system(command)
#ENDIF

print "Saving manifest as a dictionary ..."
df.to_csv(out_file,
          sep = '\t',
          index = False
)
