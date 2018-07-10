import os
import pandas as pd
import numpy as np

DATA_DIR   = '/home/jrca253/DATA/IlluminaSupportFiles/'
OUT_DIR    = '/home/jrca253/DATA/EA11101_2011-09-28_FinalReport_1_to_16/PLINK_FILES/'
locus_file = DATA_DIR + 'HumanOmni2.5-8v1_C_LocusReport.txt'
out_file   = OUT_DIR + "locus_dict.txt"


##### LOAD LOCUS REPORT #####
print "Loading locus report at " + locus_file

cols = [
    'Name',
    'Chr',
    'Position',
    'Call Freq'
]

df = pd.read_table(
    locus_file,
    sep = '\t',
    header = 0,
    usecols = cols,
    dtype = {
        'Name': object, 
        'Chr': object, 
        'Position': int, 
        'Call Freq': float
        }
)

print "Starting SNPs: " + str(df.shape[0])



##### DROP NON-SNPs FROM NAME COLUMN #####
print "Dropping NON-SNPs from name column ..."
df.drop( df.loc[ df['Name'] == 'Extension'].index, inplace=True )
df.drop( df.loc[ df['Name'] == 'Hybridization'].index, inplace=True )
df.drop( df.loc[ df['Name'] == 'Non-Polymorphic'].index, inplace=True )
df.drop( df.loc[ df['Name'] == 'Non-Specific Binding'].index, inplace=True )
df.drop( df.loc[ df['Name'] == 'Staining'].index, inplace=True )
df.drop( df.loc[ df['Name'] == 'Stringency'].index, inplace=True )
print "SNP count: " + str(df.shape[0])



##### DROP DUPLICATES AND NANS #####
print "Dropping NaNs ..."
df.dropna( inplace=True )
print "SNP count: " + str(df.shape[0])



##### DROP ROWS WHERE CHR=0, or BP=0 #####
print "Dropping SNPs where Chr = 0 and/or BP = 0 ..."
df.drop( df.loc[ df['Chr'] == '0'].index, inplace=True )
df.drop( df.loc[ df['Position'] == 0].index, inplace=True )
print "SNP count: " + str(df.shape[0])



##### FIND STILL DUPLICATED SNPs #####
print "Checking to see if there are still duplicate SNP names ... "
df_dup = df.loc[ df.duplicated(subset='Name', keep=False) ].copy()
df_dup.sort_values('Name', inplace=True)

print "Saving list ..."
df_dup.to_csv(
    "locus_report_duplicates.txt",
    sep = '\t',
    index = False
)


##### SAVE RELEVANT COLUMNS AS A DICTIONARY #####                                                                                                                                                                                             
if os.path.isfile( out_file ):
    command = 'rm ' + out_file
    print 'Removing old manifest_dict.txt ...'
    os.system(command)
#ENDIF                                                                                                                                                                                                                                        

print "Saving locus report as a dictionary ..."
df.to_csv(out_file,
          sep = '\t',
          index = False
)
