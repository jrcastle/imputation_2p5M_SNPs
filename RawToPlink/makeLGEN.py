#!/opt/anaconda/anaconda2/bin/python
import os
import pandas as pd
import numpy as np
import time
import pickle

time_start = time.time()

TEST = 0

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

rsConverter = DATA_DIR + 'PLINK_FILES/rsConverterDict.pkl'
race_file   = DATA_DIR + 'PLINK_FILES/race_dict.csv'
EUR_file    = DATA_DIR + 'PLINK_FILES/EUR_EA11101_2011-09-28.lgen'
AFR_file    = DATA_DIR + 'PLINK_FILES/AFR_EA11101_2011-09-28.lgen'
ASN_file    = DATA_DIR + 'PLINK_FILES/ASN_EA11101_2011-09-28.lgen'
OTR_file    = DATA_DIR + 'PLINK_FILES/OTR_EA11101_2011-09-28.lgen'


if TEST:
    files    = ['test.txt']
    EUR_file = 'EUR_EA11101_2011-09-28.lgen'
    AFR_file = 'AFR_EA11101_2011-09-28.lgen'
    ASN_file = 'ASN_EA11101_2011-09-28.lgen'
    OTR_file = 'OTR_EA11101_2011-09-28.lgen'
#ENDIF

if os.path.isfile( EUR_file ):
    os.system('rm ' + EUR_file)
#ENDIF 

if os.path.isfile( AFR_file ):
    os.system('rm ' + AFR_file)
#ENDIF

if os.path.isfile( ASN_file ):
    os.system('rm ' + ASN_file)
#ENDIF

if os.path.isfile( OTR_file ):
    os.system('rm ' + OTR_file)
#ENDIF 


##### GENCALL SCORE CUTOFF
# SOURCE: https://support.illumina.com/content/dam/illumina-marketing/documents/services/technote_infinium_genotyping_data_analysis.pdf 
GC_thresh = 0.15 

##### DESIRED COLUMNS #####
cols = [
    'SNP Name',
    'Sample ID',
    'Allele1 - Top',
    'Allele2 - Top',
    'GC Score',
    'Sample Name'
]



##### LOAD RS CONVERTER DICTIONARY #####
print "Loading rs converter dictionary at " + rsConverter
pkl_file = open(rsConverter, 'rb')
rs_converter_dict = pickle.load(pkl_file)



##### LOAD RACE DICTIONARY #####
print "Loading race dictionary at " + race_file
df_race = pd.read_table(race_file,
                        sep = ',',
                        header = 0,
                        dtype = {'Sample ID': object, 'Race': object}
)

race_dict = df_race.set_index('Sample ID').T.to_dict('list')
del df_race

for key in race_dict.keys():
    race_dict[key] = race_dict[key][0]
#ENDFOR



##### LOOP OVER FILES ##### 
for f in files:

    ##### LOAD SNP DATA #####
    print "Loading " + str(f) + " ..." 
    df = pd.read_table(
        f,
        sep = '\t',
        skiprows = range(0,10),
        header = 0,
        usecols = cols
    )


    ##### ADD RACE COLUMN #####
    print "Adding Race column ..."
    df['Race'] = df['Sample ID'].map(race_dict)
    print "A race could not be found for the following Sample IDs:"
    print df['Sample ID'].loc[ df['Race'].isna() ].unique().tolist()
    
    
    ##### SELECT DATA WHERE GC SCORE GREATER THAN THRESHOLD #####
    rows_before_GCCut = df.shape[0]
    df.drop(df[ df['GC Score'] < GC_thresh ].index, inplace = True)
    rows_after_GCCut = df.shape[0]
    reduction_fraction = 100. * (1. - float(rows_after_GCCut) / float(rows_before_GCCut))

    print "GC CUT:\t" + str(GC_thresh)
    print "ROWS BEFORE GC CUT:\t" + str( rows_before_GCCut )
    print "ROWS AFTER GC CUT:\t" + str( rows_after_GCCut )
    print 'DATA REDUCED BY:\t%.1f%%' % reduction_fraction
    


    ##### AFTER CUTS. DROP GC SCORE COLUMN #####
    df.drop( "GC Score",
             axis = 1,
             inplace = True
    )


    ##### CONVERT NON-RS IDs to RS IDs #####
    print "Converting non-rsIDs to rsIDs ..."
    df['SNP Name'] = df['SNP Name'].map(rs_converter_dict).fillna(df['SNP Name'])


    ##### DIVIDE DATAFRAME BY RACE #####
    print "Splitting dataframe by race ..."
    df_EUR = df.loc[ df['Race'] == 'White' ]
    df_AFR = df.loc[ df['Race'] == 'Black or African American' ]
    df_ASN = df.loc[ df['Race'] == 'Asian' ]
    df_OTR = df.loc[ (df['Race'] != 'Black or African American') & (df['Race'] != 'White') & (df['Race'] != 'Asian') ]
    del df

    
    ##### DROP THE RACE COLUMN #####
    df_EUR.drop("Race",
                axis = 1,
                inplace = True
    )

    df_AFR.drop("Race",
                axis = 1,
                inplace = True
    )

    df_ASN.drop("Race",
		axis = 1,
                inplace = True
    )

    df_OTR.drop("Race",
		axis = 1,
                inplace = True
    )
    
    ##### REORDER COLUMNS FOR LGEN FORMAT #####
    reordered_cols = [
        'Sample Name',
        'Sample ID',
        'SNP Name',
        'Allele1 - Top',
        'Allele2 - Top'
    ]
    
    df_EUR = df_EUR[ reordered_cols ]
    df_AFR = df_AFR[ reordered_cols ]
    df_ASN = df_ASN[ reordered_cols ]
    df_OTR = df_OTR[ reordered_cols ]
        
    ##### WRITE THE LGEN FILE #####
    print "Appending " + str(f) + " to EUR file ..."
    df_EUR.to_csv(EUR_file,
                  mode = 'a',
                  sep = '\t',
                  header = False,
                  index = False
    )

    print "Appending " + str(f) + " to AFR file ..."
    df_AFR.to_csv(AFR_file,
                  mode = 'a',
                  sep = '\t',
                  header = False,
                  index = False
    )

    print "Appending " + str(f) + " to ASN file ..."
    df_ASN.to_csv(ASN_file,
                  mode = 'a',
                  sep = '\t',
                  header = False,
                  index = False
    )

    print "Appending " + str(f) + " to OTR file ..."
    df_OTR.to_csv(OTR_file,
                  mode = 'a',
                  sep = '\t',
                  header = False,
                  index = False
    )

    ##### DELETE DATAFRAME, FREE MEMORY #####
    del df_EUR
    del df_AFR
    del df_ASN
    del df_OTR

#ENDFOR
        
time_stop = time.time()
duration = (time_stop - time_start) / 60 / 60
print "LGEN CREATED!"
print "Process completed in %.1f hours" % (duration)
