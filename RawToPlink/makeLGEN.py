#!/opt/anaconda/anaconda2/bin/python
import os
import pandas as pd
import numpy as np
import time
import pickle

time_start = time.time()

TEST = False

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
out_file    = DATA_DIR + 'PLINK_FILES/EA11101_2011-09-28.lgen'

if TEST:
    files = ['test.txt']
    out_file = 'EA11101_2011-09-28.lgen'


if os.path.isfile( out_file ):
    os.system('rm ' + out_file)
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



    ##### REORDER COLUMNS FOR LGEN FORMAT #####
    reordered_cols = [
        'Sample Name',
        'Sample ID',
        'SNP Name',
        'Allele1 - Top',
        'Allele2 - Top'
    ]
    df = df[ reordered_cols ]


    ##### CONVERT NON-RS IDs to RS IDs #####
    df['SNP Name'] = df['SNP Name'].map(rs_converter_dict).fillna(df['SNP Name'])


    
    ##### WRITE THE LGEN FILE #####
    print "Appending " + str(f) + " to lgen file ..."
    df.to_csv(
        out_file,
        mode = 'a',
        sep = '\t',
        header = False,
        index = False
    )



    ##### DELETE DATAFRAME, FREE MEMORY #####
    del df


#ENDFOR
        
time_stop = time.time()
duration = (time_stop - time_start) / 60 / 60
print "LGEN CREATED!"
print "Process completed in %.1f hours" % (duration)
