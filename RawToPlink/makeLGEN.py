import pandas as pd
import numpy as np

DATA_DIR='/mnt/DATA/EA11101_2011-09-28/EA11101_2011-09-28/EA11101_2011-09-28_FinalReport_1_to_16/'
test_file='test.txt'

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

##### LOAD SNP DATA #####
df = pd.read_table(
    test_file,
    sep = '\t',
    skiprows = range(0,10),
    header = 0,
    usecols = cols
)



##### SELECT DATA WHERE GC SCORE GREATER THAN THRESHOLD #####
rows_before_GCCut = df.shape[0]
df = df.drop(df[ df['GC Score'] < GC_thresh ].index)
rows_after_GCCut = df.shape[0]
reduction_fraction = 100. * (1. - float(rows_after_GCCut) / float(rows_before_GCCut))

print "ROWS BEFORE GC CUT:\t" + str( rows_before_GCCut )
print "ROWS AFTER GC CUT:\t" + str( rows_after_GCCut )
print 'DATA REDUCED BY:\t%.0f%%' % reduction_fraction



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



##### WRITE THE LGEN FILE #####
df.to_csv(
    'test.lgen',
    sep = '\t',
    header = False,
    index = False
)



unique_SNPs = df['SNP Name'].unique()
print unique_SNPs[1]
#print df.head() 
