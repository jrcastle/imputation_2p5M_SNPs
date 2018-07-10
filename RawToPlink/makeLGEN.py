#!/home/jrca253/anaconda2/bin/python
import os
import pandas as pd
import numpy as np
import time

time_start = time.time()

TEST = 0

DATA_DIR='/home/jrca253/DATA/EA11101_2011-09-28_FinalReport_1_to_16/'
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

race_file            = DATA_DIR + 'PLINK_FILES/race_dict.csv'
low_call_sample_file = DATA_DIR + 'PLINK_FILES/low_call_rate_samples.txt'

snp_file_EUR = DATA_DIR + 'PLINK_FILES/snps_to_be_removed_EUR.txt'
snp_file_AFR = DATA_DIR + 'PLINK_FILES/snps_to_be_removed_AFR.txt'
snp_file_ASN = DATA_DIR + 'PLINK_FILES/snps_to_be_removed_ASN.txt'
snp_file_OTR = DATA_DIR + 'PLINK_FILES/snps_to_be_removed_OTR.txt'

EUR_file       = DATA_DIR + 'PLINK_FILES/EUR_EA11101_2011-09-28.lgen'
AFR_file       = DATA_DIR + 'PLINK_FILES/AFR_EA11101_2011-09-28.lgen'
ASN_file       = DATA_DIR + 'PLINK_FILES/ASN_EA11101_2011-09-28.lgen'
OTR_file       = DATA_DIR + 'PLINK_FILES/OTR_EA11101_2011-09-28.lgen'


if TEST:
    files    = ['test.txt']

    snp_file_EUR = 'snps_to_be_removed_EUR.txt'
    snp_file_AFR = 'snps_to_be_removed_AFR.txt'
    snp_file_ASN = 'snps_to_be_removed_ASN.txt'
    snp_file_OTR = 'snps_to_be_removed_OTR.txt'

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
    'Allele1 - Forward',
    'Allele2 - Forward',
    'GC Score',
    'Sample Name'
]



#### LOAD LOW CALL SAMPLE FILE #####
print "Loading low-call rate sample file ..."
df_low_call = pd.read_table(
    low_call_sample_file,
    sep = '\t',
    header = None,
    names = ['Sample ID']
)



##### LOAD RACE DICTIONARY #####
print "Loading race dictionary at " + race_file
df_race = pd.read_table(
    race_file,
    sep = ',',
    header = 0,
    dtype = {'Sample ID': object, 'Race': object}
)

race_dict = df_race.set_index('Sample ID').T.to_dict('list')
del df_race

for key in race_dict.keys():
    race_dict[key] = race_dict[key][0]
#ENDFOR



##### LOAD SNPs TO BE REMOVED #####
print "Loading list of EUR SNPs to be removed..."
unmatched_df_EUR = pd.read_table(
    snp_file_EUR,
    sep = '\t',
    header = None,
    names = ['SNP']
)

print "Loading list of AFR SNPs to be removed..."
unmatched_df_AFR = pd.read_table(
    snp_file_AFR,
    sep = '\t',
    header = None,
    names = ['SNP']
)

print "Loading list of ASN SNPs to be removed..."
unmatched_df_ASN = pd.read_table(
    snp_file_ASN,
    sep = '\t',
    header = None,
    names = ['SNP']
)

print "Loading list of OTR SNPs to be removed..."
unmatched_df_OTR = pd.read_table(
    snp_file_OTR,
    sep = '\t',
    header = None,
    names = ['SNP']
)


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



    ##### DROP SAMPLES WITH LOW CALL RATE #####
    print "Dropping samples w/ call rate < 90% ..." 
    df.drop( df.loc[ df['Sample ID'].isin( df_low_call['Sample ID'] ) ].index, inplace=True )



    ##### ADD RACE COLUMN #####
    print "Adding Race column ..."
    df['Race'] = df['Sample ID'].map(race_dict)
    print "A race could not be found for the following Sample IDs:"
    print df['Sample ID'].loc[ df['Race'].isna() ].unique().tolist()
    

    
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
    


    ##### DROP SNPs THAT FAILED QC IN THE MAP STEP #####
    print "Dropping SNPs that failed QC in the map step ... "
    df_EUR.drop( df_EUR.loc[ df_EUR['SNP Name'].isin( unmatched_df_EUR['SNP'] ) ].index, inplace=True )
    df_AFR.drop( df_AFR.loc[ df_AFR['SNP Name'].isin( unmatched_df_AFR['SNP'] ) ].index, inplace=True )
    df_ASN.drop( df_ASN.loc[ df_ASN['SNP Name'].isin( unmatched_df_ASN['SNP'] ) ].index, inplace=True )
    df_OTR.drop( df_OTR.loc[ df_OTR['SNP Name'].isin( unmatched_df_OTR['SNP'] ) ].index, inplace=True )



    ##### DROP SNPs THAT ONLY HAVE 1 ALLELE #####
    print "Dropping EUR SNPs that only have 1 allele ..."
    snps_before = df_EUR.shape[0]
    df_EUR.drop( df_EUR.loc[ df_EUR['Allele1 - Forward'] == "-" ].index, inplace=True )
    df_EUR.drop( df_EUR.loc[ df_EUR['Allele2 - Forward'] == "-" ].index, inplace=True )
    snps_after = df_EUR.shape[0]
    print "%i SNPs dropped" % (snps_before - snps_after)

    print "Dropping AFR SNPs that only have 1 allele ..."
    snps_before = df_AFR.shape[0]
    df_AFR.drop( df_AFR.loc[ df_AFR['Allele1 - Forward'] == "-" ].index, inplace=True )
    df_AFR.drop( df_AFR.loc[ df_AFR['Allele2 - Forward'] == "-" ].index, inplace=True )
    snps_after = df_AFR.shape[0]
    print "%i SNPs dropped" % (snps_before - snps_after)

    print "Dropping ASN SNPs that only have 1 allele ..."
    snps_before = df_ASN.shape[0]
    df_ASN.drop( df_ASN.loc[ df_ASN['Allele1 - Forward'] == "-" ].index, inplace=True )
    df_ASN.drop( df_ASN.loc[ df_ASN['Allele2 - Forward'] == "-" ].index, inplace=True )
    snps_after = df_ASN.shape[0]
    print "%i SNPs dropped" % (snps_before - snps_after)

    print "Dropping OTR SNPs that only have 1 allele ..."
    snps_before = df_OTR.shape[0]
    df_OTR.drop( df_OTR.loc[ df_OTR['Allele1 - Forward'] == "-" ].index, inplace=True )
    df_OTR.drop( df_OTR.loc[ df_OTR['Allele2 - Forward'] == "-" ].index, inplace=True )
    snps_after = df_OTR.shape[0]
    print "%i SNPs dropped" % (snps_before - snps_after)



    ##### SELECT DATA WHERE GC SCORE GREATER THAN THRESHOLD #####
    print "Dropping EUR SNPs that are below the GC threshold ..."
    snps_before = df_EUR.shape[0]
    df_EUR.drop(df_EUR[ df_EUR['GC Score'] < GC_thresh ].index, inplace = True)
    snps_after = df_EUR.shape[0]
    print "%i SNPs dropped" % (snps_before - snps_after)

    print "Dropping AFR SNPs that are below the GC threshold ..."
    snps_before = df_AFR.shape[0]
    df_AFR.drop(df_AFR[ df_AFR['GC Score'] < GC_thresh ].index, inplace = True)
    snps_after = df_AFR.shape[0]
    print "%i SNPs dropped" % (snps_before - snps_after)

    print "Dropping ASN SNPs that are below the GC threshold ..."
    snps_before = df_ASN.shape[0]
    df_ASN.drop(df_ASN[ df_ASN['GC Score'] < GC_thresh ].index, inplace = True)
    snps_after = df_ASN.shape[0]
    print "%i SNPs dropped" % (snps_before - snps_after)

    print "Dropping OTR SNPs that are below the GC threshold ..."
    snps_before = df_OTR.shape[0]
    df_OTR.drop(df_OTR[ df_OTR['GC Score'] < GC_thresh ].index, inplace = True)
    snps_after = df_OTR.shape[0]
    print "%i SNPs dropped" % (snps_before - snps_after)



    ##### AFTER CUTS. DROP GC SCORE COLUMN #####
    df_EUR.drop("GC Score", axis = 1, inplace = True)
    df_AFR.drop("GC Score", axis = 1, inplace = True)
    df_ASN.drop("GC Score", axis = 1, inplace = True)
    df_OTR.drop("GC Score", axis = 1, inplace = True)



    ##### REORDER COLUMNS FOR LGEN FORMAT #####
    reordered_cols = [
        'Sample Name',
        'Sample ID',
        'SNP Name',
        'Allele1 - Forward',
        'Allele2 - Forward'
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
