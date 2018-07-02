#!/opt/anaconda/anaconda2/bin/python
import pandas as pd

sample_QC_file = '/home/jrca253/DATA/EA11101_2011-09-28_ArrayQC.csv'
DATA_DIR       = '/home/jrca253/DATA/EA11101_2011-09-28_FinalReport_1_to_16/PLINK_FILES/'
race_file      = DATA_DIR + 'race_dict.csv'
EUR_file       = DATA_DIR + 'EUR_EA11101_2011-09-28.fam'
AFR_file       = DATA_DIR + 'AFR_EA11101_2011-09-28.fam'
ASN_file       = DATA_DIR + 'ASN_EA11101_2011-09-28.fam'
OTR_file       = DATA_DIR + 'OTR_EA11101_2011-09-28.fam'


##### LOAD THE QC FILE #####
print "Loading Array QC file ..."
cols=['Sample ID', 'Gender']
df = pd.read_table(
    sample_QC_file,
    sep = ',',
    header = 0,
    usecols = cols,
)



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



##### ADD RACE COLUMN #####                                                                                                                                                                                                               
print "Adding Race column ..."
df['Race'] = df['Sample ID'].map(race_dict)
print "A race could not be found for the following Sample IDs:"
print df['Sample ID'].loc[ df['Race'].isna() ].unique().tolist()



##### TRANSFORM GENDER INTO SEX CODES
df.loc[(df['Gender'] != 'Female') & (df['Gender'] != 'Male'), 'Gender'] = 0
df.loc[df['Gender'] == 'Male', 'Gender'] = 1
df.loc[df['Gender'] == 'Female', 'Gender'] = 2



##### ADD MISSING COLUMNS WITH DEFAULT VALUES #####
df['Family ID']  = df['Sample ID']
df['FID_Father'] = 0
df['FID_Mother'] = 0
df['Phenotype']  = 0


##### SPLIT DATAFRAME BY RACE #####
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



##### REORDERCOLUMNS FOR FAM FORMAT #####
reordered_cols = [
    'Family ID',
    'Sample ID',
    'FID_Father',
    'FID_Mother',
    'Gender',
    'Phenotype'
]

df_EUR = df_EUR[ reordered_cols ]
df_AFR = df_AFR[ reordered_cols ]
df_ASN = df_ASN[ reordered_cols ]
df_OTR = df_OTR[ reordered_cols ]



##### WRITE FAM FILE #####
df_EUR.to_csv(
    EUR_file,
    sep = '\t',
    header = False,
    index = False
)

df_AFR.to_csv(
    AFR_file,
    sep = '\t',
    header = False,
    index = False
)

df_ASN.to_csv(
    ASN_file,
    sep = '\t',
    header = False,
    index = False
)

df_OTR.to_csv(
    OTR_file,
    sep = '\t',
    header = False,
    index = False
)
