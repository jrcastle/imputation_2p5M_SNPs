import pandas as pd

sample_QC_file = '/mnt/DATA/EA11101_2011-09-28/EA11101_2011-09-28_ArrayQC.csv'

##### LOAD THE QC FILE #####
print "Loading Array QC file ..."
cols=['Sample ID', 'Gender']
df = pd.read_table(
    sample_QC_file,
    sep = ',',
    header = 0,
    usecols = cols,
)



##### TRANSFORM GENDER INTO SEX CODES
df.loc[(df['Gender'] != 'Female') & (df['Gender'] != 'Male'), 'Gender'] = 0
df.loc[df['Gender'] == 'Male', 'Gender'] = 1
df.loc[df['Gender'] == 'Female', 'Gender'] = 2



##### ADD MISSING COLUMNS WITH DEFAULT VALUES #####
df['Family ID']  = df['Sample ID']
df['FID_Father'] = 0
df['FID_Mother'] = 0
df['Phenotype']  = 0



##### REORDER COLUMNS FOR FAM FORMAT #####
reordered_cols = [
    'Family ID',
    'Sample ID',
    'FID_Father',
    'FID_Mother',
    'Gender',
    'Phenotype'
]
df = df[ reordered_cols ]



##### WRITE FAM FILE #####
df.to_csv(
    'EA11101_2011-09-28.fam',
    sep = '\t',
    header = False,
    index = False
)
