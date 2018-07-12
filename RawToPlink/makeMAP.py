#!/opt/anaconda/anaconda2/bin/python
import os
import pandas as pd
import numpy as np
import time

start_time = time.time()

TEST = 0
DATA_DIR = '/home/jrca253/DATA/EA11101_2011-09-28_FinalReport_1_to_16/'
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

##### DICTIONARY FILES #####
race_file       = DATA_DIR + 'PLINK_FILES/race_dict.csv'
annotation_file = DATA_DIR + 'PLINK_FILES/annotation_dict.txt'



##### FILES NOTING SNPs TO BE REMOVED #####
snp_file_EUR = DATA_DIR + 'PLINK_FILES/snps_to_be_removed_EUR.txt'
snp_file_AFR = DATA_DIR + 'PLINK_FILES/snps_to_be_removed_AFR.txt'
snp_file_ASN = DATA_DIR + 'PLINK_FILES/snps_to_be_removed_ASN.txt'
snp_file_OTR = DATA_DIR + 'PLINK_FILES/snps_to_be_removed_OTR.txt'



##### FINAL MAP FILES #####
EUR_file = DATA_DIR + 'PLINK_FILES/EUR_EA11101_2011-09-28.map'
AFR_file = DATA_DIR + 'PLINK_FILES/AFR_EA11101_2011-09-28.map'
ASN_file = DATA_DIR + 'PLINK_FILES/ASN_EA11101_2011-09-28.map'
OTR_file = DATA_DIR + 'PLINK_FILES/OTR_EA11101_2011-09-28.map'



if TEST:
    files          = ['test.txt']

    snp_file_EUR = 'snps_to_be_removed_EUR.txt'
    snp_file_AFR = 'snps_to_be_removed_AFR.txt'
    snp_file_ASN = 'snps_to_be_removed_ASN.txt'
    snp_file_OTR = 'snps_to_be_removed_OTR.txt'

    EUR_file     = 'EUR_EA11101_2011-09-28.map'
    AFR_file     = 'AFR_EA11101_2011-09-28.map'
    ASN_file     = 'ASN_EA11101_2011-09-28.map'
    OTR_file     = 'OTR_EA11101_2011-09-28.map'
#ENDIF



##### REMOVE OLD FILES #####
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

if os.path.isfile( snp_file_EUR ):
    os.system('rm ' + snp_file_EUR)
#ENDIF 

if os.path.isfile( snp_file_AFR ):
    os.system('rm ' + snp_file_AFR)
#ENDIF

if os.path.isfile( snp_file_ASN ):
    os.system('rm ' + snp_file_ASN)
#ENDIF  

if os.path.isfile( snp_file_OTR ):
    os.system('rm ' + snp_file_OTR)
#ENDIF  



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



##### LOAD ANNOTATION FILE DICTIONARY #####
print "Loading annotation file at " + annotation_file
df_annotation = pd.read_table(
    annotation_file,
    sep = '\t',
    header = 0,
    dtype = {
        'Name': object,
        'Chr': object,
        'MapInfo': int
        }
)

annotation_dict = df_annotation.set_index('Name').T.to_dict('list')
del df_annotation


##### INITIALIZE THE SNP ARRAY #####
all_snps_array_EUR = np.array( [] )
all_snps_array_AFR = np.array( [] )
all_snps_array_ASN = np.array( [] )
all_snps_array_OTR = np.array( [] )



##### DESIRED COLUMNS #####
cols = [
    'SNP Name',
    'Sample ID'
]



for f in files:

    ##### LOAD DATA FILES #####
    print "Loading " + f + " ..."
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



    ##### DIVIDE DATAFRAME BY RACE #####
    print "Splitting dataframe by race ..."
    df_EUR = df.loc[ df['Race'] == 'White' ]
    df_AFR = df.loc[ df['Race'] == 'Black or African American' ]
    df_ASN = df.loc[ df['Race'] == 'Asian' ]
    df_OTR = df.loc[ (df['Race'] != 'Black or African American') & (df['Race'] != 'White') & (df['Race'] != 'Asian') ]
    del df



    ##### GET UNIQUE SNPS #####
    print "Getting list of unique SNPs ..."
    snp_list_EUR = df_EUR['SNP Name'].unique()
    snp_list_AFR = df_AFR['SNP Name'].unique()
    snp_list_ASN = df_ASN['SNP Name'].unique()
    snp_list_OTR = df_OTR['SNP Name'].unique()

    all_snps_array_EUR = np.union1d(all_snps_array_EUR, snp_list_EUR)
    all_snps_array_AFR = np.union1d(all_snps_array_AFR, snp_list_AFR)
    all_snps_array_ASN = np.union1d(all_snps_array_ASN, snp_list_ASN)
    all_snps_array_OTR = np.union1d(all_snps_array_OTR, snp_list_OTR)



    ##### FREE MEMORY #####
    del df_EUR
    del df_AFR
    del df_ASN
    del df_OTR

    del snp_list_EUR
    del snp_list_AFR
    del snp_list_ASN
    del snp_list_OTR
#ENDFOR



##### CREATE NEW DATAFRAME WITH UNIQUE SNPS
df_snps_EUR = pd.DataFrame(all_snps_array_EUR, columns = ["SNP"])
df_snps_AFR = pd.DataFrame(all_snps_array_AFR, columns = ["SNP"])
df_snps_ASN = pd.DataFrame(all_snps_array_ASN, columns = ["SNP"])
df_snps_OTR = pd.DataFrame(all_snps_array_OTR, columns = ["SNP"])

del all_snps_array_EUR
del all_snps_array_AFR
del all_snps_array_ASN
del all_snps_array_OTR



##### COUNT STARTING SNPs #####
start_snps_EUR = df_snps_EUR.shape[0]
start_snps_AFR = df_snps_AFR.shape[0]
start_snps_ASN = df_snps_ASN.shape[0]
start_snps_OTR = df_snps_OTR.shape[0]



##### MATCH SNP TO CHROMOSOME, BP POSITION #####
print "Matching SNP to chromosome and base-pair position ..."
df_snps_EUR['tmp'] = df_snps_EUR['SNP'].map(annotation_dict)
df_snps_AFR['tmp'] = df_snps_AFR['SNP'].map(annotation_dict)
df_snps_ASN['tmp'] = df_snps_ASN['SNP'].map(annotation_dict)
df_snps_OTR['tmp'] = df_snps_OTR['SNP'].map(annotation_dict)



##### FIND UNMATCHED SNPS #####
print "Finding unmatched EUR SNPs"
df_tmp = df_snps_EUR.loc[ df_snps_EUR['tmp'].isnull(), 'SNP' ]
unmatched_EUR = df_tmp.shape[0]
df_snps_EUR.dropna(inplace=True)

print "List of Unmatched EUR SNPs appended to " + snp_file_EUR
df_tmp.to_csv(
    snp_file_EUR,
    mode = 'a',
    sep = '\t',
    header = False,
    index = False
)


print "Finding unmatched AFR SNPs"
df_tmp = df_snps_AFR.loc[ df_snps_AFR['tmp'].isnull(), 'SNP' ]
unmatched_AFR = df_tmp.shape[0]
df_snps_AFR.dropna(inplace=True)

print "List of Unmatched AFR SNPs appended to " + snp_file_AFR
df_tmp.to_csv(
    snp_file_AFR,
    mode = 'a',
    sep = '\t',
    header = False,
    index = False
)


print "Finding unmatched ASN SNPs"
df_tmp = df_snps_ASN.loc[ df_snps_ASN['tmp'].isnull(), 'SNP' ]
unmatched_ASN = df_tmp.shape[0]
df_snps_ASN.dropna(inplace=True)

print "List of Unmatched ASN SNPs appended to " + snp_file_ASN
df_tmp.to_csv(
    snp_file_ASN,
    mode = 'a',
    sep = '\t',
    header = False,
    index = False
)


print "Finding unmatched OTR SNPs"
df_tmp = df_snps_OTR.loc[ df_snps_OTR['tmp'].isnull(), 'SNP' ]
unmatched_OTR = df_tmp.shape[0]
df_snps_OTR.dropna(inplace=True)

print "List of Unmatched OTR SNPs appended to " + snp_file_OTR
df_tmp.to_csv(
    snp_file_OTR,
    mode = 'a',
    sep = '\t',
    header = False,
    index = False
)

del df_tmp



##### SPLIT THE DICTIONARY MAPPING INTO MULTIPLE COLUMNS #####
print "Splitting the dictionary mapping into multiple columns ..." 
df_snps_EUR[['Chromosome','BasePairPos']] = pd.DataFrame(df_snps_EUR.tmp.values.tolist(), index = df_snps_EUR.index)
df_snps_AFR[['Chromosome','BasePairPos']] = pd.DataFrame(df_snps_AFR.tmp.values.tolist(), index = df_snps_AFR.index)
df_snps_ASN[['Chromosome','BasePairPos']] = pd.DataFrame(df_snps_ASN.tmp.values.tolist(), index = df_snps_ASN.index)
df_snps_OTR[['Chromosome','BasePairPos']] = pd.DataFrame(df_snps_OTR.tmp.values.tolist(), index = df_snps_OTR.index)

df_snps_EUR.drop("tmp", axis = 1, inplace = True)
df_snps_AFR.drop("tmp", axis = 1, inplace = True)
df_snps_ASN.drop("tmp", axis = 1, inplace = True)
df_snps_OTR.drop("tmp", axis = 1, inplace = True)



##### DROP SNPs NOT ON CHR 1-22 or X #####
df_tmp = df_snps_EUR.loc[ (df_snps_EUR['Chromosome'] == 'Y') | (df_snps_EUR['Chromosome'] == 'XY') | (df_snps_EUR['Chromosome'] == 'MT') ].copy()
print "List of EUR SNPs not on Chr 1-22 or X appended to " + snp_file_EUR
df_tmp['SNP'].to_csv(
    snp_file_EUR,
    mode = 'a',
    sep = '\t',
    header = False,
    index = False
)

df_tmp = df_snps_AFR.loc[ (df_snps_AFR['Chromosome'] == 'Y') | (df_snps_AFR['Chromosome'] == 'XY') | (df_snps_AFR['Chromosome'] == 'MT') ].copy()
print "List of AFR SNPs not on Chr 1-22 or X appended to " + snp_file_AFR
df_tmp['SNP'].to_csv(
    snp_file_AFR,
    mode = 'a',
    sep = '\t',
    header = False,
    index = False
)

df_tmp = df_snps_ASN.loc[ (df_snps_ASN['Chromosome'] == 'Y') | (df_snps_ASN['Chromosome'] == 'XY') | (df_snps_ASN['Chromosome'] == 'MT') ].copy()
print "List of ASN SNPs not on Chr 1-22 or X appended to " + snp_file_ASN
df_tmp['SNP'].to_csv(
    snp_file_ASN,
    mode = 'a',
    sep = '\t',
    header = False,
    index = False
)

df_tmp = df_snps_OTR.loc[ (df_snps_OTR['Chromosome'] == 'Y') | (df_snps_OTR['Chromosome'] == 'XY') | (df_snps_OTR['Chromosome'] == 'MT') ].copy()
print "List of OTR SNPs not on Chr 1-22 or X appended to " + snp_file_OTR
df_tmp['SNP'].to_csv(
    snp_file_OTR,
    mode = 'a',
    sep = '\t',
    header = False,
    index = False
)

del df_tmp

print "Dropping SNPs not on Chr 1-22 or X " 
chr_snps_EUR = df_snps_EUR.loc[ (df_snps_EUR['Chromosome'] == 'Y') | (df_snps_EUR['Chromosome'] == 'XY') | (df_snps_EUR['Chromosome'] == 'MT') ].shape[0]
chr_snps_AFR = df_snps_AFR.loc[ (df_snps_AFR['Chromosome'] == 'Y') | (df_snps_AFR['Chromosome'] == 'XY') | (df_snps_AFR['Chromosome'] == 'MT') ].shape[0]
chr_snps_ASN = df_snps_ASN.loc[ (df_snps_ASN['Chromosome'] == 'Y') | (df_snps_ASN['Chromosome'] == 'XY') | (df_snps_ASN['Chromosome'] == 'MT') ].shape[0]
chr_snps_OTR = df_snps_OTR.loc[ (df_snps_OTR['Chromosome'] == 'Y') | (df_snps_OTR['Chromosome'] == 'XY') | (df_snps_OTR['Chromosome'] == 'MT') ].shape[0]

df_snps_EUR.drop( df_snps_EUR.loc[ (df_snps_EUR['Chromosome'] == 'Y') | (df_snps_EUR['Chromosome'] == 'XY') | (df_snps_EUR['Chromosome'] == 'MT') ].index, inplace = True )
df_snps_AFR.drop( df_snps_AFR.loc[ (df_snps_AFR['Chromosome'] == 'Y') | (df_snps_AFR['Chromosome'] == 'XY') | (df_snps_AFR['Chromosome'] == 'MT') ].index, inplace = True )
df_snps_ASN.drop( df_snps_ASN.loc[ (df_snps_ASN['Chromosome'] == 'Y') | (df_snps_ASN['Chromosome'] == 'XY') | (df_snps_ASN['Chromosome'] == 'MT') ].index, inplace = True )
df_snps_OTR.drop( df_snps_OTR.loc[ (df_snps_OTR['Chromosome'] == 'Y') | (df_snps_OTR['Chromosome'] == 'XY') | (df_snps_OTR['Chromosome'] == 'MT') ].index, inplace = True )



##### ADD GENETIC DISTANCE (MORGANS) AND SET TO 0 #####
print "Adding genetic distance column ..."
df_snps_EUR['GeneticDist'] = 0
df_snps_AFR['GeneticDist'] = 0
df_snps_ASN['GeneticDist'] = 0
df_snps_OTR['GeneticDist'] = 0



##### REORGANIZE COLUMNS AND WRITE #####
print "Reordering columns ..." 
reordered_cols = [
    'Chromosome',
    'SNP',
    'GeneticDist',
    'BasePairPos'
]

df_snps_EUR = df_snps_EUR[ reordered_cols ]
df_snps_AFR = df_snps_AFR[ reordered_cols ]
df_snps_ASN = df_snps_ASN[ reordered_cols ]
df_snps_OTR = df_snps_OTR[ reordered_cols ]

end_snps_EUR = df_snps_EUR.shape[0]
end_snps_AFR = df_snps_AFR.shape[0]
end_snps_ASN = df_snps_ASN.shape[0]
end_snps_OTR = df_snps_OTR.shape[0]



##### WRITE MAP FILES #####
print "Writing EUR MAP file to " + EUR_file
df_snps_EUR.to_csv(
    EUR_file,
    sep = '\t',
    header = False,
    index = False
)

print "Writing AFR MAP file to " + AFR_file
df_snps_AFR.to_csv(
    AFR_file,
    sep = '\t',
    header = False,
    index = False
)

print "Writing ASN MAP file to " + ASN_file
df_snps_ASN.to_csv(
    ASN_file,
    sep = '\t',
    header = False,
    index = False
)

print "Writing OTR MAP file to " + OTR_file
df_snps_OTR.to_csv(
    OTR_file,
    sep = '\t',
    header = False,
    index = False
)



##### SUMMARY FILE #####
print "################### EUR QC REPORT ###################"
print "Starting SNPs:              \t%.0f" % (start_snps_EUR)
print "SNPs not in annotation file:\t%.0f" % (unmatched_EUR)
print "SNPs not on Chr 1-22 or X:  \t%.0f" % (chr_snps_EUR)
print "SNPs after QC:              \t%.0f" % (end_snps_EUR)
print "\n\n\n"

print "################### AFR QC REPORT ###################"
print "Starting SNPs:              \t%.0f" % (start_snps_AFR)
print "SNPs not in annotation file:\t%.0f" % (unmatched_AFR)
print "SNPs not on Chr 1-22 or X:  \t%.0f" % (chr_snps_AFR)
print "SNPs after QC:              \t%.0f" % (end_snps_AFR)
print "\n\n\n"

print "################### ASN QC REPORT ###################"
print "Starting SNPs:              \t%.0f" % (start_snps_ASN)
print "SNPs not in annotation file:\t%.0f" % (unmatched_ASN)
print "SNPs not on Chr 1-22 or X:  \t%.0f" % (chr_snps_ASN)
print "SNPs after QC:              \t%.0f" % (end_snps_ASN)
print "\n\n\n"

print "################### OTR QC REPORT ###################"
print "Starting SNPs:              \t%.0f" % (start_snps_OTR)
print "SNPs not in annotation file:\t%.0f" % (unmatched_OTR)
print "SNPs not on Chr 1-22 or X:  \t%.0f" % (chr_snps_OTR)
print "SNPs after QC:              \t%.0f" % (end_snps_OTR)
print "\n"


end_time = time.time()
run_time = (end_time - start_time) / 60 / 60
print "Process completed in %.1f hours" % (run_time)
