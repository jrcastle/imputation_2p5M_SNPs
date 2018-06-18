import os
import pandas as pd
import pickle

LGEN_file = 'EA11101_2011-09-28.lgen'
manifest_file = 'manifest_dict.pkl'

##### LOOK FOR LGEN AND MANIFEST DICTIONARY #####
if not os.path.isfile( LGEN_file ):
    print 'LGEN not found. Make LGEN file first and then run this script.'
    exit()
#ENDIF

if not os.path.isfile( LGEN_file ):
    print 'Manifest dictionary not found. Make the manifest dictionary file first and then run this script.'
    exit()
#ENDIF 



##### LOAD LGEN FILE #####
print "Loading lgen data ..."
df = pd.read_table(
    LGEN_file,
    sep = '\t',
    header = None,
    usecols = [2]
)
df.columns = ['SNP']



##### DROP DUPLICATES #####
print "Dropping duplicate SNPs ..."
df.drop_duplicates('SNP',
                   keep = 'first',
                   inplace = True
                   )




##### LOAD MANIFEST DICTIONARY #####
print "Loading manifest dictionary ... "
pkl_file = open('manifest_dict.pkl', 'rb')
manifest_dict = pickle.load(pkl_file) 



##### MATCH SNP TO CHROMOSOME AND BP POSITION #####
print "Matching SNP to chromosome and base-pair position ..."
df['tmp'] = df['SNP'].map(manifest_dict)
df[['Chromosome','BasePairPos']] = pd.DataFrame(df['tmp'].values.tolist(), index= df.index)
#df['Chromosome'].replace(manifest_dict, inplace=True)

print df.head()
