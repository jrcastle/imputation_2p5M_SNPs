#!/opt/anaconda/anaconda2/bin/python
import os
import pandas as pd
import numpy as np
import time

time_start = time.time()

TEST = False

sample_list = '/mnt/DATA/EA11101_2011-09-28/EA11101_2011-09-28_ArrayQC.csv'
race_map    = '/mnt/DATA/EA11101_2011-09-28/Race_information_for_SNP_genotyping_data.csv'
out_file    = '/mnt/DATA/EA11101_2011-09-28/EA11101_2011-09-28/EA11101_2011-09-28_FinalReport_1_to_16/PLINK_FILES/race_dict.csv' 



##### LOAD THE SAMPLE AND RACE CSVs INTO DATAFRAMES
print "Loading " +  sample_list + " ..."
df_sample = pd.read_table(sample_list,
                          sep = ',',
                          header = 0,
                          usecols = [1]
)

print "Loading " + race_map + " ..."
df_race = pd.read_table(race_map,
                        sep = ',',
                        header = 0,
	                usecols = [2,5]
)



##### KEEP ONLY THE STRING TO THE RIGHT OF THE "_" IN THE SAMPLE FILE #####
df_sample[['tmp','Decode ID']] = pd.DataFrame(df_sample['Sample ID'].str.split("_").values.tolist(), index= df_sample.index)

df_sample.drop("tmp",
               axis = 1,
               inplace = True
)



##### MAKE RACE DICTIONARY #####
race_dict = df_race.set_index('Decode ID Converted').T.to_dict('list')



##### CONVERT DICTIONARY ELEMENTS FROM LIST TO STRING #####
for key in race_dict.keys():
    race_dict[key] = race_dict[key][0]
#ENDFOR

##### MAP RACE DICTIONARY TO SAMPLES #####
df_sample['Race'] = df_sample['Decode ID'].map(race_dict)

df_sample.drop("Decode ID",
               axis = 1,
               inplace = True
)


##### SAVE DICTIONARY AS A CSV #####
df_sample.to_csv(out_file,
                 sep = ',',
                 index = False
                 )

