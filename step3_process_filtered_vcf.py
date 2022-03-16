import pandas as pd
import os
import pysam
import sys
import configparser


vcf_str = sys.argv[1]
path = sys.argv[2]

config = configparser.ConfigParser()
config.read(path)
maindir = config["Directories"]["maindir"]
processed_data = os.path.join(maindir,config["Directories"]["data_out"])
met = os.path.join(maindir,config["Files"]["metadata"])
indir = os.path.join(maindir,config["Directories"]["data_in"])

processing_list_path = os.path.join(maindir,config["Files"]["processing_list"])
list_file = open(processing_list_path, "r")
listi = list_file.readlines()
processing_list = []
for i in listi:
    i = i.replace('\n','')
    processing_list.append(i)

nucliotides = {'A', 'T', 'G', 'C'}

table = pd.DataFrame(columns=['ref','alt'])
all_vcfs = concati(processing_list,rs)
print('concati ended')
print(all_vcfs)
for i in range(all_vcfs.shape[0]):
    temp = vcf_data.iloc[i, 9]
    print(temp)
    temp_list = temp.split(':')
    ad = temp_list[1].split(',')
    ref_ad = int(ad[0])
    alt_ad = int(ad[1])
    to_append = [ref_ad,alt_ad]
    a_series = pd.Series(to_append, index=table.columns)
    table = table.append(a_series, ignore_index=True)

'''
print(table)
grp = (metadata.groupby(['ref', 'alt']).size()
       .reset_index(name='count'))
grp.to_csv(maindir+'logs/' + rs + '_ref_alt_count.tsv',sep='\t')
'''