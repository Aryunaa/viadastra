import subprocess32 as subprocess
import os
import configparser
import sys as sys
import shutil
import pandas as pd



path = sys.argv[1]

config = configparser.ConfigParser()
config.read(path)
maindir = config["Directories"]["maindir"]

vcf_calls = config["Directories"]['final_data_out']
vcf_filtered = os.path.join(maindir,config["Directories"]['vcf_filtered'])
rssnps = os.path.join(maindir,config["Directories"]['rssnps'])


met = os.path.join(maindir,config["Files"]["metadata"])
metadata = pd.read_csv(met,sep='\t')

##############
processed_data = rssnps
testdir = config["Directories"]['temp_data_out']
tmp_path = os.path.join(testdir,'babachi')
if (os.path.exists(tmp_path)):
    pass
else:
    os.mkdir(tmp_path)

lst = metadata.BADgroup.unique()
lst[0]

for i in lst:
    pulltmp = metadata[metadata['BADgroup']==i]

    intersect = list(pulltmp['ID'])
    paths_rs = []
    for my_id in intersect:
        if(os.path.exists(os.path.join(processed_data, my_id+'.snps.bed'))):
            paths_rs.append(os.path.join(processed_data, my_id+'.snps.bed'))
    print(paths_rs)

    #читаем чипсеки, смотрим распределение
    #chip_list = []
    header_list = ['#CHROM', 'POS1','POS2', 'ID', 'REF', 'ALT', 'REF_COUNT', 'ALT_COUNT','SAMPLEID']
    with open(os.path.join(tmp_path,'bedshapes'), "w") as log:
        for my_id in list(metadata['ID']):
            if (os.path.exists(os.path.join(processed_data, my_id + '.snps.bed'))):
                tempdf = pd.read_csv(os.path.join(processed_data, my_id + '.snps.bed'), sep='\t', names=header_list)
                log.write(my_id +'\t'+ str(int(tempdf.shape[0])) +'\n')


    pulledtmps_rs = pd.concat([pd.read_csv(f,sep='\t',names=header_list) for f in paths_rs])
    pulledtmps_rs.to_csv(os.path.join(tmp_path,'ppulled'+i+'.tsv'),mode='w', header=False,index=False,sep='\t')
    print(i+ ' pulled')
    process = subprocess.run(['bedtools', 'sort','-i',
                                os.path.join(tmp_path,'ppulled'+i+'.tsv')],
                               stdout=open(os.path.join(tmp_path,'pulled'+i+'.tsv'), "w"),
                               stderr=subprocess.PIPE,
                               universal_newlines=True
                               )
    print(i+ ' sorted')
    print(pulledtmps_rs.shape[0])


#_____________________________________________________________________________
