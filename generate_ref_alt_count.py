import pandas as pd
import os
import pysam
import sys
import configparser

def vcf_filter_nucli_getero(my_id,case):
    print(my_id)
    file = open(os.path.join(processed_data, my_id + '/' + my_id + '.vcf'), "r")
    line = file.readline()
    n = 0
    while line.startswith("##"):
        n += 1
        line = file.readline()
    file.close()

    vcf_data = pd.read_csv(os.path.join(processed_data, my_id + '/' + my_id + '.vcf'), sep='\t', skiprows=n)
    vcf_filtrated = pd.DataFrame(columns=vcf_data.columns)
    if (case == 'rs'):
        for i in range(vcf_data.shape[0]):
            temp = vcf_data.iloc[i, 9]
            #print(temp)
            temp_list = temp.split(':')
            #print(temp_list)
            if ((temp_list[0] == '0/1') and
                    (vcf_data.iloc[i, 3] in nucliotides) and (vcf_data.iloc[i, 4] in nucliotides) and (vcf_data.iloc[i, 2] != '.')):
                vcf_filtrated = vcf_filtrated.append(vcf_data.iloc[i, :], ignore_index=False)
        vcf_filtrated.to_csv(os.path.join(processed_data, my_id + '/' + my_id + 'rs_nucli_getero_filtrated.vcf'),
                             sep='\t')
    else:
        for i in range(vcf_data.shape[0]):
            temp = vcf_data.iloc[i, 9]
            print(temp)
            temp_list = temp.split(':')
            print(temp_list)
            if ((temp_list[0] == '0/1') and
                    (vcf_data.iloc[i, 3] in nucliotides) and (vcf_data.iloc[i, 4] in nucliotides)):
                vcf_filtrated = vcf_filtrated.append(vcf_data.iloc[i, :], ignore_index=False)
        vcf_filtrated.to_csv(os.path.join(processed_data, my_id + '/' + my_id + '_nucli_getero_filtrated.vcf'), sep='\t')
    return (vcf_filtrated)

def concati(processing_list,case):
    all_df = pd.DataFrame(columns=['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT',
       '20'])
    for my_id in processing_list:
        tempdf = vcf_filter_nucli_getero(my_id,case)
        print(my_id)
        all_df = pd.concat([all_df, tempdf],axis=0, ignore_index=True)
        #dfs.append(tempdf)

    #all_df = pd.concat(dfs,axis=0,ignore_index=True)
    print('concatenated')
    return(all_df)
rs = sys.argv[1]
path = sys.argv[2]

#path = sys.argv[4]
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


print(table)
grp = (metadata.groupby(['ref', 'alt']).size()
       .reset_index(name='count'))
grp.to_csv(maindir+'logs/' + rs + '_ref_alt_count.tsv',sep='\t')
