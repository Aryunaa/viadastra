import subprocess32 as subprocess
import os
import configparser
import sys as sys
import vcf
import pandas as pd

print('start process')

path = sys.argv[1]
config = configparser.ConfigParser()
config.read(path)
maindir = config["Directories"]["maindir"]
processed_data = os.path.join(maindir,config["Directories"]["data_out"])
met = os.path.join(maindir,config["Files"]["metadata"])
processing_list_path = os.path.join(maindir,config["Files"]["processing_list"])
indir = os.path.join(maindir,config["Directories"]["data_in"])


list_file = open(processing_list_path, "r")
listi = list_file.readlines()
processing_list = []
for i in listi:
    i = i.replace('\n','')
    processing_list.append(i)

metadata = pd.read_csv(met,sep='\t')
#badgroup1 = sys.argv[1]
#badgroup2 = sys.argv[2]
grp = (metadata.groupby(['BADgroup']).size()
       .reset_index(name='count'))
grp = grp[grp.BADgroup!='.']
bad_list = []
for i in grp.iloc[:,0]:
    bad_list.append(i)
print(bad_list)

for i in bad_list:
    pull = metadata[metadata['BADgroup'] == i]
    # pull_chip = metadata[metadata['BADgroup']==badgroup1]
    intersect = list(filter(lambda x: x in list(pull['ID']), processing_list))
    paths_rs = []
    for my_id in intersect:
        paths_rs.append(os.path.join(processed_data, my_id + '/' + my_id + '_rs_nucli_getero_filtrated.vcf'))

    header_list = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT',
                   '20']
    pulled_rs = pd.concat([pd.read_csv(f, sep='\t', comment='#', names=header_list) for f in paths_rs])
    vcf_rreader = vcf.Reader(open(paths_rs[0], 'r'))
    vcf_writer = vcf.Writer(open(os.path.join(processed_data, 'pulled_'+i+'_rs_getero_filtrated.vcf'), 'w'),
                            vcf_rreader)
    vcf_writer.close()
    pulled_rs.to_csv(os.path.join(processed_data, 'pulled_'+i+'_rs_getero_filtrated.vcf'), mode='a',
                           header=False, index=False, sep='\t')
    print(i+' pulled')
    process = subprocess.run(['bcftools', 'sort',
                              os.path.join(processed_data, 'pulled_'+i+'_rs_getero_filtrated.vcf'), '--output-file',
                              os.path.join(processed_data, 'pulled_sorted_'+i+'_rs_getero_filtrated.vcf')],
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             universal_newlines=True)
    print(i+' sorted')


print('done')