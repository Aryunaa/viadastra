import pandas as pd
import os
import sys
import configparser
import vcf
import subprocess32 as subprocess
#path = '/media/ElissarDisk/ADASTRA/parameters/CONFIG.cfg'
path = sys.argv[1]
config = configparser.ConfigParser()
config.read(path)
maindir = config["Directories"]["maindir"]
processed_data = os.path.join(maindir,config["Directories"]["data_out"])
met = os.path.join(maindir,config["Files"]["metadata"])
processing_list_path = os.path.join(maindir,config["Files"]["processing_list"])
indir = os.path.join(maindir,config["Directories"]["data_in"])
babachi = os.path.join(maindir,config["Directories"]["babachi"])

metadata = pd.read_csv(met,sep='\t')
grp = (metadata.groupby(['BADgroup']).size()
       .reset_index(name='count'))
grp = grp[grp.BADgroup!='.']
bad_list = []
for i in grp.iloc[:,0]:
    bad_list.append(i)
print(bad_list)

for i in bad_list:
    if (not os.path.exists(os.path.join(processed_data, 'pulled_'+i+'_tobabachi.tsv'))):
        vcf_reader = vcf.Reader(open(os.path.join(processed_data, 'pulled_sorted_'+i+'_rs_getero_filtrated.vcf'), 'r'))
        vcf_writer = open(os.path.join(processed_data, 'pulled_'+i+'_tobabachi.tsv'), "w")
        for record in vcf_reader:
            vcf_writer.write(record.CHROM)
            vcf_writer.write('\t')
            vcf_writer.write(str(record.POS))
            vcf_writer.write('\t')
            vcf_writer.write(record.ID)
            vcf_writer.write('\t')
            vcf_writer.write(record.REF)
            vcf_writer.write('\t')
            vcf_writer.write(str(record.ALT[0]))
            vcf_writer.write('\t')
            vcf_writer.write(str(record.genotype('20')['AD'][0]))
            vcf_writer.write('\t')
            vcf_writer.write(str(record.genotype('20')['AD'][1]))
            vcf_writer.write('\n')
        vcf_writer.close()

    print(i+' done')

for i in bad_list:
    if (not os.path.isdir(babachi)):
        os.mkdir(babachi)
    process = subprocess.run(['babachi', os.path.join(processed_data, 'pulled_'+i+'_tobabachi.tsv'),
                              '-O', babachi,
                              '--visualize'
                              ],
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             universal_newlines=True)
    #babachi /media/ElissarDisk/ADASTRA/processed_data/pulled_chipseq_tobabachi.vcf --visualize