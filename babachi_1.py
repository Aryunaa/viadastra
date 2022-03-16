import subprocess32 as subprocess
import os
import configparser
import sys as sys
import vcf
import pandas as pd
'''
def process_vcf(my_id):
    print('start process')

    inp = 'data/BAM00022_tobabachi.vcf'
    outp = 'data/BAM00022.bed'
    process = subprocess.Popen(['babachi', inp,
                                '--output', outp,'--visualize'],
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               universal_newlines=True
                               )
    stderr, stdout = process.communicate()
    loggi(tmp_log, tmp_err, stdout, stderr, 'a')
    print('done')
'''
#my_id = sys.argv[1]
def tobabachi(my_id):
    vcf_reader = vcf.Reader(open('/media/ElissarDisk/ADASTRA/processed_data/'+my_id+'/'+my_id+ '_rs_nucli_getero_filtrated.vcf', 'r'))
    vcf_writer = open('/media/ElissarDisk/ADASTRA/processed_data/'+my_id+'/'+my_id+ '_tobabachi.vcf',"w")
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

#!babachi data/BAM00022_tobabachi.vcf --output data/BAM00022.bed --visualize
#babachi /media/ElissarDisk/ADASTRA/processed_data/BAM00022/BAM00022_tobabachi.vcf --output /media/ElissarDisk/ADASTRA/processed_data/BAM00022/BAM00022.bed --visualize
print('start process')
'''
inp = '/media/ElissarDisk/ADASTRA/processed_data/BAM00028/BAM00028_tobabachi.vcf'
outp = '/media/ElissarDisk/ADASTRA/processed_data/BAM00028/BAM00028.bed'

process = subprocess.Popen(['babachi', inp,
                            '--output', outp, '--visualize'],
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE,
                           universal_newlines=True
                           )
'''

path = sys.argv[1]
#path = 'CONFIG.cfg'
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
pull_chip = metadata[metadata['BADgroup']=='chipseq']

intersect = list(filter(lambda x:x in list(pull_chip['ID']),processing_list))
paths_rs = []
for my_id in intersect:
    paths_rs.append(os.path.join(processed_data, my_id + '/' + my_id + '_rs_nucli_getero_filtrated.vcf'))

header_list = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT',
       '20']
pulled_chips_rs = pd.concat([pd.read_csv(f,sep='\t',comment='#',names=header_list) for f in paths_rs])
vcf_rreader = vcf.Reader(open(paths_rs[0], 'r'))
vcf_writer = vcf.Writer(open(os.path.join(processed_data,'pulled_chipseq_rs_getero_filtrated.vcf'), 'w'), vcf_rreader)
vcf_writer.close()
pulled_chips_rs.to_csv(os.path.join(processed_data,'pulled_chipseq_rs_getero_filtrated.vcf'),mode='a', header=False,index=False,sep='\t')
print('chipseq pulled')
process = subprocess.run(['bcftools', 'sort',
                            os.path.join(processed_data,'pulled_chipseq_rs_getero_filtrated.vcf'), '--output-file',
                            os.path.join(processed_data,'pulled_sorted_chipseq_rs_getero_filtrated.vcf')],
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE,
                           universal_newlines=True
                           )
print('chipseq sorted')


pull_atac = metadata[metadata['BADgroup']=='atacseq']
intersect = list(filter(lambda x:x in list(pull_atac['ID']),processing_list))
paths_rs = []
for my_id in intersect:
    paths_rs.append(os.path.join(processed_data, my_id + '/' + my_id + '_rs_nucli_getero_filtrated.vcf'))

header_list = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT',
       '20']
pulled_atac_rs = pd.concat([pd.read_csv(f,sep='\t',comment='#',names=header_list) for f in paths_rs])
vcf_rreader = vcf.Reader(open(paths_rs[0], 'r'))
vcf_writer = vcf.Writer(open(os.path.join(processed_data,'pulled_atacseq_rs_getero_filtrated.vcf'), 'w'), vcf_rreader)
vcf_writer.close()
pulled_atac_rs.to_csv(os.path.join(processed_data,'pulled_atacseq_rs_getero_filtrated.vcf'),mode='a', header=False,index=False,sep='\t')
print('atacseq pulled')
process = subprocess.run(['bcftools', 'sort',
                            os.path.join(processed_data,'pulled_atacseq_rs_getero_filtrated.vcf'), '--output-file',
                            os.path.join(processed_data,'pulled_sorted_atacseq_rs_getero_filtrated.vcf')],
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE,
                           universal_newlines=True
                           )
print('atacseq sorted')





'''
print('start to make available for babachi')
vcf_read = vcf.Reader(open(os.path.join(processed_data,'pulled_sorted_chipseq_rs_getero_filtrated.vcf'), 'r'))
vcf_writer = open(os.path.join(processed_data,'pulled_chipseq_tobabachi.vcf'), "w")
print(vcf_read.infos)
for record in vcf_read:
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
'''

'''
print('start to make available for babachi')
subprocess.Popen(["python", os.path.join(maindir,'scripts/babachi_2.py'),path])
'''
print('done')