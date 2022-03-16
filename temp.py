import vcf
import sys
import os
import configparser
import pandas as pd


rs = sys.argv[1]
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

for my_id in processing_list:
    print(my_id)
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

