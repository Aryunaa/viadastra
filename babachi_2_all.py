import subprocess32 as subprocess
import pandas as pd
import os
import pysam
import sys
import configparser
import vcf

#path = '/media/ElissarDisk/ADASTRA/parameters/CONFIG.cfg'
path = sys.argv[1]
config = configparser.ConfigParser()
config.read(path)
maindir = config["Directories"]["maindir"]
processed_data = os.path.join(maindir,config["Directories"]["data_out"])
met = os.path.join(maindir,config["Files"]["metadata"])
processing_list_path = os.path.join(maindir,config["Files"]["processing_list"])
indir = os.path.join(maindir,config["Directories"]["data_in"])

process = subprocess.run(['bcftools', 'sort',
                            os.path.join(processed_data,'pulled_all_rs_getero_filtrated.vcf'), '--output-file',
                            os.path.join(processed_data,'pulled_sorted_all_rs_getero_filtrated.vcf')],
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE,
                           universal_newlines=True
                           )
print('sorted')



vcf_reader = vcf.Reader(open(os.path.join(processed_data,'pulled_sorted_all_rs_getero_filtrated.vcf'), 'r'))
vcf_writer = open(os.path.join(processed_data,'pulled_all_tobabachi.vcf'), "w")
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


print('done')