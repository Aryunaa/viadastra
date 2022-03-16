
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



vcf_reader = vcf.Reader(open(os.path.join(processed_data,'pulled_sorted_chipseq_rs_getero_filtrated.vcf'), 'r'))
vcf_writer = open(os.path.join(processed_data,'pulled_chipseq_tobabachi.vcf'), "w")
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


print('chipseq done')

vcf_reader_atac = vcf.Reader(open(os.path.join(processed_data,'pulled_sorted_atacseq_rs_getero_filtrated.vcf'), 'r'))
vcf_writer_atac = open(os.path.join(processed_data,'pulled_atacseq_tobabachi.vcf'), "w")
for record in vcf_reader_atac:
    vcf_writer_atac.write(record.CHROM)
    vcf_writer_atac.write('\t')
    vcf_writer_atac.write(str(record.POS))
    vcf_writer_atac.write('\t')
    vcf_writer_atac.write(record.ID)
    vcf_writer_atac.write('\t')
    vcf_writer_atac.write(record.REF)
    vcf_writer_atac.write('\t')
    vcf_writer_atac.write(str(record.ALT[0]))
    vcf_writer_atac.write('\t')
    vcf_writer_atac.write(str(record.genotype('20')['AD'][0]))
    vcf_writer_atac.write('\t')
    vcf_writer_atac.write(str(record.genotype('20')['AD'][1]))
    vcf_writer_atac.write('\n')
vcf_writer_atac.close()
print('atacseq done')