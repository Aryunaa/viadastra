import pandas as pd
import os
import pysam
import sys
import configparser
import pybedtools
from pybedtools import BedTool
import subprocess32 as subprocess
#path = '/media/ElissarDisk/ADASTRA/parameters/CONFIG.cfg'
threshold = int(sys.argv[1])
path = sys.argv[2]
config = configparser.ConfigParser()
config.read(path)
maindir = config["Directories"]["maindir"]
processed_data = os.path.join(maindir,config["Directories"]["data_out"])
met = os.path.join(maindir,config["Files"]["metadata"])
processing_list_path = os.path.join(maindir,config["Files"]["processing_list"])
indir = os.path.join(maindir,config["Directories"]["data_in"])

######################################################
def group_by_bad(path_vcf, path_bed, out_path,threshold):
    header_list = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'ref', 'alt']
    #vcf_data_atac = pd.read_csv(os.path.join(processed_data,'pulled_atacseq_tobabachi.vcf'),sep='\t', names = header_list)
    vcf_data_atac = pd.read_csv(os.path.join(processed_data, path_vcf), sep='\t',
                                names=header_list)
    vcf_data_atac['POS2'] = vcf_data_atac['POS']
    vcf_data_atac = vcf_data_atac[['#CHROM', 'POS', 'POS2','ID', 'REF', 'ALT', 'ref', 'alt']]
    vcf_list = vcf_data_atac.values.tolist()
    test = BedTool(vcf_list)

    #bed_data = pd.read_csv(os.path.join(processed_data,'pulled_atacseq_tobabachi.bed'),sep='\t')
    bed_data = pd.read_csv(os.path.join(processed_data, path_bed), sep='\t')
    bed_data = bed_data[['#chr','start','end','BAD']]
    bed_list = bed_data.values.tolist()
    annotations = BedTool(bed_list)
    i = test.intersect(annotations, wb=True)
    df = i.to_dataframe()
    df.columns = ['#CHROM', 'POS','POS2','ID', 'REF', 'ALT', 'ref', 'alt','chr','start','end','bad']
    annotated_vcf = df[['#CHROM', 'POS','ID', 'REF', 'ALT', 'ref', 'alt','bad']]
    #annotated_vcf = annotated_vcf[(annotated_vcf.ref >= threshold) & (annotated_vcf.alt >=threshold)]
    annotated_vcf = annotated_vcf[((annotated_vcf.ref + annotated_vcf.alt) >= threshold)]
    vcf_bad1 = annotated_vcf[annotated_vcf.bad == 1]
    vcf_bad1 = vcf_bad1.drop(['bad'], axis=1)
    vcf_bad2 = annotated_vcf[annotated_vcf.bad == 2]
    vcf_bad2 = vcf_bad2.drop(['bad'], axis=1)
    vcf_bad3 = annotated_vcf[annotated_vcf.bad == 3]
    vcf_bad3 = vcf_bad3.drop(['bad'], axis=1)
    vcf_bad4 = annotated_vcf[annotated_vcf.bad == 4]
    vcf_bad4 = vcf_bad4.drop(['bad'], axis=1)
    vcf_bad5 = annotated_vcf[annotated_vcf.bad == 5]
    vcf_bad5 = vcf_bad5.drop(['bad'], axis=1)
    vcf_bad6 = annotated_vcf[annotated_vcf.bad == 6]
    vcf_bad6 = vcf_bad6.drop(['bad'], axis=1)

    vcf_bad1.to_csv(os.path.join(processed_data, out_path + '_table_BAD1.tsv'), header=True, index=False, sep='\t')
    vcf_bad2.to_csv(os.path.join(processed_data, out_path + '_table_BAD2.tsv'), header=True, index=False, sep='\t')
    vcf_bad3.to_csv(os.path.join(processed_data, out_path + '_table_BAD3.tsv'), header=True, index=False, sep='\t')
    vcf_bad4.to_csv(os.path.join(processed_data, out_path + '_table_BAD4.tsv'), header=True, index=False, sep='\t')
    vcf_bad5.to_csv(os.path.join(processed_data, out_path + '_table_BAD5.tsv'), header=True, index=False, sep='\t')
    vcf_bad6.to_csv(os.path.join(processed_data, out_path + '_table_BAD6.tsv'), header=True, index=False, sep='\t')

    grp_bad1 = (vcf_bad1.groupby(['ref','alt']).size()
           .reset_index(name='counts'))
    grp_bad2 = (vcf_bad2.groupby(['ref','alt']).size()
           .reset_index(name='counts'))
    grp_bad3 = (vcf_bad3.groupby(['ref','alt']).size()
           .reset_index(name='counts'))
    grp_bad4 = (vcf_bad4.groupby(['ref','alt']).size()
           .reset_index(name='counts'))
    grp_bad5 = (vcf_bad5.groupby(['ref','alt']).size()
           .reset_index(name='counts'))
    grp_bad6 = (vcf_bad6.groupby(['ref','alt']).size()
           .reset_index(name='counts'))

    grp_bad1.to_csv(os.path.join(processed_data,out_path + '_BAD1.tsv'),header=True,index=False,sep='\t')
    grp_bad2.to_csv(os.path.join(processed_data, out_path + '_BAD2.tsv'),header=True,index=False,sep='\t')
    grp_bad3.to_csv(os.path.join(processed_data,out_path + '_BAD3.tsv'),header=True,index=False,sep='\t')
    grp_bad4.to_csv(os.path.join(processed_data, out_path + '_BAD4.tsv'),header=True,index=False,sep='\t')
    grp_bad5.to_csv(os.path.join(processed_data, out_path + '_BAD5.tsv'),header=True,index=False,sep='\t')
    grp_bad6.to_csv(os.path.join(processed_data,out_path + '_BAD6.tsv'),header=True,index=False,sep='\t')

def negbinfit(out_path,badn,badt, threshold):
    #badt = '_BAD1.tsv'
    print(os.path.join(processed_data,out_path + badt))
    subprocess.run(['negbin_fit', os.path.join(processed_data,out_path + badt),
                                '--bad',str(badn), '--allele-reads-tr', str(threshold),
                                '-O',os.path.join(processed_data,'grouped_bads/fitted_negbin'), '--visualize', '-r'
                              ],
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               universal_newlines=True
                               )
    subprocess.run(['negbin_fit', os.path.join(processed_data,out_path + badt),
                                '--bad',str(badn), '--allele-reads-tr', str(threshold),
                                '-O',os.path.join(processed_data,'grouped_bads/fitted_negbin_lin'), '--visualize', '-r',
                              ],
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               universal_newlines=True
                               )
    subprocess.run(['negbin_fit', os.path.join(processed_data,out_path + badt),
                                '--bad',str(badn), '--allele-reads-tr', str(threshold),
                                '-O',os.path.join(processed_data,'grouped_bads/fitted_negbin_lin'), '--visualize', '-r','-l'
                              ],
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               universal_newlines=True
                               )


if(not os.path.isdir(os.path.join(processed_data,'grouped_bads'))):
    os.mkdir(os.path.join(processed_data,'grouped_bads'))
print('atac_start')
group_by_bad('pulled_atacseq_tobabachi.vcf','pulled_atacseq_tobabachi.bed', 'grouped_bads/atacseq_altref_counts',threshold)
print('atac_end')
print('chip start')
group_by_bad('pulled_chipseq_tobabachi.vcf','pulled_chipseq_tobabachi.bed', 'grouped_bads/chipseq_altref_counts',threshold)
print('chip_end')

######################################################
if(not os.path.isdir(os.path.join(processed_data,'grouped_bads/fitted_negbin'))):
    os.mkdir(os.path.join(processed_data,'grouped_bads/fitted_negbin'))
if(not os.path.isdir(os.path.join(processed_data,'grouped_bads/fitted_negbin_lin'))):
    os.mkdir(os.path.join(processed_data,'grouped_bads/fitted_negbin_lin'))

for i in range(6):
    print(str(i+1)+'_BAD.tsv')
    negbinfit('grouped_bads/atacseq_altref_counts',i+1,'_BAD'+str(i+1)+'.tsv',threshold)
    negbinfit('grouped_bads/chipseq_altref_counts', i + 1, '_BAD' + str(i + 1) + '.tsv', threshold)

#vcf_reader_atac = vcf.Reader(open(os.path.join(processed_data,'pulled_atacseq_tobabachi.vcf'), 'r'))
#vcf_writer_atac = open(os.path.join(processed_data,'pulled_atacseq_bad_annotated.vcf'), "w")

