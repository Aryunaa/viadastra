import pandas as pd
import os
import pysam
import sys
import pybedtools
from pybedtools import BedTool
import subprocess32 as subprocess

processed_data = '/home/ariuna/rafaello/bedfiles/babachi'

######################################################
def group_by_bad(path_tsv, path_badmap, out_path):
    #path_bed - badmap path
    #path_vcf - before babachi tsv
    header_list = ['#CHROM', 'POS','POS2', 'ID', 'REF', 'ALT', 'ref', 'alt']
    #vcf_data_atac = pd.read_csv(os.path.join(processed_data,'pulled_atacseq_tobabachi.vcf'),sep='\t', names = header_list)
    vcf_data = pd.read_csv(os.path.join(processed_data, path_tsv), sep='\t',
                                names=header_list)
    vcf_data['POS2'] = vcf_data['POS']
    vcf_data = vcf_data[['#CHROM', 'POS', 'POS2','ID', 'REF', 'ALT', 'ref', 'alt']]
    vcf_list = vcf_data.values.tolist()
    test = BedTool(vcf_list)

    #bed_data = pd.read_csv(os.path.join(processed_data,'pulled_atacseq_tobabachi.bed'),sep='\t')
    bed_data = pd.read_csv(os.path.join(processed_data, path_badmap), sep='\t')
    bed_data = bed_data[['#chr','start','end','BAD']]
    bed_data = bed_data[bed_data.end>=bed_data.start]
    bed_list = bed_data.values.tolist()

    annotations = BedTool(bed_list)
    i = test.intersect(annotations, wb=True)
    df = i.to_dataframe()
    df.columns = ['#CHROM', 'POS','POS2','ID', 'REF', 'ALT', 'ref', 'alt','chr','start','end','BAD']

    #annotated_vcf = annotated_vcf[(annotated_vcf.ref >= threshold) & (annotated_vcf.alt >=threshold)]
    #annotated_vcf = annotated_vcf[((annotated_vcf.ref + annotated_vcf.alt) >= threshold)]
    if (os.path.exists(os.path.join(processed_data, out_path)) == False):
        os.mkdir(os.path.join(processed_data, out_path))
    df['POS2'] = df['POS']+1
    annotated_vcf = df[['#CHROM', 'POS','POS2', 'ID', 'REF', 'ALT', 'ref', 'alt', 'BAD']]
    annotated_vcf.to_csv(os.path.join(processed_data, out_path+ '/annotated.tsv'),header=True,
                        index=False, sep='\t')

    df_list = []
    for i in range(6):
        temp_bad = annotated_vcf[annotated_vcf.BAD == (i + 1)]
        temp_bad = temp_bad.drop(['BAD'],axis=1)
        df_list.append(temp_bad)
        if (os.path.exists(os.path.join(processed_data, out_path + '/BAD' + str(i + 1))) == False):
            os.mkdir(os.path.join(processed_data, out_path + '/BAD' + str(i + 1)))
        temp_bad.to_csv(os.path.join(processed_data, out_path + '/BAD' + str(i + 1) + '/table.tsv'), header=True,
                        index=False, sep='\t')
        temp_grp_bad = (temp_bad.groupby(['ref','alt']).size().reset_index(name='counts'))
        temp_grp_bad.to_csv(os.path.join(processed_data, out_path + '/BAD' + str(i + 1) + '/stats.tsv'), header=False,
                        index=False, sep='\t')

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
group_by_bad('pulled_atacseq.tsv','pulled_atacseq.badmap.bed', 'grouped_bads/atacseq')
print('atac_end')
print('chip start')
group_by_bad('pulled_chipseq.tsv','pulled_chipseq.badmap.bed', 'grouped_bads/chipseq')
print('chip_end')

######################################################

'''

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
'''
