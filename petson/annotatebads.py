import pandas as pd
import os
import sys
import pybedtools
from pybedtools import BedTool
import subprocess32 as subprocess

processed_data = '/home/ariuna/rafaello/bedfiles/babachi'
bedfiles = '/home/ariuna/rafaello/bedfiles'
metapath = '/home/ariuna/rafaello/viadastra/additional/metadata_v5.tsv'

######################################################
def group_by_bad(path_tsv, path_badmap, out_path):
    #path_bed - badmap path
    #path_vcf - before babachi tsv
    try:
        header_list = ['#CHROM', 'POS','POS2', 'ID', 'REF', 'ALT', 'ref', 'alt']
        #vcf_data_atac = pd.read_csv(os.path.join(processed_data,'pulled_atacseq_tobabachi.vcf'),sep='\t', names = header_list)
        vcf_data = pd.read_csv(os.path.join(bedfiles, path_tsv), sep='\t',
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
        #print(df.iloc[0,:])
        df.columns = ['#CHROM', 'POS','POS2','ID', 'REF', 'ALT', 'REF_COUNTS', 'ALT_COUNTS','chr','start','end','BAD']

        #annotated_vcf = annotated_vcf[(annotated_vcf.ref >= threshold) & (annotated_vcf.alt >=threshold)]
        #annotated_vcf = annotated_vcf[((annotated_vcf.ref + annotated_vcf.alt) >= threshold)]
        df['POS2'] = df['POS']+1
        annotated_vcf = df[['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'REF_COUNTS', 'ALT_COUNTS', 'BAD']]
        annotated_vcf.to_csv(os.path.join(bedfiles, out_path+ '_annotated.tsv'),header=True,
                            index=False, sep='\t')
    except Exception:
        with open('scorefiles/logs', "a") as log:
            log.write('failed ' + i + '\t' + Exception +'\n')
        print(Exception)


if(not os.path.isdir(os.path.join(bedfiles,'scorefiles'))):
    os.mkdir(os.path.join(bedfiles,'scorefiles'))
print('start')
listi = os.listdir(bedfiles)
bedlist = [k for k in listi if '.bed' in k]
lis=[k.split('.')[0] for k in bedlist]


metadata = pd.read_csv(metapath, sep = '\t')

with open('scorefiles/logs', "w") as log:
    log.write('start' + '\n')

for i in bedlist:
    j = i.split('.')[0]
    row = metadata[metadata.ID == j]
    a = row.values.tolist()

    if (a[0][4]=='ChIPseq'):
        group_by_bad(i, 'pulled_chipseq.badmap.bed', 'scorefiles/chipmap_'+ j)
    elif (a[0][4] == 'ATACseq'):
        group_by_bad(i, 'pulled_atacseq.badmap.bed', 'scorefiles/atacmap_' + j)
    else:
        group_by_bad(i, 'pulled_chipseq.badmap.bed', 'scorefiles/chipmap_' + j)
        group_by_bad(i, 'pulled_atacseq.badmap.bed', 'scorefiles/atacmap_' + j)




print('end')