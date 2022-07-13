import pandas as pd
import os
import numpy as np
import sys
import pybedtools
from pybedtools import BedTool
import subprocess32 as subprocess

processed_data = '/home/ariuna/rafaello/bedfiles/babachi'
bedfiles = '/home/ariuna/rafaello/bedfiles'
metapath = '/home/ariuna/rafaello/viadastra/additional/metadata_yes_no.tsv'

######################################################
def annotate_by_bad(path_tsv, path_badmap, out_path):
    #path_bed - badmap path
    #path_vcf - before babachi tsv


    try:
        header_list = ['#CHROM', 'POS', 'POS2', 'ID', 'REF', 'ALT', 'ref', 'alt']
        # vcf_data_atac = pd.read_csv(os.path.join(processed_data,'pulled_atacseq_tobabachi.vcf'),sep='\t', names = header_list)
        vcf_data = pd.read_csv(os.path.join(bedfiles, path_tsv), sep='\t',
                               names=header_list)

        #vcf_data = vcf_data[['#CHROM', 'POS', 'POS2', 'ID', 'REF', 'ALT', 'ref', 'alt']]
        vcf_list = vcf_data.values.tolist()
        test = BedTool(vcf_list)
        print(vcf_data.shape[0])
        # bed_data = pd.read_csv(os.path.join(processed_data,'pulled_atacseq_tobabachi.bed'),sep='\t')
        bed_data_old = pd.read_csv(os.path.join(processed_data, path_badmap), sep='\t')
        bed_data = bed_data_old[['#chr', 'start', 'end']]
        bed_data = bed_data[bed_data.end >= bed_data.start]
        bed_list = bed_data.values.tolist()
        print(bed_data.shape[0])
        annotations = BedTool(bed_list)
        i = test.intersect(annotations, wb=True)
        df = i.to_dataframe()
        # print(df.iloc[0,:])
        df.columns = ['#CHROM', 'POS', 'POS2', 'ID', 'REF', 'ALT', 'REF_COUNTS', 'ALT_COUNTS', '#chr', 'start', 'end']
        df = pd.merge(df, bed_data_old[['#chr', 'start', 'end', 'BAD']], how='left', on=['#chr', 'start', 'end'])
        print(df.shape[0])
        # annotated_vcf = annotated_vcf[(annotated_vcf.ref >= threshold) & (annotated_vcf.alt >=threshold)]
        # annotated_vcf = annotated_vcf[((annotated_vcf.ref + annotated_vcf.alt) >= threshold)]
        #df['POS2'] = df['POS'] + 1
        annotated_vcf = df[['#CHROM', 'POS2', 'ID', 'REF', 'ALT', 'REF_COUNTS', 'ALT_COUNTS', 'BAD']]
        print('shape= ' + str(annotated_vcf.shape[0]))
        annotated_vcf.to_csv(os.path.join(bedfiles, out_path + '_annotated.tsv'), header=True,
                             index=False, sep='\t')
    except Exception:
        with open(os.path.join(bedfiles, 'scorefiles/logs'), "a") as log:
            log.write('failed ' + path_tsv + '\t' + str(Exception) +'\n')
        print(Exception)


if(not os.path.isdir(os.path.join(bedfiles,'scorefiles'))):
    os.mkdir(os.path.join(bedfiles,'scorefiles'))
print('start')
listi = os.listdir(bedfiles)
bedlist = [k for k in listi if '.bed' in k]
lis=[k.split('.')[0] for k in bedlist]


metadata = pd.read_csv(metapath, sep = '\t')
metadata['ChipTFrepair'] = metadata['ChipTFrepair'].replace(np.nan, 0)
with open(os.path.join(bedfiles,'scorefiles/logs'), "w") as log:
    log.write('start' + '\n')

for i in bedlist:
    j = i.split('.')[0]
    row = metadata[metadata.ID == j]
    a = row.values.tolist()

    print(row)
    print(i)
    print(j)
    if (a[0][4] == 'ChIPseq'):
        annotate_by_bad(i, 'pulled_chipseq.badmap.bed', 'scorefiles/chipmap_' + j)
    elif (a[0][4] == 'ATACseq'):
        annotate_by_bad(i, 'pulled_atacseq.badmap.bed', 'scorefiles/atacmap_' + j)
    else:
        annotate_by_bad(i, 'pulled_chipseq.badmap.bed', 'scorefiles/chipmap_' + j)
        annotate_by_bad(i, 'pulled_atacseq.badmap.bed', 'scorefiles/atacmap_' + j)

    #by condition
    if(a[0][8]=='yes'):
        print(a[0][8])
        print(row)
        print(i)
        print(j)
        if (a[0][4]=='ChIPseq'):
            annotate_by_bad(i, 'pulled_chipseq.badmap.bed', 'scorefiles/chipmap_yes'+ j)
        elif (a[0][4] == 'ATACseq'):
            annotate_by_bad(i, 'pulled_atacseq.badmap.bed', 'scorefiles/atacmap_yes' + j)
        else:
            annotate_by_bad(i, 'pulled_chipseq.badmap.bed', 'scorefiles/chipmap_yes' + j)
            annotate_by_bad(i, 'pulled_atacseq.badmap.bed', 'scorefiles/atacmap_yes' + j)




print('end')