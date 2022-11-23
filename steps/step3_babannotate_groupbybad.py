import subprocess32 as subprocess
import os
import configparser
import sys as sys
import pathlib
import pandas as pd
from pybedtools import BedTool
import shlex

def babachi(configpath):
    config = configparser.ConfigParser()
    config.read(configpath)
    maindir = config["Directories"]["maindir"]
    print(maindir)
    filt_bed_rs = os.path.join(maindir, config["Directories"]["rssnps"])
    mainlogs = os.path.join(maindir, config["Directories"]["mainlogs"])
    all_log = os.path.join(mainlogs, 'babachilogs')
    met = os.path.join(maindir, config["Files"]["metadata"])
    metadata = pd.read_csv(met, sep='\t')
    tmp_path = os.path.join(maindir, config["Directories"]["babachi"])
    with open(all_log, "w") as log:
        log.write('babachi' + '\n')

    lst = metadata.BADgroup.unique()
    lst[0]
    header_list = ['#CHROM', 'POS1', 'POS2', 'ID', 'REF', 'ALT', 'REF_COUNT', 'ALT_COUNT', 'SAMPLEID']
    # list of bads
    for i in lst:
        process = subprocess.Popen(
            shlex.split('babachi '+os.path.join(tmp_path,'pulled'+i+'.tsv')+' -j 8 --visualize -e png -O '+tmp_path),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True)
        stderr, stdout = process.communicate()
        with open(all_log, "a") as log:
            log.write(stdout)
            log.write(stderr)
    print('babachi done')

def annotate_by_bads(path_tsv, path_badmap, out_path):
    config = configparser.ConfigParser()
    config.read(configpath)
    maindir = config["Directories"]["maindir"]
    print(maindir)
    filt_bed_rs = os.path.join(maindir, config["Directories"]["rssnps"])
    met = os.path.join(maindir, config["Files"]["metadata"])
    metadata = pd.read_csv(met, sep='\t')
    tmp_path = os.path.join(maindir, config["Directories"]["babachi"])

    try:
        header_list = ['#CHROM', 'POS1', 'POS2', 'ID', 'REF', 'ALT', 'ref', 'alt','sampleid']
        # vcf_data_atac = pd.read_csv(os.path.join(processed_data,'pulled_atacseq_tobabachi.vcf'),sep='\t', names = header_list)
        vcf_data = pd.read_csv(os.path.join(filt_bed_rs, path_tsv), sep='\t',
                               names=header_list)

        #vcf_data = vcf_data[['#CHROM', 'POS', 'POS2', 'ID', 'REF', 'ALT', 'ref', 'alt']]
        vcf_list = vcf_data.values.tolist()
        test = BedTool(vcf_list)
        print(vcf_data.shape[0])
        # bed_data = pd.read_csv(os.path.join(processed_data,'pulled_atacseq_tobabachi.bed'),sep='\t')
        bed_data_old = pd.read_csv(os.path.join(tmp_path, path_badmap), sep='\t')
        bed_data = bed_data_old[['#chr', 'start', 'end']]
        bed_data = bed_data[bed_data.end >= bed_data.start]
        bed_list = bed_data.values.tolist()
        print(bed_data.shape[0])
        annotations = BedTool(bed_list)
        i = test.intersect(annotations, wb=True)
        df = i.to_dataframe()
        # print(df.iloc[0,:])
        df.columns = ['#CHROM', 'POS1', 'POS2', 'ID', 'REF', 'ALT', 'REF_COUNTS', 'ALT_COUNTS', '#chr', 'start', 'end']
        df = pd.merge(df, bed_data_old[['#chr', 'start', 'end', 'BAD']], how='left', on=['#chr', 'start', 'end'])
        print(df.shape[0])
        # annotated_vcf = annotated_vcf[(annotated_vcf.ref >= threshold) & (annotated_vcf.alt >=threshold)]
        # annotated_vcf = annotated_vcf[((annotated_vcf.ref + annotated_vcf.alt) >= threshold)]
        #df['POS2'] = df['POS'] + 1
        df['POS']=df['POS2']
        annotated_vcf = df[['#CHROM','POS', 'ID', 'REF', 'ALT', 'REF_COUNTS', 'ALT_COUNTS', 'BAD']]
        print('shape= ' + str(annotated_vcf.shape[0]))
        annotated_vcf.to_csv(os.path.join(filt_bed_rs, out_path + '_annotated.tsv'), header=True,
                             index=False, sep='\t')
        print('annotate')
    except Exception:
        with open(os.path.join(filt_bed_rs, 'scorefiles/logs'), "a") as log:
            log.write('failed ' + path_tsv + '\t' + str(Exception) +'\n')
        print(Exception)


def group_by_bads():
    print('group')