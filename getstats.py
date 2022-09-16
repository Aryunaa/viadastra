import subprocess32 as subprocess
import os
import configparser
import sys as sys
import shutil
import pandas as pd




path = sys.argv[1]

config = configparser.ConfigParser()
config.read(path)
maindir = config["Directories"]["maindir"]

vcf_calls = config["Directories"]['final_data_out']
vcf_filtered = os.path.join(maindir,config["Directories"]['vcf_filtered'])
rssnps = os.path.join(maindir,config["Directories"]['rssnps'])


met = os.path.join(maindir,config["Files"]["metadata"])
metadata = pd.read_csv(met,sep='\t')
metadata['Extra1'] = metadata['Extra1'].apply(lambda x: x.lower())
outdir = config["Directories"]['temp_data_out']
tmp_path = os.path.join(outdir,'stats')
statfile = os.path.join(tmp_path, 'stats.tsv')
metadata['starting_snps']=0
metadata['filtrated_snps']=0
metadata['rssnps']=0
for i in range(metadata.shape[0]):

    strf = os.path.join(vcf_calls,i+'.vcf')
    file_rs = open(strf, "r")
    line_rs = file_rs.readline()
    n = 0
    while line_rs.startswith("##"):
        n += 1
        line_rs = file_rs.readline()
    file_rs.close()

    vcf = pd.read_csv(strf,
                      sep='\t', skiprows=n)


    strfb = os.path.join(vcf_filtered, i + '.snps.bed')
    bedf = pd.read_csv(strfb,sep='\t')

    rss = os.path.join(rssnps, i + '.snps.bed')
    bedrs = pd.read_csv(rss,sep='\t')

    metadata.iloc[i, 6] = vcf.shape[0]
    metadata.iloc[i, 7] = bedf.shape[0]
    metadata.iloc[i, 8] = bedrs.shape[0]

metadata.to_csv(statfile, index=False, sep ='\t')

