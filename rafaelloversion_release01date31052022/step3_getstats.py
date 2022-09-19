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
metadata['bamsize MB']=0
myids=list(metadata['ID'])

#get shapes of dataframes
for i in myids:

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
    bedf['sample_id']==i
    bedf.to_csv(strf, sep ='\t',index=False)
    header_list = ['#chr','start','end','ID','ref','alt','ref_counts','alt_counts','sample_id']
    if (os.stat(os.path.join(rssnps, i + '.snps.bed')).st_size != 0):
        rss = os.path.join(rssnps, i + '.snps.bed')
        bedrs = pd.read_csv(rss,sep='\t', names=header_list)
        bedrs['sample_id']==i
        bedrs.to_csv(rss, sep='\t', index=False)
        a = list(metadata.index[metadata['ID'] == i])
        loc = a[0]
        metadata.iloc[loc, 6] = vcf.shape[0]
        metadata.iloc[loc, 7] = bedf.shape[0]
        metadata.iloc[loc, 8] = bedrs.shape[0]
    else:
        a = list(metadata.index[metadata['ID'] == i])
        loc = a[0]

        metadata.iloc[loc, 6] = vcf.shape[0]
        metadata.iloc[loc, 7] = bedf.shape[0]

#############get bam size#####iloc[i,9]
source = os.path.join(maindir,config["Directories"]["bam"])
bamsize = []
for i in range(metadata.shape[0]):
    print(source)
    print(metadata.iloc[i,0])
    pathfile = os.path.join(source,metadata.iloc[i,0])
    metadata.iloc[i,9]= round(os.path.getsize(pathfile)/(1024*1024), 2)




metadata.to_csv(statfile, index=False, sep ='\t')

