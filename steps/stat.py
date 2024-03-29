import subprocess32 as subprocess
import os
import configparser
import sys as sys
import shutil
import pandas as pd
import pysam

path = sys.argv[1]

config = configparser.ConfigParser()
config.read(path)
maindir = config["Directories"]["maindir"]

vcf_calls = config["Directories"]['final_data_out']
vcf_filtered = os.path.join(maindir,config["Directories"]['vcf_filtered'])
rssnps = os.path.join(maindir,config["Directories"]['rssnps'])


met = os.path.join(maindir,config["Files"]["metadata"])
metadata = pd.read_csv('/media/ElissarDisk/ADASTRA/oct23/logs/stats.tsv',sep='\t')
#metadata['Extra1'] = metadata['Extra1'].apply(lambda x: x.lower())
outdir = config["Directories"]['temp_data_out']
tmp_path = os.path.join(outdir,'stats')
statfile = os.path.join(tmp_path, 'stats.tsv')
#metadata['starting_snps']=0
#metadata['filtrated_snps']=0
#metadata['rssnps']=0
metadata['bamsize MB']=0
metadata['readsnum'] = 0
myids=list(metadata['ID'])

pathloc = metadata.columns.get_loc("path")
startloc = metadata.columns.get_loc("starting_snps")
filloc = metadata.columns.get_loc("filtrated_snps")
rsloc = metadata.columns.get_loc("rssnps")
sizeloc = metadata.columns.get_loc("bamsize MB")
nreadsloc = metadata.columns.get_loc("readsnum")


def readsnum_func(BAM):
    readsnum = pysam.view("-c","-F","260", BAM)

    # for row in samfile:
    # print(row)

    # print(samfile)
    return (readsnum)


#get shapes of dataframes
'''
for i in myids:

    startvcf = os.path.join(vcf_calls,i+'.vcf')
    if os.path.exists(startvcf):
        file_rs = open(startvcf, "r")
        line_rs = file_rs.readline()
        n = 0
        while line_rs.startswith("##"):
            n += 1
            line_rs = file_rs.readline()
        file_rs.close()

        vcf = pd.read_csv(startvcf,
                          sep='\t', skiprows=n)


    header_list = ['#chr', 'start', 'end', 'ID', 'ref', 'alt', 'ref_counts', 'alt_counts', 'sample_id']
    strfb = os.path.join(vcf_filtered, i + '.snps.bed')
    if os.path.exists(strfb):
        bedf = pd.read_csv(strfb,sep='\t', names=header_list)
        bedf['sample_id']=i
        print(i)
        bedf.to_csv(strfb, sep ='\t',index=False, header=None)

    if os.path.exists(os.path.join(rssnps, i + '.snps.bed')):
        if (os.stat(os.path.join(rssnps, i + '.snps.bed')).st_size != 0):
            rss = os.path.join(rssnps, i + '.snps.bed')
            bedrs = pd.read_csv(rss,sep='\t', names=header_list)
            bedrs['sample_id']=i
            bedrs.to_csv(rss, sep='\t', index=False,header=None)
            a = list(metadata.index[metadata['ID'] == i])
            loc = a[0]
            metadata.iloc[loc, startloc] = vcf.shape[0]
            metadata.iloc[loc, filloc] = bedf.shape[0]
            metadata.iloc[loc, rsloc] = bedrs.shape[0]
        else:
            a = list(metadata.index[metadata['ID'] == i])
            loc = a[0]

            metadata.iloc[loc, startloc] = vcf.shape[0]
            metadata.iloc[loc, filloc] = bedf.shape[0]
'''
#############get bam size#####iloc[i,9]
source = os.path.join(maindir,config["Directories"]["bam"])



for i in range(metadata.shape[0]):
    try:
        pathfile = os.path.join(source,metadata.iloc[i,0])
        print(pathfile)
        metadata.iloc[i,sizeloc]= round(os.path.getsize(pathfile)/(1024*1024), 2)
        print(metadata.iloc[i,sizeloc])

        metadata.iloc[i,nreadsloc] = readsnum_func(pathfile)
        print(metadata.iloc[i,nreadsloc])
    except:
        pass




metadata.to_csv(statfile, index=False, sep ='\t')

