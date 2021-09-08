import subprocess32 as subprocess
import os
import pandas as pd
import configparser
import sys as sys

def readConfig_SNP(path, dir):
    config = configparser.ConfigParser()
    config.read(path)

    """
    maindir = '/media/ElissarDisk/ADASTRA/'
    indir = maindir + 'data/'
    processed_ref = maindir + 'processed_ref/genome-norm.fasta'
    ref_vcf = maindir + 'reference/00-common_all.vcf.gz'

    outdir = maindir + 'processed_data/'
    logdir = maindir + 'logs/data_processing/'
    javapars = '-Xmx12G -XX:ParallelGCThreads=4'
    """
    # Читаем некоторые
    #    значения из конфиг. файла.
    maindir = config.get("Directories", "maindir")
    indir = config.get("Directories", "data_in")
    processed_ref = config.get("Files", "ref_out1")
    ref_vcf = config.get("Files", "ref_vcf")
    metadata = config.get("Files","metadata")
    data_out = config.get("Directories", "data_out")
    data_log = config.get("Directories", "data_log")

    javapars = config.get("Parameters", "JavaParameters")

    if (dir == 'maindir'):
        return (maindir)
    elif (dir == 'indir'):
        return (maindir + indir)
    elif (dir == 'processed_ref'):
        return (maindir + processed_ref)
    elif (dir == 'ref_vcf'):
        return (maindir + ref_vcf)
    elif (dir == 'outdir'):
        return (maindir + data_out)
    elif (dir == 'logdir'):
        return (maindir + data_log)
    elif (dir == 'javapars'):
        return (javapars)
    elif (dir == 'metadata'):
        return (metadata)


path = sys.argv[1]
maindir = readConfig_SNP(path,'maindir')
indir = readConfig_SNP(path,'indir')
processed_ref = readConfig_SNP(path,'processed_ref')
ref_vcf = readConfig_SNP(path,'ref_vcf')
outdir = readConfig_SNP(path,'outdir')
logdir = readConfig_SNP(path,'logdir')
javapars = readConfig_SNP(path,'javapars')
met = readConfig_SNP(path,'metadata')

met = maindir + 'parameters/metadata.tsv'
metadata = pd.read_csv(met,sep='\t')
norna = metadata[metadata["Extra1"] != 'RNA-seq']
idid = norna['ID']
idid = list(idid)
with open(maindir + 'parameters/idid', "w") as outfile:
    outfile.write("\n".join(idid))

subprocess.run(['parallel', '-j', '4','python', 'step2_snp_calling.py',path,'::::',maindir + 'parameters/idid'],
               stdout=subprocess.PIPE, stderr=subprocess.PIPE)