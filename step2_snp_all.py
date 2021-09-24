import subprocess32 as subprocess
import os
import pandas as pd
import configparser
import sys as sys
from config import readConfig_SNP
import pathlib


jobs = sys.argv[1]
#jobs = 4
path = sys.argv[2]
#path ='C:/Users/Aryuna/Desktop/IB/viadastra_pretty/config.cfg'
# reading config --------------------------
dicti = readConfig_SNP(path)

maindir = dicti['maindir']
print(maindir)
indir = dicti['indir']
print(indir)
processed_ref = dicti['processed_ref']
ref_vcf = dicti['ref_vcf']
outdir = dicti['outdir']
logdir = dicti['logdir']
javapars = dicti['javapars']
met = dicti['metadata']
processing_list = dicti['processing_list']

dir = pathlib.Path(__file__).parent.absolute()
script = os.path.join(dir,'step2_snp_calling.py')

# reading metadata, filtrating ---------------
'''
metadata = pd.read_csv(met,sep='\t')
norna = metadata[metadata["Extra1"] != 'RNA-seq']

idid = norna['ID']
idid = list(idid)
with open(maindir + 'parameters/idid', "w") as outfile:
    outfile.write("\n".join(idid))


subprocess.run(['parallel', '-j', jobs,'python', script,path,'::::',processing_list],
               stdout=subprocess.PIPE, stderr=subprocess.PIPE)
'''

process = subprocess.Popen(['parallel', '-j', jobs,'python', script,path,'::::',processing_list],
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               universal_newlines=True
                               )
stderr, stdout = process.communicate()
tmp_log = os.path.join(maindir,'logs/whole_log')
with open(tmp_log, "w") as log:
    log.write(stdout)

    