import subprocess32 as subprocess
import os
import pandas as pd
import configparser
import sys as sys

import pathlib


jobs = sys.argv[1]
#jobs = 4
path = sys.argv[2]
#path ='C:/Users/Aryuna/Desktop/IB/viadastra_pretty/config.cfg'
# reading config --------------------------
config = configparser.ConfigParser()
config.read(path)
maindir = config["Directories"]["maindir"]
print(maindir)
indir = os.path.join(maindir,config["Directories"]["data_in"])
print(indir)
processed_ref = os.path.join(maindir,config["Files"]["ref_out1"])
ref_vcf = os.path.join(maindir,config["Files"]["ref_vcf"])
print(ref_vcf)
#outdir = os.path.join(maindir,config["Directories"]["data_out"])
outdir = config["Directories"]["temp_data_out"]
print(outdir)
final_outdir = config["Directories"]["final_data_out"]
print(final_outdir)
logdir = os.path.join(maindir,config["Directories"]["data_log"])
javapars = config["Parameters"]["javaparameters"]
met = os.path.join(maindir,config["Files"]["metadata"])
processing_list = os.path.join(maindir, config["Files"]["processing_list"])

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
all_log = os.path.join(maindir, 'logs/whole_log')
with open(all_log, "w") as log:
    log.write('STARTING! all'+ '\n')
try:
    process = subprocess.Popen(['parallel', '--memfree','40G','--retry-failed','--joblog',
                                os.path.join(maindir, 'logs/parallel_log'),'-j', jobs,'python', script ,path ,'::::',processing_list],
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE,
                                   universal_newlines=True
                                   )
    stderr, stdout = process.communicate()
    with open(all_log, "a") as log:
        log.write(stdout)
    with open(all_log, "a") as err:
        err.write(stderr)
except KeyboardInterrupt:
    with open(all_log, "a") as log:
        log.write('exception KeyboardInterrupt emerged'+ '\n')
    sys.exit(10)
except:
    with open(all_log, "a") as log:
        log.write('Some exception emerged'+ '\n')
    sys.exit(10)


#--memfree 100G --retries 5