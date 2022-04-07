import subprocess32 as subprocess
import os
#from config import readConfig_SNP
import configparser
import sys as sys
#now
import shutil
import pandas as pd
path = sys.argv[1]
'''
dicti = readConfig_SNP(path)
maindir = dicti['maindir']
indir = dicti['indir']

processed_ref = dicti['processed_ref']
ref_vcf = dicti['ref_vcf']
outdir = dicti['outdir']
logdir = dicti['logdir']
javapars = dicti['javapars']
met = dicti['metadata']
'''
config = configparser.ConfigParser()
config.read(path)
maindir = config["Directories"]["maindir"]
print(maindir)
indir = os.path.join(maindir,config["Directories"]["data_in"])
print(indir)
processed_ref = os.path.join(maindir,config["Files"]["ref_out1"])
ref_vcf = os.path.join(maindir,config["Files"]["ref_vcf"])
outdir = config["Directories"]["temp_data_out"]
print(outdir)
final_outdir = config["Directories"]["final_data_out"]
print(final_outdir)
logdir = os.path.join(maindir,config["Directories"]["data_log"])
javapars = config["Parameters"]["javaparameters"]
met = os.path.join(maindir,config["Files"]["metadata"])

# pipe ----------------------------------------------------
my_id = sys.argv[2]


process = subprocess.Popen(['java', javapars, '-jar', '$PICARD', 'SortSam', 'I=' + os.path.join(outdir,my_id) + '/' + my_id + '_chop.bam',
                                    'O=' + os.path.join(outdir,my_id) + '/' + my_id + '_sorted.bam',
                                    'SORT_ORDER=coordinate','VALIDATION_STRINGENCY=LENIENT'],
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE,
                                   universal_newlines=True)
stderr, stdout = process.communicate()
print('done')
rc = process.returncode
if (rc == 0):
    print('picard SortSam done with ' + my_id + '\n')
else:
    print('picard SortSam failed with ' + my_id + '\n')
    sys.exit(4)