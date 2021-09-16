import subprocess32 as subprocess
import sys
import os
import configparser
from config import readConfig_ref
'''
maindir = '/media/ElissarDisk/ADASTRA/'
outdir = maindir + 'processed_ref/'
logdir = maindir + 'logs/reference_processing/'

indir = maindir + 'reference/hg19.fa'
out1 = outdir +'genome-norm.fasta'
out2 = outdir + 'genome-norm.dict'
loglog = logdir + "stdout"
logerr = logdir + "stderr"
'''

path = sys.argv[1]
#path = "/media/ElissarDisk/ADASTRA/parameters/CONFIG.cfg"
#path ='C:/Users/Aryuna/Desktop/IB/viadastra_pretty/config.cfg'
dicti = readConfig_ref(path)

maindir = dicti['maindir']
print(maindir)
indir = dicti['indir']
print(indir)
logdir = dicti['logdir']
print(logdir)
out1 = dicti['out1']
print(out1)
out2 = dicti['out2']
print(out2)
loglog = dicti['loglog']
print(loglog)
logerr = dicti['logerr']
print(logerr)


#___________
process = subprocess.Popen(['picard', 'NormalizeFasta', 'I='+indir, 'O=' + out1],
                     stdout=subprocess.PIPE,
                     stderr=subprocess.PIPE,
                     universal_newlines=True)
stderr, stdout = process.communicate()
stdout.split('\n')
stderr.split('\n')


with open(loglog, "w") as log:
    log.write(stdout)
with open(logerr, "w") as err:
    err.write(stderr)

#----------
process = subprocess.Popen(['samtools', 'faidx', out1],
                     stdout=subprocess.PIPE,
                     stderr=subprocess.PIPE,
                     universal_newlines=True)
stderr, stdout = process.communicate()
stdout.split('\n')
stderr.split('\n')

with open(loglog, "a") as log:
    log.write(stdout)
with open(logerr, "a") as err:
    err.write(stderr)
#___________
process = subprocess.Popen(['picard', 'CreateSequenceDictionary', 'R='+out1, 'O='+out2],
                     stdout=subprocess.PIPE,
                     stderr=subprocess.PIPE,
                     universal_newlines=True)
stderr, stdout = process.communicate()
stdout.split('\n')
stderr.split('\n')

with open(loglog, "a") as log:
    log.write(stdout)
with open(logerr, "a") as err:
    err.write(stderr)
#___________

